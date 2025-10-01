"""Utilities for ranking variants via an LLM."""
from __future__ import annotations

import json
import pickle
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional
from pdb import set_trace


DEFAULT_DATA_PATH = Path(__file__).resolve().parent / "data" / "aggregated_variants_object.pkl"


@dataclass
class VariantEvidence:
    variant_id: str
    coordinate: Optional[Dict[str, Any]] = None
    clinvar: Optional[Dict[str, Any]] = None
    literature_summary: Optional[str] = None
    evo2_scores: List[Dict[str, Any]] = field(default_factory=list)
    alpha_predictions: List[Dict[str, Any]] = field(default_factory=list)

    def to_prompt_block(self) -> str:
        lines = [f"Variant: {self.variant_id}"]

        if self.coordinate:
            chrom = self.coordinate.get("chrom")
            pos = self.coordinate.get("pos")
            ref = self.coordinate.get("ref")
            alt = self.coordinate.get("alt")
            lines.append(f" Genomic position: chr{chrom}:{pos} {ref}>{alt}")

        if self.clinvar:
            description = self.clinvar.get("description") or "unknown"
            review = self.clinvar.get("review_status") or "not provided"
            traits = ", ".join(self.clinvar.get("trait_names") or []) or "unspecified traits"
            lines.append(f" ClinVar significance: {description} (review: {review})")
            lines.append(f" Reported conditions: {traits}")

        if self.evo2_scores:
            for score in self.evo2_scores:
                delta = score.get("Evo2_deltaScore")
                if delta is not None:
                    lines.append(f" Evo2 delta score: {delta}")
                    break

        if self.alpha_predictions:
            alpha_summary = _summarize_alpha(self.alpha_predictions[0])
            if alpha_summary:
                lines.append(f" AlphaGenome signal summary: {alpha_summary}")

        if self.literature_summary:
            lines.append(" Key literature notes: " + _truncate(self.literature_summary))

        return "\n".join(lines)


def _truncate(text: str, limit: int = 500) -> str:
    cleaned = " ".join(text.split())
    if len(cleaned) <= limit:
        return cleaned
    return cleaned[: limit - 3] + "..."


def _summarize_alpha(record: Dict[str, Any], limit: int = 5) -> str:
    metrics = {
        key: value
        for key, value in record.items()
        if isinstance(value, (int, float)) and key != "_id"
    }
    if not metrics:
        return ""
    ranked = sorted(metrics.items(), key=lambda kv: abs(kv[1]), reverse=True)[:limit]
    return ", ".join(f"{k}={v:.3f}" for k, v in ranked)


def _normalise_variant_id(raw_id: Optional[str]) -> Optional[str]:
    if not raw_id:
        return None
    if raw_id.startswith("chr"):
        return raw_id
    if raw_id[0].isdigit() or raw_id.startswith("X") or raw_id.startswith("Y"):
        return f"chr{raw_id}"
    return raw_id


def _load_aggregated_data(pickle_path: Path) -> Dict[str, Any]:
    with pickle_path.open("rb") as handle:
        return pickle.load(handle)


def _gather_variant_evidence(aggregated: Dict[str, Any]) -> List[VariantEvidence]:
    by_variant: Dict[str, VariantEvidence] = {}

    for entry in aggregated.get("VariantGetterBioMCPAgent") or []:
        if not isinstance(entry, dict):
            continue
        coordinate = entry.get("coordinate") or {}
        chrom = coordinate.get("chrom")
        pos = coordinate.get("pos")
        ref = coordinate.get("ref")
        alt = coordinate.get("alt")
        if not all([chrom, pos, ref, alt]):
            continue
        variant_id = f"chr{chrom}:g.{pos}{ref}>{alt}"
        evidence = by_variant.setdefault(variant_id, VariantEvidence(variant_id=variant_id))
        evidence.coordinate = coordinate
        evidence.clinvar = entry.get("clinvar_data") or evidence.clinvar
        summary = entry.get("lit_article_summary")
        if summary:
            evidence.literature_summary = summary

    for entry in aggregated.get("Evo2Agent") or []:
        if not isinstance(entry, dict):
            continue
        variant_id = _normalise_variant_id(entry.get("_id"))
        if not variant_id:
            continue
        evidence = by_variant.setdefault(variant_id, VariantEvidence(variant_id=variant_id))
        evidence.evo2_scores.append(entry)

    alpha_entries: Iterable[Any]
    alpha_raw = aggregated.get("AlphaGenomeAgent")
    if isinstance(alpha_raw, list):
        alpha_entries = alpha_raw
    elif isinstance(alpha_raw, dict):
        alpha_entries = [alpha_raw]
    else:
        alpha_entries = []

    for entry in alpha_entries:
        if not isinstance(entry, dict):
            continue
        variant_id = _normalise_variant_id(entry.get("_id"))
        if not variant_id:
            continue
        evidence = by_variant.setdefault(variant_id, VariantEvidence(variant_id=variant_id))
        evidence.alpha_predictions.append(entry)

    return sorted(by_variant.values(), key=lambda item: item.variant_id)


def _build_prompt(evidence_blocks: List[VariantEvidence]) -> str:
    header = (
        "You are a clinical genomics expert. Review the supplied variant evidence "
        "and return a JSON array that ranks the variants from highest to lowest "
        "clinical concern (1 = highest priority for follow-up). Use conservative "
        "language and justify each ranking."
    )

    expectations = (
        "JSON schema: [ {\"variant_id\": str, \"rank\": int, "
        "\"priority_reason\": str, \"overall_assessment\": str} ]. "
        "The array must be sorted by ascending rank."
    )

    blocks = "\n\n".join(item.to_prompt_block() for item in evidence_blocks)

    instructions = (
        "When forming the rationale consider: ClinVar annotations, literature "
        "notes, Evo2 delta scores (more negative may imply larger effect), and "
        "AlphaGenome signals. If evidence indicates low risk, explain that."
    )

    return f"{header}\n\n{expectations}\n\nEvidence:\n{blocks}\n\n{instructions}"


def _parse_llm_rankings(raw_response: str) -> Optional[List[Dict[str, Any]]]:
    try:
        return json.loads(raw_response)
    except json.JSONDecodeError:
        pass

    start = raw_response.find("[")
    end = raw_response.rfind("]")
    if start != -1 and end != -1 and start < end:
        snippet = raw_response[start : end + 1]
        try:
            return json.loads(snippet)
        except json.JSONDecodeError:
            return None
    return None


def rank_variants(
    pickle_path: Path = DEFAULT_DATA_PATH,
    *,
    model: Optional[str] = None,
    temperature: Optional[float] = None,
) -> Dict[str, Any]:
    aggregated = _load_aggregated_data(pickle_path)
    # set_trace()
    evidence_blocks = _gather_variant_evidence(aggregated)
    if not evidence_blocks:
        raise ValueError("No variant evidence found in aggregated data.")

    prompt = _build_prompt(evidence_blocks)

    llm_query = _import_llm_query()

    llm_kwargs: Dict[str, Any] = {}
    if model is not None:
        llm_kwargs["model"] = model
    if temperature is not None:
        llm_kwargs["temperature"] = temperature

    if llm_kwargs:
        response = llm_query(prompt, **llm_kwargs)
    else:
        response = llm_query(prompt)

    parsed = _parse_llm_rankings(response)

    return {
        "prompt": prompt,
        "raw_response": response,
        "parsed_rankings": parsed,
    }


def _import_llm_query() -> Callable[..., str]:
    try:  # prefer llm_query when available
        from llm_utils import llm_query as imported  # type: ignore
    except ImportError:
        from llm_utils import query_llm as imported  # type: ignore
    return imported


if __name__ == "__main__":  # manual usage
    result = rank_variants()
    print(json.dumps(result["parsed_rankings"], indent=2) if result["parsed_rankings"] else result["raw_response"])
