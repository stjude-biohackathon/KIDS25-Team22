from llm_utils import query_llm
import prompts
import os
import subprocess
import utils
from config import LLM_CONFIG
import re
import requests
from duckduckgo_search import DDGS
import xml.etree.ElementTree as ET
import argparse
import json
from bs4 import BeautifulSoup
from pdb import set_trace

# add persistent context memory

# Principal Investigator agent
class OrchestratorAgent:
    def __init__(
        self,
        verbose=False,
        data_retrieval_agent=None,
        variant_getter_biomcp_agent=None,
        variant_getter_clinvar_agent =None,
        article_getter_agent=None,
        alpha_genome_agent=None,
        evo2_agent=None,
        gpn_agent=None,
        ranking_agent=None,
        reporter_agent=None,
        variant_aggregation_agent=None,
        conda_env=None,
    ):
        self.verbose = verbose
        self.data_retrieval_agent = data_retrieval_agent
        self.variant_getter_agent = variant_getter_biomcp_agent
        self.variant_getter_clinvar_agent = variant_getter_clinvar_agent
        self.article_getter_agent = article_getter_agent
        self.alpha_genome_agent = alpha_genome_agent
        self.evo2_agent = evo2_agent
        self.gpn_agent = gpn_agent
        self.ranking_agent = ranking_agent
        self.reporter_agent = reporter_agent
        self.variant_aggregation_agent = variant_aggregation_agent
        self.conda_env = conda_env
        self.last_variant_results = None

    def coordinate(self, file_path, phenotype):
        if self.variant_getter_agent is None:
            raise ValueError("Variant getter agent is not configured.")

        coordinates = self._prepare_variant_coordinates(file_path, phenotype)
        if self.variant_aggregation_agent:
            result = self.variant_aggregation_agent.gather_variants(
                variant_getter=self.variant_getter_agent,
                alpha_genome_agent=self.alpha_genome_agent,
                evo2_agent=self.evo2_agent,
                gpn_agent=self.gpn_agent,
                coordinates=coordinates,
                phenotype=phenotype,
            )
            self.last_variant_results = result
            return result

        result = {
            "VariantGetterBioMCPAgent": self.variant_getter_agent.search_variants(coordinates)
        }
        self.last_variant_results = result
        return result

    def _prepare_variant_coordinates(self, file_path, phenotype):
        self.verify_inputs(file_path, phenotype)
        coordinates = utils.extract_coordinates(file_path)
        if not coordinates:
            raise ValueError("No variant coordinates could be extracted from the provided file.")
        return coordinates

    def verify_inputs(self, file_path, phenotype):
        if not file_path:
            raise ValueError("A variant file path must be provided.")

        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Variant file not found: {file_path}")

        if not file_path.lower().endswith(".vcf"):
            raise ValueError("Currently only VCF files are supported for variant input.")

        if phenotype is None or not str(phenotype).strip():
            raise ValueError("A phenotype description must be provided.")

        with open(file_path, "r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                columns = line.rstrip().split("\t")
                if len(columns) < 5:
                    continue

                chrom, pos, _id, ref, alt = columns[:5]
                if not re.match(r"^(chr)?[0-9XYM]+$", chrom, re.IGNORECASE):
                    continue
                if not pos.isdigit():
                    continue
                if not re.fullmatch(r"[ACGTN]+", ref, re.IGNORECASE):
                    continue
                if not re.fullmatch(r"[ACGTN,]+", alt, re.IGNORECASE):
                    continue
                return True

        raise ValueError("Input VCF does not contain valid variant records.")




    
class VariantGetterBioMCPAgent:
    def __init__(self, verbose=True, assembly="GRCh38"):
        self.verbose = verbose
        self.assembly = assembly

    def search_variants(self, coordinates):
        """Search for variant details using genomic coordinates."""
        variant_data = {}

        for record in coordinates[:1]:  # Limit to first 10 for testing
            chrom = record["chrom"]
            pos = record["pos"]
            ref = record["ref"]
            alt = record["alt"]
            var_id = chrom + ':g.' + str(pos) + ref + '>' + alt
            # var_id = 'chr11:g.32392032G>A'
            var_id = 'chr11:g.32413578G>A'
            print(f"var_id: {var_id}")
            command = [
                "biomcp",
                "variant",
                "get",
                "-j",
                var_id,
            ]

            try:
                if self.verbose:
                    print(f"Searching variants for {var_id}")

                result = subprocess.run(
                    command,
                    capture_output=True,
                    text=True,
                    timeout=60,
                )

                # set_trace()

                if result.returncode == 0:
                    formatted_output = utils.format_biomcp_variant_output(result.stdout)
                    if self.verbose:
                        print(formatted_output)
                    variant_data = {
                            "coordinate": record,
                            "output": result.stdout.strip(),
                            "output_table": formatted_output,
                        }
                    
                else:
                    variant_data = {
                            "coordinate": record,
                            "error": result.stderr.strip() or "biomcp returned a non-zero exit code",
                        }
                # set_trace()
            except Exception as exc:
                variant_data.append(
                    {
                        "coordinate": record,
                        "error": str(exc),
                    }
                )

            # Send the Clinvar IDs to the ClinVar agent to fetch more data
            # clinvar_id = query_llm(f"What is the ClinVar ID (Variant Id in section Clinvar) fetched for the variant {var_id} from {result}? If there are multiple, list them all separated by commas. If there is none, respond with 'None'.")
            # print(f"ClinVar ID from LLM: {clinvar_id}")
            # if clinvar_id.lower() == 'none':
            #     continue
            # In the result.stdout, find the line that starts with "variant_id" and extract the Variant Id
            # Extract variant_id from result.stdout
            variant_id_match = re.search(r'"variant_id":\s*(\d+)', result.stdout)
            variant_id = variant_id_match.group(1) if variant_id_match else None
            # set_trace()
            variant_id = int(variant_id)
            # if variant_id:
            #     variant_data[-1]["variant_id"] = variant_id


            # set_trace()
            if variant_id is None:
                return {"error": "No ClinVar ID found in BioMCP output."}
            else:
                clinvar_agent = VariantGetterClinVarAgent(verbose=self.verbose)
                clinvar_results = clinvar_agent.fetch_clinvar_data(variant_id)

            variant_data.update({"clinvar_data": clinvar_results})
            # set_trace()
            # if clinvar_results:
            #     variant_data[-1]["clinvar_data"] = clinvar_results[0]

        return variant_data
      
# Gets output from biomcp variant get to fetch clinvar ids and obtains clinvar data from ncbi    
class VariantGetterClinVarAgent:
    def __init__(self, verbose=False):
        self.verbose = verbose

    def fetch_clinvar_data_old(self, coordinates):
        clinvar_data = []
        for record in coordinates:
            chrom = record["chrom"]
            pos = record["pos"]
            ref = record["ref"]
            alt = record["alt"]
            url = (
                f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
                f"db=clinvar&term={chrom}[CHR]+AND+{pos}[CHRPOS]+AND+{ref}[REF]+AND+{alt}[ALT]"
            )
            try:
                if self.verbose:
                    print(f"Fetching ClinVar data for {chrom}:{pos} {ref}>{alt}")

                response = requests.get(url, timeout=30)
                response.raise_for_status()

                root = ET.fromstring(response.content)
                id_list = root.find("IdList")
                ids = [id_elem.text for id_elem in id_list.findall("Id")] if id_list is not None else []

                clinvar_data.append(
                    {
                        "coordinate": record,
                        "clinvar_ids": ids,
                    }
                )
            except Exception as exc:
                clinvar_data.append(
                    {
                        "coordinate": record,
                        "error": str(exc),
                    }
                )

        return clinvar_data
    
    def fetch_clinvar_data(self, variant_id):
        clinvar_data = {}
        print("In fetch_clinvar_data")
        # for record in variant_id:
        # chrom = record["chrom"]
        # pos = record["pos"]
        # ref = record["ref"]
        # alt = record["alt"]
        url = (
            # f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
            # f"db=clinvar&term={chrom}[CHR]+AND+{pos}[CHRPOS]+AND+{ref}[REF]+AND+{alt}[ALT]"
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={variant_id}&retmode=json"
        )
        try:
            if self.verbose:
                print(f"Fetching ClinVar data for {variant_id}")

            response = requests.get(url, timeout=30)
            # response.raise_for_status()
            # print(f"ClinVar response: {response.text}")
            response = response.json()
            # from the json format, extract the "description" field from the key "germline_classification"
            description = response["result"][str(variant_id)]["germline_classification"]["description"]
            review_status = response["result"][str(variant_id)]["germline_classification"]["review_status"]
            # trait_name = response["result"][str(variant_id)]["germline_classification"]["trait_set"][-1]['trait_name']
            trait_names = [trait['trait_name'] for trait in response["result"][str(variant_id)]["germline_classification"]["trait_set"]]
            clinvar_data = {
                "variant_id": variant_id,
                "description": description,
                "review_status": review_status,
                "trait_names": trait_names,
            }
            print(f"ClinVar data: {clinvar_data}"
            )

            # set_trace()

        #     root = ET.fromstring(response.content)
        #     id_list = root.find("IdList")
        #     ids = [id_elem.text for id_elem in id_list.findall("Id")] if id_list is not None else []

        #     detailed_records = []
        #     set_trace()
        #     for clinvar_id in ids:
        #         fetch_url = (
        #             f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        #             f"db=clinvar&id={clinvar_id}&retmode=xml"
        #         )
        #         fetch_response = requests.get(fetch_url, timeout=30)
        #         fetch_response.raise_for_status()
        #         detailed_records.append(fetch_response.text)

        #     clinvar_data.append(
        #         {
        #             "variant_id": variant_id,
        #             "clinvar_ids": ids,
        #             "detailed_records": detailed_records,
        #         }
        #     )
        except Exception as exc:
            clinvar_data.append(
                {
                    "variant_id": variant_id,
                    "error": str(exc),
                }
            )

        return clinvar_data

class DataRetrievalAgent:
    def __init__(self, verbose=False):
        self.verbose = verbose

    def retrieve(self, *args, **kwargs):
        return {}


class ArticleGetterBioMCPAgent:
    def __init__(self, verbose=False):
        self.verbose = verbose

    def fetch_articles(self, *args, **kwargs):
        return []


class AlphaGenomeBioMCPAgent:
    def __init__(self, verbose=False):
        self.verbose = verbose

    def analyze(self, *args, **kwargs):
        return None


class Evo2Agent:
    def __init__(self, verbose=False):
        self.verbose = verbose

    def run(self, *args, **kwargs):
        return None


class GPNAgent:
    def __init__(self, verbose=False):
        self.verbose = verbose

    def score(self, *args, **kwargs):
        return None


class RankingAgent:
    def __init__(self, verbose=False):
        self.verbose = verbose

    def rank(self, *args, **kwargs):
        return []


class ReporterAgent:
    def __init__(self, verbose=False):
        self.verbose = verbose

    def create_report(self, *args, **kwargs):
        return ""


class VariantAggregationAgent:
    def __init__(self, verbose=False):
        self.verbose = verbose
        self.latest_results = None

    def gather_variants(
        self,
        *,
        variant_getter=None,
        alpha_genome_agent=None,
        evo2_agent=None,
        gpn_agent=None,
        coordinates=None,
        phenotype=None,
    ):
        if coordinates is None:
            raise ValueError("Coordinates must be provided for aggregation.")

        aggregated = {}

        if variant_getter is not None:
            aggregated["VariantGetterBioMCPAgent"] = variant_getter.search_variants(coordinates)

        if alpha_genome_agent is not None and hasattr(alpha_genome_agent, "analyze"):
            aggregated["AlphaGenomeBioMCPAgent"] = alpha_genome_agent.analyze(
                coordinates=coordinates,
                phenotype=phenotype,
            )

        if evo2_agent is not None and hasattr(evo2_agent, "run"):
            aggregated["Evo2Agent"] = evo2_agent.run(
                coordinates=coordinates,
                phenotype=phenotype,
            )

        if gpn_agent is not None and hasattr(gpn_agent, "score"):
            aggregated["GPNAgent"] = gpn_agent.score(
                coordinates=coordinates,
                phenotype=phenotype,
            )

        self.latest_results = aggregated
        return aggregated
