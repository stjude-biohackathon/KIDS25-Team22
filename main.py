# on laptop: conda activate /mnt/c/Users/tnandi/Downloads/ai_codes/ai_py3p12_env OR conda activate llm_env
# on Polaris: module load conda;conda activate /lus/grand/projects/GeomicVar/tarak/ai_codes/ai_py3p12_env

# Now using deepseek and qwen reasoning models hosted on Sophia/Polaris using Ollama. Will move to the ALCF inference endpoints when they make these models available
# the deepseek 70b and qwen 32B models work fine, but the 671b model throws error related to the number of experts being used is more than that allowed by the ollama llama.cpp installation

# to check if the ollama server is active try: curl --noproxy localhost http://localhost:11434/api/generate -d '{"model": "deepseek-r1:70b", "prompt": "Explain polygenic risk scores.", "temperature": 0.3}'
# OR for the qwen model: curl --noproxy localhost http://localhost:11434/api/generate -d '{"model": "qwq:latest", "prompt": "Explain polygenic risk scores.", "temperature": 0.3}'
# curl --noproxy localhost http://localhost:11434/api/generate -d '{"model": "gpt-oss:20b", "prompt": "Explain polygenic risk scores.", "temperature": 0.3}'

# curl --noproxy localhost http://localhost:11434/api/generate -d '{"model": "codellama:latest", "prompt": "Write optimized matrix vector multiplication CUDA code without using cuBLAS", "temperature": 0.3}' -o output.jsonl

# to list all models available (but may not be currently active): curl http://localhost:11434/api/tags | jq '.models[] | {name, parameter_size: .details.parameter_size, quant: .details.quantization_level}'
# 

from agents import (
    OrchestratorAgent,
    DataRetrievalAgent,
    VariantGetterBioMCPAgent,
    VariantGetterClinVarAgent,
    ArticleGetterBioMCPAgent,
    AlphaGenomeBioMCPAgent,
    Evo2Agent,
    GPNAgent,
    RankingAgent,
    ReporterAgent,
    VariantAggregationAgent,

)
import config
import utils
import argparse
import os

# Argument parser for topic input
parser = argparse.ArgumentParser(description="Carry out variant prioritization from VCF.")
parser.add_argument(
    "--vcf_file",
    type=str,
    default="/Users/tnandi/Downloads/stjudes_hackathon/pediatric_cancer_variants.vcf",
    help="Path to the VCF containing variants to analyze.",
)
parser.add_argument(
    "--phenotype",
    type=str,
    default="breast cancer",
    help="Phenotype or disease term associated with the sample.",
)
parser.add_argument(
    "--conda_env",
    type=str,
    default="/Users/tnandi/Downloads/agents/agentic_lab/agentic_lab_env",
    help="Path to conda environment for code execution (e.g., /path/to/env)",
)


def main(file_path=None, phenotype=None):
    args = parser.parse_args()

    vcf_file = file_path or args.vcf_file
    phenotype_term = phenotype or args.phenotype

    # # Process files directory if provided
    # files_dir_content = ""
    # if args.files_dir:
    #     print(f"Exploring files directory: {args.files_dir}")
    #     files_dir_content = utils.explore_files_directory(args.files_dir)
    #     if files_dir_content:
    #         print(f"Successfully explored files directory")
    #     else:
    #         print("Warning: Could not explore files directory")
    
    # Initialize agents
    data_retrieval_agent = DataRetrievalAgent(verbose=True)
    variant_getter_biomcp_agent = VariantGetterBioMCPAgent(verbose=True)
    variant_getter_clinvar_agent = VariantGetterClinVarAgent(verbose=True)
    article_getter_agent = ArticleGetterBioMCPAgent(verbose=True)
    alpha_genome_agent = AlphaGenomeBioMCPAgent(verbose=True)
    evo2_agent = Evo2Agent(verbose=True)
    gpn_agent = GPNAgent(verbose=True)
    ranking_agent = RankingAgent(verbose=True)
    reporter_agent = ReporterAgent(verbose=True)
    variant_aggregation_agent = VariantAggregationAgent(verbose=True)

    
    
    orchestrator_agent = OrchestratorAgent(
        verbose=True,
        data_retrieval_agent=data_retrieval_agent,
        variant_getter_biomcp_agent=variant_getter_biomcp_agent,
        variant_getter_clinvar_agent=variant_getter_clinvar_agent,
        article_getter_agent=article_getter_agent,
        alpha_genome_agent=alpha_genome_agent,
        evo2_agent=evo2_agent,
        gpn_agent=gpn_agent,
        ranking_agent=ranking_agent,
        reporter_agent=reporter_agent,
        variant_aggregation_agent=variant_aggregation_agent,
        # files_dir_content=files_dir_content,
        # mode=args.mode,
        conda_env=args.conda_env,
    )
    
    if vcf_file is None:
        raise ValueError("A VCF file must be provided via argument or function call.")

    if phenotype_term is None:
        raise ValueError("A phenotype must be provided via argument or function call.")

    return orchestrator_agent.coordinate(vcf_file, phenotype_term)

if __name__ == "__main__":
    main()
