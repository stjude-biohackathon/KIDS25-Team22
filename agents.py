from llm_utils import query_llm
#import prompts
import os, sys
import subprocess
import utils
from config import LLM_CONFIG
import re
import requests
#from duckduckgo_search import DDGS
import xml.etree.ElementTree as ET
import argparse
import json
from bs4 import BeautifulSoup
from pdb import set_trace

from Bio import SeqIO
import gzip
import base64
from io import BytesIO
#import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#import seaborn as sns
#from sklearn.metrics import roc_auc_score
import requests
from dotenv import load_dotenv

from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome.models import variant_scorers

from liftover import get_lifter

load_dotenv()
alpha_api_key = os.getenv("ALPHAGENOME_API_KEY")


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
        alpha_genome_biomcp_agent=None,
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
        set_trace()
        result = {
            "VariantGetterBioMCPAgent": self.variant_getter_agent.search_variants(coordinates, phenotype),
            "AlphaGenomeAgent": self.alpha_genome_agent.predict_variants_effects(coordinates),
            "Evo2Agent": self.evo2_agent.getEvo2score(coordinates)
        }
        set_trace()
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

    def search_variants(self, coordinates, phenotype):
        """Search for variant details using genomic coordinates."""
        
        variant_data_list = []
        for record in coordinates[:2]:  # Limit to first 10 for testing
            variant_data = {}
            chrom = record["chrom"]
            pos = record["pos"]
            ref = record["ref"]
            alt = record["alt"]
            var_id = 'chr' + chrom + ':g.' + str(pos) + ref + '>' + alt
            # var_id = 'chr11:g.32392032G>A'
            # var_id = 'chr11:g.32413578G>A'
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
                print("Completed biomcp variant get call")
                if result.returncode == 0:
                    formatted_output = utils.format_biomcp_variant_output(result.stdout)
                    # if self.verbose:
                    #     print(formatted_output)
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
            if variant_id is None:
                return {"error": "No ClinVar ID found in BioMCP output."}
            else:
                print("Fetching clinVar data")
                variant_id = int(variant_id)
                clinvar_agent = VariantGetterClinVarAgent(verbose=self.verbose)
                clinvar_results = clinvar_agent.fetch_clinvar_data(variant_id)

            variant_data.update({"clinvar_data": clinvar_results})
            # set_trace()

            # Extract genename from variant_data['output'] string (e.g., WT1 from "genename": "WT1")
            gene_match = re.search(r'"genename":\s*"([^"]+)"', variant_data['output'])
            gene_name = gene_match.group(1) if gene_match else "Unknown"
            # variant_data.update({"gene_name": gene_name})

            # Create an instance of the ArticleGetterBioMCPAgent; pass the gene name and phenotype for it to retrieve articles and return in json format for ArticleProcessorAgent to process
            print("fetching articles using ArticleGetterBioMCPAgent.fetch_articles()")
            article_getter = ArticleGetterBioMCPAgent(verbose=self.verbose)
            articles = article_getter.fetch_articles(gene=gene_name, phenotype=phenotype)

            # Pass articles to ArticleProcessorAgent to process and return a summary
            article_processor = ArticleProcessorAgent(verbose=self.verbose)
            # set_trace()
            # articles is in json format
            # Process articles using the process function in the ArticleProcessorAgent class
            article_summaries = article_processor.process(articles=articles, gene=gene_name, phenotype=phenotype)

            # Collect article summary for all variants 
            variant_data.update({"lit_article_summary": article_summaries})
            # Save as pickl for Franz (his agent accepts a lit of dicts as input)
            
            variant_data_list.append(variant_data)
        #set_trace()    
        #with open('my_object.pkl', 'wb') as file:
        #    pickle.dump(variant_data_list, file)
        return variant_data_list
      
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
            # print(f"ClinVar data: {clinvar_data}")

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

    def fetch_clinvar_data_batch(self, variant_id):
        pass

class DataRetrievalAgent:
    def __init__(self, verbose=False):
        self.verbose = verbose

    def retrieve(self, *args, **kwargs):
        return {}

class ArticleProcessorAgent:
    def __init__(self, verbose=False):
        self.verbose = verbose
    # Takes as input the articles in json format fetched by ArticleGetterBioMCPAgent and processes them to return a summary
    def process(self, *args, **kwargs):
        print("In ArticleProcessorAgent process")
        articles = kwargs.get("articles")
        gene = kwargs.get("gene")
        phenotype = kwargs.get("phenotype")
        lit_articles = articles['articles']
        # Get the number of articles 
        num_lit_articles = len(lit_articles)
        print(f"Number of articles from the literature: {num_lit_articles}")

        if lit_articles is None or "error" in lit_articles:
            return {"error": "No lit_articles to process or error in fetching lit_articles."}
        if len(lit_articles) == 0:
            return {"error": "No lit_articles found to process."}
        summaries = []
        article_dict = []
        for lit_article in lit_articles[:30]:
            title = lit_article.get("title", "No title")
            # print(f"Title: {title}")
            abstract = lit_article.get("abstract")
            # print(f"Abstract: {abstract}")
            # if abstract == "" or abstract is None:
            #     abstract = "No abstract"
            # else:
            #     summary_prompt = f"Summarize the following article titled '{title}' solely using the abstract: {abstract} in 2-3 sentences focusing on its relevance to the gene and phenotype provided. You need have to access the full text of the article."
            #     print(f"Querying LLM with prompt: {summary_prompt}")
            #     # set_trace()
            #     summary = query_llm(summary_prompt, model=LLM_CONFIG["default_model"], temperature=0.3)
            #     summaries.append({
            #         "title": title,
            #         "summary": summary
            #     })
            # Create summary of all articles together
            # Create a dictionary of abstract and title
            article_dict.append({
                "title": title,
                "abstract": abstract
            })
        # print(f"Article dict: {article_dict}")
        # Pass the article_dict to the LLM to get a summary for all articles to be output in json format with title and summary as keys and values
        #summary_prompt = f"Summarize the articles in {article_dict} in 2-3 sentences focusing on its relevance to the gene and phenotype provided. You need have to access the full text of the article if abstract is not available. Return summary for all articles in json format with title and summary as keys and values"
        # set_trace()

        summary_prompt = (
            f"Summarize the articles in {article_dict} focusing on their relevance to the {gene} and {phenotype} provided. "
            "You need have to access the full text of the article if abstract is not available. "
            "Comment on whether the articles suggest a role for the gene in the phenotype."
            "If you use citations, create a bibliography near the end of the response"
        )
        print(f"Querying LLM with given prompt")
        summary = query_llm(summary_prompt, model=LLM_CONFIG["default_model"], temperature=0.3)
        # set_trace()
            
        
        return summary
        


class ArticleGetterBioMCPAgent:
    def __init__(self, verbose=False):
        self.verbose = verbose

    def fetch_articles(self, *args, **kwargs):
        # Execute biomcp article search --gene BRAF --limit 5 and return a json file containing all the articles found
        gene = kwargs.get("gene")
        phenotype = kwargs.get("phenotype")
        if gene is None:
            raise ValueError("Gene must be provided to fetch articles.")
        command = [
            "biomcp",
            "article",
            "search",
            "--gene",
            gene,
            "--disease",
            phenotype,
            "-j",
        ]
        try:
            if self.verbose:
                print(f"Fetching articles for gene {gene} and phenotype {phenotype}")

            result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                timeout=60,
            )
            # set_trace()
            if result.returncode == 0:
                articles = json.loads(result.stdout)
                if self.verbose:
                    print(f"Fetched {len(articles)} articles for gene {gene}")
                # set_trace()
                return articles #json
                
            else:
                if self.verbose:
                    print(f"Error fetching articles: {result.stderr.strip()}")
                return {"error": result.stderr.strip() or "biomcp returned a non-zero exit code"}
        except Exception as exc:
            if self.verbose:
                print(f"Exception fetching articles: {str(exc)}")
            return {"error": str(exc)}  
                
        return []


class AlphaGenomeBioMCPAgent:
    def __init__(self, verbose=False):
        self.verbose = verbose

    def analyze(self, *args, **kwargs):
        return None


class AlphaGenomeAgent:
    def __init__(self, verbose=False):
        self.verbose = verbose
        alpha_api_key = os.getenv("ALPHAGENOME_API_KEY")
        self.alphaGenomeModel = dna_client.create(alpha_api_key)

    def predict_variants_effects(
            self, coordinates, phenotype=None,
    ):
        """Search for variant details using genomic coordinates."""
        variant_data = {}
        converter = get_lifter('hg19', 'hg38', one_based=True)
        for record in coordinates[:2]:  # Limit to first 10 for testing
            chrom = record["chrom"]
            pos = record["pos"]
            ref = record["ref"]
            alt = record["alt"]
            variant_id = chrom + ':' + str(pos) + ":" + ref + '>' + alt
            
            newCoords = converter[chrom][pos][0]
            chrom_hg38 = newCoords[0]
            pos_hg38 = newCoords[1]            
            variant = genome.Variant(
                chromosome=chrom_hg38,
                position=pos_hg38,
                reference_bases=ref,  # Can differ from the true reference genome base.
                alternate_bases=alt,
            )

            interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_16KB)

            try:
                if self.verbose:
                    print(f"Searching variant effects for {variant_id} using AlphaGenome")
                    variant_scores = self.alphaGenomeModel.score_variant(
                        interval=interval, variant=variant, 
                        variant_scorers=list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())
                    )

                    df_scores = variant_scorers.tidy_scores(variant_scores)
                    top_scores = df_scores.groupby(['output_type']).agg('first')
                    top_scores_dict = dict(zip(top_scores.index, top_scores['quantile_score']))
                    variant_data[variant_id] = top_scores_dict

            except Exception as exc:
                score_assays = [
                    'ATAC', 'CAGE', 'CHIP_HISTONE', 'CHIP_TF', 'CONTACT_MAPS', 'DNASE', 'PROCAP', 
                    'RNA_SEQ', 'SPLICE_JUNCTIONS', 'SPLICE_SITES', 'SPLICE_SITE_USAGE'
                ]

                variant_data [variant_id] = {assay: 0.0 for assay in score_assays}
            
        print(variant_data)
        return variant_data
            



class Evo2Agent:
    def __init__(self, verbose=False):
        self.verbose = verbose

    def calculateSequenceScore(self,input_sequence, base64_data) -> float:
        """
        Calculates the total logit score for a DNA sequence from the Evo2 /forward API response.

        Args:
            input_sequence: The original DNA sequence provided to the API (e.g., "ATGC").
            base64_data: The Base64 encoded string from the 'data' field of the API response.

        Returns:
            The total logit score of the sequence as a float.
        """

        # Step 1: Decode the Base64 string and load the NPZ data
        try:
            decoded_data = base64.b64decode(base64_data)
            npz_file = np.load(BytesIO(decoded_data))
            #print(list(npz_file.keys()))
            logits = npz_file['output_layer.output']
        except Exception as e:
            print(f"Error processing input data: {e}")
            return 0.0

        # Step 2: Define the mapping from DNA characters to their logit indices
        # From the NVIDIA NIM docs: A=65, C=67, G=71, T=84
        char_to_index = {'A': 65, 'C': 67, 'G': 71, 'T': 84}

        # Step 3: Iterate through the sequence to sum the relevant logit scores
        total_score = 0.0

        # We loop from the second character (index 1) to the end of the sequence.
        # The score for each character is found in the logit predictions from the *previous* step.
        for i in range(1, len(input_sequence)):
            # The character we are currently scoring (e.g., 'T' in "ATGC")
            current_char = input_sequence[i]

            # Ensure the character is valid before proceeding
            if current_char not in char_to_index:
                #print(f"Warning: Skipping unrecognized character '{current_char}' in sequence.")
                continue

            logit_index = char_to_index[current_char]

            # Logits at position i-1 contain the predictions for the character at position i.
            # We use np.squeeze() to remove the unnecessary middle dimension (batch size of 1).
            previous_step_logits = np.squeeze(logits[i - 1])

            # Get the specific logit score for our character from the 512-length vocabulary array.
            char_score = previous_step_logits[logit_index]

            total_score += char_score
        return total_score


    def runEvo2(self,dna_sequence):
        load_dotenv()

        MODEL_ENDPOINT = "https://evo2-40b-predictor-austaadmin-stju-b700e7ae.ai-application.stjude.org"
        #there are 2 different endpoints, one is for sequence generation, hence the name "generate"
        # the other is for forward passing for likelihood calculation, extracting embeddings, etc, hence the name "forward"
        #API_PATH = "/biology/arc/evo2/generate"
        API_PATH = "/biology/arc/evo2/forward"
        
        AUTH_TOKEN = os.getenv("PCAI_EVO2_TOKEN")
        #print(AUTH_TOKEN)
        # Your BODY_DATA dictionary as a Python object
        #the following input is example to generate dna sequence given a prompt as a sequence
        BODY_DATA = {
            "sequence": dna_sequence,
            "output_layers": ['output_layer']
        }

        headers = {
            #"Content-Type": "application/json",
            "Authorization": f"Bearer {AUTH_TOKEN}",
        }

        try:
            #print(json.dumps(BODY_DATA))
            response = requests.post(f"{MODEL_ENDPOINT}{API_PATH}",
                                    headers=headers,
                                    json=BODY_DATA,
                                    verify=False)

            response.raise_for_status()

            # Print the full response
            # Note than the following lines are for using forward API to get likelihood score,
            # If you want to use the generate API, you need to adjust accordingly, i.e. no need to calculate the score again
            print("Status Code:", response.status_code)
        except requests.exceptions.RequestException as e:
            print(f"An error occurred: {e}")
            print(f"Server's detailed error message: {response.text}")

        return response.json()['data']

    def getEvo2score(self, coordinates):
        evo2_results = {}

        variant_data = {}
        for record in coordinates[:1]:  # Limit to first 10 for testing
            chrom = record["chrom"]
            coordinates = record["pos"]
            ref = record["ref"]
            alt = record["alt"]
            variant_id = chrom + ':' + str(coordinates) + ":" + ref + '>' + alt
            
            #print(chr,coordinates,ref,alt)
            WINDOW_SIZE = 8192

            # Read the reference genome sequence of chromosome 17
            with gzip.open(os.path.join('hg19','chr'+str(chrom)+'.fa.gz'), "rt") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    seq_chr = str(record.seq)
                    break

            def parse_sequences(pos, ref, alt):
                """
                Parse reference and variant sequences from the reference genome sequence.
                """
                p = pos-1  # Convert to 0-indexed position
                full_seq = seq_chr

                ref_seq_start = max(0, p - WINDOW_SIZE//2)
                ref_seq_end = min(len(full_seq), p + WINDOW_SIZE//2)
                ref_seq = seq_chr[ref_seq_start:ref_seq_end]
                snv_pos_in_ref = min(WINDOW_SIZE//2, p)
                #print(snv_pos_in_ref)
                var_seq = ref_seq[:snv_pos_in_ref] + alt + ref_seq[snv_pos_in_ref+1:]
                #print(var_seq)
                # Sanity checks
                ##assert len(var_seq) == len(ref_seq)
                #print(ref_seq[snv_pos_in_ref-3:snv_pos_in_ref+3],ref)
                ##assert ref_seq[snv_pos_in_ref] == ref
                ##assert var_seq[snv_pos_in_ref] == alt

                return ref_seq, var_seq

            ref_seq, var_seq = parse_sequences(coordinates,ref,alt)
            # set_trace()

            varScore = self.calculateSequenceScore(var_seq, self.runEvo2(var_seq))
            refScore = self.calculateSequenceScore(ref_seq, self.runEvo2(ref_seq))
            print("refScore",refScore)
            print("varScore",varScore)
            finalScore={"Variant":str(chr)+":"+str(coordinates)+ref+">"+alt, "Evo2_deltaScore":float(varScore-refScore)}
            print(finalScore)
            variant_data[variant_id] = finalScore

        return variant_data



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
            aggregated["VariantGetterBioMCPAgent"] = variant_getter.search_variants(coordinates, phenotype)

        if alpha_genome_agent is not None and hasattr(alpha_genome_agent, "predict_variants_effects"):
            aggregated["AlphaGenomeAgent"] = alpha_genome_agent.predict_variants_effects(
                coordinates=coordinates,
                phenotype=phenotype,
            )

        if evo2_agent is not None:
            aggregated["Evo2Agent"] = evo2_agent.getEvo2score(coordinates=coordinates)


        if gpn_agent is not None and hasattr(gpn_agent, "score"):
            aggregated["GPNAgent"] = gpn_agent.score(
                coordinates=coordinates,
                phenotype=phenotype,
            )

        self.latest_results = aggregated
        return aggregated
