# Required imports
from Bio import SeqIO
import gzip
import base64
from io import BytesIO
#import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os,sys
#import seaborn as sns
#from sklearn.metrics import roc_auc_score
import requests
import json
from dotenv import load_dotenv


def calculateSequenceScore(input_sequence, base64_data) -> float:
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


def runEvo2(dna_sequence):
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





human_chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']
if sys.argv[1] in human_chromosomes:
        chr = sys.argv[1]
else:
        print(f"'{sys.argv[1]}' is not a valid human chromosome.")
coordinates=int(sys.argv[2])
ref=sys.argv[3]
alt=sys.argv[4]
#print(chr,coordinates,ref,alt)

WINDOW_SIZE = 8192

# Read the reference genome sequence of chromosome 17
with gzip.open(os.path.join('hg19','chr'+str(chr)+'.fa.gz'), "rt") as handle:
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
#print(row)
#print('--')
#print(f'Reference, SNV 0: ...{ref_seq[4082:4112]}...')
#print(f'Variant, SNV 0:   ...{var_seq[4082:4112]}...')




varScore = calculateSequenceScore(var_seq, runEvo2(var_seq))
refScore = calculateSequenceScore(ref_seq, runEvo2(ref_seq))
print("refScore",refScore)
print("varScore",varScore)
finalScore={"Variant":str(chr)+":"+str(coordinates)+ref+">"+alt, "Evo2_deltaScore":float(varScore-refScore)}
print(finalScore)
