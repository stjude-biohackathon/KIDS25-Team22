import requests
import os
from config import LLM_CONFIG

total_tokens_used = 0
output_log = []
os.environ["NO_PROXY"] = "localhost"

def query_llm(prompt, model=LLM_CONFIG["default_model"], temperature=0):
    global total_tokens_used
    url = "http://localhost:11434/api/generate"
    payload = {
        "model": model,
        "prompt": prompt,
        "temperature": temperature,
        "stream": False,
    }
    response = requests.post(url, json=payload)
    if response.status_code == 200:
        response_data = response.json()
        total_tokens_used += response_data.get("eval_count", 0)  # track tokens
        output_log.append(
            {
                "model": model,
                "prompt": prompt,
                "response": response_data.get("response", "").strip(),
                "tokens_used": response_data.get("eval_count", 0),
            }
        )
        return response.json().get("response", "").strip()
    else:
        raise Exception(f"Error: {response.status_code}, {response.text}")
