import os
import re
from datetime import datetime
from docx import Document
from duckduckgo_search import DDGS
from bs4 import BeautifulSoup
import requests
from llm_utils import query_llm
import prompts
import io

def variant_getter_biomcp_prompt():
    """Create a comprehensive plan prompt for the Principal Investigator"""
    
    variant_getter_biomcp_prompt = f"""
    As a Principal Investigator, analyze the following sources and create a detailed plan for the topic: '{topic}'
    
    Sources:
    {sources}
    
    Mode: {mode}
    
    Create a detailed plan by THINKING STEP BY STEP that includes:
    1. Key insights from the sources
    2. Analysis of any files found in the directory (if applicable)
    3. Specific tasks for each agent:
       - Research Agent: What aspects to focus on
       - Code Writer Agent: What code to implement, what packages to import
       - Code Executor Agent: How to execute the code, what packages to install
       - Code Reviewer Agent: What to review in the code and the execution result or error messages
       - Critic Agent: What to evaluate in the report and the code
    
    
    Provide a clear, actionable plan that all agents can follow. You do not need to provide time estimates for tasks. 
    ABSOLUTELY DO NOT plan for tasks that haven't been asked for, if you do, it will destroy the pipeline.
    I REPEAT, DO NOT PLAN FOR TASKS THAT HAVEN'T BEEN ASKED FOR.
    """

    # if changes:
    #     plan_prompt += f"""
    # User requested the following changes to the plan:
    # {changes}

    # Reasoning about the changes and incorporating them into the new plan.
    # """
    
    return variant_getter_biomcp_prompt

