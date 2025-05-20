from django.shortcuts import render
from .forms import LiteratureForm
from .utils import fetch_pubmed_abstracts, summarize_with_llm
from collections import defaultdict
import datetime
from itertools import combinations
from collections import Counter
import networkx as nx
import requests

def pubtator_annotate(request):
    query = request.GET.get("query", "")
    annotated = []

    if query:
        results = pubmed_query(query)
        for res in results[:10]:  # limit to 10 abstracts
            pmid = res["pmid"]
            r = requests.get(f"https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/pubtator?pmids={pmid}")
            if r.status_code == 200:
                annotations = parse_pubtator(r.text)
                highlighted = highlight_entities(res["abstract"], annotations)
                annotated.append({"title": res["title"], "text": highlighted})

    return render(request, "literature_summarizer/annotated_text.html", {"results": annotated})

def keyword_network(request):
    query = request.GET.get("query", "")
    co_occurrences = Counter()

    if query:
        abstracts = [r["abstract"] for r in pubmed_query(query)]
        for abs_text in abstracts:
            keywords = extract_keywords(abs_text)  # You can use RAKE, TF-IDF, etc.
            for combo in combinations(set(keywords), 2):
                co_occurrences[tuple(sorted(combo))] += 1

    G = nx.Graph()
    for (k1, k2), weight in co_occurrences.items():
        if weight >= 2:
            G.add_edge(k1, k2, weight=weight)

    nodes = [{"data": {"id": n}} for n in G.nodes()]
    edges = [{"data": {"source": u, "target": v, "weight": d['weight']}} for u, v, d in G.edges(data=True)]

    return render(request, "literature_summarizer/cooccurrence.html", {"elements": nodes + edges})

def trend_analysis(request):
    keyword = request.GET.get("keyword", "")
    year_freq = defaultdict(int)

    if keyword:
        # Simulate or query PubMed API, extract years from articles
        results = pubmed_query(keyword)
        for result in results:
            year = extract_year(result.get("pubDate", ""))
            if year:
                year_freq[year] += 1

    data = sorted(year_freq.items())  # e.g., [(2021, 3), (2022, 5), ...]
    return render(request, "literature_summarizer/trends.html", {"data": data, "keyword": keyword})

def summarize_view(request):
    summary = ""
    if request.method == "POST":
        form = LiteratureForm(request.POST)
        if form.is_valid():
            query = form.cleaned_data["query"]
            num = form.cleaned_data["num_results"]
            abstracts = fetch_pubmed_abstracts(query, max_results=num)
            summary = summarize_with_llm(abstracts)
    else:
        form = LiteratureForm()
    return render(request, "literature_summarizer/summarize.html", {"form": form, "summary": summary})
