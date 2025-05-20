from django.shortcuts import render
from .forms import LiteratureForm
from .utils import fetch_pubmed_abstracts, summarize_with_llm

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
