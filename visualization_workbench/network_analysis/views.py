from django.shortcuts import render
from django.http import JsonResponse
from .utils.network_utils import get_sample_network, detect_modules

def network_home(request):
    return render(request, 'network_analysis/network.html')

def fetch_network_data(request):
    network_data = get_sample_network()  # This would eventually be dynamic
    return JsonResponse(network_data)

def fetch_modules(request):
    module_data = detect_modules()
    return JsonResponse(module_data)
