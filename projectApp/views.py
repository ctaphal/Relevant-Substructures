from django.shortcuts import render
from django.http import HttpResponse
from subunitIdent import subunitIdent

def home(request):
    return render(request, "projectApp/home.html")

def oneMolecule(request):
    return render(request, "projectApp/oneMolecule.html")

def twoMolecules(request):
    return render(request, "projectApp/TwoMolecules.html")

def collectInputOne(request):
    return render(request, 'projectApp/oneMolecule.html')

def collectInputTwo(request):
    return render(request, 'projectApp/twoMolecules.html')

def processData(request):
    input = request.POST.get("user-input")
    output = subunitIdent(input)
    return render(request, "projectApp/oneMolecules.html", {"molecular_data": output.molecular_data})
