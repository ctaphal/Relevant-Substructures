from django.shortcuts import render
from django.http import HttpResponse
#from .forms import oneMolecule

def home(request):
    return render(request, "projectApp/home.html")

def inputMolecules(request):
    molecule = request.POST.get("molecule")
    print(molecule)
    return HttpResponse("Recieved molecule!")

def inputOneMolecule(request):
    return render(request, "projectApp/oneMolecule.html")

def inputTwoMolecules(request):
    return render(request, "projectApp/twoMolecules.html")