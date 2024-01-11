from django.shortcuts import render
from django.http import HttpResponse
from .forms import oneMolecule

def getMolecules(request):
    return render(request, "projectApp/get-molecules.html")

def inputMolecules(request):
    molecule = request.POST.get("molecule")
    print(molecule)
    return HttpResponse("Recieved molecule!")

"""

def inputOneMolecule(request):
    form = oneMolecules()
    return render(request, 'projectApp/oneMolecule.html', {'form': form})
"""