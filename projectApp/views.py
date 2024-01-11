from django.shortcuts import render
from django.http import HttpResponse

def getMolecules(request):
    return render(request, "projectApp/get-molecules.html")

def inputMolecules(request):
    molecule = request.POST.get("molecule")
    print(molecule)
    return HttpResponse("Recieved molecule!")

