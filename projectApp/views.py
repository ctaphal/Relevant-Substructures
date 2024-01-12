from django.shortcuts import render
from django.http import HttpResponse
from main.graphutils import convertMolecule
from main.graphutils import iterate_submolecules
from main.graphutils import draw_submolecules

def home(request):
    return render(request, "projectApp/home.html")

def oneMolecule(request):
    return render(request, "projectApp/oneMolecule.html")

def twoMolecules(request):
    return render(request, "projectApp/twoMolecules.html")

def collectInputTwo(request):
    return render(request, 'projectApp/twoMolecules.html')

def displayImages(smile):
    images = draw_submolecules(*iterate_submolecules(smile))

    indexed = images[-2:-1]
    return render(smile, "projectApp/processedOne.html", {"images": indexed.images})

def processDataOne(request):
    input = request.POST.get("user-input")
    print(input)
    output = convertMolecule(input)

    displayImages(output)
    return render(request, "projectApp/processedOne.html", {"molecular_data": output.molecular_data})
