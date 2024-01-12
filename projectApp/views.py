from django.shortcuts import render
from django.http import HttpResponse
from main.graphutils import average_tuples
from main.graphutils import flatten
from main.graphutils import find_best_repeated_submolecule 
from main.graphutils import nth_smallest_submolecule 
from main.graphutils import convertMolecule
from main.graphutils import iterate_submolecules
from main.graphutils import draw_submolecules
import base64



def home(request):
    return render(request, "projectApp/home.html")

def oneMolecule(request):
    return render(request, "projectApp/oneMolecule.html")

def twoMolecules(request):
    return render(request, "projectApp/twoMolecules.html")

def collectInputTwo(request):
    return render(request, 'projectApp/twoMolecules.html')

def displayImages(request, mol):
    smallest = nth_smallest_submolecule(1, *iterate_submolecules(mol))
    second_smallest = nth_smallest_submolecule(2, *iterate_submolecules(mol))
    if second_smallest:
        return render(request, "projectApp/processedOne.html", {
            "smile_1": smallest[0],
            "submol_base64_1": base64.b64encode(smallest[1]).decode(),
            "annotated_molecule_base64_1": base64.b64encode(smallest[2]).decode(),
            "smile_2": second_smallest[0],
            "submol_base64_2": base64.b64encode(second_smallest[1]).decode(),
            "annotated_molecule_base64_2": base64.b64encode(second_smallest[2]).decode(),
        })
    else:
        return render(request, "projectApp/processedOne.html", {
            "smile_1": smallest[0],
            "submol_base64_1": base64.b64encode(smallest[1]).decode(),
            "annotated_molecule_base64_1": base64.b64encode(smallest[2]).decode(),
        })

def processDataOne(request):
    input = request.POST.get("user-input")
    print(input)
    
    output = convertMolecule(input)
    
    return displayImages(request, output)
    
    # return render(request, "projectApp/processedOne.html", {"molecular_data": output.molecular_data})

