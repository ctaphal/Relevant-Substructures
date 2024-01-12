from django.shortcuts import render
from django.http import HttpResponse
from main.graphutils import average_tuples
from main.graphutils import flatten
from main.graphutils import find_best_repeated_submolecule 
from main.graphutils import nth_smallest_submolecule 
from main.graphutils import convertMolecule
from main.graphutils import iterate_submolecules
from main.graphutils import draw_submolecules
from main.graphutils import two_molecules
import base64


def start(request):
    return render(request, "projectApp/home.html")

def home(request):
    return render(request, "projectApp/home.html")  

def oneMolecule(request):
    return render(request, "projectApp/oneMolecule.html")

def twoMolecule(request):
    mol1 =  request.POST.get("firstMoleculeInput")
    mol2 =  request.POST.get("secondMoleculeInput")
    # print(mol1, mol2)
    ms = []
    ms.append(mol1)
    ms.append(mol2)
    
    all_mols = two_molecules(ms) 
    # print(all_mols)
    return render(request, "projectApp/processedTwo.html", {
        "mols": base64.b64encode(all_mols[0]).decode(),
        "common_mol": base64.b64encode(all_mols[1]).decode(),
        })

def collectInputTwo(request):
    return render(request, 'projectApp/twoMolecule.html')

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

