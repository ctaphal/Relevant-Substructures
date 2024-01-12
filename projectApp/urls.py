from django.urls import path
from . import views

urlpatterns = [
    path("home/", views.home, name="home"),
    path("inputMolecules/", views.inputMolecules, name="inputMolecules"),
    path("inputOneMolecule/", views.inputOneMolecule, name="inputOneMolecule"),
    path("inputTwoMolecules/", views.inputTwoMolecules, name="inputTwoMolecules")
]