from django.urls import path
from . import views

urlpatterns = [
    path("getMolecules/", views.getMolecules, name="getMolecules"),
    path("inputMolecules/", views.inputMolecules, name="inputMolecules")
]