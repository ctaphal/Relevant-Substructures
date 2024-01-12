from django.urls import path
from django.conf import settings
from django.conf.urls.static import static
from . import views

urlpatterns = [
    path("", views.home, name="home"),
    path("oneMolecule/", views.oneMolecule, name="oneMolecule"),
    path("views.twoMolecules/", views.twoMolecules, name="twoMolecules"),
    path("collectDataTwo/", views.collectInputTwo, name="collectInputTwo"),
    path("processDataOne/", views.processDataOne, name="processDataOne")
] + static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
