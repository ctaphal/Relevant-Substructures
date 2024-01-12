from django.urls import path
from django.conf import settings
from django.conf.urls.static import static
from . import views

urlpatterns = [
    path("home/", views.home, name="home"),
    path("oneMolecule/", views.oneMolecule, name="oneMolecule"),
    path("views.twoMolecules/", views.twoMolecules, name="twoMolecules"),
    path("collectDataOne/", views.collectInputOne, name="collectInputOne"),
    path("collectDataTwo/", views.collectInputTwo, name="collectInputTwo"),
    path("processData/", views.processData, name="processData")
] + static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)