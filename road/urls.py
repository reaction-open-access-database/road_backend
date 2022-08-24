from django.urls import path, include
from rest_framework import routers
from . import views

router = routers.DefaultRouter()
router.register(r'molecules', views.MoleculeViewSet)

urlpatterns = [
    path('', include(router.urls)),
    path('session-auth/', include('rest_framework.urls', namespace='rest_framework')),
    path('token-auth', include('knox.urls')),
]
