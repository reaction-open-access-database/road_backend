from django.urls import path, include
from rest_framework import routers
from .views import MoleculeViewSet
from dj_rest_auth.registration.views import VerifyEmailView, RegisterView, ConfirmEmailView
from dj_rest_auth.views import LoginView, LogoutView

router = routers.DefaultRouter()
router.register('molecules', MoleculeViewSet)

urlpatterns = [
    path('', include(router.urls)),
    path('register/', RegisterView.as_view()),
    path('login/', LoginView.as_view()),
    path('logout/', LogoutView.as_view()),
    path('account-confirm-email/', VerifyEmailView.as_view(), name='account_email_verification_sent'),
    path('account-confirm-email/<key>/', ConfirmEmailView.as_view(), name='account_confirm_email'),
]
