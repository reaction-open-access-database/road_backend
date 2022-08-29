from django.urls import path, include
from rest_framework import routers
from .views import MoleculeViewSet, UserViewSet
from dj_rest_auth.registration.views import VerifyEmailView, RegisterView, ConfirmEmailView, ResendEmailVerificationView
from dj_rest_auth.views import LoginView, LogoutView
from rest_framework_simplejwt.views import TokenVerifyView
from dj_rest_auth.jwt_auth import get_refresh_view

router = routers.DefaultRouter()
router.register('molecules', MoleculeViewSet)
router.register('user-profiles', UserViewSet)

urlpatterns = [
    path('', include(router.urls)),
    path('register/', RegisterView.as_view()),
    path('login/', LoginView.as_view()),
    path('logout/', LogoutView.as_view()),
    path('resend-email-verification/', ResendEmailVerificationView.as_view()),
    path('token/verify/', TokenVerifyView.as_view(), name='token_verify'),
    path('token/refresh/', get_refresh_view().as_view(), name='token_refresh'),
    path('account-confirm-email/', VerifyEmailView.as_view(), name='account_email_verification_sent'),
    path('account-confirm-email/<key>/', ConfirmEmailView.as_view(), name='account_confirm_email'),
]
