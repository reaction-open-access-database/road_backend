# Load secret key, database, etc from .env file
from dotenv import load_dotenv, find_dotenv
load_dotenv(find_dotenv())

# Base settings
from .base import *

# Development settings
DEBUG = True

ALLOWED_HOSTS = []
