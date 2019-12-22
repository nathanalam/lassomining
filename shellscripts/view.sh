#!/bin/bash

source bin/activate

django-admin runserver --pythonpath="." --settings=Server 50.116.48.39:8080