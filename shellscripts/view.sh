#!/bin/bash

source bin/activate

django-admin runserver --pythonpath="." --settings=Server 5000