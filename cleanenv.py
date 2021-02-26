import numpy as np
import re

envfile = open("xvsy_environment_final.yml", "w")
with open('environment.yml', 'r') as f:
    all_lines = f.readlines()
    for line in all_lines:
        pkg = line.strip()
        pkg = "=".join(pkg.split("=")[:2])
        envfile.write(pkg + '\n')
envfile.close()