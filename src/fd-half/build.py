#
# Скрипт сборки проекта
#
import os
import sys
import subprocess
import platform

def main():
    build_dir = os.getcwd() + "/build"
    if not os.path.exists(build_dir):
        os.makedirs(build_dir)
    os.chdir(build_dir)
    # TODO: function run_subprocess_sync
    if platform.system() == "Windows":
        subprocess.call(["cmake","-DCMAKE_GENERATOR_PLATFORM=x64","-DCMAKE_BUILD_TYPE=Release",".."])
    else:
        subprocess.call(["cmake","-DCMAKE_BUILD_TYPE=Release",".."])
    subprocess.call(["cmake", "--build", "."])

main()