from globals import source_url
import requests
import subprocess
import os
from pathlib import Path

def retrieve_and_rename_source(file_name: str = "saltwatereos-1.7.0.zip"):
    """
    Retrieve the zip File
    """
    response = requests.get(source_url)
    with open(file_name, 'wb') as file:
        file.write(response.content)

def unpack_source(file_name, folder_name: str = 'extracted_source'):
    """
    Unpack source
    """
    subprocess.run(['mkdir', '-p', folder_name], check=True)
    subprocess.run(['unzip', '-o', file_name, '-d', folder_name], check=True)

def run_command(command, cwd=None):
    """
    Runs a shell command in the specified directory and prints the output.
    """
    result = subprocess.run(command, shell=True, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print(result.stdout.decode())
    if result.returncode != 0:
        print(result.stderr.decode())
        raise Exception(f"Command failed with return code {result.returncode}")

def build_library(folder_name: str = 'extracted_source'):
    path = Path(folder_name)
    library_folder = path / 'saltwatereos-1.7.0' / "Library"
    library_dir = str(library_folder.absolute())
    # Define directories
    os.path.abspath(library_dir)  # Directory containing the C++ source and CMakeLists.txt
    # Run CMake to configure the project
    run_command("sh ./build.sh", cwd=library_dir)

def compile_cli_source(folder_name: str = 'extracted_source'):

    path = Path(folder_name)
    build_folder = path / 'saltwatereos-1.7.0' / "commandline" / "build"
    build_folder.mkdir(parents=True, exist_ok=True)

    build_dir = str(build_folder.absolute())
    # Define directories
    os.path.abspath(build_dir)

    # Create the build directory
    os.makedirs(build_dir, exist_ok=True)

    # Run CMake to configure the project
    run_command('cmake -DCMAKE_CXX_FLAGS="-I/usr/local/include -I/opt/local/include" ..', cwd=build_dir)

    # Run make to build the project
    run_command("make", cwd=build_dir)

def compile_library_and_cli_source(file_name: str = "saltwatereos-1.7.0.zip", destination_folder: str = 'extracted_source'):

    # Step 1: Download the Zip File
    retrieve_and_rename_source(file_name)

    # Step 2: Unpack the Zip File
    unpack_source(file_name, destination_folder)

    # Step 3: build and install the library
    build_library(destination_folder)

    # Step 4: compile the cli application
    compile_cli_source(destination_folder)