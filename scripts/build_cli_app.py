from globals import source_name
from generator.source_handler import compile_library_and_cli_source
from generator.source_handler import move_cli

def build_cli_application(file_name: str = source_name, destination_folder: str = 'extracted_source'):
    try:
        # Compile library and cli source for data generation
        compile_library_and_cli_source(file_name, destination_folder)
    except Exception as e:
        print(f"Error: {e} - An unexpected error occurred.")
    finally:

        # move cli
        origin_folder: str = destination_folder
        move_cli(origin_folder)
        print("CLI is built.")

