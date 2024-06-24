from globals import source_name
from scripts.build_cli_app import build_cli_application
from generator.source_handler import run_command

try:
    # # build cli application
    destination_folder: str = "extracted_source"
    build_cli_application(destination_folder=destination_folder)
except Exception as e:
    print(f"Error: {e} - An unexpected error occurred.")
finally:
    # invoke the cli application
    prefix = "./cli_" + source_name + "/"
    run_command(prefix + "swEOS -help")
