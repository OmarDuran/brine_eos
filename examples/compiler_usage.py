from generator.source_handler import retrieve_and_rename_source
from generator.source_handler import compile_library_and_cli_source


file_name: str = "saltwatereos-1.7.0.zip"
destination_folder: str = "extracted_source"

# compile library and cli source for data generation
compile_library_and_cli_source(file_name,destination_folder)
