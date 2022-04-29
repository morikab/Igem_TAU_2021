import os
import requests
from pathlib import Path

supported_formats = ("FASTA", "GFF")


def download_genome_files(genome, download_url):
    download_folder = os.path.join("mgnify_files", genome)
    Path(download_folder).mkdir(parents=True, exist_ok=True)
    downloads = requests.get(download_url)
    downloads_json = downloads.json()
    files_list = downloads_json["data"]
    for file in files_list:
        file_format = file["attributes"]["file-format"]["name"]
        if file_format not in supported_formats:
            continue
        file_url = file["links"]["self"]
        file_response = requests.get(file_url, allow_redirects=True)
        file_name = file["id"]
        with open(os.path.join(download_folder, file_name), "wb") as fd:
            fd.write(file_response.content)


def download_files():
    url = "https://www.ebi.ac.uk/metagenomics/api/v1/genomes"
    page_count = 0
    while url is not None:
        print(F"Starting page no. {page_count}")
        genomes_response = requests.get(url)
        genomes_json = genomes_response.json()

        genomes_list = genomes_json["data"]
        for genome in genomes_list:
            print(F"Iterating id {genome['id']}")
            genome_type = genome["attributes"]["type"]
            if genome_type != "MAG":
                print(F"Skip id {genome['id']}. Type: {genome_type}")
                continue
            download_link = genome["relationships"]["downloads"]["links"]["related"]
            download_genome_files(genome=genome["id"], download_url=download_link)

        page_count += 1

        url = genomes_json["links"]["next"]


if __name__ == "__main__":
    print("Start")
    download_files()
    print("The end")

