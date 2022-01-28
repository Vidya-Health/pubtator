"""Makes requests to the Pubtator api.

This file contains a single function that can query the pubtator API for annotated PMID abstracts or annotated full-text
PMCIDs.
"""
import logging
from typing import Union
from enum import Enum

import omni
import requests


class PubtatorDatabase(Enum):
    pubmed: str = "pmids"
    pmc: str = "pmcids"


def request_annotations(
    ids: list[Union[int, str]],
    database: PubtatorDatabase = PubtatorDatabase.pubmed,
    concepts: tuple[str, ...] = ("gene", "disease"),
    pubtator: bool = True,
) -> str:
    """Requests annotated documents from the pubtator central api.

    Follows the details laid out in:
        https://www.ncbi.nlm.nih.gov/research/pubtator/api.html

    Makes a request to the url:
        'https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export'
    with the specified options for requesting the documents and concept annotations according to the specified format.

    Arguments:
        ids: A list of pmids to request information on.
        database: A PubtatorDatabase enum value specifying the database to request information on.
        concepts: A list of concepts from (gene, disease, chemical, species, mutation and/or cellline).
        pubtator: Set True to return in pubtator format, False will instead result in biocjson.

    Returns:
        A single string containing the request in pubtator format (unless pubtator is set to False in which case it is
        in biocjson format).
    """
    # Base url
    assert not (pubtator and database == PubtatorDatabase.pmc), "PMC must have pubtator set to False."
    fmt = "pubtator" if pubtator else "biocxml"
    url = f"https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/{fmt}"

    # id and concept params
    params = {
        database.value: ",".join([str(x) for x in ids]),
        "concepts": ",".join(concepts),
    }

    # Make the get request
    response = requests.get(url, params=params)

    # Handle error
    if not response:
        logging.warning(
            "Response returned status code: {}\nResponse: {}\nurl: {}".format(
                response.status_code, response, response.url
            )
        )

    if database == PubtatorDatabase.pmc:
        if len(response.text) < 500:
            raise FileNotFoundError("No PMC articles found.")

    return response.text
