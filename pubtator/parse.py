"""Parsing functions to handle the output of the Pubtator requests.

"""
import logging
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np


@dataclass
class AnnotatedPubtatorAbstract:
    """Holds text and annotation information for a specified PMID.

    Attributes:
        pmid: An integer containing the id of the data.
        text: A string containing the abstract.
        string_representations: A dictionary that contains the representations of all the concepts in the text. This is
            best seen via an example:
                {
                    "Gene": {
                        "123": ["TSH"],
                        "345": ["Crp, "C-reactive protein]
                    },
                    "Disease": {...},
                    ...
                }
            This indicates the text contains references to genes with ids 123 and 345 and the list denotes the names
            that each id takes in the string.
        concept_ids: A dictionary with concepts as keys and a list of ids for that concept as the values.
        associated_ids: A list of all associated ids.
    """

    pmid: int
    text: str
    string_representations: Dict[str, Dict[str, List[str]]]

    def __post_init__(self) -> None:
        self.string_representations = {
            concept: {id.strip("MESH:"): s_reprs for id, s_reprs in id_dict.items()}
            for concept, id_dict in self.string_representations.items()
        }
        self.concepts = list(self.string_representations.keys())

    @property
    def concept_ids(self) -> Dict[str, List[str]]:
        return {concept: list(id_dict.keys()) for concept, id_dict in self.string_representations.items()}

    @property
    def associated_ids(self) -> List[str]:
        return [x for ls in self.concept_ids.values() for x in ls]


def parse_pubtator_request_string(
    request_string: str, required_concepts: Tuple[str, ...] = ("Gene", "Disease")
) -> List[AnnotatedPubtatorAbstract]:
    """Handles the parsing of multiple annotated abstracts returned from a pubtator request.

    This converts the text returned from a request to the Pubtator API onto a list of annotated abstracts that are
    significantly easier to deal with.

    Arguments:
        request_string: The return text object from a request to the pubtator-api. This will usually just be the return
            of the request_annotations function.
        required_concepts: Concepts that must have at least one annotation for the element to be considered.

    Return:
        A list of AnnotatedPubtatorAbstracts which provide a more user friendly version of the data.
    """
    # Results are separated by two lines
    annotated_abstracts = request_string.split("\n\n")
    if annotated_abstracts[-1] == "":
        annotated_abstracts.pop(-1)

    # Parse each pmid into a dictionary
    annotated_pubtator_abstracts = []
    for annotated_abstract in annotated_abstracts:
        if all(["\t{}\t".format(s) in annotated_abstract for s in required_concepts]):
            try:
                annotated_pubtator_abstracts.append(_parse_single_annotated_abstract(annotated_abstract))
            except Exception as e:
                logging.warning("Abstract could not be parsed with error {}.".format(e))

    return annotated_pubtator_abstracts


def _parse_single_annotated_abstract(
    annotated_abstract: str,
) -> AnnotatedPubtatorAbstract:
    """Handles the parsing of a single abstract with annotation in pubtator format.

    Pubtator returns results from multiple PMIDs split by two lines (\n\n). This converts one of the elements of the
    request string after they have been split on \n\n.

    Arguments:
        annotated_abstract: Annotated abstract string as returned by the pubtator API.

    Returns:
        A populated AnnotatedPubtatorAbstract object.
    """
    split_text = annotated_abstract.split("\n")
    title_line, abstract_line = split_text[0], split_text[1]

    # Make title
    split_title = title_line.split("|")
    pmid, title = split_title[0], split_title[2]
    abstract = abstract_line.split("|")[-1]
    text = title + " " + abstract

    # Handle the annotations
    string_representations = _parse_annotation_lines(split_text[2:])

    # Add into the dataclass
    return AnnotatedPubtatorAbstract(pmid=int(pmid), text=text, string_representations=string_representations)


def _parse_annotation_lines(annotations: List[str]) -> Dict[str, Dict[str, List[str]]]:
    """Parses the important information from the annotation lines into dicts keyed with concepts and ids.

    We are primarily interested in each concept id that occurs in the text, and what strings were used to represent that
    concept. To this end, we parse the annotation component of the request sting into a dictionary that can be keyed
    with [concept][concept_id] to return the string representations of the concept_id in the text. If we then wish to
    replace said annotations we can do this easily by replacing the string representations.

    Arguments:
        annotations: The annotations split by line in a list format.

    Returns:
        A multilevel dictionary with first level being the concept and second level the id. The values are the string
        representations that the ids can take in the text.
    """
    # Only take lines of length 6 else there is a parsing fail
    csv_lines = [x.split("\t") for x in annotations]
    annotation_array = np.array([x for x in csv_lines if len(x) == 6])

    # Now build gene/disease dicts that map ids to word occurrences according to annotation locations
    concepts = np.unique(annotation_array[:, -2])

    concept_representations = dict.fromkeys(concepts)
    for concept in concepts:
        concept_representations[concept] = {}
        for id in np.unique(annotation_array[annotation_array[:, -2] == concept][:, -1]):
            sub_annotation = annotation_array[annotation_array[:, -1] == id]
            string_representations = np.unique(sub_annotation[:, 3]).tolist()
            concept_representations[concept][id] = string_representations

    return concept_representations
