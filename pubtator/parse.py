"""Parsing functions to handle the output of the Pubtator requests.

"""
import collections
import logging
from dataclasses import dataclass
from typing import Any
from typing import Dict, List, Tuple, Optional

import numpy as np
import omni
import pandas as pd
import xmltodict


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


def _filter_passages(passages: list[dict[Any, Any]]) -> list[dict[Any, Any]]:
    keep_passages = []
    for passage in passages:
        # We care about it if it doesnt contain annotations
        if "annotation" in passage:
            passage["infon"]


def get_section_type(passage: dict[Any, Any]) -> Optional[str]:
    """Gets the section type information from the passage"""
    section_type = None
    section_type_dict = [p for p in passage["infon"] if p["@key"] == "section_type"]
    if len(section_type_dict) > 0:
        section_type = section_type_dict[0].get("#text")
    return section_type


def parse_annotations(annotations: list[collections.OrderedDict[str, str]]) -> list[list[str]]:
    """Pull out and neaten the annotations."""
    data_list = []
    for annotation in annotations:
        # Get the text
        data = {"text": annotation.get("text")}

        # New get the identifier and the type
        infon = annotation.get("infon")
        if isinstance(infon, list):
            for ann in annotation.get("infon"):
                key = ann.get("@key")
                if key in ("identifier", "type"):
                    data[key] = ann["#text"]

        # Provided the data is all there, add to the full data list
        if all(x in data for x in ["identifier", "type"]):
            data_list.append([data["identifier"].replace("MESH:", ""), data["type"], data["text"]])

    return data_list


def list_mentions_to_dict(mentions: list[list[str]]) -> dict[str, dict[str, list[str]]]:
    """Converts a list of [id, type, text] into a dictionary of {type: {id: [text1,text2, ...]}}"""
    # Create an initial dict containing the types
    mentions_frame = pd.DataFrame(mentions, columns=["id", "type", "mentions"])
    mentions_frame = mentions_frame.groupby(["type", "id"], as_index=True)["mentions"].apply(
        lambda x: list(set(x.tolist()))
    )

    mentions_dict = {}
    for tp in mentions_frame.index.get_level_values(0).unique():
        mentions_dict[tp] = mentions_frame.loc[tp].to_dict()

    return mentions_dict


def append_text(old_text: str, new_text: str) -> str:
    """Some basic append rules with punctuation additions."""
    if len(old_text) == 0:
        text = new_text
    elif old_text[-1] in (".", "?", "!"):
        text = old_text + " " + new_text
    else:
        text = old_text + ". " + new_text
    return text


def parse_pubtator_pmc_request_string(
    request_string: str, required_concepts: Optional[tuple[str]] = ("Disease",), pmc_to_pmid: dict[str, int] = None
) -> list[list[AnnotatedPubtatorAbstract]]:
    """Handles the parsing of multiple annotated abstracts returned from a pubtator request.

    This splits the string into sections as indexed by the passages in the xml. Each is made its own
    AnnotatedPubtatorAbstract and given an additional section_count field to keep track. Text is only kept if disease
    and chemical or gene are present.

    Args:
        request_string: The output of the pmc xml request.
        required_concepts: The required concepts to keep the text.
        pmc_to_pmid: An optional dictionary that maps the pmc back to the pmid for the id.

    Returns:
        A list of lists of AnnotatedPubtatorAbstracts where each inner list contains the abstract for one single pmcid.
    """
    # Make dict and get the only bits that matter
    xml_dict = xmltodict.parse(request_string, process_namespaces=True, xml_attribs=True)
    docs = xml_dict["collection"]["document"]
    docs = docs if isinstance(docs, list) else [docs]

    # Iterate through each doc and get the pmcid and concatenated text
    outputs = []
    for doc in docs:
        # Setup id
        pmcid = "PMC" + doc["id"]
        if pmc_to_pmid:
            pmcid = pmc_to_pmid[pmcid]

        section_count = 1
        pmcid_outputs = []

        for passage in doc["passage"]:
            # Require multiple annotations, i.e. a list
            if not isinstance(passage, collections.OrderedDict):
                continue

            # Get annotations, skip if they do not exist
            annotations = passage.get("annotation")
            annotations = [annotations] if isinstance(annotations, collections.OrderedDict) else annotations
            if not isinstance(annotations, list):
                continue

            # Find the annotations
            mentions = list_mentions_to_dict(parse_annotations(annotations))

            # Check if the passage has the required concepts
            if required_concepts:
                if not all(s in mentions for s in required_concepts):
                    continue

            # Get the text
            text = passage.get("text")
            abstract_obj = AnnotatedPubtatorAbstract(pmcid, text, mentions)
            abstract_obj.section_count = section_count

            # Add to outputs and update
            pmcid_outputs.append(abstract_obj)
            section_count += 1

        # Add if we had a passage added
        if len(pmcid_outputs) > 0:
            outputs.append(pmcid_outputs)

    return outputs


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
