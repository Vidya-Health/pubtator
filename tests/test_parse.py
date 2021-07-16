import pubtator


def test_parse_single_annotated_abstract() -> None:
    # Test with a single request.
    annotated_pubtator_request = "\n".join(
        [
            "34223537|t|Enteric Nervous System Remodeling in a Rat Model of Spinal Cord...",
            "34223537|a|The physiopathology of digestive disorders in patients with spinal cord injury (SCI) remains "
            "largely...",
            "34223537	52	70	Spinal Cord Injury	Disease	MESH:D013119",
            "34223537	147	165	spinal cord injury	Disease	MESH:D013119",
            "34223537	825	850	choline acetyltransferase	Gene	290567",
            "34223537	852	859	Colonic	Disease	MESH:D003110",
            "34223537	1613	1646	intercellular adhesion molecule-1	Gene	25464",
            "34223537	1648	1654	ICAM-1	Gene	25464",
        ]
    )

    # Create true dataclass
    should_return = pubtator.AnnotatedPubtatorAbstract(
        pmid=34223537,
        text="Enteric Nervous System Remodeling in a Rat Model of Spinal Cord... The physiopathology of digestive "
        "disorders in patients with spinal cord injury (SCI) remains largely...",
        string_representations={
            "Gene": {
                "25464": ["ICAM-1", "intercellular adhesion " "molecule-1"],
                "290567": ["choline acetyltransferase"],
            },
            "Disease": {
                "D003110": ["Colonic"],
                "D013119": ["Spinal Cord Injury", "spinal cord injury"],
            },
        },
    )

    # Check its right
    parsed = pubtator.parse_pubtator_request_string(annotated_pubtator_request)
    assert parsed[0] == should_return
