import pytest
import pubtator


def test_request_annotations() -> None:
    # Check pmids return fine
    annotations = pubtator.request_annotations([28483577, 28483578, 28483579])
    assert all(
        [
            s in annotations
            for s in (
                "Elastography for Assessing Liver Fibrosis",
                "myofibroblast-like phenotype",
                "Chronic Hepatitis B	Disease	MESH:D019694",
            )
        ]
    )


def test_request_pmcid() -> None:
    # Check we can get pmcid information
    fail_ids = ["PMC1666769", "PMC1666842", "PMC1666801"]
    with pytest.raises(FileNotFoundError):
        _ = pubtator.request_annotations(
            fail_ids,
            database=pubtator.request.PubtatorDatabase.pmc,
            concepts=("gene", "disease", "chemical"),
            pubtator=False,
        )

    pass_ids = ["PMC7571963", "PMC4623938"]
    annotations = pubtator.request_annotations(
        pass_ids,
        database=pubtator.request.PubtatorDatabase.pmc,
        concepts=("gene", "disease", "chemical"),
        pubtator=False,
    )
    assert (
        isinstance(p, pubtator.AnnotatedPubtatorAbstract)
        for p in pubtator.parse.parse_pubtator_pmc_request_string(annotations)
    )
