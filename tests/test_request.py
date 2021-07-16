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
