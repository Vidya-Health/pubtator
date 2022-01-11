from Bio import Entrez


def get_pubmed_ids_in_daterange(email: str, start_date: int, end_date: int = 3000) -> list[int]:
    """Gets a list of PubMed IDs in a given date range.

    Code taken from: https://stackoverflow.com/questions/37539333/download-a-list-of-all-pubmed-ids-by-date-from-to

    Args:
        email: Email address for Entrez.
        start_date: Start date for the date range, e.g. '2014/12/20'.
        end_date: End year for the date range, defaults to 3000 i.e. all papers from the start_date.

    Returns:
        A list of ints each representing an existing PubMed ID.
    """
    Entrez.email = email

    # Get the PubMed IDs for the given date range.
    handle = Entrez.esearch(
        db="pubmed",
        term='("%s"[Date - Publication] : "%s"[Date - Publication]) ' % (start_date, end_date),
        retmax=100000000,
    )
    records = Entrez.read(handle)

    # Concat and return the list of PubMed IDs.
    pmids = sorted([int(record) for record in records["IdList"]])

    return pmids
