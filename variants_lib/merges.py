from typing import List, Dict
import os
from variants_lib.utils import get_ddb_client

TABLE_NAME = os.environ.get("MERGES_TABLE", "merged-rsids")


def canonical_rsids(rsids: List[str]) -> Dict[str, str]:
    """Return a dictionary mapping each provided rsid to the
    corresponding most recent rsid."""
    client = get_ddb_client("merged-rsids")
    payload = {TABLE_NAME: {"Keys": [{"rsid": {"S": rsid}} for rsid in rsids]}}
    resp = client.batch_get_item(RequestItems=payload)
    result = {rsid: rsid for rsid in rsids}
    try:
        items = resp["Responses"][TABLE_NAME]
    except KeyError:
        # None of the rsids was found in the table (= no merge).
        items = []
    for item in items:
        if "merged_into" in item:
            result[item["rsid"]["S"]] = item["merged_into"]["S"]
    return result
