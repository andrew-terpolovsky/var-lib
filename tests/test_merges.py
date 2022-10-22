from pytest import fixture
from variants_lib import merges
from .utils import MockDdbClient


class TestCanonicalRsid:
    @fixture(autouse=True)
    def mock_merges_table(self, monkeypatch):
        monkeypatch.setattr(
            merges,
            "get_ddb_client",
            lambda role: MockDdbClient(
                table_name="merged-rsids",
                items=[
                    {"rsid": {"S": "rs123"}, "merged_into": {"S": "rs456"}},
                    {"rsid": {"S": "rs321"}},
                ],
                key_name="rsid",
            )
            if role == "merged-rsids"
            else None,
        )

    def test_canonical_rsids(self):
        assert merges.canonical_rsids(["rs123", "rs321", "rs456", "rs789"]) == {
            "rs123": "rs456",
            "rs321": "rs321",
            "rs456": "rs456",
            "rs789": "rs789",
        }
