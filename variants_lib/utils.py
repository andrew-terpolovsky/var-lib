from typing import Literal
import boto3

ROLE_MAPPING = {
    "build38": "dynamodb_crossaccount_readonlyaccess_role",
    "merged-rsids": "merged_rsids_crossaccount_readonlyaccess_role",
}


def get_ddb_client(role: Literal["build38", "merged-rsids"]):

    sts_client = boto3.client("sts")
    sts_session = sts_client.assume_role(
        RoleArn=f"arn:aws:iam::308889746780:role/{ROLE_MAPPING[role]}",  # pylint: disable=line-too-long
        RoleSessionName="variants-lib-assumed-role-shared-ddb",
    )
    credentials = sts_session["Credentials"]
    client = boto3.client(
        "dynamodb",
        region_name="us-east-1",
        aws_access_key_id=credentials["AccessKeyId"],
        aws_secret_access_key=credentials["SecretAccessKey"],
        aws_session_token=credentials["SessionToken"],
    )
    return client
