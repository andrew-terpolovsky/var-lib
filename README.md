# Permissions

## Dynamodb tables

The resource using the package must be able to access the `merged-rsids` and the `genome_reference_build38` dynamodb tables from the shared-environment, as well as the Athena metadata table from the environment holding the desired user genome files. The former two permissions can be hardcoded because they do not depend on the execution environment:

```
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Sid": "VisualEditor0",
            "Effect": "Allow",
            "Action": "sts:AssumeRole",
            "Resource": "arn:aws:iam::308889746780:role/merged_rsids_crossaccount_readonlyaccess_role"
        }, 
        {
            "Sid": "VisualEditor0",
            "Effect": "Allow",
            "Action": "sts:AssumeRole",
            "Resource": "arn:aws:iam::308889746780:role/dynamodb_crossaccount_readonlyaccess_role"
        }
    ]
}
```

For the metadata table, the `dynamodb:GetItem` action must be allowed, over the dynamodb exported by `sd-platform-{STAGE}` as `sd-platform-{STAGE}:AthenaSnpTableWorkflowDdbTable`.

## S3 Bucket

The resource must be able to read (`s3:ListBucket`, `s3:GetObject`) from the bucket where Athena stores the parquet files, which is exported by `sd-platform-{STAGE}` as `sd-platform-{STAGE}:AthenaDataS3Bucket`.


# Environment variables

The `ATHENA_TABLE_METADATA` should be equal to the Athena metadata dynamodb table's name, see above for the exported value.