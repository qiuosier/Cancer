# Accessing FASTQ files on Illumina BaseSpace
The `basespace` package includes Python API for accessing [Illumina BaseSpace](https://developer.basespace.illumina.com/docs/content/documentation/getting-started/overview). Sequencing data from Illumina instruments, including MiSeq and HiSeq, are automatically uploaded to BaseSpace over the Internet. This API wraps the [BaseSpace V1 REST API](https://developer.basespace.illumina.com/docs/content/documentation/rest-api/v1-api-reference) and focus on downloading and transferring the FASTQ files.

See also:
* [BaseSpace V1 API Reference](https://developer.basespace.illumina.com/docs/content/documentation/rest-api/v1-api-reference)
* [BaseSpace Python SDK](https://developer.basespace.illumina.com/docs/content/documentation/sdk-samples/python-sdk-overview)

## Access Token
An access token is required for accessing BaseSpace. This program reads the access token from the `BASESPACE_CREDENTIALS` environment variable. You can also store the access token in a JSON file and set the `BASESPACE_CREDENTIALS` variable to point to your json file.
```
{
    "access_token": "YOUR_ACCESS_TOKEN"
}
```

## Command
```
python -m Cancer.run basespace [-h] [--properties PROPERTIES [PROPERTIES ...]]
                        [--url [URL]]
                        [collection] [basespace_id]

positional arguments:
  collection            "projects", "runs", "samples" or "appsessions"
  basespace_id          ID of BaseSpace object

optional arguments:
  -h, --help            show this help message and exit
  --properties PROPERTIES [PROPERTIES ...]
                        Property name of BaseSpace object
  --url [URL]           Endpoint URL
```

For example, to list all projects:
```
BASESPACE_CREDENTIALS=YOUR_ACCESS_TOKEN python -m Cancer.run basespace projects
```
