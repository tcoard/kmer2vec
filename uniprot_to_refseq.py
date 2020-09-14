"""Get more information about the sequences, such as the organism's name, from uniprot"""
####################################################################################################
# TODO this file is still a work in progress
####################################################################################################
import json
import requests

requestURL = (
    "https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=100&accession=P21802"
)

r = requests.get(requestURL, headers={"Accept": "application/json"})

if not r.ok:
    r.raise_for_status()
    sys.exit()

payload = json.loads(r.text)
if len(payload) > 1:
    print("warning")
print(payload[0].get("dbReferences"))
