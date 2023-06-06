import requests
import json
import argparse

# Define the autolink reference pattern because GitHub's API doesn't
# support keys being a substring of another key e.g. DO- and DOC- overlap
key_prefixes = ['AB-1a23v', 'ABC-1a23v', "AC-1a23v"]

# Set up argument parser
parser = argparse.ArgumentParser(
    description='Update repos with GitHub token and organization name.')
parser.add_argument('--token', required=True, help='Your GitHub token')
parser.add_argument('--org', required=True,
                    help='Your organization name on GitHub')
parser.add_argument('--repo', help='A specific repo to update')

args = parser.parse_args()

# Get your tokens and organization name
GITHUB_TOKEN = args.token
GITHUB_ORG_NAME = args.org

# Set the headers
github_headers = {
    'Authorization': f'token {GITHUB_TOKEN}',
    'Accept': 'application/vnd.github.v3+json'
}

# List of repositories to update
repos_to_update = []

# If a specific repo is specified, update only this one
if args.repo:
    repos_to_update.append(args.repo)
else:
    # Get the repositories of the organization
    page = 1
    while True:
        response = requests.get(
            f'https://api.github.com/orgs/{GITHUB_ORG_NAME}/repos?page={page}&per_page=30', headers=github_headers)

        # Ensure we got a successful response
        if response.status_code != 200:
            raise Exception(
                f'Failed to fetch repos from GitHub: {response.text}')

        # Get the JSON data from the response
        data = response.json()

        # Get the names of the repositories
        repos = [repo['name'] for repo in data]

        # If there are no more repositories, break the loop
        if not repos:
            break

        repos_to_update.extend(repos)
        page += 1

# reverse sort prefixes to prevent overlaps (e.g. DO- and DOC-)
key_prefixes.sort(key=len, reverse=True)

# Define the autolink reference pattern
# Append '-1a23v' to each prefix
for repo_name in repos_to_update:
    for prefix in key_prefixes:
        # Create the payload for the GitHub request
        autolink_payload = {
            'key_prefix': prefix,  # The key_prefix now includes '-1a23v'
            # Remove '1a23v' in the URL template
            'url_template': f"https://bitgoinc.atlassian.net/browse/{prefix.rstrip('1a23v')}<num>"
        }

        # Create the autolink reference in the repository
        autolink_response = requests.post(f'https://api.github.com/repos/{GITHUB_ORG_NAME}/{repo_name}/autolinks',
                                          headers=github_headers, data=json.dumps(autolink_payload))

        # Check if we got a successful response
        if autolink_response.status_code == 201:
            print(
                f'Successfully created autolink reference in repo "{repo_name}" with prefix "{prefix}".')
        elif autolink_response.status_code == 422 and "already_exists" in autolink_response.text:
            print(
                f'Skipping prefix "{prefix}" in repo "{repo_name}", it already exists.')
        else:
            print(
                f'Failed to create autolink reference in repo "{repo_name}": {autolink_response.text}')

print('Autolink references processing finished.')
