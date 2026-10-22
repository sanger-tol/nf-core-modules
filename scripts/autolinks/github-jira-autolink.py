import argparse
import json
import subprocess
from urllib.error import URLError
from urllib.request import urlopen

# Define the autolink reference pattern because GitHub's API doesn't
# support keys being a substring of another key e.g. DO- and DOC- overlap
key_prefixes = ["TOLIT-", "TOLP-", "TOLA-", "TOLSD-"]

jira_server = "https://jira.sanger.ac.uk/browse/"

# Set up argument parser
parser = argparse.ArgumentParser(description="Update Jira autolinks using an authenticated GitHub CLI.")
parser.add_argument("--gh", "--gh-path", default="gh", help="Path to the GitHub CLI executable (default: gh from PATH)")
parser.add_argument(
    "--repo", "-R", metavar="OWNER/REPO", help="A specific repo to update; defaults to all ToL pipelines"
)
parser.add_argument("--dry-run", action="store_true", help="Report missing required autolinks without making changes")

args = parser.parse_args()


def gh_autolink(repo_name, *options):
    try:
        return subprocess.check_output(
            [args.gh, "repo", "autolink", *options, "-R", repo_name],
            stderr=subprocess.PIPE,
            text=True,
        )
    except OSError as error:
        parser.exit(1, f'Failed to run GitHub CLI "{args.gh}": {error}\n')


if args.repo:
    if len(args.repo.split("/")) != 2 or not all(args.repo.split("/")):
        parser.error("--repo must be in owner/repo format")
    repos_to_update = [args.repo]
else:
    try:
        with urlopen("https://pipelines.tol.sanger.ac.uk/pipelines.json", timeout=30) as response:
            data = json.load(response)
        repos_to_update = [repo["full_name"] for repo in data["remote_workflows"]]
    except (URLError, OSError, ValueError, KeyError, TypeError) as error:
        parser.exit(1, f"Failed to fetch pipeline repositories: {error}\n")

# reverse sort prefixes to prevent overlaps (e.g. DO- and DOC-)
key_prefixes.sort(key=len, reverse=True)

failures = 0
for repo_name in repos_to_update:
    try:
        autolink_output = gh_autolink(repo_name, "list", "--json", "id,keyPrefix,urlTemplate,isAlphanumeric")
    except subprocess.CalledProcessError as error:
        print(f'Failed to fetch autolinks in repo "{repo_name}": {error.output}{error.stderr}')
        failures += 1
        continue

    try:
        autolinks = json.loads(autolink_output.strip() or "[]")
    except json.JSONDecodeError as error:
        print(f'Failed to parse autolinks in repo "{repo_name}": {error}')
        failures += 1
        continue

    existing_autolinks = {autolink["keyPrefix"]: autolink for autolink in autolinks}
    for prefix in key_prefixes:
        required_autolink = {
            "keyPrefix": prefix,
            "urlTemplate": f"{jira_server}{prefix}<num>",
            "isAlphanumeric": False,
        }

        existing = existing_autolinks.get(prefix)
        if existing is not None and all(existing.get(field) == value for field, value in required_autolink.items()):
            continue

        if args.dry_run:
            replacement = " (existing rule needs replacement)" if existing else ""
            print(
                f'Missing required autolink in repo "{repo_name}": '
                f"{prefix} -> {required_autolink['urlTemplate']} "
                f"(alphanumeric: false){replacement}"
            )
            continue

        if existing is not None:
            try:
                gh_autolink(repo_name, "delete", str(existing["id"]), "--yes")
            except subprocess.CalledProcessError as error:
                print(f'Failed to replace prefix "{prefix}" in repo "{repo_name}": {error.output}{error.stderr}')
                failures += 1
                continue

        try:
            gh_autolink(
                repo_name,
                "create",
                prefix,
                required_autolink["urlTemplate"],
                "--numeric",
            )
        except subprocess.CalledProcessError as error:
            print(f'Failed to create prefix "{prefix}" in repo "{repo_name}": {error.output}{error.stderr}')
            failures += 1
        else:
            print(f'Successfully created autolink reference in repo "{repo_name}" with prefix "{prefix}".')

print(f"Autolink references {'dry run' if args.dry_run else 'processing'} finished.")
if failures:
    parser.exit(1)
