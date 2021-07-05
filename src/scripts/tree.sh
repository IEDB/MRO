#!/usr/bin/env bash
#
# This simple CGI script helps create a tree browser for ONTIE

cd ../..

URL="http://example.com?${QUERY_STRING}"
ID=$(urlp --query --query_field=id "${URL}")
PROJECT=$(urlp --query --query_field=project-name "${URL}")
BRANCH=$(urlp --query --query_field=branch-name "${URL}")

# Check that the sqlite database exists
if ! [[ -s build/mro.db ]]; then
	rm build/mro.db > /dev/null 2>&1
	make build/mro.db > /dev/null 2>&1
fi

if [[ ${ID} ]]; then
	python3 -m gizmos.tree build/mro.db ${ID}
else
	python3 -m gizmos.tree build/mro.db
fi

echo "Content-Type: text/html"
echo ""
echo "<a href=\"/${PROJECT}/branches/${BRANCH}\"><b>Return Home</b></a>"