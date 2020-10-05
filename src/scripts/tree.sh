#!/usr/bin/env bash
#
# This simple CGI script helps create a tree browser for ONTIE

cd ../..

URL="http://example.com?${QUERY_STRING}"
ID=$(urlp --query --query_field=id "${URL}")
PROJECT=$(urlp --query --query_field=project-name "${URL}")
BRANCH=$(urlp --query --query_field=branch-name "${URL}")

if [[ ${ID} ]]; then
	python3 -m gizmos.tree build/mro.db ${ID}
else
	python3 -m gizmos.tree build/mro.db
fi

echo "<a href=\"/${PROJECT}/branches/${BRANCH}\"><b>Return Home</b></a>"