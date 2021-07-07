#!/usr/bin/env bash
#
# This simple CGI script helps create a tree browser for ONTIE

cd ../..

URL="http://example.com?${QUERY_STRING}"
ID=$(urlp --query --query_field=id "${URL}")
TEXT=$(urlp --query --query_field=text "${URL}")
PROJECT=$(urlp --query --query_field=project-name "${URL}")
BRANCH=$(urlp --query --query_field=branch-name "${URL}")

# Check that the sqlite database exists
if ! [[ -s build/mro.db ]]; then
	rm build/mro.db > /dev/null 2>&1
	make build/mro.db > /dev/null 2>&1
fi

if [[ ${TEXT} ]]; then
	echo "Content-Type: application/json"
	echo ""
	python3 -m gizmos.search build/mro.db "${TEXT}"
else
	echo "Content-Type: text/html"
	echo ""
	echo "<a href=\"/${PROJECT}/branches/${BRANCH}\"><b>Return Home</b></a>"

	if [[ ${ID} ]]; then
		python3 -m gizmos.tree build/mro.db --include-search ${ID}
	else
		python3 -m gizmos.tree build/mro.db --include-search
	fi
fi
