#!/usr/bin/env bash
#
# This simple CGI script helps create a new Google Sheet
# and share it with the user.

echo "Content-Type: text/html"
echo ""

cd ../..

URL="http://example.com?${QUERY_STRING}" 
EMAIL=$(urlp --query --query_field=email "${URL}")
PROJECT=$(urlp --query --query_field=project-name "${URL}")
BRANCH=$(urlp --query --query_field=branch-name "${URL}")

if [[ ${EMAIL} ]]; then
  echo "<p>Creating a Google Sheet and sharing it with ${EMAIL} ...</p>"
  TITLE="MRO ${BRANCH}"
  cogs init -t "${TITLE}" -u ${EMAIL} -r writer
  make load push
fi

if [ -d .cogs ]; then
  LINK=$(cogs open)
  echo "<ul><li><a href='${LINK}' target='_blank'>Open Google Sheet</a><br></li>"
  echo "<li><a href='/${PROJECT}/branches/${BRANCH}'>Return to ${BRANCH}</a></li></ul>"
  echo "<script type='text/javascript'>window.open('${LINK}');</script>"
  # echo "<meta http-equiv='refresh' content='0; /${PROJECT}/branches/${BRANCH}'/>"
  exit 0
fi

echo '<p>Please provide your Google Account name, so we can create and share a Google Sheet with you.</p>'
echo '<form action="">'
echo 'Google account (email):'
echo '<input name="email" type="text"/>'
echo '<input type="submit"/>'
echo '</form>'
