#!/usr/bin/env bash
#
# This simple CGI script helps create a new Google Sheet
# and share it with the user.

echo "Content-Type: text/html"
echo ""
echo '<head>'
echo '<link href="https://cdn.jsdelivr.net/npm/bootstrap@5.0.0-beta2/dist/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-BmbxuPwQa2lc/FVzBcNJ7UAyJxM6wuqIj61tLrc4wSX0szH/Ev+nYRRuWlolflfl" crossorigin="anonymous">'
echo '</head>'
echo '<body>'
echo '<div class="container" style="padding-top:20px;">'

cd ../..

EMAIL_PAT="^[a-z0-9!#\$%&'*+/=?^_\`{|}~-]+(\.[a-z0-9!#$%&'*+/=?^_\`{|}~-]+)*@([a-z0-9]([a-z0-9-]*[a-z0-9])?\.)+[a-z0-9]([a-z0-9-]*[a-z0-9])?\$"
URL="http://example.com?${QUERY_STRING}" 
EMAIL=$(urlp --query --query_field=email "${URL}")
BRANCH=$(git branch --show-current)
INVALID=$(urlp --query --query_field=invalid "${URL}")

if [ -d .cogs ]; then
  LINK=$(cogs open)
  echo "<meta http-equiv='refresh' content='0; ${LINK}'/>"
  exit 0
fi

if [[ ${EMAIL} ]]; then
  TITLE="ONTIE ${BRANCH}"
  if [[ ${EMAIL} =~ ${EMAIL_PAT} ]]; then
    cogs init -c credentials.json -t "${TITLE}" -u ${EMAIL} -r writer || exit 1
    make load push || exit 1
  else
    echo '<meta http-equiv="refresh" content="0; ?invalid=true"/>'
  fi
fi

if [[ ${INVALID} ]]; then
  echo '<div class="alert alert-danger">Invalid email entered. Please try again.</div>'
fi
echo '<p>Please provide your Google Account name, so we can create and share a Google Sheet with you.</p>'
echo '<form action="">'
echo 'Google account (email):'
echo '<input name="email" type="text"/>'
echo '<input type="submit"/>'
echo '</form>'
echo '</div>'
echo '</body>'
