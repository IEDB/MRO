FROM obolibrary/odkfull

COPY requirements.txt /tools/mro-requirements.txt
RUN pip install -r /tools/mro-requirements.txt

RUN apt-get install sqlite3
