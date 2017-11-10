module load python3/3.6.2
if [ ! -d "env" ]; then
  python3 -m venv env
fi

source env/bin/activate
pip3 install -r requirements.txt
deactivate

module load agdc-py3-prod
export PYTHONPATH="${PYTHONPATH}":"$(pwd -P)/env/lib/python3.6/site-packages"
