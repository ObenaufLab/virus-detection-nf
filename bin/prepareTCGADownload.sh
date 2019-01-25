mkdir -p $NAME
mkdir -p $NAME/raw

echo "/groups/obenauf/Software/gdc-client download -m gdc_manifest.txt -t /groups/obenauf/Tobias_Neumann/TCGA/token/gdc-user-token.2019-01-17T18_55_35.738Z.txt -n 100 --log-file gdc.log --no-related-files --no-segment-md5sums --debug" >> download.sh
