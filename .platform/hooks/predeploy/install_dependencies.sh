sudo apt install openbabel
obabel_path=$(which openbabel)
sudo apt-get install jq
jsonStr=$(cat ./autogrow4-4.0.2/sample_submit_autogrow.json)
jq 'del(.mgltools_directory)' <<<"$jsonStr"
temp='. + { "obabel_path": "'
temp+=obabel_path
temp+='" }'
jq temp <<<"$jsonStr" > tmp.$$.json && mv tmp.$$.json ./autogrow4-4.0.2/sample_submit_autogrow.json