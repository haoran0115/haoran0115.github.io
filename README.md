# Academic page
Forked from https://github.com/academicpages/academicpages.github.io

## Install Bundler
Using conda to install bundler.
```
conda install -c conda-forge ruby
```

## Directory
```
./_config.yml   # general configs, including: profile photo & introduction
```

## Run locally
```
# run '. ~/.proxy' if needed
# run 'bundle add webrick' if jekyll not working
bundle exec jekyll liveserve --config _config.yml,_config.dev.yml --port=10000
```

Script
```
# preview.sh
#!/usr/bin/bash
rm -r ./_site
bundle exec jekyll liveserve --config _config.yml,_config.dev.yml --port=10000
```

## Notes
```
# navigation settings
_data/navigation.yml
# absolute directories: add{{ site.baseurl}}
Test link: [pdb]({{ site.baseurl}}/files/pdb/6j83.pdb).
```
