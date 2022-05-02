# Academic page
Forked from https://github.com/academicpages/academicpages.github.io

## Run locally
```
bundle exec jekyll liveserve --config _config.yml,_config.dev.yml
```

Script
```
# preview.sh
#!/usr/bin/bash
rm -r ./_site
bundle exec jekyll liveserve --config _config.yml,_config.dev.yml
```

## Notes
```
# navigation settings
_data/navigation.yml
# absolute directories: add{{ site.baseurl}}
Test link: [pdb]({{ site.baseurl}}/files/pdb/6j83.pdb).
```
