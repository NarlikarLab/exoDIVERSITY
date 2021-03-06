#!/bin/sh
bedver=$(bedtools --version);
successBEd=$([[ $bedver =~ \s" "v2\.2[0-9]\.[0-9] ]] && echo "1" || echo "Please keep bedtools version >=2.25 in your PATH");

pyver=$(python -c "import sys;t='{v[0]}.{v[1]}'.format(v=list(sys.version_info[:2]));sys.stdout.write(t)");
successpython=$([[ pyver=~2\.[6-9]* ]] && echo "1" || echo "Incompatible python version present. Python should be >= 2.6");
