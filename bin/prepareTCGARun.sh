#!/usr/bin/env bash

# MIT License
#
# Copyright (c) 2018 Tobias Neumann
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

NAME=$1
REPODIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Prepare directory structure

mkdir -p $NAME
mkdir -p $NAME/analysis
mkdir -p $NAME/virus-integration/
mkdir -p $NAME/virus-integration/centrifuge
mkdir -p $NAME/virus-integration/salmon
mkdir -p $NAME/virus-integration/TCGAconversion

# s3 upload

cp ${REPODIR}/s3templates/upload.s3 $NAME/s3.sh
sed -i "s/SAMPLE/$NAME/g" $NAME/s3.sh

# Nextflow AWS runs

cp ${REPODIR}/NFtemplates/conversionTemplate.nf $NAME/virus-integration/TCGAconversion/callNextflow.sh
sed -i "s/SAMPLE/$NAME/g" $NAME/virus-integration/TCGAconversion/callNextflow.sh

cp ${REPODIR}/NFtemplates/centrifugeTemplate.nf $NAME/virus-integration/centrifuge/callNextflow.sh
sed -i "s/SAMPLE/$NAME/g" $NAME/virus-integration/centrifuge/callNextflow.sh

cp ${REPODIR}/NFtemplates/salmonTemplate.nf $NAME/virus-integration/salmon/callNextflow.sh
sed -i "s/SAMPLE/$NAME/g" $NAME/virus-integration/salmon/callNextflow.sh

# s3 fetch results

cp ${REPODIR}/s3templates/downloadConversion.s3 $NAME/virus-integration/TCGAconversion/s3.sh
sed -i "s/SAMPLE/$NAME/g" $NAME/virus-integration/TCGAconversion/s3.sh

cp ${REPODIR}/s3templates/downloadCentrifuge.s3 $NAME/virus-integration/centrifuge/s3.sh
sed -i "s/SAMPLE/$NAME/g" $NAME/virus-integration/centrifuge/s3.sh

cp ${REPODIR}/s3templates/downloadSalmon.s3 $NAME/virus-integration/salmon/s3.sh
sed -i "s/SAMPLE/$NAME/g" $NAME/virus-integration/salmon/s3.sh
