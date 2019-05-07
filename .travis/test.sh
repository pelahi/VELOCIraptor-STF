#!/bin/bash
#
# Travis CI test script
#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2018
# Copyright by UWA (in the framework of the ICRAR)
# All rights reserved
#
# Contributed by Pascal Elahi
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA 02111-1307  USA
#

try() {
	"$@"
	status=$?
	if [ $status -ne 0 ]; then
		echo "Command exited with status $status, aborting build now: $@" 1>&2
		exit 1
	fi
}

upload_to_dropbox() {
	# Upload and create a shareable link
	dropbox_dir=/pr${TRAVIS_PULL_REQUEST}
	dropbox_path=$dropbox_dir/`date +%s`_$1
	curl -X POST https://api.dropboxapi.com/2/files/create_folder_v2 \
	    -H "Authorization: Bearer $DROPBOX_TOKEN" \
	    -H 'Content-Type: application/json' \
	    --data "{\"path\": \"$dropbox_dir\", \"autorename\": false}" > /dev/null
	try curl -X POST https://content.dropboxapi.com/2/files/upload \
	    -H "Authorization: Bearer $DROPBOX_TOKEN" \
	    -H 'Content-Type: application/octet-stream' \
	    -H "Dropbox-API-Arg: {\"path\": \"$dropbox_path\"}" \
	    --data-binary @$1 > /dev/null
	try curl -X POST https://api.dropboxapi.com/2/sharing/create_shared_link_with_settings \
	    -H "Authorization: Bearer $DROPBOX_TOKEN" \
	    -H "Content-Type: application/json" \
	    --data "{\"path\": \"$dropbox_path\",\"settings\": {\"requested_visibility\": \"public\"}}" \
	    -o output.json > /dev/null
	url=`sed -n 's/.*"url": "\([^"]*\)?dl=0".*/\1/p' output.json`?raw=1
	echo $url
}

comment_on_github() {
	echo "Uploading comment to GitHub: $1"
	try curl \
	    -H "Authorization: token ${GITHUB_TOKEN}" \
	    -X POST \
	    -d "{\"body\": \"$1\"}" \
	    "https://api.github.com/repos/${TRAVIS_REPO_SLUG}/issues/${TRAVIS_PULL_REQUEST}/comments"
}

cd ${TRAVIS_BUILD_DIR}/build

# Run unit tests first
# Need to update test so as to run VR in full on an input simulation.


# If this is a pull request, run VR against our testing data, produce a couple of standard plots
# and post the results back to the pull request
# We only have travis-ci.org configured, so let's check that the environment variables
# are set to distinguish when and when not to actually run this
if [ "$TRAVIS_PULL_REQUEST" != false -a "$MAKE_PLOTS" = true -a -n "$VR_CONFIG_URL" ]; then
	try wget --no-verbose -O vr.conf "$VR_CONFIG_URL"
	try wget --no-verbose -O EAGLE-L0012N0188-z0.tar.gz "$EAGLE_DATA_URL"
	try tar xf EAGLE-L0012N0188-z0.tar.gz
	OMP_NUM_THREADS=2 try ./stf -C vr.conf -i EAGLE-L0012N0188-z0/snap_028_z000p000 -o small-run -I 2 -s 16

	# Produce a small plot
	echo "backend: Agg" >> matplotlibrc
	try python <<EOF
import h5py
from matplotlib import pyplot
f = h5py.File('small-run.properties.0')
m200, r200 = f['/Mass_200mean'], f['/R_200mean']
pyplot.xscale('log'); pyplot.yscale('log')
pyplot.plot(r200, m200, marker='.', linestyle=' ')
pyplot.savefig('mass_radius.png')
EOF
	mass_radius_url=`upload_to_dropbox mass_radius.png`
	comment='Mass Radius relationship: ![mass-radius-rel]('"$mass_radius_url"' \"Mass Radius relationship\")'
	comment_on_github "$comment"
fi
