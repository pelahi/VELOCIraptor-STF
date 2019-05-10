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

create_dropbox_folder() {
	dropbox_dir="$1"
	curl -X POST https://api.dropboxapi.com/2/files/create_folder_v2 \
	    -H "Authorization: Bearer $DROPBOX_TOKEN" \
	    -H 'Content-Type: application/json' \
	    --data "{\"path\": \"$dropbox_dir\", \"autorename\": false}" > /dev/null
}

upload_to_dropbox() {
	# Upload and create a shareable link
	dropbox_path=$1/$2
	flavour=$3
	try curl -X POST https://content.dropboxapi.com/2/files/upload \
	    -H "Authorization: Bearer $DROPBOX_TOKEN" \
	    -H 'Content-Type: application/octet-stream' \
	    -H "Dropbox-API-Arg: {\"path\": \"$dropbox_path\"}" \
	    --data-binary @$2 > /dev/null
	try curl -X POST https://api.dropboxapi.com/2/sharing/create_shared_link_with_settings \
	    -H "Authorization: Bearer $DROPBOX_TOKEN" \
	    -H "Content-Type: application/json" \
	    --data "{\"path\": \"$dropbox_path\",\"settings\": {\"requested_visibility\": \"public\"}}" \
	    -o output.json > /dev/null
	url=`sed -n 's/.*"url": "\([^"]*\)?dl=0".*/\1/p' output.json`
	if [ $flavour == raw ]; then
		url+=?raw=1
	elif [ $flavour == download ]; then
		url+=?dl=1
	else
		url+=?dl=0
	fi
	echo $url
}

post_comment() {
	# Transform real newlines into escaped newlines
	comment="`echo "$1" | sed ':a; N; $!ba; s/\n/\\\\n/g'`"
	echo "Uploading comment to GitHub: $comment"
	try curl \
	    -H "Authorization: token ${GITHUB_TOKEN}" \
	    -X POST \
	    --data '{"body": "'"$comment"'"}' \
	    "https://api.github.com/repos/${TRAVIS_REPO_SLUG}/issues/${TRAVIS_PULL_REQUEST}/comments"
}

config_param_as_row() {
	sed -n "s/$2=\([^# ]*\).*/$2 | \\1/p" "$1"
}

run_vr() {
	config_file_url="$1"
	run_name=$2
	dropbox_dir="$3"

	title="`echo $run_name | sed '{h; p; x; s/./=/g}'`"

	# Get config file and run VR
	try wget --no-verbose -O $run_name.conf "$config_file_url"
	OMP_NUM_THREADS=2 try ./stf -C $run_name.conf -i EAGLE-L0012N0188-z0/snap_028_z000p000 -o $run_name -I 2 -s 16

	# Get relevant configuration parameters for display in comment
	config_table="Relevant configuration parameters:\n\n"
	config_table+="Parameter | Value\n--- | ---\n"
	config_table+="`config_param_as_row $run_name.conf Bound_halos`\n"
	config_table+="`config_param_as_row $run_name.conf Physical_linking_length`\n"
	config_table+="`config_param_as_row $run_name.conf FoF_Field_search_type`\n"
	config_table+="`config_param_as_row $run_name.conf Halo_6D_linking_length_factor`\n"
	config_table+="`config_param_as_row $run_name.conf Halo_6D_vel_linking_length_factor`\n"

	# M200 v/s R200
	try python ${TRAVIS_BUILD_DIR}/.travis/plots/m200_r200.py ${run_name}.properties.0 ${run_name}_m200_r200.png
	m200_r200_url=`upload_to_dropbox $dropbox_dir ${run_name}_m200_r200.png raw`
	m200_r200_comment='M200 v/s R200: ![m200-r200]('"$m200_r200_url"' \"M200 v/s R200\")'

	comment+="\n\n$title\n\n$config_table\n\n$m200_r200_comment"
}

cd ${TRAVIS_BUILD_DIR}/build

# Run unit tests first
# Need to update test so as to run VR in full on an input simulation.


# If this is a pull request, run VR against our testing data, produce a couple of standard plots
# and post the results back to the pull request
# We only have travis-ci.org configured, so let's check that the environment variables
# are set to distinguish when and when not to actually run this
if [ "$TRAVIS_PULL_REQUEST" != false -a "$MAKE_PLOTS" = true -a -n "$EAGLE_DATA_URL" ]; then

	# Get input data, prepare for plotting and upload
	try wget --no-verbose -O EAGLE-L0012N0188-z0.tar.gz "$EAGLE_DATA_URL"
	try tar xf EAGLE-L0012N0188-z0.tar.gz
	echo "backend: Agg" >> matplotlibrc
	dropbox_dir=/pr${TRAVIS_PULL_REQUEST}/${TRAVIS_BUILD_ID}
	create_dropbox_folder $dropbox_dir

	# Run VR using 3D and 6D configuration files. This implicitly
	# augments the $comment variable
	comment=''
	run_vr "$VR_3DCONFIG_URL" 3d-run "$dropbox_dir"
	run_vr "$VR_6DCONFIG_URL" 6d-run "$dropbox_dir"

	post_comment "$comment"
fi
