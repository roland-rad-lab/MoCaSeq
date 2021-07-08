#!/usr/bin/env bash

function moc_mktemp_file {
    local base=${1:-/tmp}
    local suff=""
    if [[ ! -z "${2}" ]]; then
        suff=".${2}"
    fi
    if [[ $(uname) = Darwin ]]; then mktemp $base/moc.XXXXXXXXXX${suff}
    else TMPDIR="$base" mktemp -t moc.XXXXXXXXXX${suff}
    fi
}

function extract_if_zip {
	local -n result=$2
	if [[ $1 =~ \.gz$ ]]; then
		zcat ${1} > $3
		result=$3
	else
		result=$1
	fi
}

