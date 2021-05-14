
process bash_expand_path {

	input:
		val (path_string)

	output:
		stdout

	script:
	"""#!/usr/bin/env bash

for f in ${path_string};
do
	echo \${f}
done
	"""
}

