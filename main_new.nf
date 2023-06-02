#!/usr/bin/env nextflow

process sayHello {
    """
    echo 'Hello world!' > file
    env >> file
    """
}

workflow {
	main:
		sayHello()
}




