#!/bin/bash

################################################################################
# This checks if USERID and GRPID were specified as environment variables      #
# (docker run -e USERID -e GRPID ....) when the container was started; if that #
# is the case, we create a user and group 'docker' with these IDs, give        #
# ownership of the output directories to them and use the 'gosu' command to    #
# start the main script as user 'docker' and pass to it all arguments          #
################################################################################

if [ ${USERID:-0} -ne 0 ] && [ ${GRPID:-0} -ne 0 ]; then \
    groupadd -g ${GRPID} docker &&\
    useradd -l -u ${USERID} -d /home/docker -m -g docker docker &&\
    chmod 775 /var/pipeline &&\
    chown --silent --no-dereference --recursive \
          --from=0:0 ${USERID}:${GRPID} \
        /var/pipeline &&\
    exec gosu docker ${PACKAGE_DIR}/DNA/MoCaSeq.sh "$@" \
;fi

exec ${PACKAGE_DIR}/DNA/MoCaSeq.sh "$@"
