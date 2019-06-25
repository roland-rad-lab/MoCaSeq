#!/bin/sh
datetime=$(date);
timestamp=$(date +%s);
mem=`free -mh | awk '{print $3}'|sed -n 2p|grep -Eo '[+-]?[0-9]+([.][0-9]+)?'`
cpu=`iostat | grep -A 2 'avg-cpu'|sed -n 2p| awk '{print $1}'`
echo "$datetime \t $timestamp \t $mem \t $cpu"