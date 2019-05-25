#!/bin/sh
datetime=$(date);
timestamp=$(date +%s);
mem=$(cat /sys/fs/cgroup/memory/memory.usage_in_bytes);
cpu=$(cat /sys/fs/cgroup/cpuacct/cpuacct.usage_all | awk 'BEGIN{FS=" "};{a=a + $2};END{printf "%.0f", a; print "\t" NR-1}');
echo "$datetime \t $timestamp \t $mem \t $cpu"
