#!/usr/bin/env bash

watch -n 1 "sacct --format='JobID,JobName%30,Partition%20,State,Elapsed,User%10,NodeList'"
