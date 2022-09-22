#!/bin/bash
cd ../../submodules/pytorch_ssd_server
while sleep 1; do
	python3 server.py
done
