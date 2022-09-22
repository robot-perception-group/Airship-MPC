#!/bin/bash
while sleep 1; do
	./LibrePilot/build/firmware/fw_simposix/fw_simposix.elf $1
done
