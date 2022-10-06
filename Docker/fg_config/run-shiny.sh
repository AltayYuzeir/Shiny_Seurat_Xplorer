#!/bin/bash

echo "configuring session $1"
sed -i -e "s/{sessionsId}/$1/g" /etc/shiny-server/shiny-server.conf

echo "starting shiny-server"
/usr/bin/shiny-server
