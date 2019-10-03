Links to lab shiny apps
https://leeoesterreich.org/resources


README for lab shiny app maintenance 

#Links for setting up rstudio shiny server on digital ocean

https://www.digitalocean.com/community/tutorials/initial-server-setup-with-ubuntu-16-04
https://www.digitalocean.com/docs/droplets/how-to/connect-with-ssh/
https://www.digitalocean.com/community/tutorials/how-to-set-up-shiny-server-on-ubuntu-16-04#step-1-—-installing-shiny
https://deanattali.com/2015/05/09/setup-rstudio-shiny-server-digital-ocean/#user-libraries
https://www.r-bloggers.com/deploying-your-very-own-shiny-server/

#Update shiny apps

Edit repository to include new directory with ui.R, server.R, and any accompanying code

#ssh to droplet

cd /srv/shiny-server

git pull https://github.com/leeoesterreich/shiny-server/

#ensure necessary libraries are installed - can be performed after signing into the server by modifying example command below

sudo su - -c "R -e \"install.packages('data.table', repos='http://cran.rstudio.com/')\""
    

