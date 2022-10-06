# Docker Base Image for R Shiny Apps on FASTGenomics

This is the [FASTGenomics](https://beta.fastgenomics.org/) base image for [R Shiny](https://shiny.rstudio.com/) analyses.
Clone this repository if you want to bring your own Shiny App to FASTGenomics.
You can find an example Shiny analyses [on FASTGenomics](https://beta.fastgenomics.org/analyses/detail-analysis-d16bad01d96b4a08b25b5b68e718da49#Run), and a basic Shiny tutorial [here](https://shiny.rstudio.com/tutorial/).

See also our FASTGenomics documentation on [Working with Custom Docker Images](https://beta.fastgenomics.org/docs/analyses.html#working-with-custom-docker-images).

## Structure of this Repository

### Directories

- `APP`  
  All your shiny app code and runtime data has to go here.
  At least an `app.R` file has to be present.
- `fg_config`  
  **Don't modify this directory**  
  The files in this directory are required to correctly integrate your app into FASTGenomics.

### Files

- `Dockerfile`  
  The instructions on how to build your docker container.
  You do not _need_ to modify it, unless you want to
  - Use another base image
  - Install additional system packages
  - Copy additional data for the use during image building
- `install.R`  
  This `R`file is executed during image build.
  Usually it is used to install `R`packages.

# Build Your Own Shiny App and Publish it on FASTGenomics

1.  **Prerequisites.**
    You need to have [Docker Desktop](https://www.docker.com/) installed on your computer as well as an user account at Docker.

1.  **Clone this repository.**

1.  **Add your code and data.**
    All your code and data goes to the `APP` directory.
    The whole folder and it's contents are available at runtime.
    It is mandatory that there is at least an `app.R` file which will get executed at container runtime (replace the existing one).

1.  **Install additional packages.**
    If your app requires additional packages there are usually two places to install them:

    - `Dockerfile`  
       Add system libraries or packages via `apt-get install`. Just append the existing list.  
       You can also select another base image by changing the `FROM` at the top of the `Dockerfile`.
      However, the [`rocker/shiny-verse`](https://hub.docker.com/r/rocker/shiny-verse) is usually a good start as it contains shiny and all `tidy-verse` packages.
      Currently, we use it with `R` version 4.
    - `install.R`  
       In this file you add all your code to install required `R` packages.
      Common commands to use are `install.packages()` or `devtools::install_github()`.
      The script will be executed during the build of the Docker image.

1.  **Build your docker image.**
    To [build](https://docs.docker.com/engine/reference/commandline/build/) your Docker image just run the following command in the root directory of this repository:

    ```bash
    docker build -t <IMAGE_NAME> .
    ```

    `<IMAGE_NAME>` should be in the form `<DOCKER_HUB_USERNAME>/<APP_NAME>:<TAG>`.
    So if you want to build your app `fancy_shiny` as user `fastgenomics` for the first time, your tag could be version `1.0`.
    Then, the command would be:

    ```bash
    docker build -t fastgenomics/fancy_shiny:1.0 .
    ```

1.  **Test your docker image locally.**
    Once your image is build, you should test it locally.
    To start it just run:

    ```bash
    docker run -p 3838:3838 <IMAGE_NAME>
    ```

    The first port number after `-p` flag defines the port on your host machine.
    You can change it to something else, just make sure to use it when starting the app.
    The second port number is the default internal port of the shiny app and is set to `3838`.
    This must be the case for use in FASTGenomics.

    To see your app in action head over to your browser and type

    ```
    http://localhost:3838/proxy/shiny
    ```

    If anything goes wrong you can access the logs from within the docker container in the folder `/var/log/shiny-server`.

1.  **Push your image to Docker Hub.**
    The final step is to [push](https://docs.docker.com/engine/reference/commandline/push/) your image to Docker Hub by using:

    ```bash
    docker push <IMAGE_NAME>
    ```

    In our example this would be:

    ```bash
    docker push fastgenomics/fancy_shiny:1.0
    ```

    Please note that depending on your internet connection and the image size this can take a while.
    Alternatively, after testing you can also use [_automated builds_](https://docs.docker.com/docker-hub/builds/) on Docker Hub directly.

1.  **Create an analysis using your new image on FASTGenomics.**
    To use your new Shiny app in FASTGenomics, create a new analysis with your own image as described [here](https://beta.fastgenomics.org/docs/analyses.html#create-new-analyses).

# Further Notes

- Try to keep the Images size as small as possible.
  In the free-to-use version of FASTGenomics, the images has top be pulled from [Docker Hub](https://hub.docker.com/) every time it is started.
  The bigger your image, the longer it takes.
  So if you modify the `Dockerfile` always try to clean temporary files in each layer.
- FASTGenomics can only use images that are publicly available on [Docker Hub](https://hub.docker.com/).
  If you need a more private solution, please contact us via [e-mail](mailto:contact@fastgenomics.org) or [Slack](https://join.slack.com/t/fastgenomics/shared_invite/enQtNjU2ODk0OTk5MTA3LTkwZTgxN2EzYzAyMmExZTJiYmYxMjRhYjM2ODBiMWIwYmQ3MzZhYmIzZDkxZTI4OGFhYjQ4ODIzMTU3OWQ2NTc)
- If you need help with Docker, a good start are the official [docker docs](https://docs.docker.com/).
- You are responsible for the content of your custom images.
