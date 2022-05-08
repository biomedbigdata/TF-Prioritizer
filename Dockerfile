FROM ubuntu:20.04

RUN apt-get update && apt-get -y install sudo

RUN addgroup --gid 1000 docker && adduser --ingroup docker --disabled-password --uid 1000 docker

ARG DEBIAN_FRONTEND=noninteractive

RUN for i in \
    /srv/dependencies/ext \
    /srv/dependencies/scripts \
    /srv/dependencies/config_templates \
    /src/install \
    /srv/app \
    /srv/wd \
    ; do mkdir -p $i && chown -R 1000:1000 $i; done

COPY --chown=docker:docker install/ /srv/install
RUN bash -c "/srv/install/install.sh -d"

COPY --chown=docker:docker ext/ /srv/dependencies/ext
COPY --chown=docker:docker scripts /srv/dependencies/scripts
COPY --chown=docker:docker config_templates/defaultConfigs.json /srv/dependencies/config_templates/defaultConfigs.json
COPY --chown=docker:docker build/TFPRIO.jar /srv/app

USER docker

CMD ["java", "-jar", "/srv/app/TFPRIO.jar", "-c", "/srv/wd/input/configs.json", "-w", "/srv/wd", "-p", \
"/srv/dependencies"]
