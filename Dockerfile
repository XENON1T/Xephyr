FROM rootproject/root-fedora

# RUN python3 -m pip install --upgrade pip # This causes error, skipping

RUN python3 -m pip install --ignore-installed tornado==4.4.3
RUN python3 -m pip install jupyter # install jupyter
# RUN python3 -m pip install prompt-toolkit==1.0.15 # deal with prompt-toolkit issue (downgrade)
RUN python3 -m pip install --upgrade ipykernel

# Create the xephyrian user environemnt
RUN useradd -ms /bin/bash xephyrian
# copy xephyr files
ADD ./ /home/xephyrian/Xephyr
WORKDIR /home/xephyrian
RUN mkdir xephyr_projects
RUN chown -R xephyrian /home/xephyrian
ENV XEPHYR_DIR /home/xephyrian/
USER xephyrian

# compile xephyr libraries
RUN root -b "/home/xephyrian/Xephyr/loadXephyr.C"

EXPOSE 8080
# Run jupyter
CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--port=8080", "--allow-root"]
