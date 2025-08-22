# Extension of Bioconductor Docker Image

This folder contains the `Dockerfile` used to extend the official [`bioconductor/bioconductor_docker:RELEASE_3_21`](https://hub.docker.com/r/bioconductor/bioconductor_docker) image with additional R and Bioconductor packages required by this pipeline.


## Build the image

From this directory, run:

```bash
docker build -t bioliners/bioconductor_r3_21-vc:latest .
docker login
docker push bioliners/bioconductor_r3_21-vc:latest
```