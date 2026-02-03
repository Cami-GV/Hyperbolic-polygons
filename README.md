# Masters-Code
This repository constains the Sage code I've implemented during my Master's Degree. It returns a graphic representing a compact Riemann surface endowed with the action of a given finite group, with a given signature.

### Running Sage code: beginner's guide

You can follow [Sage installation guide](https://doc.sagemath.org/html/en/installation/index.html#) and lauching it following [this instructions](https://doc.sagemath.org/html/en/installation/launching.html).
However, I have installed Sage through a Docker image, and I have run it into a jupyter notebook in VS Code. Next, I will detail the process I followed.

1. [Install Docker](https://docs.docker.com/engine/install/) and launch it.

2. Look in [Docker Hub](https://hub.docker.com/) a Sagemath image. I recommend the versions in sagemath/sagemath-dev.

3. Select a tag and run the pull command.
  ```sh
  docker pull sagemath/sagemath-dev:<version>
  ```

4. Create and run a container. The -it option stands for "interactive", and -p is for setting up a port
  ```sh
  docker run -it --name <container_name> -p <host_port>:<container_port> <image_name>
  ```

5. To launch a jupyter notebook, run
```sh
sage -n jupyter --ip 0.0.0.0 --port <port_number> --no-browser --allow-root --NotebookApp.password=''
```
and then copy the shown url.

6. Open a jupyter notebook in VS Code, click "Select Kernel" and paste the url.

7. Done! When you finish your work session, use
```sh
docker stop
```
to stop the container. To initialize it, use 
```sh
docker start -i <container_name>
```
and run the command in step 5.

I refer you to the `Usage.ipynb` file to consult an example of my code's usage.

## References

> Behn A., Rodíguez R.E, & Rojas A.M. (2013).
> Adapted hyperbolic polygons and symplectic representations for gtoup actions on Riemann surfaces.
> Journal of Pure and Applied Algebra, 217(3), 409-426.
