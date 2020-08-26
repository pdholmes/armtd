# ARMTD
Trying to clean up the rotatotope portions of the [arm_planning](https://github.com/ramvasudevan/arm_planning) repo.
To run this code, you'll need MATLAB R2018b or newer. Also, you will need [CORA_2018](https://tumcps.github.io/CORA/)... there's a newer version of CORA, but I haven't seen if it's compatible yet.

The file `fetch_example_new_code` is the same as `arm_example_11_rotatotope_RTD_planning_3D_fetch.m` in the other repo, and is a good place to start.
I updated the JRS computation (`precompute_joint_reachable_sets`), the rotatotope class (see `rotatotope_v2`), the FRS classes (now a more general class called `robot_rotatotope_FRS`) and the planner (more general, see `robot_rotatotope_RTD_planner`).

There are still some updates I want to make (like multiple k's per joint) that will come in the future.
I'm confident everything works the same as the [arm_planning](https://github.com/ramvasudevan/arm_planning) repo in 3D, but feedback about the 2D case or any bugs you find would be appreciated!
Check also the scripts `test_fetch_rotatotope_FRS_obstacle` and `test_fetch_rotatotope_FRS_self_intersection` to see some self-contained examples checking that the constraints are reasonable.

Also, if there're portions that need more explanation let me know!
I tried to add comments that reference the correct portions of the paper.

When running the code in this repo, you'll need to have the `rotatotopes` folder from this repo and the `simulator` and `simulator_files` from the [arm_planning](https://github.com/ramvasudevan/arm_planning) repo on your path!
