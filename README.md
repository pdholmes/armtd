# ARMTD
Just trying to clean up the rotatotope portions of the arm_planning repo.
I updated the JRS computation (`precompute_joint_reachable_sets`), the rotatotope class (see `rotatotope_v2`), the FRS classes (now a more general class called `robot_rotatotope_FRS`) and the planner (more general, see `robot_rotatotope_RTD_planner`).

I haven't generalized the code to multiple k's per joint but I figure we (hi Jon and Shannon) can work together on that!
Also, I'm confident everything works the same as before in 3D, but feedback about the 2D case would be appreciated!
Check also the scripts `test_fetch_rotatotope_FRS_obstacle` and `test_fetch_rotatotope_FRS_self_intersection` to see some self-contained examples checking that the constraints are reasonable.

Also, if there're portions that need more explanation let me know!
I tried to add comments that reference the correct portions of the paper.
Feel free to update this repo, then we can push to the arm_planning repo when we feel good.