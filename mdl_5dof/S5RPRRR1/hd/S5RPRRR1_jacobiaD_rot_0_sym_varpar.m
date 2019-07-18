% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 0 (0=Basis) von
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:26
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPRRR1_jacobiaD_rot_0_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_jacobiaD_rot_0_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_jacobiaD_rot_0_sym_varpar: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_jacobiaD_rot_0_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_0_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:26:25
% EndTime: 2019-07-18 13:26:25
% DurationCPUTime: 0.01s
% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
JaD_rot  = t1;
