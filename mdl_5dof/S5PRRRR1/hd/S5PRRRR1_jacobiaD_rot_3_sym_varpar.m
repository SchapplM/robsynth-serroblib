% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S5PRRRR1
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
%
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PRRRR1_jacobiaD_rot_3_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_jacobiaD_rot_3_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_jacobiaD_rot_3_sym_varpar: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_jacobiaD_rot_3_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:29:14
% EndTime: 2019-07-18 13:29:14
% DurationCPUTime: 0.04s
% Computational Cost: add. (30->4), mult. (82->10), div. (42->6), fcn. (80->4), ass. (0->8)
t50 = sin(qJ(2));
t52 = cos(qJ(2));
t61 = t50 ^ 2 / t52 ^ 2;
t64 = 0.1e1 + t61;
t38 = cos(atan2(0, -t52));
t36 = 0.1e1 / t38 ^ 2;
t32 = t36 * t61 + 0.1e1;
t1 = [0, 0, 0, 0, 0; 0, 0.2e1 * (-0.1e1 / t32 + t64 / t32 ^ 2 * t36) * t64 / t52 * t50 * qJD(2) / t38, 0, 0, 0; 0, 0, 0, 0, 0;];
JaD_rot  = t1;
