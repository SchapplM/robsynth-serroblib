% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
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

function JaD_rot = S5PRRRR1_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_jacobiaD_rot_4_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_jacobiaD_rot_4_sym_varpar: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_jacobiaD_rot_4_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:29:14
% EndTime: 2019-07-18 13:29:14
% DurationCPUTime: 0.07s
% Computational Cost: add. (95->4), mult. (119->10), div. (57->6), fcn. (102->4), ass. (0->8)
t69 = sin(qJ(2));
t70 = cos(qJ(2));
t79 = t69 ^ 2 / t70 ^ 2;
t82 = 0.1e1 + t79;
t54 = cos(atan2(0, -t70));
t53 = 0.1e1 / t54 ^ 2;
t48 = t53 * t79 + 0.1e1;
t1 = [0, 0, 0, 0, 0; 0, 0.2e1 * (-0.1e1 / t48 + t82 / t48 ^ 2 * t53) * t82 / t70 * t69 * qJD(2) / t54, 0, 0, 0; 0, 0, 0, 0, 0;];
JaD_rot  = t1;
