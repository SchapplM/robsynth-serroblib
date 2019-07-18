% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4RRPR2
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
%
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S4RRPR2_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_jacobiaD_rot_4_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_jacobiaD_rot_4_sym_varpar: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_jacobiaD_rot_4_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:41
% EndTime: 2019-07-18 18:16:41
% DurationCPUTime: 0.10s
% Computational Cost: add. (514->14), mult. (492->30), div. (62->4), fcn. (544->4), ass. (0->23)
t79 = qJ(1) + qJ(2);
t68 = cos(t79);
t70 = sin(qJ(4));
t76 = sin(t79);
t84 = cos(qJ(4));
t74 = t68 * t84 + t76 * t70;
t58 = 0.1e1 / t74 ^ 2;
t86 = t58 * t74;
t85 = qJD(4) - qJD(1) - qJD(2);
t57 = 0.1e1 / t74;
t73 = -t68 * t70 + t76 * t84;
t52 = t85 * t73;
t56 = t73 ^ 2;
t81 = t56 * t58;
t55 = 0.1e1 + t81;
t82 = t85 * t86;
t78 = t73 * t82;
t59 = t57 * t58;
t80 = t56 * t59;
t83 = (-t52 * t80 - t78) / t55 ^ 2;
t53 = 0.1e1 / t55;
t49 = 0.2e1 * (t57 * t74 + t81) * t83 + (0.2e1 * t78 + (-t57 + 0.2e1 * t80 + t86) * t52) * t53;
t1 = [0, 0, 0, 0; 0, 0, 0, 0; t49, t49, 0, -0.2e1 * t83 - 0.2e1 * (t53 * t82 - (-t52 * t53 * t59 - t58 * t83) * t73) * t73;];
JaD_rot  = t1;
