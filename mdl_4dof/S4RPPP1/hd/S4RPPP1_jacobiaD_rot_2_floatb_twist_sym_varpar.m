% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 2 (0=Basis) von
% S4RPPP1
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
%
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:45
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_rot = S4RPPP1_jacobiaD_rot_2_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_jacobiaD_rot_2_floatb_twist_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_jacobiaD_rot_2_floatb_twist_sym_varpar: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_jacobiaD_rot_2_floatb_twist_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_2_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:47
% EndTime: 2018-11-14 13:45:47
% DurationCPUTime: 0.24s
% Computational Cost: add. (275->34), mult. (660->98), div. (108->12), fcn. (792->13), ass. (0->50)
t88 = sin(pkin(4));
t81 = t88 ^ 2;
t90 = cos(pkin(4));
t83 = 0.1e1 / t90 ^ 2;
t92 = cos(qJ(1));
t86 = t92 ^ 2;
t77 = t81 * t83 * t86 + 0.1e1;
t91 = sin(qJ(1));
t85 = t91 ^ 2;
t107 = 0.1e1 / t77 ^ 2 * t85;
t111 = t107 * t83;
t105 = t92 * t88;
t76 = atan2(t105, t90);
t72 = sin(t76);
t73 = cos(t76);
t59 = t72 * t105 + t73 * t90;
t56 = 0.1e1 / t59;
t78 = pkin(4) + pkin(6);
t79 = pkin(4) - pkin(6);
t70 = sin(t78) / 0.2e1 - sin(t79) / 0.2e1;
t89 = cos(pkin(6));
t101 = t91 * t70 - t89 * t92;
t63 = 0.1e1 / t101;
t82 = 0.1e1 / t90;
t57 = 0.1e1 / t59 ^ 2;
t64 = 0.1e1 / t101 ^ 2;
t110 = t57 * t91;
t71 = cos(t79) / 0.2e1 + cos(t78) / 0.2e1;
t87 = sin(pkin(6));
t68 = t91 * t71 + t87 * t92;
t109 = t64 * t68;
t67 = -t70 * t92 - t91 * t89;
t108 = t67 * t68;
t106 = t81 * t82;
t104 = qJD(1) * t92;
t74 = 0.1e1 / t77;
t103 = (t74 - 0.1e1) * t88;
t102 = -0.2e1 * t82 * t111;
t66 = t71 * t92 - t91 * t87;
t49 = (-t73 * t74 * t92 * t106 + t72 * t103) * t91;
t80 = t88 * t81;
t65 = t63 * t64;
t62 = t68 ^ 2;
t61 = t67 * qJD(1);
t60 = t66 * qJD(1);
t58 = t56 * t57;
t55 = t57 * t81 * t85 + 0.1e1;
t52 = t62 * t64 + 0.1e1;
t48 = qJD(1) * t49;
t1 = [(-t74 * t82 * t88 + t80 * t102) * t104, 0, 0, 0; (0.2e1 * (t49 * t110 - t56 * t92) / t55 ^ 2 * (-t48 * t58 * t85 + t104 * t110) * t81 + ((0.2e1 * t49 * t58 * t91 - t57 * t92) * t48 + (-t91 * t56 + ((-t49 + (-t80 * t111 - t103) * t91 * t72) * t92 - (t86 * t81 ^ 2 * t102 + (-t107 + (0.2e1 * t85 - t86) * t74) * t106) * t91 * t73) * t57) * qJD(1)) / t55) * t88, 0, 0, 0; 0.2e1 * (t64 * t108 + t63 * t66) / t52 ^ 2 * (t61 * t62 * t65 + t60 * t109) + (-t67 * t60 * t64 + (-0.2e1 * t65 * t108 - t66 * t64) * t61 + (-t101 * t109 + t68 * t63) * qJD(1)) / t52, 0, 0, 0;];
JaD_rot  = t1;
