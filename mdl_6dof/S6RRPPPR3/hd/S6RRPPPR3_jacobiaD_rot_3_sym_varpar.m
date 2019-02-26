% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPPR3_jacobiaD_rot_3_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_jacobiaD_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_jacobiaD_rot_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_jacobiaD_rot_3_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:23:03
% EndTime: 2019-02-26 21:23:04
% DurationCPUTime: 0.50s
% Computational Cost: add. (776->71), mult. (1835->159), div. (470->13), fcn. (2177->7), ass. (0->71)
t82 = sin(qJ(1));
t113 = qJD(1) * t82;
t133 = 0.2e1 * t82;
t73 = t82 ^ 2;
t84 = cos(qJ(1));
t77 = t84 ^ 2;
t78 = 0.1e1 / t84;
t131 = (t73 / t77 + 0.1e1) * t78 * t113;
t81 = sin(qJ(2));
t114 = t82 * t81;
t83 = cos(qJ(2));
t63 = atan2(-t114, -t83);
t61 = sin(t63);
t105 = t61 * t114;
t62 = cos(t63);
t57 = -t62 * t83 - t105;
t54 = 0.1e1 / t57;
t74 = 0.1e1 / t83;
t55 = 0.1e1 / t57 ^ 2;
t75 = 0.1e1 / t83 ^ 2;
t130 = -0.2e1 * t81;
t71 = t81 ^ 2;
t118 = t71 * t75;
t68 = t73 * t118 + 0.1e1;
t64 = 0.1e1 / t68;
t129 = t64 - 0.1e1;
t112 = qJD(1) * t84;
t103 = t81 * t112;
t111 = qJD(2) * t82;
t120 = t62 * t81;
t110 = qJD(2) * t83;
t50 = (-(-t82 * t110 - t103) * t74 + t111 * t118) * t64;
t46 = (-t50 * t82 + qJD(2)) * t120 + (-t103 + (t50 - t111) * t83) * t61;
t128 = t46 * t54 * t55;
t127 = t50 * t61;
t126 = t50 * t81;
t125 = t55 * t81;
t124 = t55 * t84;
t115 = t74 * t81;
t70 = t81 * t71;
t76 = t74 * t75;
t92 = qJD(2) * (t70 * t76 + t115);
t96 = t71 * t82 * t112;
t123 = (t73 * t92 + t75 * t96) / t68 ^ 2;
t102 = 0.1e1 + t118;
t60 = t102 * t82 * t64;
t122 = t60 * t82;
t121 = t61 * t83;
t119 = t71 * t74;
t117 = t71 * t77;
t116 = t73 / t84 ^ 2;
t53 = t55 * t117 + 0.1e1;
t109 = 0.2e1 * (-t117 * t128 + (t77 * t81 * t110 - t96) * t55) / t53 ^ 2;
t108 = 0.2e1 * t128;
t69 = t75 * t116 + 0.1e1;
t107 = 0.2e1 * (t76 * qJD(2) * t81 * t116 + t75 * t131) / t69 ^ 2;
t106 = t81 * t124;
t104 = t64 * t119;
t101 = 0.1e1 + t116;
t100 = t81 * t109;
t99 = t123 * t130;
t98 = t123 * t133;
t97 = t82 * t104;
t95 = t102 * t84;
t93 = t101 * t81 * t75;
t66 = 0.1e1 / t69;
t51 = 0.1e1 / t53;
t49 = (t129 * t81 * t61 - t62 * t97) * t84;
t48 = -t82 * t121 + t120 + (-t62 * t114 + t121) * t60;
t47 = -t102 * t98 + (qJD(1) * t95 + t92 * t133) * t64;
t1 = [t84 * t74 * t99 + (qJD(2) * t95 - t113 * t115) * t64, t47, 0, 0, 0, 0; (t54 * t100 + (-t54 * t110 + (qJD(1) * t49 + t46) * t125) * t51) * t82 + (t55 * t100 * t49 + (-((t129 * t110 + t50 * t97 + t99) * t61 + (t98 * t119 - t126 + (t126 + (-t70 * t75 + t130) * t111) * t64) * t62) * t106 + (t81 * t108 - t55 * t110) * t49 + (-t54 + ((-t73 + t77) * t62 * t104 + t129 * t105) * t55) * t81 * qJD(1)) * t51) * t84 (t48 * t125 - t54 * t83) * t84 * t109 + ((-t54 * t113 + (-qJD(2) * t48 - t46) * t124) * t83 + (-t84 * qJD(2) * t54 - (-t47 * t62 * t82 + t61 * t111 + t122 * t127 - t127 + (-qJD(2) * t61 - t112 * t62) * t60) * t106 + (t84 * t108 + t55 * t113) * t48 - ((t47 - t112) * t61 + ((0.1e1 - t122) * qJD(2) + (t60 - t82) * t50) * t62) * t83 * t124) * t81) * t51, 0, 0, 0, 0; t101 * t74 * t107 + (-qJD(2) * t93 - 0.2e1 * t74 * t131) * t66, t78 * t75 * t107 * t114 + ((-0.2e1 * t71 * t76 - t74) * t78 * t111 - qJD(1) * t93) * t66, 0, 0, 0, 0;];
JaD_rot  = t1;
