% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR1_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_jacobiaD_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:28:06
% EndTime: 2019-02-26 21:28:07
% DurationCPUTime: 0.54s
% Computational Cost: add. (1893->71), mult. (1835->159), div. (470->13), fcn. (2177->7), ass. (0->71)
t92 = sin(qJ(1));
t122 = qJD(1) * t92;
t141 = 0.2e1 * t92;
t119 = qJD(2) * t92;
t93 = cos(qJ(1));
t121 = qJD(1) * t93;
t85 = qJ(2) + pkin(10);
t83 = sin(t85);
t112 = t83 * t121;
t79 = t83 ^ 2;
t84 = cos(t85);
t81 = 0.1e1 / t84 ^ 2;
t127 = t79 * t81;
t87 = t92 ^ 2;
t74 = t87 * t127 + 0.1e1;
t72 = 0.1e1 / t74;
t80 = 0.1e1 / t84;
t58 = (-(-t84 * t119 - t112) * t80 + t119 * t127) * t72;
t139 = t58 - t119;
t88 = t93 ^ 2;
t89 = 0.1e1 / t93;
t138 = (t87 / t88 + 0.1e1) * t89 * t122;
t123 = t92 * t83;
t71 = atan2(-t123, -t84);
t69 = sin(t71);
t114 = t69 * t123;
t70 = cos(t71);
t65 = -t70 * t84 - t114;
t62 = 0.1e1 / t65;
t63 = 0.1e1 / t65 ^ 2;
t137 = -0.2e1 * t83;
t136 = t72 - 0.1e1;
t129 = t70 * t83;
t54 = (-t58 * t92 + qJD(2)) * t129 + (t139 * t84 - t112) * t69;
t135 = t54 * t62 * t63;
t134 = t58 * t83;
t133 = t63 * t83;
t132 = t63 * t93;
t125 = t80 * t83;
t78 = t83 * t79;
t82 = t80 * t81;
t101 = qJD(2) * (t78 * t82 + t125);
t105 = t79 * t92 * t121;
t131 = (t87 * t101 + t81 * t105) / t74 ^ 2;
t130 = t69 * t92;
t128 = t79 * t80;
t126 = t79 * t88;
t124 = t87 / t93 ^ 2;
t120 = qJD(2) * t84;
t61 = t63 * t126 + 0.1e1;
t118 = 0.2e1 * (-t126 * t135 + (t83 * t88 * t120 - t105) * t63) / t61 ^ 2;
t117 = 0.2e1 * t135;
t77 = t81 * t124 + 0.1e1;
t116 = 0.2e1 * (t82 * qJD(2) * t83 * t124 + t138 * t81) / t77 ^ 2;
t115 = t83 * t132;
t113 = t72 * t128;
t111 = 0.1e1 + t127;
t110 = 0.1e1 + t124;
t109 = t83 * t118;
t108 = t131 * t137;
t107 = t131 * t141;
t106 = t92 * t113;
t104 = t111 * t93;
t102 = t81 * t83 * t110;
t75 = 0.1e1 / t77;
t67 = t111 * t92 * t72;
t59 = 0.1e1 / t61;
t57 = (t136 * t83 * t69 - t70 * t106) * t93;
t56 = -t84 * t130 + t129 + (-t70 * t123 + t69 * t84) * t67;
t55 = -t111 * t107 + (qJD(1) * t104 + t101 * t141) * t72;
t1 = [t93 * t80 * t108 + (qJD(2) * t104 - t122 * t125) * t72, t55, 0, 0, 0, 0; (t62 * t109 + (-t62 * t120 + (qJD(1) * t57 + t54) * t133) * t59) * t92 + (t63 * t109 * t57 + (-((t58 * t106 + t136 * t120 + t108) * t69 + (t107 * t128 - t134 + (t134 + (-t78 * t81 + t137) * t119) * t72) * t70) * t115 + (t83 * t117 - t63 * t120) * t57 + (-t62 + ((-t87 + t88) * t70 * t113 + t136 * t114) * t63) * t83 * qJD(1)) * t59) * t93 (t56 * t133 - t62 * t84) * t93 * t118 + ((-t62 * t122 + (-qJD(2) * t56 - t54) * t132) * t84 + (-t93 * qJD(2) * t62 - (-t55 * t70 * t92 - t139 * t69 + (-qJD(2) * t69 - t121 * t70 + t130 * t58) * t67) * t115 + (t93 * t117 + t63 * t122) * t56 - ((t55 - t121) * t69 + ((-t67 * t92 + 0.1e1) * qJD(2) + (t67 - t92) * t58) * t70) * t84 * t132) * t83) * t59, 0, 0, 0, 0; t110 * t80 * t116 + (-qJD(2) * t102 - 0.2e1 * t138 * t80) * t75, t89 * t81 * t116 * t123 + ((-0.2e1 * t79 * t82 - t80) * t89 * t119 - qJD(1) * t102) * t75, 0, 0, 0, 0;];
JaD_rot  = t1;
