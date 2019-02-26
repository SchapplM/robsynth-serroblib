% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRPR8_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_jacobiaD_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:29:44
% EndTime: 2019-02-26 20:29:45
% DurationCPUTime: 0.48s
% Computational Cost: add. (1694->69), mult. (1839->158), div. (436->14), fcn. (2165->7), ass. (0->74)
t93 = cos(qJ(1));
t141 = 0.2e1 * t93;
t89 = 0.1e1 / t93;
t92 = sin(qJ(1));
t123 = t89 * t92;
t87 = t92 ^ 2;
t88 = t93 ^ 2;
t101 = qJD(1) * (t87 / t88 + 0.1e1) * t123;
t85 = pkin(9) + qJ(4);
t83 = sin(t85);
t119 = qJD(4) * t83;
t84 = cos(t85);
t104 = t84 * t87 * t119;
t90 = 0.1e1 / t93 ^ 2;
t124 = t87 * t90;
t77 = t83 ^ 2;
t76 = t124 * t77 + 0.1e1;
t140 = -0.2e1 * (t101 * t77 + t104 * t90) / t76 ^ 2;
t117 = qJD(4) * t93;
t121 = qJD(1) * t92;
t111 = t84 * t121;
t79 = 0.1e1 / t83 ^ 2;
t82 = t84 ^ 2;
t126 = t79 * t82;
t73 = t126 * t88 + 0.1e1;
t71 = 0.1e1 / t73;
t78 = 0.1e1 / t83;
t57 = ((t117 * t83 + t111) * t78 + t117 * t126) * t71;
t139 = -t57 + t117;
t122 = t93 * t84;
t70 = atan2(-t122, t83);
t68 = sin(t70);
t69 = cos(t70);
t64 = -t122 * t68 + t69 * t83;
t61 = 0.1e1 / t64;
t62 = 0.1e1 / t64 ^ 2;
t138 = t71 - 0.1e1;
t120 = qJD(1) * t93;
t105 = t82 * t92 * t120;
t125 = t82 * t87;
t129 = t69 * t84;
t53 = (-t57 * t93 + qJD(4)) * t129 + (t139 * t83 + t111) * t68;
t136 = t53 * t61 * t62;
t60 = t125 * t62 + 0.1e1;
t137 = (-t125 * t136 + (-t104 + t105) * t62) / t60 ^ 2;
t135 = t57 * t84;
t134 = t62 * t84;
t133 = t62 * t92;
t127 = t78 * t84;
t81 = t84 * t82;
t100 = qJD(4) * (-t78 / t77 * t81 - t127);
t132 = (t100 * t88 - t105 * t79) / t73 ^ 2;
t130 = t68 * t93;
t128 = t78 * t82;
t118 = qJD(4) * t92;
t116 = -0.2e1 * t136;
t115 = t84 * t137;
t114 = t84 * t133;
t113 = t84 * t132;
t112 = t71 * t128;
t110 = 0.1e1 + t126;
t109 = 0.1e1 + t124;
t108 = t132 * t141;
t107 = t93 * t112;
t106 = t138 * t84 * t68;
t103 = t110 * t92;
t102 = t109 * t84;
t74 = 0.1e1 / t76;
t66 = t110 * t93 * t71;
t58 = 0.1e1 / t60;
t56 = (-t107 * t69 - t106) * t92;
t55 = t83 * t130 + t129 + (-t122 * t69 - t68 * t83) * t66;
t54 = -t110 * t108 + (-qJD(1) * t103 + t100 * t141) * t71;
t1 = [-0.2e1 * t92 * t78 * t113 + (-qJD(4) * t103 + t120 * t127) * t71, 0, 0, t54, 0, 0; (0.2e1 * t61 * t115 + (t61 * t119 + (qJD(1) * t56 + t53) * t134) * t58) * t93 + (-0.2e1 * t62 * t115 * t56 + (((t107 * t57 + t119 * t138 + 0.2e1 * t113) * t68 + (t108 * t128 + t135 + (-t135 + (t79 * t81 + 0.2e1 * t84) * t117) * t71) * t69) * t114 + (t116 * t84 - t119 * t62) * t56 + (t61 + ((t87 - t88) * t69 * t112 - t93 * t106) * t62) * t84 * qJD(1)) * t58) * t92, 0, 0, 0.2e1 * (-t134 * t55 - t61 * t83) * t92 * t137 + ((t61 * t120 + (-qJD(4) * t55 - t53) * t133) * t83 + (t61 * t118 + (-t54 * t69 * t93 + t139 * t68 + (-qJD(4) * t68 + t121 * t69 + t130 * t57) * t66) * t114 + (t116 * t92 + t120 * t62) * t55 + ((-t54 - t121) * t68 + ((t66 * t93 - 0.1e1) * qJD(4) + (-t66 + t93) * t57) * t69) * t83 * t133) * t84) * t58, 0, 0; t109 * t83 * t140 + (qJD(4) * t102 + 0.2e1 * t101 * t83) * t74, 0, 0, t84 * t123 * t140 + (-t118 * t83 * t89 + qJD(1) * t102) * t74, 0, 0;];
JaD_rot  = t1;
