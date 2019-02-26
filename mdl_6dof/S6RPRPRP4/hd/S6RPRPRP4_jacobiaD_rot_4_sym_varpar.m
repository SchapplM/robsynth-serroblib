% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRP4_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_jacobiaD_rot_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:45:18
% EndTime: 2019-02-26 20:45:18
% DurationCPUTime: 0.51s
% Computational Cost: add. (1452->71), mult. (1839->155), div. (436->14), fcn. (2165->7), ass. (0->76)
t85 = qJ(1) + pkin(9);
t83 = sin(t85);
t143 = 0.2e1 * t83;
t78 = 0.1e1 / t83;
t84 = cos(t85);
t130 = t78 * t84;
t77 = t83 ^ 2;
t82 = t84 ^ 2;
t142 = qJD(1) * (0.1e1 / t77 * t82 + 0.1e1) * t130;
t92 = sin(qJ(3));
t127 = t83 * t92;
t93 = cos(qJ(3));
t73 = atan2(-t127, -t93);
t68 = sin(t73);
t114 = t68 * t127;
t69 = cos(t73);
t65 = -t69 * t93 - t114;
t62 = 0.1e1 / t65;
t89 = 0.1e1 / t93;
t63 = 0.1e1 / t65 ^ 2;
t141 = -0.2e1 * t92;
t87 = t92 ^ 2;
t90 = 0.1e1 / t93 ^ 2;
t124 = t87 * t90;
t76 = t77 * t124 + 0.1e1;
t74 = 0.1e1 / t76;
t140 = t74 - 0.1e1;
t121 = qJD(1) * t92;
t112 = t84 * t121;
t120 = qJD(3) * t83;
t131 = t69 * t92;
t119 = qJD(3) * t93;
t57 = (-(-t83 * t119 - t112) * t89 + t120 * t124) * t74;
t53 = (-t57 * t83 + qJD(3)) * t131 + (-t112 + (t57 - t120) * t93) * t68;
t139 = t53 * t62 * t63;
t138 = t57 * t68;
t137 = t57 * t92;
t136 = t63 * t84;
t135 = t63 * t92;
t86 = t92 * t87;
t88 = t93 ^ 2;
t100 = qJD(3) * (t86 / t88 + t92) * t89;
t122 = qJD(1) * t84;
t105 = t83 * t87 * t122;
t134 = (t77 * t100 + t90 * t105) / t76 ^ 2;
t110 = 0.1e1 + t124;
t67 = t110 * t83 * t74;
t133 = t67 * t83;
t132 = t68 * t93;
t79 = 0.1e1 / t83 ^ 2;
t129 = t79 * t82;
t128 = t82 * t87;
t126 = t84 * t92;
t125 = t87 * t89;
t123 = qJD(1) * t83;
t104 = t82 * t92 * t119;
t60 = t63 * t128 + 0.1e1;
t118 = 0.2e1 * (-t128 * t139 + (t104 - t105) * t63) / t60 ^ 2;
t117 = 0.2e1 * t139;
t72 = t88 * t129 + 0.1e1;
t116 = 0.2e1 * (-t79 * t104 - t88 * t142) / t72 ^ 2;
t115 = t63 * t126;
t113 = t74 * t125;
t111 = 0.1e1 + t129;
t109 = t92 * t118;
t108 = t134 * t143;
t107 = t134 * t141;
t106 = t83 * t113;
t103 = t111 * t92;
t102 = t110 * t84;
t70 = 0.1e1 / t72;
t58 = 0.1e1 / t60;
t56 = (t140 * t92 * t68 - t69 * t106) * t84;
t55 = -t83 * t132 + t131 + (-t69 * t127 + t132) * t67;
t54 = -t110 * t108 + (qJD(1) * t102 + t100 * t143) * t74;
t1 = [t84 * t89 * t107 + (-t83 * t89 * t121 + qJD(3) * t102) * t74, 0, t54, 0, 0, 0; (t62 * t109 + (-t62 * t119 + (qJD(1) * t56 + t53) * t135) * t58) * t83 + (t63 * t109 * t56 + (-((t57 * t106 + t140 * t119 + t107) * t68 + (t108 * t125 - t137 + (t137 + (-t86 * t90 + t141) * t120) * t74) * t69) * t115 + (t92 * t117 - t63 * t119) * t56 + (-t62 + ((-t77 + t82) * t69 * t113 + t140 * t114) * t63) * t121) * t58) * t84, 0 (t55 * t135 - t62 * t93) * t84 * t118 + ((-t62 * t123 + (-qJD(3) * t55 - t53) * t136) * t93 + (-t84 * qJD(3) * t62 - (-t54 * t69 * t83 + t68 * t120 + t133 * t138 - t138 + (-qJD(3) * t68 - t122 * t69) * t67) * t115 + (t84 * t117 + t63 * t123) * t55 - ((t54 - t122) * t68 + ((0.1e1 - t133) * qJD(3) + (t67 - t83) * t57) * t69) * t93 * t136) * t92) * t58, 0, 0, 0; t111 * t93 * t116 + (qJD(3) * t103 + 0.2e1 * t93 * t142) * t70, 0, t78 * t116 * t126 + (qJD(1) * t103 - t119 * t130) * t70, 0, 0, 0;];
JaD_rot  = t1;
