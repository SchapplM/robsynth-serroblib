% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPPR2_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_jacobiaD_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:22:26
% EndTime: 2019-02-26 21:22:27
% DurationCPUTime: 0.50s
% Computational Cost: add. (1883->70), mult. (1839->155), div. (436->14), fcn. (2165->7), ass. (0->75)
t95 = sin(qJ(1));
t145 = 0.2e1 * t95;
t90 = 0.1e1 / t95;
t96 = cos(qJ(1));
t129 = t90 * t96;
t123 = qJD(2) * t95;
t125 = qJD(1) * t96;
t88 = qJ(2) + pkin(9);
t86 = sin(t88);
t115 = t86 * t125;
t81 = t86 ^ 2;
t87 = cos(t88);
t84 = 0.1e1 / t87 ^ 2;
t132 = t81 * t84;
t89 = t95 ^ 2;
t76 = t89 * t132 + 0.1e1;
t74 = 0.1e1 / t76;
t83 = 0.1e1 / t87;
t60 = (-(-t87 * t123 - t115) * t83 + t123 * t132) * t74;
t144 = t60 - t123;
t94 = t96 ^ 2;
t143 = qJD(1) * (0.1e1 / t89 * t94 + 0.1e1) * t129;
t127 = t95 * t86;
t73 = atan2(-t127, -t87);
t71 = sin(t73);
t117 = t71 * t127;
t72 = cos(t73);
t67 = -t72 * t87 - t117;
t64 = 0.1e1 / t67;
t65 = 0.1e1 / t67 ^ 2;
t142 = -0.2e1 * t86;
t141 = t74 - 0.1e1;
t134 = t72 * t86;
t56 = (-t60 * t95 + qJD(2)) * t134 + (t144 * t87 - t115) * t71;
t140 = t56 * t64 * t65;
t139 = t60 * t86;
t138 = t65 * t86;
t137 = t65 * t96;
t130 = t83 * t86;
t80 = t86 * t81;
t82 = t87 ^ 2;
t103 = qJD(2) * (t80 * t83 / t82 + t130);
t108 = t81 * t95 * t125;
t136 = (t89 * t103 + t84 * t108) / t76 ^ 2;
t135 = t71 * t95;
t133 = t81 * t83;
t131 = t81 * t94;
t91 = 0.1e1 / t95 ^ 2;
t128 = t91 * t94;
t126 = qJD(1) * t95;
t124 = qJD(2) * t87;
t122 = qJD(2) * t96;
t107 = t86 * t94 * t124;
t63 = t65 * t131 + 0.1e1;
t121 = 0.2e1 * (-t131 * t140 + (t107 - t108) * t65) / t63 ^ 2;
t120 = 0.2e1 * t140;
t79 = t82 * t128 + 0.1e1;
t119 = 0.2e1 * (-t91 * t107 - t82 * t143) / t79 ^ 2;
t118 = t86 * t137;
t116 = t74 * t133;
t114 = 0.1e1 + t132;
t113 = 0.1e1 + t128;
t112 = t86 * t121;
t111 = t136 * t142;
t110 = t136 * t145;
t109 = t95 * t116;
t106 = t114 * t96;
t105 = t113 * t86;
t77 = 0.1e1 / t79;
t69 = t114 * t95 * t74;
t61 = 0.1e1 / t63;
t59 = (t141 * t86 * t71 - t72 * t109) * t96;
t58 = -t87 * t135 + t134 + (-t72 * t127 + t71 * t87) * t69;
t57 = -t114 * t110 + (qJD(1) * t106 + t103 * t145) * t74;
t1 = [t96 * t83 * t111 + (qJD(2) * t106 - t126 * t130) * t74, t57, 0, 0, 0, 0; (t64 * t112 + (-t64 * t124 + (qJD(1) * t59 + t56) * t138) * t61) * t95 + (t65 * t112 * t59 + (-((t60 * t109 + t141 * t124 + t111) * t71 + (t110 * t133 - t139 + (t139 + (-t80 * t84 + t142) * t123) * t74) * t72) * t118 + (t86 * t120 - t65 * t124) * t59 + (-t64 + ((-t89 + t94) * t72 * t116 + t141 * t117) * t65) * t86 * qJD(1)) * t61) * t96 (t58 * t138 - t64 * t87) * t96 * t121 + ((-t64 * t126 + (-qJD(2) * t58 - t56) * t137) * t87 + (-t64 * t122 - (-t57 * t72 * t95 - t144 * t71 + (-qJD(2) * t71 - t125 * t72 + t135 * t60) * t69) * t118 + (t96 * t120 + t65 * t126) * t58 - ((t57 - t125) * t71 + ((-t69 * t95 + 0.1e1) * qJD(2) + (t69 - t95) * t60) * t72) * t87 * t137) * t86) * t61, 0, 0, 0, 0; t113 * t87 * t119 + (qJD(2) * t105 + 0.2e1 * t87 * t143) * t77, t86 * t119 * t129 + (-t87 * t90 * t122 + qJD(1) * t105) * t77, 0, 0, 0, 0;];
JaD_rot  = t1;
