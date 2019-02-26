% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRPR2_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:26:16
% EndTime: 2019-02-26 20:26:17
% DurationCPUTime: 0.54s
% Computational Cost: add. (2561->72), mult. (1839->155), div. (436->14), fcn. (2165->7), ass. (0->77)
t99 = qJ(1) + pkin(9);
t97 = cos(t99);
t128 = qJD(1) * t97;
t95 = sin(t99);
t150 = 0.2e1 * t95;
t84 = t95 ^ 2;
t85 = 0.1e1 / t95;
t93 = t97 ^ 2;
t148 = (0.1e1 + 0.1e1 / t84 * t93) * t85 * t128;
t98 = pkin(10) + qJ(4);
t94 = sin(t98);
t130 = t95 * t94;
t96 = cos(t98);
t75 = atan2(-t130, -t96);
t73 = sin(t75);
t120 = t73 * t130;
t74 = cos(t75);
t69 = -t74 * t96 - t120;
t66 = 0.1e1 / t69;
t89 = 0.1e1 / t96;
t147 = -0.2e1 * t94;
t67 = 0.1e1 / t69 ^ 2;
t83 = t94 ^ 2;
t90 = 0.1e1 / t96 ^ 2;
t135 = t83 * t90;
t80 = t84 * t135 + 0.1e1;
t76 = 0.1e1 / t80;
t146 = t76 - 0.1e1;
t118 = t94 * t128;
t127 = qJD(4) * t95;
t137 = t74 * t94;
t126 = qJD(4) * t96;
t62 = (-(-t95 * t126 - t118) * t89 + t127 * t135) * t76;
t58 = (-t62 * t95 + qJD(4)) * t137 + (-t118 + (t62 - t127) * t96) * t73;
t145 = t58 * t66 * t67;
t144 = t62 * t73;
t143 = t62 * t94;
t142 = t67 * t94;
t141 = t67 * t97;
t132 = t89 * t94;
t82 = t94 * t83;
t88 = t96 ^ 2;
t106 = qJD(4) * (t82 * t89 / t88 + t132);
t111 = t83 * t95 * t128;
t140 = (t84 * t106 + t90 * t111) / t80 ^ 2;
t117 = 0.1e1 + t135;
t72 = t117 * t95 * t76;
t139 = t72 * t95;
t138 = t73 * t96;
t136 = t83 * t89;
t134 = t83 * t93;
t86 = 0.1e1 / t95 ^ 2;
t133 = t86 * t93;
t131 = t94 * t97;
t129 = qJD(1) * t95;
t125 = qJD(4) * t97;
t110 = t93 * t94 * t126;
t65 = t67 * t134 + 0.1e1;
t124 = 0.2e1 * (-t134 * t145 + (t110 - t111) * t67) / t65 ^ 2;
t123 = 0.2e1 * t145;
t81 = t88 * t133 + 0.1e1;
t122 = 0.2e1 * (-t86 * t110 - t88 * t148) / t81 ^ 2;
t121 = t67 * t131;
t119 = t76 * t136;
t116 = 0.1e1 + t133;
t115 = t94 * t124;
t114 = t140 * t147;
t113 = t140 * t150;
t112 = t95 * t119;
t109 = t117 * t97;
t108 = t116 * t94;
t78 = 0.1e1 / t81;
t63 = 0.1e1 / t65;
t61 = (t146 * t94 * t73 - t74 * t112) * t97;
t60 = -t95 * t138 + t137 + (-t74 * t130 + t138) * t72;
t59 = -t117 * t113 + (qJD(1) * t109 + t106 * t150) * t76;
t1 = [t97 * t89 * t114 + (qJD(4) * t109 - t129 * t132) * t76, 0, 0, t59, 0, 0; (t66 * t115 + (-t66 * t126 + (qJD(1) * t61 + t58) * t142) * t63) * t95 + (t67 * t115 * t61 + (-((t62 * t112 + t146 * t126 + t114) * t73 + (t113 * t136 - t143 + (t143 + (-t82 * t90 + t147) * t127) * t76) * t74) * t121 + (t94 * t123 - t67 * t126) * t61 + (-t66 + ((-t84 + t93) * t74 * t119 + t146 * t120) * t67) * t94 * qJD(1)) * t63) * t97, 0, 0 (t60 * t142 - t66 * t96) * t97 * t124 + ((-t66 * t129 + (-qJD(4) * t60 - t58) * t141) * t96 + (-t66 * t125 - (-t59 * t74 * t95 + t73 * t127 + t139 * t144 - t144 + (-qJD(4) * t73 - t128 * t74) * t72) * t121 + (t97 * t123 + t67 * t129) * t60 - ((t59 - t128) * t73 + ((0.1e1 - t139) * qJD(4) + (t72 - t95) * t62) * t74) * t96 * t141) * t94) * t63, 0, 0; t116 * t96 * t122 + (qJD(4) * t108 + 0.2e1 * t96 * t148) * t78, 0, 0, t85 * t122 * t131 + (-t85 * t96 * t125 + qJD(1) * t108) * t78, 0, 0;];
JaD_rot  = t1;
