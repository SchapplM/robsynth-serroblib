% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPPR7_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_jacobiaD_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:42:33
% EndTime: 2019-02-26 20:42:33
% DurationCPUTime: 0.50s
% Computational Cost: add. (1694->69), mult. (1839->158), div. (436->14), fcn. (2165->7), ass. (0->74)
t96 = cos(qJ(1));
t144 = 0.2e1 * t96;
t92 = 0.1e1 / t96;
t95 = sin(qJ(1));
t126 = t92 * t95;
t90 = t95 ^ 2;
t91 = t96 ^ 2;
t104 = qJD(1) * (t90 / t91 + 0.1e1) * t126;
t88 = qJ(3) + pkin(9);
t86 = sin(t88);
t122 = qJD(3) * t86;
t87 = cos(t88);
t107 = t87 * t90 * t122;
t93 = 0.1e1 / t96 ^ 2;
t127 = t90 * t93;
t80 = t86 ^ 2;
t79 = t80 * t127 + 0.1e1;
t143 = -0.2e1 * (t80 * t104 + t93 * t107) / t79 ^ 2;
t120 = qJD(3) * t96;
t124 = qJD(1) * t95;
t114 = t87 * t124;
t82 = 0.1e1 / t86 ^ 2;
t85 = t87 ^ 2;
t129 = t82 * t85;
t76 = t91 * t129 + 0.1e1;
t74 = 0.1e1 / t76;
t81 = 0.1e1 / t86;
t60 = ((t86 * t120 + t114) * t81 + t120 * t129) * t74;
t142 = -t60 + t120;
t125 = t96 * t87;
t73 = atan2(-t125, t86);
t71 = sin(t73);
t72 = cos(t73);
t67 = -t71 * t125 + t72 * t86;
t64 = 0.1e1 / t67;
t65 = 0.1e1 / t67 ^ 2;
t141 = t74 - 0.1e1;
t123 = qJD(1) * t96;
t108 = t85 * t95 * t123;
t128 = t85 * t90;
t132 = t72 * t87;
t56 = (-t60 * t96 + qJD(3)) * t132 + (t142 * t86 + t114) * t71;
t139 = t56 * t64 * t65;
t63 = t65 * t128 + 0.1e1;
t140 = (-t128 * t139 + (-t107 + t108) * t65) / t63 ^ 2;
t138 = t60 * t87;
t137 = t65 * t87;
t136 = t65 * t95;
t130 = t81 * t87;
t84 = t87 * t85;
t103 = qJD(3) * (-t81 / t80 * t84 - t130);
t135 = (t91 * t103 - t82 * t108) / t76 ^ 2;
t133 = t71 * t96;
t131 = t81 * t85;
t121 = qJD(3) * t95;
t119 = -0.2e1 * t139;
t118 = t87 * t140;
t117 = t87 * t136;
t116 = t87 * t135;
t115 = t74 * t131;
t113 = 0.1e1 + t129;
t112 = 0.1e1 + t127;
t111 = t135 * t144;
t110 = t96 * t115;
t109 = t141 * t87 * t71;
t106 = t113 * t95;
t105 = t112 * t87;
t77 = 0.1e1 / t79;
t69 = t113 * t96 * t74;
t61 = 0.1e1 / t63;
t59 = (-t72 * t110 - t109) * t95;
t58 = t86 * t133 + t132 + (-t72 * t125 - t71 * t86) * t69;
t57 = -t113 * t111 + (-qJD(1) * t106 + t103 * t144) * t74;
t1 = [-0.2e1 * t95 * t81 * t116 + (-qJD(3) * t106 + t123 * t130) * t74, 0, t57, 0, 0, 0; (0.2e1 * t64 * t118 + (t64 * t122 + (qJD(1) * t59 + t56) * t137) * t61) * t96 + (-0.2e1 * t65 * t118 * t59 + (((t60 * t110 + t141 * t122 + 0.2e1 * t116) * t71 + (t111 * t131 + t138 + (-t138 + (t82 * t84 + 0.2e1 * t87) * t120) * t74) * t72) * t117 + (t87 * t119 - t65 * t122) * t59 + (t64 + ((t90 - t91) * t72 * t115 - t96 * t109) * t65) * t87 * qJD(1)) * t61) * t95, 0, 0.2e1 * (-t58 * t137 - t64 * t86) * t95 * t140 + ((t64 * t123 + (-qJD(3) * t58 - t56) * t136) * t86 + (t64 * t121 + (-t57 * t72 * t96 + t142 * t71 + (-qJD(3) * t71 + t124 * t72 + t133 * t60) * t69) * t117 + (t95 * t119 + t65 * t123) * t58 + ((-t57 - t124) * t71 + ((t69 * t96 - 0.1e1) * qJD(3) + (-t69 + t96) * t60) * t72) * t86 * t136) * t87) * t61, 0, 0, 0; t112 * t86 * t143 + (qJD(3) * t105 + 0.2e1 * t104 * t86) * t77, 0, t87 * t126 * t143 + (-t86 * t92 * t121 + qJD(1) * t105) * t77, 0, 0, 0;];
JaD_rot  = t1;
