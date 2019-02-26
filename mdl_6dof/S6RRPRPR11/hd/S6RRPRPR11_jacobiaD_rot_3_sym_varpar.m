% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPRPR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR11_jacobiaD_rot_3_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_jacobiaD_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR11_jacobiaD_rot_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_jacobiaD_rot_3_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:43:46
% EndTime: 2019-02-26 21:43:46
% DurationCPUTime: 0.48s
% Computational Cost: add. (774->70), mult. (1839->155), div. (436->14), fcn. (2165->7), ass. (0->75)
t87 = cos(qJ(1));
t116 = qJD(1) * t87;
t85 = sin(qJ(1));
t138 = 0.2e1 * t85;
t74 = t85 ^ 2;
t75 = 0.1e1 / t85;
t83 = t87 ^ 2;
t136 = (0.1e1 + 0.1e1 / t74 * t83) * t75 * t116;
t84 = sin(qJ(2));
t118 = t85 * t84;
t86 = cos(qJ(2));
t65 = atan2(-t118, -t86);
t63 = sin(t65);
t108 = t63 * t118;
t64 = cos(t65);
t59 = -t64 * t86 - t108;
t56 = 0.1e1 / t59;
t79 = 0.1e1 / t86;
t57 = 0.1e1 / t59 ^ 2;
t135 = -0.2e1 * t84;
t73 = t84 ^ 2;
t80 = 0.1e1 / t86 ^ 2;
t123 = t73 * t80;
t70 = t74 * t123 + 0.1e1;
t66 = 0.1e1 / t70;
t134 = t66 - 0.1e1;
t106 = t84 * t116;
t115 = qJD(2) * t85;
t125 = t64 * t84;
t114 = qJD(2) * t86;
t52 = (-(-t85 * t114 - t106) * t79 + t115 * t123) * t66;
t48 = (-t52 * t85 + qJD(2)) * t125 + (-t106 + (t52 - t115) * t86) * t63;
t133 = t48 * t56 * t57;
t132 = t52 * t63;
t131 = t52 * t84;
t130 = t57 * t84;
t129 = t57 * t87;
t120 = t79 * t84;
t72 = t84 * t73;
t78 = t86 ^ 2;
t94 = qJD(2) * (t72 * t79 / t78 + t120);
t99 = t73 * t85 * t116;
t128 = (t74 * t94 + t80 * t99) / t70 ^ 2;
t105 = 0.1e1 + t123;
t62 = t105 * t85 * t66;
t127 = t62 * t85;
t126 = t63 * t86;
t124 = t73 * t79;
t122 = t73 * t83;
t76 = 0.1e1 / t85 ^ 2;
t121 = t76 * t83;
t119 = t84 * t87;
t117 = qJD(1) * t85;
t113 = qJD(2) * t87;
t55 = t57 * t122 + 0.1e1;
t98 = t83 * t84 * t114;
t112 = 0.2e1 * (-t122 * t133 + (t98 - t99) * t57) / t55 ^ 2;
t111 = 0.2e1 * t133;
t71 = t78 * t121 + 0.1e1;
t110 = 0.2e1 * (-t78 * t136 - t76 * t98) / t71 ^ 2;
t109 = t57 * t119;
t107 = t66 * t124;
t104 = 0.1e1 + t121;
t103 = t84 * t112;
t102 = t128 * t135;
t101 = t128 * t138;
t100 = t85 * t107;
t97 = t105 * t87;
t96 = t104 * t84;
t68 = 0.1e1 / t71;
t53 = 0.1e1 / t55;
t51 = (t134 * t84 * t63 - t64 * t100) * t87;
t50 = -t85 * t126 + t125 + (-t64 * t118 + t126) * t62;
t49 = -t105 * t101 + (qJD(1) * t97 + t94 * t138) * t66;
t1 = [t87 * t79 * t102 + (qJD(2) * t97 - t117 * t120) * t66, t49, 0, 0, 0, 0; (t56 * t103 + (-t56 * t114 + (qJD(1) * t51 + t48) * t130) * t53) * t85 + (t57 * t103 * t51 + (-((t52 * t100 + t134 * t114 + t102) * t63 + (t101 * t124 - t131 + (t131 + (-t72 * t80 + t135) * t115) * t66) * t64) * t109 + (t84 * t111 - t57 * t114) * t51 + (-t56 + ((-t74 + t83) * t64 * t107 + t134 * t108) * t57) * t84 * qJD(1)) * t53) * t87 (t50 * t130 - t56 * t86) * t87 * t112 + ((-t56 * t117 + (-qJD(2) * t50 - t48) * t129) * t86 + (-t56 * t113 - (-t49 * t64 * t85 + t63 * t115 + t127 * t132 - t132 + (-qJD(2) * t63 - t116 * t64) * t62) * t109 + (t87 * t111 + t57 * t117) * t50 - ((t49 - t116) * t63 + ((0.1e1 - t127) * qJD(2) + (t62 - t85) * t52) * t64) * t86 * t129) * t84) * t53, 0, 0, 0, 0; t104 * t86 * t110 + (qJD(2) * t96 + 0.2e1 * t86 * t136) * t68, t75 * t110 * t119 + (-t75 * t86 * t113 + qJD(1) * t96) * t68, 0, 0, 0, 0;];
JaD_rot  = t1;
