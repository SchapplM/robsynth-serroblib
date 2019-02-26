% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPPR2_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:39:45
% EndTime: 2019-02-26 20:39:45
% DurationCPUTime: 0.52s
% Computational Cost: add. (2561->72), mult. (1839->154), div. (436->14), fcn. (2165->7), ass. (0->77)
t102 = qJ(1) + pkin(9);
t100 = cos(t102);
t129 = qJD(1) * t100;
t98 = sin(t102);
t153 = 0.2e1 * t98;
t87 = t98 ^ 2;
t88 = 0.1e1 / t98;
t96 = t100 ^ 2;
t151 = (0.1e1 + 0.1e1 / t87 * t96) * t88 * t129;
t101 = qJ(3) + pkin(10);
t97 = sin(t101);
t136 = t98 * t97;
t99 = cos(t101);
t78 = atan2(-t136, -t99);
t76 = sin(t78);
t124 = t76 * t136;
t77 = cos(t78);
t72 = -t77 * t99 - t124;
t69 = 0.1e1 / t72;
t92 = 0.1e1 / t99;
t150 = -0.2e1 * t97;
t70 = 0.1e1 / t72 ^ 2;
t86 = t97 ^ 2;
t93 = 0.1e1 / t99 ^ 2;
t139 = t86 * t93;
t83 = t87 * t139 + 0.1e1;
t79 = 0.1e1 / t83;
t149 = t79 - 0.1e1;
t119 = t97 * t129;
t131 = qJD(3) * t98;
t141 = t77 * t97;
t130 = qJD(3) * t99;
t65 = (-(-t98 * t130 - t119) * t92 + t131 * t139) * t79;
t61 = (-t65 * t98 + qJD(3)) * t141 + (-t119 + (t65 - t131) * t99) * t76;
t148 = t61 * t69 * t70;
t147 = t65 * t76;
t146 = t65 * t97;
t145 = t70 * t97;
t85 = t97 * t86;
t91 = t99 ^ 2;
t109 = qJD(3) * (t85 / t91 + t97) * t92;
t113 = t86 * t98 * t129;
t144 = (t87 * t109 + t93 * t113) / t83 ^ 2;
t121 = 0.1e1 + t139;
t112 = t121 * t79;
t75 = t98 * t112;
t143 = t75 * t98;
t142 = t76 * t99;
t140 = t86 * t92;
t138 = t86 * t96;
t89 = 0.1e1 / t98 ^ 2;
t137 = t89 * t96;
t135 = t100 * t70;
t134 = t100 * t97;
t133 = qJD(1) * t97;
t132 = qJD(1) * t98;
t128 = qJD(3) * t100;
t114 = t96 * t97 * t130;
t68 = t70 * t138 + 0.1e1;
t127 = 0.2e1 * (-t138 * t148 + (-t113 + t114) * t70) / t68 ^ 2;
t126 = 0.2e1 * t148;
t84 = t91 * t137 + 0.1e1;
t125 = 0.2e1 * (-t89 * t114 - t91 * t151) / t84 ^ 2;
t123 = t79 * t92 * t98;
t122 = t70 * t134;
t120 = 0.1e1 + t137;
t118 = t97 * t127;
t117 = t144 * t150;
t116 = t144 * t153;
t115 = t86 * t123;
t111 = t120 * t97;
t81 = 0.1e1 / t84;
t66 = 0.1e1 / t68;
t64 = (t149 * t97 * t76 - t77 * t115) * t100;
t63 = -t98 * t142 + t141 + (-t77 * t136 + t142) * t75;
t62 = -t121 * t116 + (t109 * t153 + t121 * t129) * t79;
t1 = [-t123 * t133 + (qJD(3) * t112 + t92 * t117) * t100, 0, t62, 0, 0, 0; (t69 * t118 + (-t69 * t130 + (qJD(1) * t64 + t61) * t145) * t66) * t98 + (t70 * t118 * t64 + (-((t65 * t115 + t149 * t130 + t117) * t76 + (t116 * t140 - t146 + (t146 + (-t85 * t93 + t150) * t131) * t79) * t77) * t122 + (t97 * t126 - t70 * t130) * t64 + (-t69 + ((-t87 + t96) * t79 * t77 * t140 + t149 * t124) * t70) * t133) * t66) * t100, 0 (t63 * t145 - t69 * t99) * t100 * t127 + ((-t69 * t132 + (-qJD(3) * t63 - t61) * t135) * t99 + (-t69 * t128 - (-t62 * t77 * t98 + t76 * t131 + t143 * t147 - t147 + (-qJD(3) * t76 - t129 * t77) * t75) * t122 + (t100 * t126 + t70 * t132) * t63 - ((t62 - t129) * t76 + ((0.1e1 - t143) * qJD(3) + (t75 - t98) * t65) * t77) * t99 * t135) * t97) * t66, 0, 0, 0; t120 * t99 * t125 + (qJD(3) * t111 + 0.2e1 * t99 * t151) * t81, 0, t88 * t125 * t134 + (-t88 * t99 * t128 + qJD(1) * t111) * t81, 0, 0, 0;];
JaD_rot  = t1;
