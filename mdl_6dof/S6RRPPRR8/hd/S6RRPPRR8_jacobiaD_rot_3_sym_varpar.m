% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPPRR8
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
% Datum: 2019-02-26 21:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR8_jacobiaD_rot_3_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_jacobiaD_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR8_jacobiaD_rot_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_jacobiaD_rot_3_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:32:39
% EndTime: 2019-02-26 21:32:40
% DurationCPUTime: 0.62s
% Computational Cost: add. (892->82), mult. (2191->191), div. (456->12), fcn. (2616->9), ass. (0->85)
t100 = sin(qJ(2));
t92 = t100 ^ 2;
t102 = cos(qJ(2));
t95 = 0.1e1 / t102 ^ 2;
t141 = t92 * t95;
t101 = sin(qJ(1));
t122 = 0.1e1 + t141;
t93 = t101 ^ 2;
t90 = t93 * t141 + 0.1e1;
t88 = 0.1e1 / t90;
t113 = t122 * t88;
t74 = t101 * t113;
t157 = t101 * t74 - 0.1e1;
t103 = cos(qJ(1));
t127 = qJD(2) * t103;
t118 = t100 * t127;
t134 = t101 * t102;
t98 = sin(pkin(10));
t99 = cos(pkin(10));
t82 = t103 * t98 - t99 * t134;
t76 = t82 * qJD(1) - t99 * t118;
t133 = t102 * t103;
t84 = t101 * t98 + t99 * t133;
t78 = 0.1e1 / t84;
t79 = 0.1e1 / t84 ^ 2;
t80 = t78 * t79;
t145 = t76 * t80;
t81 = -t103 * t99 - t98 * t134;
t75 = t81 * qJD(1) - t98 * t118;
t146 = t75 * t79;
t83 = -t101 * t99 + t98 * t133;
t77 = t83 ^ 2;
t72 = t77 * t79 + 0.1e1;
t156 = (-t77 * t145 + t83 * t146) / t72 ^ 2;
t155 = t100 * t141;
t143 = t83 * t99;
t112 = t79 * t143 - t78 * t98;
t70 = 0.1e1 / t72;
t154 = t112 * t70;
t135 = t101 * t100;
t87 = atan2(-t135, -t102);
t85 = sin(t87);
t123 = t85 * t135;
t86 = cos(t87);
t69 = -t102 * t86 - t123;
t66 = 0.1e1 / t69;
t94 = 0.1e1 / t102;
t67 = 0.1e1 / t69 ^ 2;
t153 = 0.2e1 * t100;
t152 = t88 - 0.1e1;
t130 = qJD(1) * t103;
t114 = t101 * t92 * t130;
t128 = qJD(2) * t102;
t97 = t103 ^ 2;
t140 = t92 * t97;
t129 = qJD(2) * t101;
t138 = t102 * t85;
t62 = (-(-t100 * t130 - t101 * t128) * t94 + t129 * t141) * t88;
t57 = (t62 - t129) * t138 + (-t85 * t130 + (-t101 * t62 + qJD(2)) * t86) * t100;
t150 = t57 * t66 * t67;
t65 = t67 * t140 + 0.1e1;
t151 = (-t140 * t150 + (t100 * t97 * t128 - t114) * t67) / t65 ^ 2;
t63 = 0.1e1 / t65;
t148 = t63 * t67;
t111 = qJD(2) * (t100 + t155) * t94;
t147 = (t93 * t111 + t95 * t114) / t90 ^ 2;
t144 = t82 * t83;
t142 = t88 * t94;
t137 = t103 * t67;
t136 = qJD(2) * t74;
t132 = qJD(1) * t100;
t131 = qJD(1) * t101;
t126 = 0.2e1 * t150;
t125 = t66 * t151;
t124 = t101 * t142;
t121 = t100 * t152;
t120 = t63 * t128;
t119 = t100 * t129;
t117 = 0.2e1 * t67 * t151;
t116 = -0.2e1 * t94 * t147;
t115 = t92 * t124;
t61 = (-t86 * t115 + t85 * t121) * t103;
t59 = (-t101 + t74) * t138 - t157 * t86 * t100;
t58 = t113 * t130 + 0.2e1 * (t111 * t88 - t122 * t147) * t101;
t1 = [-t124 * t132 + (qJD(2) * t113 + t100 * t116) * t103, t58, 0, 0, 0, 0; (-t66 * t120 + (0.2e1 * t125 + (qJD(1) * t61 + t57) * t148) * t100) * t101 + (-t61 * t67 * t120 + (t61 * t117 + (t61 * t126 + ((-t62 * t115 - t152 * t128 + t147 * t153) * t85 + (-t62 * t121 + (t92 * t116 + (t153 + t155) * t88 * qJD(2)) * t101) * t86) * t137) * t63) * t100 + (-t66 + (-(t93 - t97) * t92 * t86 * t142 + t152 * t123) * t67) * t63 * t132) * t103 (-t66 * t63 * t131 + (-0.2e1 * t125 + (-qJD(2) * t59 - t57) * t148) * t103) * t102 + (t59 * t103 * t117 + (-t66 * t127 - ((-t101 * t58 - t130 * t74) * t86 + (t157 * t62 + t129 - t136) * t85) * t100 * t137 + (t103 * t126 + t67 * t131) * t59) * t63 - ((t58 - t130) * t85 + (t62 * t74 + qJD(2) + (-t62 - t136) * t101) * t86) * t133 * t148) * t100, 0, 0, 0, 0; 0.2e1 * (t79 * t144 - t78 * t81) * t156 + ((-t83 * qJD(1) + t98 * t119) * t78 + 0.2e1 * t144 * t145 + (-t81 * t76 - (-t84 * qJD(1) + t99 * t119) * t83 - t82 * t75) * t79) * t70, t102 * t127 * t154 + (-t131 * t154 + (-0.2e1 * t112 * t156 + (t99 * t146 + (-0.2e1 * t80 * t143 + t79 * t98) * t76) * t70) * t103) * t100, 0, 0, 0, 0;];
JaD_rot  = t1;
