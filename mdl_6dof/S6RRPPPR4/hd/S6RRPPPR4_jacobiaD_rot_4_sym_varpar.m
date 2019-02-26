% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPPR4_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_jacobiaD_rot_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:23:35
% EndTime: 2019-02-26 21:23:36
% DurationCPUTime: 0.68s
% Computational Cost: add. (703->81), mult. (2191->188), div. (456->12), fcn. (2616->9), ass. (0->87)
t101 = sin(qJ(2));
t93 = 0.1e1 / t101 ^ 2;
t103 = cos(qJ(2));
t97 = t103 ^ 2;
t151 = t93 * t97;
t102 = sin(qJ(1));
t126 = 0.1e1 + t151;
t95 = t102 ^ 2;
t91 = t95 * t151 + 0.1e1;
t89 = 0.1e1 / t91;
t116 = t126 * t89;
t75 = t102 * t116;
t161 = t102 * t75 - 0.1e1;
t159 = t103 * t151;
t92 = 0.1e1 / t101;
t113 = qJD(2) * (-t103 - t159) * t92;
t104 = cos(qJ(1));
t138 = qJD(1) * t104;
t117 = t102 * t97 * t138;
t160 = (t95 * t113 + t93 * t117) / t91 ^ 2;
t142 = t102 * t103;
t88 = atan2(-t142, t101);
t86 = sin(t88);
t127 = t86 * t142;
t87 = cos(t88);
t70 = t87 * t101 - t127;
t67 = 0.1e1 / t70;
t100 = cos(pkin(9));
t143 = t102 * t100;
t144 = t101 * t104;
t99 = sin(pkin(9));
t83 = t99 * t144 + t143;
t79 = 0.1e1 / t83;
t68 = 0.1e1 / t70 ^ 2;
t80 = 0.1e1 / t83 ^ 2;
t158 = t89 - 0.1e1;
t136 = qJD(2) * t102;
t149 = t101 * t86;
t63 = ((t101 * t136 - t103 * t138) * t92 + t136 * t151) * t89;
t58 = (-t63 + t136) * t149 + (-t86 * t138 + (-t102 * t63 + qJD(2)) * t87) * t103;
t157 = t58 * t67 * t68;
t98 = t104 ^ 2;
t150 = t97 * t98;
t66 = t68 * t150 + 0.1e1;
t64 = 0.1e1 / t66;
t155 = t64 * t68;
t141 = t104 * t100;
t147 = t102 * t99;
t114 = t101 * t141 - t147;
t154 = t80 * t114;
t153 = t80 * t99;
t152 = t92 * t97;
t146 = t104 * t68;
t145 = qJD(2) * t75;
t140 = qJD(1) * t102;
t139 = qJD(1) * t103;
t137 = qJD(2) * t101;
t135 = qJD(2) * t103;
t134 = qJD(2) * t104;
t133 = 0.2e1 * (-t150 * t157 + (-t101 * t98 * t135 - t117) * t68) / t66 ^ 2;
t132 = 0.2e1 * t157;
t78 = t114 ^ 2;
t73 = t78 * t80 + 0.1e1;
t122 = t103 * t134;
t84 = t101 * t143 + t104 * t99;
t76 = t84 * qJD(1) - t100 * t122;
t85 = -t101 * t147 + t141;
t77 = t85 * qJD(1) + t99 * t122;
t81 = t79 * t80;
t131 = 0.2e1 * (-t78 * t81 * t77 - t76 * t154) / t73 ^ 2;
t130 = 0.2e1 * t160;
t129 = -0.2e1 * t81 * t114;
t128 = t102 * t89 * t92;
t125 = t103 * t158;
t124 = t64 * t137;
t123 = t102 * t135;
t121 = t67 * t133;
t120 = t68 * t133;
t119 = t103 * t130;
t118 = t97 * t128;
t115 = t100 * t79 - t114 * t153;
t71 = 0.1e1 / t73;
t112 = t115 * t71;
t62 = (t87 * t118 + t86 * t125) * t104;
t60 = -t161 * t87 * t103 + (t102 - t75) * t149;
t59 = t116 * t138 + 0.2e1 * (t113 * t89 - t126 * t160) * t102;
t1 = [t128 * t139 + (qJD(2) * t116 + t92 * t119) * t104, t59, 0, 0, 0, 0; (t67 * t124 + (t121 + (qJD(1) * t62 + t58) * t155) * t103) * t102 + (t62 * t68 * t124 + (t62 * t120 + (t62 * t132 + ((t63 * t118 + t158 * t137 + t119) * t86 + (-t63 * t125 + (t130 * t152 + (0.2e1 * t103 + t159) * t89 * qJD(2)) * t102) * t87) * t146) * t64) * t103 + (-t67 + (-(-t95 + t98) * t89 * t87 * t152 + t158 * t127) * t68) * t64 * t139) * t104 (t67 * t64 * t140 + (t121 + (qJD(2) * t60 + t58) * t155) * t104) * t101 + (t60 * t104 * t120 + (-t67 * t134 - ((-t102 * t59 - t138 * t75) * t87 + (t161 * t63 + t136 - t145) * t86) * t103 * t146 + (t104 * t132 + t68 * t140) * t60) * t64 - ((-t59 + t138) * t86 + (-t63 * t75 - qJD(2) + (t63 + t145) * t102) * t87) * t144 * t155) * t103, 0, 0, 0, 0; (-t85 * t154 - t79 * t84) * t131 + ((t114 * qJD(1) + t100 * t123) * t79 + t85 * t77 * t129 + (-t84 * t77 + (-t83 * qJD(1) - t99 * t123) * t114 - t85 * t76) * t80) * t71, t101 * t112 * t134 + (t112 * t140 + (t115 * t131 + (-t76 * t153 + (t100 * t80 + t99 * t129) * t77) * t71) * t104) * t103, 0, 0, 0, 0;];
JaD_rot  = t1;
