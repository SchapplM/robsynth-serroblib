% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPPR3_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_jacobiaD_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:23:03
% EndTime: 2019-02-26 21:23:04
% DurationCPUTime: 0.62s
% Computational Cost: add. (703->82), mult. (2191->188), div. (456->12), fcn. (2616->9), ass. (0->87)
t101 = sin(qJ(2));
t93 = 0.1e1 / t101 ^ 2;
t103 = cos(qJ(2));
t97 = t103 ^ 2;
t149 = t93 * t97;
t102 = sin(qJ(1));
t125 = 0.1e1 + t149;
t95 = t102 ^ 2;
t91 = t95 * t149 + 0.1e1;
t89 = 0.1e1 / t91;
t114 = t125 * t89;
t75 = t102 * t114;
t161 = t102 * t75 - 0.1e1;
t158 = t103 * t149;
t92 = 0.1e1 / t101;
t113 = qJD(2) * (-t103 - t158) * t92;
t104 = cos(qJ(1));
t136 = qJD(1) * t104;
t115 = t102 * t97 * t136;
t160 = (t95 * t113 + t93 * t115) / t91 ^ 2;
t100 = cos(pkin(9));
t132 = qJD(2) * t104;
t121 = t103 * t132;
t141 = t102 * t100;
t99 = sin(pkin(9));
t83 = -t101 * t141 - t104 * t99;
t77 = t83 * qJD(1) + t100 * t121;
t139 = t104 * t100;
t145 = t102 * t99;
t85 = t101 * t139 - t145;
t80 = 0.1e1 / t85 ^ 2;
t159 = t77 * t80;
t140 = t102 * t103;
t88 = atan2(-t140, t101);
t86 = sin(t88);
t126 = t86 * t140;
t87 = cos(t88);
t70 = t87 * t101 - t126;
t67 = 0.1e1 / t70;
t79 = 0.1e1 / t85;
t68 = 0.1e1 / t70 ^ 2;
t157 = t89 - 0.1e1;
t134 = qJD(2) * t102;
t147 = t101 * t86;
t63 = ((t101 * t134 - t103 * t136) * t92 + t134 * t149) * t89;
t58 = (-t63 + t134) * t147 + (-t86 * t136 + (-t102 * t63 + qJD(2)) * t87) * t103;
t156 = t58 * t67 * t68;
t142 = t101 * t104;
t84 = t99 * t142 + t141;
t151 = t80 * t84;
t152 = t79 * t159;
t78 = t84 ^ 2;
t73 = t78 * t80 + 0.1e1;
t82 = -t101 * t145 + t139;
t76 = t82 * qJD(1) + t99 * t121;
t155 = (t76 * t151 - t78 * t152) / t73 ^ 2;
t98 = t104 ^ 2;
t148 = t97 * t98;
t66 = t68 * t148 + 0.1e1;
t64 = 0.1e1 / t66;
t153 = t64 * t68;
t150 = t92 * t97;
t144 = t104 * t68;
t143 = qJD(2) * t75;
t138 = qJD(1) * t102;
t137 = qJD(1) * t103;
t135 = qJD(2) * t101;
t133 = qJD(2) * t103;
t131 = 0.2e1 * (-t148 * t156 + (-t101 * t98 * t133 - t115) * t68) / t66 ^ 2;
t130 = 0.2e1 * t156;
t129 = 0.2e1 * t155;
t128 = 0.2e1 * t160;
t127 = t102 * t89 * t92;
t124 = t157 * t103;
t123 = t64 * t135;
t122 = t102 * t133;
t120 = t67 * t131;
t119 = t68 * t131;
t118 = 0.2e1 * t84 * t152;
t117 = t103 * t128;
t116 = t97 * t127;
t71 = 0.1e1 / t73;
t112 = (t100 * t151 - t79 * t99) * t71;
t62 = (t87 * t116 + t86 * t124) * t104;
t60 = -t161 * t87 * t103 + (t102 - t75) * t147;
t59 = t114 * t136 + 0.2e1 * (t113 * t89 - t125 * t160) * t102;
t1 = [t127 * t137 + (qJD(2) * t114 + t92 * t117) * t104, t59, 0, 0, 0, 0; (t67 * t123 + (t120 + (qJD(1) * t62 + t58) * t153) * t103) * t102 + (t62 * t68 * t123 + (t62 * t119 + (t62 * t130 + ((t63 * t116 + t157 * t135 + t117) * t86 + (-t63 * t124 + (t128 * t150 + (0.2e1 * t103 + t158) * t89 * qJD(2)) * t102) * t87) * t144) * t64) * t103 + (-t67 + (-(-t95 + t98) * t89 * t87 * t150 + t157 * t126) * t68) * t64 * t137) * t104 (t67 * t64 * t138 + (t120 + (qJD(2) * t60 + t58) * t153) * t104) * t101 + (t60 * t104 * t119 + (-t67 * t132 - ((-t102 * t59 - t136 * t75) * t87 + (t161 * t63 + t134 - t143) * t86) * t103 * t144 + (t104 * t130 + t68 * t138) * t60) * t64 - ((-t59 + t136) * t86 + (-t63 * t75 - qJD(2) + (t63 + t143) * t102) * t87) * t142 * t153) * t103, 0, 0, 0, 0; (t83 * t151 - t79 * t82) * t129 + ((-t84 * qJD(1) - t99 * t122) * t79 + t83 * t118 + (-t82 * t77 - (-t85 * qJD(1) - t100 * t122) * t84 - t83 * t76) * t80) * t71, t101 * t112 * t132 + (t112 * t138 + ((-0.2e1 * t79 * t155 - t71 * t159) * t99 + (t129 * t151 + (-t76 * t80 + t118) * t71) * t100) * t104) * t103, 0, 0, 0, 0;];
JaD_rot  = t1;
