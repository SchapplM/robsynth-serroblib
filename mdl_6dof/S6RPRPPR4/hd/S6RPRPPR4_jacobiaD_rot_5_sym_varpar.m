% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPPR4
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
% Datum: 2019-02-26 20:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPPR4_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:40:54
% EndTime: 2019-02-26 20:40:54
% DurationCPUTime: 0.53s
% Computational Cost: add. (2051->78), mult. (3228->189), div. (613->15), fcn. (4191->9), ass. (0->87)
t130 = cos(pkin(10));
t181 = 0.2e1 * t130;
t127 = pkin(9) + qJ(3);
t123 = sin(t127);
t120 = 0.1e1 / t123;
t124 = cos(t127);
t132 = cos(qJ(1));
t159 = t132 * t130;
t129 = sin(pkin(10));
t131 = sin(qJ(1));
t162 = t131 * t129;
t111 = t124 * t162 + t159;
t125 = 0.1e1 / t129;
t167 = t111 * t125;
t149 = t120 * t167;
t163 = t123 * t129;
t101 = atan2(-t111, t163);
t97 = sin(t101);
t98 = cos(t101);
t106 = t111 ^ 2;
t121 = 0.1e1 / t123 ^ 2;
t126 = 0.1e1 / t129 ^ 2;
t104 = t106 * t121 * t126 + 0.1e1;
t99 = 0.1e1 / t104;
t180 = (t98 * t149 + t97) * t99 - t97;
t93 = -t97 * t111 + t98 * t163;
t90 = 0.1e1 / t93;
t160 = t132 * t129;
t147 = t124 * t160;
t161 = t131 * t130;
t114 = t147 - t161;
t82 = t180 * t114;
t179 = 0.2e1 * t82;
t115 = t124 * t159 + t162;
t108 = 0.1e1 / t115;
t109 = 0.1e1 / t115 ^ 2;
t91 = 0.1e1 / t93 ^ 2;
t107 = t114 ^ 2;
t154 = qJD(3) * t132;
t144 = t123 * t154;
t94 = t111 * qJD(1) + t129 * t144;
t176 = t91 * t94;
t156 = qJD(3) * t124;
t170 = t123 * t97;
t173 = t111 * t98;
t146 = t121 * t156;
t155 = qJD(3) * t131;
t158 = qJD(1) * t131;
t96 = qJD(1) * t147 - t130 * t158 - t155 * t163;
t83 = (t111 * t146 - t120 * t96) * t99 * t125;
t80 = -t83 * t173 - t97 * t96 + (t98 * t156 - t83 * t170) * t129;
t177 = t80 * t90 * t91;
t87 = t107 * t91 + 0.1e1;
t178 = (-t107 * t177 - t114 * t176) / t87 ^ 2;
t119 = t123 ^ 2;
t122 = t120 / t119;
t175 = 0.1e1 / t104 ^ 2 * (-t106 * t122 * t156 + t111 * t121 * t96) * t126;
t113 = -t124 * t161 + t160;
t140 = t130 * t144;
t95 = t113 * qJD(1) - t140;
t174 = t108 * t109 * t95;
t172 = t114 * t98;
t171 = t120 * t99;
t169 = t97 * t114;
t164 = t121 * t124;
t89 = (t164 * t167 + t131) * t99;
t168 = t131 - t89;
t166 = t113 * t132;
t128 = t132 ^ 2;
t165 = t119 * t128;
t157 = qJD(1) * t132;
t153 = -0.2e1 * t175;
t148 = t109 * t165;
t105 = 0.1e1 + t148;
t141 = t119 * t131 * t157;
t142 = t165 * t174;
t145 = qJD(3) * t123 * t128;
t152 = 0.2e1 / t105 ^ 2 * (-t142 + (t124 * t145 - t141) * t109);
t151 = t91 * t178;
t150 = t91 * t169;
t143 = 0.2e1 * t90 * t178;
t138 = 0.2e1 * t120 * t175 + t99 * t146;
t102 = 0.1e1 / t105;
t85 = 0.1e1 / t87;
t81 = -t89 * t173 + (t124 * t98 + t168 * t170) * t129;
t79 = t131 * t153 + t99 * t157 + (t96 * t99 * t164 + (t153 * t164 + (-0.2e1 * t122 * t124 ^ 2 - t120) * t99 * qJD(3)) * t111) * t125;
t1 = [(t138 * t114 + t94 * t171) * t125, 0, t79, 0, 0, 0; t111 * t143 + (-t96 * t90 + (t111 * t80 + t82 * t94) * t91) * t85 + (t151 * t179 + (t177 * t179 - (-t83 * t99 * t149 + t153) * t150 - ((t99 - 0.1e1) * t83 + (-t138 * t111 + t96 * t171) * t125) * t91 * t172 + t180 * t176) * t85) * t114, 0, t81 * t85 * t176 + (-(-t89 * t98 * t96 + (t83 * t89 * t97 - t79 * t98) * t111) * t91 * t85 + 0.2e1 * (t85 * t177 + t151) * t81) * t114 + (t132 * t143 * t123 + ((-t90 * t154 - (t168 * qJD(3) - t83) * t150) * t124 + (t90 * t158 + (t132 * t80 - (-t79 + t157) * t169 - (t168 * t83 - qJD(3)) * t172) * t91) * t123) * t85) * t129, 0, 0, 0; (-t108 * t131 - t109 * t166) * t123 * t152 + (-0.2e1 * t123 * t166 * t174 + (t123 * t157 + t124 * t155) * t108 + (((-t95 + t140) * t131 - t115 * t157) * t123 + (-t123 * t158 + t124 * t154) * t113) * t109) * t102, 0 (t108 * t124 * t132 + t130 * t148) * t152 + (t142 * t181 + (t124 * t158 + t144) * t108 + (t141 * t181 + (-0.2e1 * t130 * t145 + t132 * t95) * t124) * t109) * t102, 0, 0, 0;];
JaD_rot  = t1;
