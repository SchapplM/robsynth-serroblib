% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR6_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR6_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_jacobiaD_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:18:48
% EndTime: 2019-02-26 22:18:48
% DurationCPUTime: 0.74s
% Computational Cost: add. (1350->93), mult. (2519->209), div. (480->12), fcn. (2968->9), ass. (0->94)
t135 = sin(qJ(1));
t129 = t135 ^ 2;
t134 = sin(qJ(2));
t128 = t134 ^ 2;
t136 = cos(qJ(2));
t131 = 0.1e1 / t136 ^ 2;
t184 = t128 * t131;
t121 = t129 * t184 + 0.1e1;
t127 = t134 * t128;
t130 = 0.1e1 / t136;
t181 = t130 * t134;
t146 = qJD(2) * (t127 * t130 * t131 + t181);
t137 = cos(qJ(1));
t173 = qJD(1) * t137;
t182 = t128 * t135;
t150 = t173 * t182;
t190 = (t129 * t146 + t131 * t150) / t121 ^ 2;
t200 = -0.2e1 * t190;
t126 = qJ(3) + pkin(11);
t125 = cos(t126);
t175 = t136 * t137;
t124 = sin(t126);
t179 = t135 * t124;
t114 = t125 * t175 + t179;
t109 = 0.1e1 / t114 ^ 2;
t178 = t135 * t125;
t113 = t124 * t175 - t178;
t187 = t113 * t125;
t108 = 0.1e1 / t114;
t189 = t108 * t124;
t148 = t109 * t187 - t189;
t107 = t113 ^ 2;
t100 = t107 * t109 + 0.1e1;
t98 = 0.1e1 / t100;
t199 = t148 * t98;
t157 = 0.1e1 + t184;
t198 = t135 * t157;
t177 = t135 * t134;
t118 = atan2(-t177, -t136);
t116 = cos(t118);
t115 = sin(t118);
t163 = t115 * t177;
t104 = -t116 * t136 - t163;
t101 = 0.1e1 / t104;
t102 = 0.1e1 / t104 ^ 2;
t119 = 0.1e1 / t121;
t197 = t119 - 0.1e1;
t133 = t137 ^ 2;
t171 = qJD(2) * t136;
t183 = t128 * t133;
t172 = qJD(2) * t135;
t159 = t131 * t172;
t160 = t134 * t173;
t94 = (-(-t135 * t171 - t160) * t130 + t128 * t159) * t119;
t154 = t94 - t172;
t155 = -t135 * t94 + qJD(2);
t186 = t116 * t134;
t88 = t155 * t186 + (t154 * t136 - t160) * t115;
t193 = t101 * t102 * t88;
t97 = t102 * t183 + 0.1e1;
t196 = (-t183 * t193 + (t133 * t134 * t171 - t150) * t102) / t97 ^ 2;
t188 = t109 * t113;
t152 = -qJD(1) * t136 + qJD(3);
t153 = qJD(3) * t136 - qJD(1);
t170 = qJD(2) * t137;
t158 = t134 * t170;
t185 = t124 * t137;
t93 = -t153 * t185 + (t152 * t135 - t158) * t125;
t192 = t108 * t109 * t93;
t176 = t135 * t136;
t147 = t124 * t176 + t125 * t137;
t92 = t147 * qJD(1) - t114 * qJD(3) + t124 * t158;
t195 = (-t107 * t192 - t92 * t188) / t100 ^ 2;
t95 = 0.1e1 / t97;
t194 = t102 * t95;
t180 = t134 * t137;
t174 = qJD(1) * t135;
t169 = 0.2e1 * t196;
t168 = -0.2e1 * t195;
t167 = 0.2e1 * t193;
t166 = t101 * t196;
t165 = t113 * t192;
t164 = t95 * t171;
t162 = t119 * t128 * t130;
t156 = t130 * t200;
t151 = t135 * t162;
t149 = t157 * t137;
t145 = t134 * t172 + t152 * t137;
t112 = -t125 * t176 + t185;
t106 = t119 * t198;
t91 = (t197 * t134 * t115 - t116 * t151) * t137;
t90 = -t115 * t176 + t186 + (t115 * t136 - t116 * t177) * t106;
t89 = t198 * t200 + (qJD(1) * t149 + 0.2e1 * t135 * t146) * t119;
t1 = [t156 * t180 + (qJD(2) * t149 - t174 * t181) * t119, t89, 0, 0, 0, 0; (-t101 * t164 + (0.2e1 * t166 + (qJD(1) * t91 + t88) * t194) * t134) * t135 + ((-t91 * t164 + (t91 * t169 + ((0.2e1 * t134 * t190 - t94 * t151 - t197 * t171) * t115 + (t156 * t182 + t134 * t94 + (t127 * t159 - (t94 - 0.2e1 * t172) * t134) * t119) * t116) * t95 * t137) * t134) * t102 + (t91 * t167 + (-t101 + ((-t129 + t133) * t116 * t162 + t197 * t163) * t102) * qJD(1)) * t134 * t95) * t137 (-t101 * t95 * t174 + (-0.2e1 * t166 + (-qJD(2) * t90 - t88) * t194) * t137) * t136 + (t90 * t137 * t102 * t169 + ((-qJD(2) * t101 + t90 * t167) * t137 + (t90 * t174 + (-(-t106 * t173 - t135 * t89) * t116 - ((t106 * t135 - 0.1e1) * t94 + (-t106 + t135) * qJD(2)) * t115) * t180) * t102) * t95 - ((t89 - t173) * t115 + (t154 * t106 + t155) * t116) * t175 * t194) * t134, 0, 0, 0, 0; 0.2e1 * (t108 * t147 + t112 * t188) * t195 + (0.2e1 * t112 * t165 - t153 * t108 * t178 + t145 * t189 + (-t153 * t113 * t179 + t112 * t92 - t145 * t187 + t147 * t93) * t109) * t98, t136 * t170 * t199 + (-t174 * t199 + (t148 * t168 + ((-qJD(3) * t108 - 0.2e1 * t165) * t125 + (-t125 * t92 + (-qJD(3) * t113 + t93) * t124) * t109) * t98) * t137) * t134, t168 + 0.2e1 * (-t109 * t92 * t98 + (-t109 * t195 - t98 * t192) * t113) * t113, 0, 0, 0;];
JaD_rot  = t1;
