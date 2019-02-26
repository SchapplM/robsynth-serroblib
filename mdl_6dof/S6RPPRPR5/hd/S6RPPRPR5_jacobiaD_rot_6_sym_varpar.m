% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRPR5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:28:01
% EndTime: 2019-02-26 20:28:01
% DurationCPUTime: 0.70s
% Computational Cost: add. (972->92), mult. (2519->207), div. (480->12), fcn. (2968->9), ass. (0->95)
t135 = sin(qJ(1));
t134 = sin(qJ(4));
t128 = 0.1e1 / t134 ^ 2;
t136 = cos(qJ(4));
t132 = t136 ^ 2;
t185 = t128 * t132;
t159 = 0.1e1 + t185;
t198 = t135 * t159;
t179 = t135 * t136;
t121 = atan2(t179, t134);
t118 = cos(t121);
t117 = sin(t121);
t165 = t117 * t179;
t104 = t118 * t134 + t165;
t101 = 0.1e1 / t104;
t126 = pkin(9) + qJ(6);
t125 = cos(t126);
t137 = cos(qJ(1));
t182 = t134 * t137;
t163 = t125 * t182;
t124 = sin(t126);
t181 = t135 * t124;
t114 = t163 - t181;
t108 = 0.1e1 / t114;
t127 = 0.1e1 / t134;
t102 = 0.1e1 / t104 ^ 2;
t109 = 0.1e1 / t114 ^ 2;
t130 = t135 ^ 2;
t122 = t130 * t185 + 0.1e1;
t119 = 0.1e1 / t122;
t197 = t119 - 0.1e1;
t180 = t135 * t125;
t113 = t124 * t182 + t180;
t107 = t113 ^ 2;
t190 = t109 * t113;
t152 = qJD(1) * t134 + qJD(6);
t153 = qJD(6) * t134 + qJD(1);
t171 = qJD(4) * t137;
t160 = t136 * t171;
t186 = t124 * t137;
t93 = -t153 * t186 + (-t152 * t135 + t160) * t125;
t193 = t108 * t109 * t93;
t175 = qJD(1) * t137;
t92 = -qJD(6) * t163 - t124 * t160 - t125 * t175 + t152 * t181;
t97 = t107 * t109 + 0.1e1;
t196 = (-t107 * t193 - t92 * t190) / t97 ^ 2;
t133 = t137 ^ 2;
t184 = t132 * t133;
t100 = t102 * t184 + 0.1e1;
t98 = 0.1e1 / t100;
t195 = t102 * t98;
t173 = qJD(4) * t135;
t161 = t128 * t173;
t162 = t136 * t175;
t94 = ((-t134 * t173 + t162) * t127 - t132 * t161) * t119;
t154 = -t94 - t173;
t155 = t135 * t94 + qJD(4);
t187 = t118 * t136;
t88 = t155 * t187 + (t154 * t134 + t162) * t117;
t194 = t101 * t102 * t88;
t131 = t136 * t132;
t148 = (t128 * t131 + t136) * t127;
t183 = t132 * t135;
t150 = t175 * t183;
t192 = (-t148 * t130 * qJD(4) + t128 * t150) / t122 ^ 2;
t191 = t108 * t124;
t189 = t113 * t125;
t188 = t117 * t134;
t178 = t136 * t137;
t177 = qJD(1) * t135;
t176 = qJD(1) * t136;
t174 = qJD(4) * t134;
t172 = qJD(4) * t136;
t170 = -0.2e1 * (-t184 * t194 + (-t133 * t134 * t172 - t150) * t102) / t100 ^ 2;
t169 = 0.2e1 * t196;
t168 = -0.2e1 * t194;
t167 = 0.2e1 * t192;
t166 = t98 * t174;
t164 = t119 * t127 * t132;
t158 = t101 * t170;
t157 = 0.2e1 * t113 * t193;
t156 = -0.2e1 * t127 * t192;
t151 = t135 * t164;
t149 = t159 * t137;
t147 = t109 * t189 - t191;
t95 = 0.1e1 / t97;
t146 = t147 * t95;
t145 = -t135 * t172 - t152 * t137;
t112 = -t134 * t180 - t186;
t111 = t125 * t137 - t134 * t181;
t106 = t119 * t198;
t91 = (-t197 * t136 * t117 + t118 * t151) * t137;
t90 = -t135 * t188 + t187 - (t118 * t179 - t188) * t106;
t89 = t167 * t198 + (-qJD(1) * t149 + 0.2e1 * t148 * t173) * t119;
t1 = [t156 * t178 + (-t127 * t135 * t176 - qJD(4) * t149) * t119, 0, 0, t89, 0, 0; (-t101 * t166 + (t158 + (-qJD(1) * t91 - t88) * t195) * t136) * t135 + (t91 * t136 * t98 * t168 + (-t91 * t166 + (t91 * t170 + ((t136 * t167 - t94 * t151 + t197 * t174) * t117 + (t156 * t183 + t136 * t94 + (-t131 * t161 + (-t94 - 0.2e1 * t173) * t136) * t119) * t118) * t98 * t137) * t136) * t102 + (t101 + ((-t130 + t133) * t118 * t164 + t197 * t165) * t102) * t98 * t176) * t137, 0, 0 (-t101 * t98 * t177 + (t158 + (-qJD(4) * t90 - t88) * t195) * t137) * t134 + (t90 * t137 * t102 * t170 + ((qJD(4) * t101 + t90 * t168) * t137 + (-t90 * t177 + ((-t106 * t175 + t135 * t89) * t118 + ((t106 * t135 - 0.1e1) * t94 + (t106 - t135) * qJD(4)) * t117) * t178) * t102) * t98 + ((-t89 - t175) * t117 + (-t154 * t106 - t155) * t118) * t182 * t195) * t136, 0, 0; (-t108 * t111 + t112 * t190) * t169 + (t112 * t157 - t153 * t108 * t180 + t145 * t191 + (-t153 * t113 * t181 - t111 * t93 + t112 * t92 - t145 * t189) * t109) * t95, 0, 0, t134 * t146 * t171 + (t146 * t177 + (t147 * t169 + ((qJD(6) * t108 + t157) * t125 + (t125 * t92 + (qJD(6) * t113 - t93) * t124) * t109) * t95) * t137) * t136, 0, -0.2e1 * t196 + 0.2e1 * (-t109 * t92 * t95 + (-t109 * t196 - t95 * t193) * t113) * t113;];
JaD_rot  = t1;
