% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRP9_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP9_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:48:11
% EndTime: 2019-02-26 20:48:12
% DurationCPUTime: 0.92s
% Computational Cost: add. (4821->119), mult. (6168->264), div. (1114->15), fcn. (7752->9), ass. (0->112)
t157 = sin(qJ(1));
t158 = cos(qJ(3));
t207 = t157 * t158;
t150 = pkin(9) + qJ(5);
t148 = sin(t150);
t156 = sin(qJ(3));
t228 = cos(qJ(1));
t189 = t228 * t156;
t149 = cos(t150);
t208 = t157 * t149;
t137 = t148 * t189 + t208;
t206 = t158 * t148;
t125 = atan2(t137, t206);
t122 = cos(t125);
t121 = sin(t125);
t219 = t121 * t137;
t116 = t122 * t206 + t219;
t113 = 0.1e1 / t116;
t136 = t228 * t148 + t156 * t208;
t131 = 0.1e1 / t136;
t145 = 0.1e1 / t148;
t153 = 0.1e1 / t158;
t114 = 0.1e1 / t116 ^ 2;
t132 = 0.1e1 / t136 ^ 2;
t146 = 0.1e1 / t148 ^ 2;
t209 = t157 * t148;
t135 = -t228 * t149 + t156 * t209;
t130 = t135 ^ 2;
t111 = t114 * t130 + 0.1e1;
t184 = qJD(1) * t228;
t202 = qJD(3) * t158;
t170 = -t156 * t184 - t157 * t202;
t167 = t228 * qJD(5) - t170;
t179 = qJD(5) * t156 + qJD(1);
t119 = t167 * t148 + t179 * t208;
t222 = t119 * t114;
t134 = t137 ^ 2;
t154 = 0.1e1 / t158 ^ 2;
t212 = t146 * t154;
t126 = t134 * t212 + 0.1e1;
t123 = 0.1e1 / t126;
t200 = qJD(5) * t158;
t204 = qJD(3) * t156;
t171 = -t148 * t204 + t149 * t200;
t191 = t137 * t212;
t188 = t228 * t158;
t176 = qJD(3) * t188;
t177 = t149 * t189;
t117 = -qJD(5) * t177 - t148 * t176 - t149 * t184 + (qJD(1) * t156 + qJD(5)) * t209;
t214 = t145 * t153;
t194 = t117 * t214;
t105 = (-t171 * t191 - t194) * t123;
t169 = -t105 * t137 - t171;
t101 = (-t105 * t206 - t117) * t121 - t169 * t122;
t115 = t113 * t114;
t226 = t101 * t115;
t227 = 0.1e1 / t111 ^ 2 * (-t130 * t226 + t135 * t222);
t147 = t145 * t146;
t152 = t158 ^ 2;
t155 = t153 / t152;
t201 = qJD(5) * t149;
t186 = t154 * t201;
t225 = (-t117 * t191 + (t146 * t155 * t204 - t147 * t186) * t134) / t126 ^ 2;
t151 = t157 ^ 2;
t211 = t151 * t152;
t193 = t132 * t211;
t129 = 0.1e1 + t193;
t168 = -t151 * t156 * t202 + t152 * t157 * t184;
t120 = t167 * t149 - t179 * t209;
t221 = t120 * t131 * t132;
t178 = t211 * t221;
t224 = (t168 * t132 - t178) / t129 ^ 2;
t223 = t114 * t135;
t220 = t121 * t135;
t218 = t121 * t158;
t217 = t122 * t135;
t216 = t122 * t137;
t215 = t122 * t156;
t213 = t146 * t149;
t210 = t156 * t157;
t205 = qJD(1) * t157;
t203 = qJD(3) * t157;
t199 = 0.2e1 * t227;
t198 = -0.2e1 * t225;
t197 = 0.2e1 * t224;
t196 = 0.2e1 * t115 * t135;
t195 = t114 * t220;
t192 = t137 * t214;
t190 = t145 * t154 * t156;
t187 = t154 * t204;
t173 = t137 * t190 + t228;
t112 = t173 * t123;
t185 = t228 - t112;
t183 = -0.2e1 * t113 * t227;
t182 = t114 * t199;
t181 = 0.2e1 * t153 * t225;
t180 = -0.2e1 * t135 * t207;
t175 = t145 * t181;
t138 = t177 - t209;
t174 = t137 * t213 - t138 * t145;
t172 = t132 * t138 * t157 - t228 * t131;
t166 = t119 * t214 - (t146 * t153 * t201 - t145 * t187) * t135;
t127 = 0.1e1 / t129;
t118 = t136 * qJD(1) + t137 * qJD(5) - t149 * t176;
t109 = 0.1e1 / t111;
t108 = t174 * t153 * t123;
t104 = (-t121 + (-t122 * t192 + t121) * t123) * t135;
t103 = t112 * t216 + (t185 * t218 - t215) * t148;
t102 = t122 * t149 * t158 + t121 * t138 - (-t121 * t206 + t216) * t108;
t100 = t173 * t198 + (-t117 * t190 - t205 + (-t146 * t156 * t186 + (0.2e1 * t155 * t156 ^ 2 + t153) * t145 * qJD(3)) * t137) * t123;
t98 = t174 * t181 + (-t174 * t187 + (t117 * t213 - t118 * t145 + (-t138 * t213 + (0.2e1 * t147 * t149 ^ 2 + t145) * t137) * qJD(5)) * t153) * t123;
t1 = [-t166 * t123 + t135 * t175, 0, t100, 0, t98, 0; t137 * t183 + (-t117 * t113 + (-t101 * t137 - t104 * t119) * t114) * t109 + (t104 * t182 + (0.2e1 * t104 * t226 + (-t119 * t123 + t119 - (t105 * t123 * t192 + t198) * t135) * t114 * t121 + (-(t137 * t175 - t105) * t223 + (-(t105 + t194) * t135 + t166 * t137) * t114 * t123) * t122) * t109) * t135, 0, t103 * t135 * t182 + (-(t100 * t216 + (-t105 * t219 - t117 * t122) * t112) * t223 + (t101 * t196 - t222) * t103 + (t113 * t207 - (-t112 * t218 + t121 * t188 - t215) * t223) * t201) * t109 + (t183 * t207 + ((-t113 * t203 - (-t185 * qJD(3) + t105) * t195) * t156 + (t113 * t184 + (-t157 * t101 - (-t100 - t205) * t220 - (t185 * t105 - qJD(3)) * t217) * t114) * t158) * t109) * t148, 0 (t102 * t223 - t113 * t136) * t199 + (-t102 * t222 + t120 * t113 + (t102 * t196 - t136 * t114) * t101 - (-t149 * t204 - t148 * t200 + t108 * t117 + t137 * t98 + (t108 * t206 + t138) * t105) * t114 * t217 - (-t118 + (-t105 * t149 - t148 * t98) * t158 - t169 * t108) * t195) * t109, 0; t172 * t158 * t197 + (t172 * t204 + ((-qJD(1) * t131 + 0.2e1 * t138 * t221) * t157 + (t118 * t157 - t228 * t120 - t138 * t184) * t132) * t158) * t127, 0 (t131 * t210 + t149 * t193) * t197 + (0.2e1 * t149 * t178 + t170 * t131 + (qJD(5) * t148 * t211 + t120 * t210 - 0.2e1 * t149 * t168) * t132) * t127, 0, t132 * t180 * t224 + (t180 * t221 + (t119 * t207 + (-t156 * t203 + t158 * t184) * t135) * t132) * t127, 0;];
JaD_rot  = t1;
