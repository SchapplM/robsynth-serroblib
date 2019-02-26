% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRP5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:27:27
% EndTime: 2019-02-26 21:27:28
% DurationCPUTime: 1.16s
% Computational Cost: add. (4821->118), mult. (6168->265), div. (1114->15), fcn. (7752->9), ass. (0->113)
t150 = pkin(9) + qJ(5);
t148 = sin(t150);
t156 = sin(qJ(2));
t159 = cos(qJ(1));
t149 = cos(t150);
t157 = sin(qJ(1));
t208 = t157 * t149;
t137 = t148 * t159 + t156 * t208;
t158 = cos(qJ(2));
t207 = t158 * t149;
t125 = atan2(t137, t207);
t121 = sin(t125);
t122 = cos(t125);
t134 = t137 ^ 2;
t146 = 0.1e1 / t149 ^ 2;
t153 = 0.1e1 / t158 ^ 2;
t213 = t146 * t153;
t126 = t134 * t213 + 0.1e1;
t123 = 0.1e1 / t126;
t145 = 0.1e1 / t149;
t186 = t145 * t153 * t156;
t173 = t137 * t186 + t157;
t112 = t173 * t123;
t205 = t112 - t157;
t232 = t205 * t121 * t158 + t122 * t156;
t209 = t156 * t159;
t136 = t148 * t209 + t208;
t132 = 0.1e1 / t136 ^ 2;
t151 = t158 ^ 2;
t155 = t159 ^ 2;
t210 = t151 * t155;
t189 = t132 * t210;
t129 = 0.1e1 + t189;
t201 = qJD(2) * t158;
t203 = qJD(1) * t159;
t169 = -t151 * t157 * t203 - t155 * t156 * t201;
t177 = qJD(1) * t156 + qJD(5);
t178 = qJD(5) * t156 + qJD(1);
t200 = qJD(2) * t159;
t182 = t158 * t200;
t211 = t149 * t159;
t120 = t178 * t211 + (-t177 * t157 + t182) * t148;
t131 = 0.1e1 / t136;
t222 = t120 * t131 * t132;
t176 = t210 * t222;
t231 = (t169 * t132 - t176) / t129 ^ 2;
t152 = 0.1e1 / t158;
t212 = t148 * t157;
t138 = -t156 * t212 + t211;
t214 = t146 * t148;
t171 = t137 * t214 + t138 * t145;
t230 = t152 * t171;
t206 = t158 * t159;
t220 = t121 * t137;
t116 = t122 * t207 + t220;
t113 = 0.1e1 / t116;
t114 = 0.1e1 / t116 ^ 2;
t185 = t149 * t209;
t135 = -t185 + t212;
t228 = 0.2e1 * t135;
t130 = t135 ^ 2;
t111 = t114 * t130 + 0.1e1;
t119 = t137 * qJD(1) + t136 * qJD(5) - t149 * t182;
t223 = t119 * t114;
t197 = qJD(5) * t158;
t202 = qJD(2) * t156;
t170 = -t148 * t197 - t149 * t202;
t187 = t137 * t213;
t174 = t178 * t157;
t183 = t157 * t201;
t198 = qJD(5) * t149;
t117 = -qJD(1) * t185 + t148 * t174 - t149 * t183 - t159 * t198;
t215 = t145 * t152;
t190 = t117 * t215;
t105 = (-t170 * t187 - t190) * t123;
t168 = -t105 * t137 - t170;
t101 = (-t105 * t207 - t117) * t121 - t168 * t122;
t115 = t113 * t114;
t226 = t101 * t115;
t227 = 0.1e1 / t111 ^ 2 * (-t130 * t226 + t135 * t223);
t147 = t145 * t146;
t154 = t152 / t151;
t199 = qJD(5) * t148;
t181 = t153 * t199;
t225 = (-t117 * t187 + (t146 * t154 * t202 + t147 * t181) * t134) / t126 ^ 2;
t224 = t114 * t135;
t221 = t121 * t135;
t218 = t122 * t135;
t217 = t122 * t137;
t204 = qJD(1) * t157;
t196 = 0.2e1 * t227;
t195 = -0.2e1 * t225;
t194 = 0.2e1 * t231;
t193 = t115 * t228;
t192 = t113 * t227;
t191 = t114 * t221;
t188 = t137 * t215;
t184 = t153 * t202;
t180 = t114 * t196;
t179 = t206 * t228;
t175 = 0.2e1 * t215 * t225;
t172 = -t132 * t138 * t159 - t131 * t157;
t167 = t119 * t215 - (-t146 * t152 * t199 - t145 * t184) * t135;
t127 = 0.1e1 / t129;
t118 = t149 * t174 + (t177 * t159 + t183) * t148;
t109 = 0.1e1 / t111;
t108 = t123 * t230;
t104 = (-t121 + (-t122 * t188 + t121) * t123) * t135;
t103 = t112 * t217 - t149 * t232;
t102 = -t122 * t148 * t158 + t121 * t138 + (-t121 * t207 + t217) * t108;
t100 = t173 * t195 + (-t117 * t186 + t203 + (t146 * t156 * t181 + (0.2e1 * t154 * t156 ^ 2 + t152) * t145 * qJD(2)) * t137) * t123;
t98 = t195 * t230 + (t171 * t184 + (-t117 * t214 - t118 * t145 + (t138 * t214 + (0.2e1 * t147 * t148 ^ 2 + t145) * t137) * qJD(5)) * t152) * t123;
t1 = [-t167 * t123 + t135 * t175, t100, 0, 0, t98, 0; -0.2e1 * t137 * t192 + (-t117 * t113 + (-t101 * t137 - t104 * t119) * t114) * t109 + (t104 * t180 + (0.2e1 * t104 * t226 + (-t119 * t123 + t119 - (t105 * t123 * t188 + t195) * t135) * t114 * t121 + (-(t137 * t175 - t105) * t224 + (-(t105 + t190) * t135 + t167 * t137) * t114 * t123) * t122) * t109) * t135, t103 * t135 * t180 + (-(t100 * t217 + (-t105 * t220 - t117 * t122) * t112) * t224 + (t101 * t193 - t223) * t103 + (t113 * t206 - t224 * t232) * t199) * t109 + (0.2e1 * t192 * t206 + ((t113 * t200 - (t205 * qJD(2) + t105) * t191) * t156 + (t113 * t204 + (t159 * t101 - (-t100 + t203) * t221 - (-t205 * t105 - qJD(2)) * t218) * t114) * t158) * t109) * t149, 0, 0 (t102 * t224 - t113 * t136) * t196 + (-t102 * t223 + t120 * t113 + (t102 * t193 - t114 * t136) * t101 - (t148 * t202 - t149 * t197 - t108 * t117 + t137 * t98 + (-t108 * t207 + t138) * t105) * t114 * t218 - (-t118 + (t105 * t148 - t149 * t98) * t158 + t168 * t108) * t191) * t109, 0; t172 * t158 * t194 + (t172 * t202 + ((qJD(1) * t131 - 0.2e1 * t138 * t222) * t159 + (-t118 * t159 + (-qJD(1) * t138 - t120) * t157) * t132) * t158) * t127 (-t131 * t209 - t148 * t189) * t194 + (-0.2e1 * t148 * t176 + (-t156 * t204 + t182) * t131 + (-t120 * t209 + 0.2e1 * t148 * t169 + t198 * t210) * t132) * t127, 0, 0, t132 * t179 * t231 + (t179 * t222 + (-t119 * t206 + (t156 * t200 + t158 * t204) * t135) * t132) * t127, 0;];
JaD_rot  = t1;
