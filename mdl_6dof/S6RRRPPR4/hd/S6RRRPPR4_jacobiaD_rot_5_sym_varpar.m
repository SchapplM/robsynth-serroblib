% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPPR4_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR4_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:05:09
% EndTime: 2019-02-26 22:05:10
% DurationCPUTime: 1.03s
% Computational Cost: add. (4821->123), mult. (6168->268), div. (1114->15), fcn. (7752->9), ass. (0->114)
t158 = qJ(3) + pkin(10);
t157 = cos(t158);
t240 = 0.2e1 * t157;
t156 = sin(t158);
t165 = cos(qJ(2));
t166 = cos(qJ(1));
t216 = t166 * t157;
t236 = sin(qJ(1));
t142 = t236 * t156 + t165 * t216;
t136 = 0.1e1 / t142 ^ 2;
t164 = sin(qJ(2));
t159 = t164 ^ 2;
t163 = t166 ^ 2;
t221 = t159 * t163;
t201 = t136 * t221;
t132 = 0.1e1 + t201;
t190 = qJD(1) * t236;
t214 = qJD(2) * t165;
t174 = t159 * t166 * t190 - t163 * t164 * t214;
t213 = qJD(2) * t166;
t194 = t164 * t213;
t177 = t165 * t190 + t194;
t189 = t236 * qJD(3);
t217 = t166 * t156;
t121 = (-qJD(3) * t165 + qJD(1)) * t217 + (t189 - t177) * t157;
t135 = 0.1e1 / t142;
t231 = t121 * t135 * t136;
t184 = t221 * t231;
t239 = (-t174 * t136 - t184) / t132 ^ 2;
t219 = t164 * t166;
t196 = t236 * t165;
t138 = t156 * t196 + t216;
t182 = t156 * t189;
t210 = qJD(3) * t166;
t192 = t157 * t210;
t120 = t138 * qJD(1) + t156 * t194 - t165 * t192 - t182;
t141 = -t236 * t157 + t165 * t217;
t153 = 0.1e1 / t156;
t154 = 0.1e1 / t156 ^ 2;
t160 = 0.1e1 / t164;
t161 = 0.1e1 / t164 ^ 2;
t195 = t161 * t214;
t212 = qJD(3) * t157;
t224 = t153 * t160;
t238 = t141 * (t154 * t160 * t212 + t153 * t195) + t120 * t224;
t220 = t164 * t156;
t128 = atan2(-t138, t220);
t125 = cos(t128);
t124 = sin(t128);
t230 = t124 * t138;
t119 = t125 * t220 - t230;
t116 = 0.1e1 / t119;
t117 = 0.1e1 / t119 ^ 2;
t237 = 0.2e1 * t141;
t133 = t138 ^ 2;
t222 = t154 * t161;
t129 = t133 * t222 + 0.1e1;
t126 = 0.1e1 / t129;
t211 = qJD(3) * t164;
t178 = t156 * t214 + t157 * t211;
t199 = t138 * t222;
t197 = t236 * t164;
t183 = qJD(2) * t197;
t215 = qJD(1) * t166;
t122 = (t189 * t165 - t190) * t157 + (t215 * t165 - t183 - t210) * t156;
t202 = t122 * t224;
t108 = (t178 * t199 - t202) * t126;
t175 = -t108 * t138 + t178;
t104 = (-t108 * t220 - t122) * t124 + t175 * t125;
t118 = t116 * t117;
t235 = t104 * t118;
t155 = t153 * t154;
t162 = t160 / t159;
t193 = t161 * t212;
t234 = (t122 * t199 + (-t154 * t162 * t214 - t155 * t193) * t133) / t129 ^ 2;
t233 = t117 * t141;
t232 = t120 * t117;
t229 = t124 * t141;
t228 = t124 * t164;
t227 = t125 * t138;
t226 = t125 * t141;
t225 = t125 * t165;
t223 = t154 * t157;
t218 = t165 * t166;
t134 = t141 ^ 2;
t114 = t117 * t134 + 0.1e1;
t209 = 0.2e1 * (-t134 * t235 - t141 * t232) / t114 ^ 2;
t208 = -0.2e1 * t234;
t207 = 0.2e1 * t239;
t206 = t118 * t237;
t205 = t160 * t234;
t204 = t117 * t229;
t200 = t138 * t224;
t198 = t153 * t161 * t165;
t180 = t138 * t198 + t236;
t115 = t180 * t126;
t191 = t236 - t115;
t188 = t116 * t209;
t187 = t117 * t209;
t186 = t219 * t237;
t185 = t153 * t205;
t140 = t157 * t196 - t217;
t181 = t138 * t223 - t140 * t153;
t179 = t136 * t140 * t166 - t236 * t135;
t130 = 0.1e1 / t132;
t123 = t142 * qJD(1) - t157 * t183 - t165 * t182 - t192;
t112 = 0.1e1 / t114;
t111 = t181 * t160 * t126;
t107 = (-t124 + (t125 * t200 + t124) * t126) * t141;
t106 = -t115 * t227 + (t191 * t228 + t225) * t156;
t105 = t125 * t157 * t164 - t124 * t140 + (-t124 * t220 - t227) * t111;
t103 = t180 * t208 + (t122 * t198 + t215 + (-t154 * t165 * t193 + (-0.2e1 * t162 * t165 ^ 2 - t160) * t153 * qJD(2)) * t138) * t126;
t101 = -0.2e1 * t181 * t205 + (-t181 * t195 + (t122 * t223 - t123 * t153 + (t140 * t223 + (-0.2e1 * t155 * t157 ^ 2 - t153) * t138) * qJD(3)) * t160) * t126;
t1 = [t238 * t126 + t185 * t237, t103, t101, 0, 0, 0; t138 * t188 + (-t122 * t116 + (t104 * t138 + t107 * t120) * t117) * t112 + (t107 * t187 + (0.2e1 * t107 * t235 + (t120 * t126 - t120 - (-t108 * t126 * t200 + t208) * t141) * t117 * t124 + (-(-0.2e1 * t138 * t185 - t108) * t233 + (-(t108 + t202) * t141 + t238 * t138) * t117 * t126) * t125) * t112) * t141, t106 * t141 * t187 + (-(-t103 * t227 + (t108 * t230 - t122 * t125) * t115) * t233 + (t104 * t206 + t232) * t106 + (-t116 * t219 - (-t115 * t228 + t124 * t197 + t225) * t233) * t212) * t112 + (t188 * t219 + ((-t116 * t213 - (qJD(2) * t191 - t108) * t204) * t165 + (t116 * t190 + (t166 * t104 - (-t103 + t215) * t229 - (t108 * t191 - qJD(2)) * t226) * t117) * t164) * t112) * t156 (t105 * t233 - t116 * t142) * t209 + (t105 * t232 + t121 * t116 + (t105 * t206 - t117 * t142) * t104 - (t157 * t214 - t156 * t211 - t101 * t138 - t111 * t122 + (-t111 * t220 - t140) * t108) * t117 * t226 - (-t123 + (-t101 * t156 - t108 * t157) * t164 - t175 * t111) * t204) * t112, 0, 0, 0; t179 * t164 * t207 + (-t179 * t214 + ((qJD(1) * t135 + 0.2e1 * t140 * t231) * t166 + (-t121 * t236 - t123 * t166 + t140 * t190) * t136) * t164) * t130 (t135 * t218 + t157 * t201) * t207 + (t184 * t240 + t177 * t135 + (qJD(3) * t156 * t221 + t121 * t218 + t174 * t240) * t136) * t130, t136 * t186 * t239 + (t186 * t231 + (t120 * t219 + (t164 * t190 - t165 * t213) * t141) * t136) * t130, 0, 0, 0;];
JaD_rot  = t1;
