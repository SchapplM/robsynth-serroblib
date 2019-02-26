% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPP5_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_jacobiaD_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:58:48
% EndTime: 2019-02-26 20:58:49
% DurationCPUTime: 0.98s
% Computational Cost: add. (3877->124), mult. (6168->272), div. (1114->15), fcn. (7752->9), ass. (0->115)
t159 = pkin(9) + qJ(3);
t158 = cos(t159);
t164 = sin(qJ(4));
t238 = sin(qJ(1));
t198 = t238 * t164;
t165 = cos(qJ(4));
t166 = cos(qJ(1));
t218 = t166 * t165;
t142 = t158 * t218 + t198;
t136 = 0.1e1 / t142 ^ 2;
t157 = sin(t159);
t153 = t157 ^ 2;
t163 = t166 ^ 2;
t226 = t153 * t163;
t203 = t136 * t226;
t131 = 0.1e1 + t203;
t190 = qJD(1) * t238;
t215 = qJD(3) * t166;
t194 = t157 * t215;
t176 = t158 * t190 + t194;
t189 = t238 * qJD(4);
t219 = t166 * t164;
t121 = (-qJD(4) * t158 + qJD(1)) * t219 + (t189 - t176) * t165;
t135 = 0.1e1 / t142;
t233 = t121 * t135 * t136;
t184 = t226 * t233;
t195 = qJD(3) * t157 * t163;
t241 = (-t184 + (-t153 * t166 * t190 + t158 * t195) * t136) / t131 ^ 2;
t221 = t157 * t166;
t138 = t158 * t198 + t218;
t181 = t164 * t189;
t212 = qJD(4) * t166;
t192 = t165 * t212;
t120 = t138 * qJD(1) - t158 * t192 + t164 * t194 - t181;
t197 = t238 * t165;
t141 = t158 * t219 - t197;
t154 = 0.1e1 / t157;
t160 = 0.1e1 / t164;
t161 = 0.1e1 / t164 ^ 2;
t213 = qJD(4) * t165;
t193 = t161 * t213;
t155 = 0.1e1 / t157 ^ 2;
t216 = qJD(3) * t158;
t196 = t155 * t216;
t225 = t154 * t160;
t240 = (t154 * t193 + t160 * t196) * t141 + t120 * t225;
t222 = t157 * t164;
t130 = atan2(-t138, t222);
t125 = cos(t130);
t124 = sin(t130);
t232 = t124 * t138;
t119 = t125 * t222 - t232;
t116 = 0.1e1 / t119;
t117 = 0.1e1 / t119 ^ 2;
t239 = 0.2e1 * t141;
t133 = t138 ^ 2;
t223 = t155 * t161;
t132 = t133 * t223 + 0.1e1;
t128 = 0.1e1 / t132;
t177 = t157 * t213 + t164 * t216;
t201 = t138 * t223;
t199 = t238 * t157;
t182 = qJD(3) * t199;
t183 = t165 * t190;
t217 = qJD(1) * t166;
t122 = t165 * t189 * t158 - t183 + (t217 * t158 - t182 - t212) * t164;
t204 = t122 * t225;
t108 = (t177 * t201 - t204) * t128;
t174 = -t108 * t138 + t177;
t104 = (-t108 * t222 - t122) * t124 + t174 * t125;
t118 = t116 * t117;
t237 = t104 * t118;
t156 = t154 / t153;
t162 = t160 * t161;
t236 = (t122 * t201 + (-t155 * t162 * t213 - t156 * t161 * t216) * t133) / t132 ^ 2;
t235 = t117 * t141;
t234 = t120 * t117;
t231 = t124 * t141;
t230 = t124 * t157;
t229 = t125 * t138;
t228 = t125 * t141;
t227 = t125 * t158;
t224 = t155 * t158;
t220 = t161 * t165;
t214 = qJD(4) * t164;
t134 = t141 ^ 2;
t114 = t117 * t134 + 0.1e1;
t211 = 0.2e1 * (-t134 * t237 - t141 * t234) / t114 ^ 2;
t210 = 0.2e1 * t241;
t209 = -0.2e1 * t236;
t208 = t118 * t239;
t207 = t154 * t236;
t206 = t117 * t231;
t202 = t138 * t225;
t200 = t160 * t224;
t179 = t138 * t200 + t238;
t115 = t179 * t128;
t191 = t238 - t115;
t188 = t116 * t211;
t187 = t117 * t211;
t186 = t221 * t239;
t185 = t160 * t207;
t140 = t158 * t197 - t219;
t180 = t138 * t220 - t140 * t160;
t178 = t136 * t140 * t166 - t238 * t135;
t126 = 0.1e1 / t131;
t123 = t142 * qJD(1) - t158 * t181 - t165 * t182 - t192;
t112 = 0.1e1 / t114;
t111 = t180 * t154 * t128;
t107 = (-t124 + (t125 * t202 + t124) * t128) * t141;
t106 = -t115 * t229 + (t191 * t230 + t227) * t164;
t105 = t125 * t157 * t165 - t124 * t140 + (-t124 * t222 - t229) * t111;
t103 = t179 * t209 + (t122 * t200 + t217 + (-t193 * t224 + (-0.2e1 * t156 * t158 ^ 2 - t154) * t160 * qJD(3)) * t138) * t128;
t101 = -0.2e1 * t180 * t207 + (-t180 * t196 + (t122 * t220 - t123 * t160 + (t140 * t220 + (-0.2e1 * t162 * t165 ^ 2 - t160) * t138) * qJD(4)) * t154) * t128;
t1 = [t240 * t128 + t185 * t239, 0, t103, t101, 0, 0; t138 * t188 + (-t122 * t116 + (t104 * t138 + t107 * t120) * t117) * t112 + (t107 * t187 + (0.2e1 * t107 * t237 + (t120 * t128 - t120 - (-t108 * t128 * t202 + t209) * t141) * t117 * t124 + (-(-0.2e1 * t138 * t185 - t108) * t235 + (-(t108 + t204) * t141 + t240 * t138) * t117 * t128) * t125) * t112) * t141, 0, t106 * t141 * t187 + (-(-t103 * t229 + (t108 * t232 - t122 * t125) * t115) * t235 + (t104 * t208 + t234) * t106 + (-t116 * t221 - (-t115 * t230 + t124 * t199 + t227) * t235) * t213) * t112 + (t188 * t221 + ((-t116 * t215 - (t191 * qJD(3) - t108) * t206) * t158 + (t116 * t190 + (t166 * t104 - (-t103 + t217) * t231 - (t191 * t108 - qJD(3)) * t228) * t117) * t157) * t112) * t164 (t105 * t235 - t116 * t142) * t211 + (t105 * t234 + t121 * t116 + (t105 * t208 - t117 * t142) * t104 - (t165 * t216 - t157 * t214 - t101 * t138 - t111 * t122 + (-t111 * t222 - t140) * t108) * t117 * t228 - (-t123 + (-t101 * t164 - t108 * t165) * t157 - t174 * t111) * t206) * t112, 0, 0; t178 * t157 * t210 + (-t178 * t216 + ((qJD(1) * t135 + 0.2e1 * t140 * t233) * t166 + (-t238 * t121 - t123 * t166 + t140 * t190) * t136) * t157) * t126, 0 (t135 * t158 * t166 + t165 * t203) * t210 + (0.2e1 * t165 * t184 + t176 * t135 + ((t121 * t166 - 0.2e1 * t165 * t195) * t158 + (t163 * t214 + 0.2e1 * t166 * t183) * t153) * t136) * t126, t136 * t186 * t241 + (t186 * t233 + (t120 * t221 + (t157 * t190 - t158 * t215) * t141) * t136) * t126, 0, 0;];
JaD_rot  = t1;
