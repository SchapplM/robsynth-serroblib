% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function JaD_rot = S6RRPPRR8_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR8_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:32:34
% EndTime: 2019-02-26 21:32:35
% DurationCPUTime: 0.82s
% Computational Cost: add. (1540->107), mult. (3805->232), div. (486->12), fcn. (4569->11), ass. (0->107)
t187 = sin(qJ(2));
t178 = t187 ^ 2;
t189 = cos(qJ(2));
t181 = 0.1e1 / t189 ^ 2;
t237 = t178 * t181;
t188 = sin(qJ(1));
t179 = t188 ^ 2;
t173 = t179 * t237 + 0.1e1;
t180 = 0.1e1 / t189;
t234 = t180 * t187;
t256 = t187 * t237;
t198 = qJD(2) * (t180 * t256 + t234);
t190 = cos(qJ(1));
t227 = qJD(1) * t190;
t235 = t178 * t188;
t202 = t227 * t235;
t240 = (t179 * t198 + t181 * t202) / t173 ^ 2;
t257 = -0.2e1 * t240;
t210 = 0.1e1 + t237;
t255 = t188 * t210;
t176 = qJD(5) + qJD(6);
t186 = cos(pkin(10));
t228 = qJD(1) * t188;
t185 = sin(pkin(10));
t230 = t190 * t185;
t254 = -t176 * t230 + t186 * t228;
t229 = t190 * t186;
t253 = t176 * t229 + t185 * t228;
t166 = -t188 * t186 + t189 * t230;
t167 = t188 * t185 + t189 * t229;
t184 = qJ(5) + qJ(6);
t174 = sin(t184);
t175 = cos(t184);
t147 = t166 * t174 + t167 * t175;
t141 = 0.1e1 / t147;
t232 = t188 * t187;
t172 = atan2(t232, t189);
t169 = cos(t172);
t168 = sin(t172);
t219 = t168 * t232;
t154 = t169 * t189 + t219;
t151 = 0.1e1 / t154;
t142 = 0.1e1 / t147 ^ 2;
t152 = 0.1e1 / t154 ^ 2;
t252 = 0.2e1 * t187;
t170 = 0.1e1 / t173;
t251 = t170 - 0.1e1;
t231 = t188 * t189;
t165 = -t186 * t231 + t230;
t224 = qJD(2) * t190;
t211 = t187 * t224;
t206 = t165 * qJD(1) + t166 * t176 - t186 * t211;
t164 = -t185 * t231 - t229;
t207 = -t164 * qJD(1) + t167 * t176 + t185 * t211;
t133 = t206 * t174 + t207 * t175;
t146 = -t166 * t175 + t167 * t174;
t140 = t146 ^ 2;
t137 = t140 * t142 + 0.1e1;
t244 = t142 * t146;
t134 = -t207 * t174 + t206 * t175;
t143 = t141 * t142;
t247 = t134 * t143;
t250 = (t133 * t244 - t140 * t247) / t137 ^ 2;
t183 = t190 ^ 2;
t236 = t178 * t183;
t150 = t152 * t236 + 0.1e1;
t225 = qJD(2) * t189;
t213 = t187 * t227;
t226 = qJD(2) * t188;
t139 = ((t188 * t225 + t213) * t180 + t226 * t237) * t170;
t238 = t169 * t187;
t130 = (t139 * t188 - qJD(2)) * t238 + (t213 + (-t139 + t226) * t189) * t168;
t248 = t130 * t151 * t152;
t249 = (-t236 * t248 + (t183 * t187 * t225 - t202) * t152) / t150 ^ 2;
t246 = t139 * t168;
t245 = t139 * t187;
t199 = -t174 * t185 - t175 * t186;
t233 = t187 * t190;
t162 = t199 * t233;
t243 = t142 * t162;
t242 = t152 * t187;
t241 = t152 * t190;
t156 = t170 * t255;
t239 = t156 * t188;
t223 = 0.2e1 * t250;
t222 = -0.2e1 * t248;
t221 = 0.2e1 * t143 * t146;
t220 = t152 * t233;
t218 = t170 * t178 * t180;
t212 = t187 * t226;
t209 = -0.2e1 * t187 * t249;
t208 = t180 * t257;
t205 = t166 * qJD(1) + t165 * t176 - t185 * t212;
t204 = -t167 * qJD(1) + t164 * t176 + t186 * t212;
t203 = t188 * t218;
t201 = t210 * t190;
t200 = -t174 * t186 + t175 * t185;
t161 = t200 * t233;
t148 = 0.1e1 / t150;
t145 = t164 * t174 + t165 * t175;
t144 = -t164 * t175 + t165 * t174;
t138 = (-t251 * t187 * t168 + t169 * t203) * t190;
t135 = 0.1e1 / t137;
t132 = t168 * t231 - t238 + (-t168 * t189 + t169 * t232) * t156;
t131 = t255 * t257 + (qJD(1) * t201 + 0.2e1 * t188 * t198) * t170;
t127 = -0.2e1 * t250 + 0.2e1 * (t133 * t142 * t135 + (-t135 * t247 - t142 * t250) * t146) * t146;
t1 = [t208 * t233 + (qJD(2) * t201 - t228 * t234) * t170, t131, 0, 0, 0, 0; (t151 * t209 + (t151 * t225 + (-qJD(1) * t138 - t130) * t242) * t148) * t188 + (t152 * t209 * t138 + (((-t139 * t203 - t251 * t225 + t240 * t252) * t168 + (t208 * t235 + t245 + (-t245 + (t252 + t256) * t226) * t170) * t169) * t220 + (t152 * t225 + t187 * t222) * t138 + (t151 + ((-t179 + t183) * t169 * t218 + t251 * t219) * t152) * t187 * qJD(1)) * t148) * t190, 0.2e1 * (-t132 * t242 + t151 * t189) * t190 * t249 + ((t151 * t228 + (qJD(2) * t132 + t130) * t241) * t189 + (t151 * t224 + (t131 * t169 * t188 - t168 * t226 - t239 * t246 + t246 + (qJD(2) * t168 + t169 * t227) * t156) * t220 + (-t152 * t228 + t190 * t222) * t132 + ((-t131 + t227) * t168 + ((-0.1e1 + t239) * qJD(2) + (-t156 + t188) * t139) * t169) * t189 * t241) * t187) * t148, 0, 0, 0, 0; (-t141 * t144 + t145 * t244) * t223 + ((t204 * t174 + t205 * t175) * t141 + t145 * t134 * t221 + (-t144 * t134 - (-t205 * t174 + t204 * t175) * t146 - t145 * t133) * t142) * t135 (-t141 * t161 + t146 * t243) * t223 + (-t133 * t243 + (-t161 * t142 + t162 * t221) * t134 + (t200 * t141 - t199 * t244) * t189 * t224 + ((-t253 * t141 - t254 * t244) * t175 + (t254 * t141 - t253 * t244) * t174) * t187) * t135, 0, 0, t127, t127;];
JaD_rot  = t1;
