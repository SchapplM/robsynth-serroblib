% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPPR8_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR8_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_jacobiaD_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:07:39
% EndTime: 2019-02-26 22:07:41
% DurationCPUTime: 1.13s
% Computational Cost: add. (3667->120), mult. (11094->248), div. (691->12), fcn. (14128->11), ass. (0->107)
t203 = sin(qJ(3));
t205 = cos(qJ(3));
t202 = cos(pkin(6));
t204 = sin(qJ(2));
t265 = sin(qJ(1));
t232 = t265 * t204;
t224 = t202 * t232;
t206 = cos(qJ(2));
t207 = cos(qJ(1));
t244 = t207 * t206;
t217 = t224 - t244;
t201 = sin(pkin(6));
t233 = t201 * t265;
t218 = t203 * t217 + t205 * t233;
t266 = -0.2e1 * t218;
t246 = t201 * t207;
t273 = -qJD(1) * t224 - qJD(2) * t232 + (qJD(2) * t202 + qJD(1)) * t244 - qJD(3) * t246;
t272 = 0.2e1 * t205;
t231 = t265 * t206;
t245 = t207 * t204;
t216 = t202 * t231 + t245;
t180 = t203 * t233 - t205 * t217;
t169 = 0.1e1 / t180;
t188 = t202 * t245 + t231;
t164 = t188 * qJD(1) + t216 * qJD(2);
t230 = qJD(1) * t246;
t151 = t218 * qJD(3) - t164 * t205 + t203 * t230;
t170 = 0.1e1 / t180 ^ 2;
t270 = t151 * t170;
t258 = t169 * t270;
t271 = -0.2e1 * t216 * t258;
t174 = t188 * t203 + t205 * t246;
t248 = t201 * t204;
t185 = -t202 * t205 + t203 * t248;
t161 = atan2(-t174, t185);
t154 = sin(t161);
t155 = cos(t161);
t149 = -t154 * t174 + t155 * t185;
t147 = 0.1e1 / t149 ^ 2;
t168 = t218 ^ 2;
t145 = t168 * t147 + 0.1e1;
t150 = t180 * qJD(3) - t164 * t203 - t205 * t230;
t259 = t150 * t147;
t167 = t174 ^ 2;
t182 = 0.1e1 / t185 ^ 2;
t160 = t167 * t182 + 0.1e1;
t156 = 0.1e1 / t160;
t223 = qJD(1) * t233;
t242 = qJD(3) * t205;
t152 = t188 * t242 + t273 * t203 - t205 * t223;
t186 = t202 * t203 + t205 * t248;
t247 = t201 * t206;
t229 = qJD(2) * t247;
t172 = t186 * qJD(3) + t203 * t229;
t181 = 0.1e1 / t185;
t251 = t174 * t182;
t221 = -t152 * t181 + t172 * t251;
t138 = t221 * t156;
t222 = -t154 * t185 - t155 * t174;
t134 = t222 * t138 - t154 * t152 + t155 * t172;
t146 = 0.1e1 / t149;
t148 = t146 * t147;
t263 = t134 * t148;
t241 = 0.2e1 * (-t168 * t263 - t218 * t259) / t145 ^ 2;
t269 = t172 * t182;
t187 = -t202 * t244 + t232;
t219 = t181 * t187 + t247 * t251;
t268 = t203 * t219;
t153 = (-qJD(3) * t188 + t223) * t203 + t273 * t205;
t267 = -0.2e1 * t174;
t184 = t216 ^ 2;
t250 = t184 * t170;
t162 = 0.1e1 + t250;
t163 = t187 * qJD(1) + t217 * qJD(2);
t234 = t184 * t258;
t254 = t170 * t216;
t262 = (-t163 * t254 - t234) / t162 ^ 2;
t253 = t181 * t269;
t261 = (t152 * t251 - t167 * t253) / t160 ^ 2;
t260 = t147 * t218;
t257 = t154 * t218;
t256 = t155 * t218;
t252 = t174 * t181;
t249 = t216 * t203;
t243 = qJD(2) * t204;
t240 = -0.2e1 * t261;
t239 = t148 * t266;
t238 = t181 * t261;
t237 = t147 * t257;
t236 = t147 * t256;
t227 = t253 * t267;
t176 = t188 * t205 - t203 * t246;
t225 = t254 * t262;
t220 = -t176 * t181 + t186 * t251;
t214 = -t154 + (t155 * t252 + t154) * t156;
t173 = -t185 * qJD(3) + t205 * t229;
t165 = t216 * qJD(1) + t188 * qJD(2);
t158 = 0.1e1 / t162;
t143 = 0.1e1 / t145;
t142 = t156 * t268;
t141 = t220 * t156;
t137 = t214 * t218;
t136 = (t154 * t187 + t155 * t247) * t203 + t222 * t142;
t135 = t222 * t141 - t154 * t176 + t155 * t186;
t133 = t220 * t240 + (t186 * t227 - t153 * t181 + (t152 * t186 + t172 * t176 + t173 * t174) * t182) * t156;
t131 = t240 * t268 + (t219 * t242 + (t227 * t247 + t165 * t181 + (-t172 * t187 + (t152 * t206 - t174 * t243) * t201) * t182) * t203) * t156;
t1 = [t238 * t266 + (-t150 * t181 - t218 * t269) * t156, t131, t133, 0, 0, 0; t174 * t146 * t241 + (-t152 * t146 + (t134 * t174 + t137 * t150) * t147) * t143 - (-t137 * t147 * t241 + (-0.2e1 * t137 * t263 + (-t138 * t156 * t252 + t240) * t237 + (t238 * t267 - t138 + (t138 - t221) * t156) * t236 - t214 * t259) * t143) * t218 (-t136 * t260 + t146 * t249) * t241 + (-t136 * t259 + (t163 * t203 - t216 * t242) * t146 + (t136 * t239 + t147 * t249) * t134 + (-t131 * t174 - t142 * t152 + (-t203 * t243 + t206 * t242) * t201 + (-t142 * t185 + t187 * t203) * t138) * t236 + (t187 * t242 - t131 * t185 - t142 * t172 + t165 * t203 + (t142 * t174 - t203 * t247) * t138) * t237) * t143 (-t135 * t260 - t146 * t180) * t241 + (t135 * t134 * t239 + t151 * t146 + (-t180 * t134 - t135 * t150 + (-t133 * t174 - t141 * t152 + t173 + (-t141 * t185 - t176) * t138) * t256 + (-t133 * t185 - t141 * t172 - t153 + (t141 * t174 - t186) * t138) * t257) * t147) * t143, 0, 0, 0; -0.2e1 * t187 * t169 * t262 + 0.2e1 * t176 * t225 + (-(-t163 * t170 + t271) * t176 - t153 * t254 + t165 * t169 - t187 * t270) * t158, 0.2e1 * (-t169 * t217 + t205 * t250) * t262 + (t234 * t272 + t164 * t169 + (qJD(3) * t184 * t203 + t163 * t216 * t272 - t151 * t217) * t170) * t158, t225 * t266 + (t218 * t271 + (-t150 * t216 - t163 * t218) * t170) * t158, 0, 0, 0;];
JaD_rot  = t1;
