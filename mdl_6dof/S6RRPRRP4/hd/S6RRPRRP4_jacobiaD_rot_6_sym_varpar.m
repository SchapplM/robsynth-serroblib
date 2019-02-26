% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRP4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:47:45
% EndTime: 2019-02-26 21:47:47
% DurationCPUTime: 1.17s
% Computational Cost: add. (9776->126), mult. (8378->269), div. (1558->15), fcn. (10537->9), ass. (0->118)
t204 = qJ(2) + pkin(10);
t200 = cos(t204);
t206 = qJ(4) + qJ(5);
t201 = sin(t206);
t278 = sin(qJ(1));
t237 = t278 * t201;
t202 = cos(t206);
t207 = cos(qJ(1));
t257 = t207 * t202;
t181 = t200 * t257 + t237;
t175 = 0.1e1 / t181 ^ 2;
t199 = sin(t204);
t192 = t199 ^ 2;
t205 = t207 ^ 2;
t266 = t192 * t205;
t245 = t175 * t266;
t171 = 0.1e1 + t245;
t203 = qJD(4) + qJD(5);
t230 = qJD(1) * t278;
t254 = qJD(2) * t207;
t232 = t199 * t254;
t217 = t200 * t230 + t232;
t235 = t278 * t203;
t258 = t207 * t201;
t160 = (-t200 * t203 + qJD(1)) * t258 + (t235 - t217) * t202;
t174 = 0.1e1 / t181;
t273 = t160 * t174 * t175;
t225 = t266 * t273;
t233 = qJD(2) * t199 * t205;
t282 = (-t225 + (-t192 * t207 * t230 + t200 * t233) * t175) / t171 ^ 2;
t197 = 0.1e1 / t201 ^ 2;
t259 = t202 * t203;
t241 = t197 * t259;
t261 = t199 * t207;
t177 = t200 * t237 + t257;
t224 = t201 * t235;
t239 = t203 * t257;
t159 = t177 * qJD(1) - t200 * t239 + t201 * t232 - t224;
t236 = t278 * t202;
t180 = t200 * t258 - t236;
t193 = 0.1e1 / t199;
t196 = 0.1e1 / t201;
t194 = 0.1e1 / t199 ^ 2;
t255 = qJD(2) * t200;
t234 = t194 * t255;
t265 = t193 * t196;
t281 = (t193 * t241 + t196 * t234) * t180 + t159 * t265;
t262 = t199 * t201;
t167 = atan2(-t177, t262);
t164 = cos(t167);
t163 = sin(t167);
t272 = t163 * t177;
t158 = t164 * t262 - t272;
t155 = 0.1e1 / t158;
t156 = 0.1e1 / t158 ^ 2;
t280 = -0.2e1 * t177;
t279 = 0.2e1 * t180;
t172 = t177 ^ 2;
t264 = t194 * t197;
t168 = t172 * t264 + 0.1e1;
t165 = 0.1e1 / t168;
t218 = t199 * t259 + t201 * t255;
t243 = t177 * t264;
t238 = t199 * t278;
t222 = qJD(2) * t238;
t223 = t202 * t230;
t256 = qJD(1) * t207;
t161 = -t201 * t222 - t203 * t258 - t223 + (t201 * t256 + t202 * t235) * t200;
t246 = t161 * t265;
t147 = (t218 * t243 - t246) * t165;
t215 = -t147 * t177 + t218;
t142 = (-t147 * t262 - t161) * t163 + t215 * t164;
t157 = t155 * t156;
t277 = t142 * t157;
t195 = t193 / t192;
t240 = t196 * t241;
t276 = (t161 * t243 + (-t195 * t197 * t255 - t194 * t240) * t172) / t168 ^ 2;
t275 = t156 * t180;
t274 = t159 * t156;
t271 = t163 * t180;
t270 = t163 * t199;
t269 = t164 * t177;
t268 = t164 * t180;
t267 = t164 * t200;
t263 = t194 * t200;
t260 = t201 * t203;
t173 = t180 ^ 2;
t153 = t156 * t173 + 0.1e1;
t253 = 0.2e1 * (-t173 * t277 - t180 * t274) / t153 ^ 2;
t252 = -0.2e1 * t276;
t251 = 0.2e1 * t282;
t250 = t157 * t279;
t249 = t193 * t276;
t248 = t156 * t271;
t244 = t177 * t265;
t242 = t196 * t263;
t220 = t177 * t242 + t278;
t154 = t220 * t165;
t231 = t278 - t154;
t229 = t155 * t253;
t228 = t156 * t253;
t227 = t261 * t279;
t226 = t196 * t249;
t179 = t200 * t236 - t258;
t221 = t177 * t197 * t202 - t179 * t196;
t219 = t175 * t179 * t207 - t278 * t174;
t169 = 0.1e1 / t171;
t162 = t181 * qJD(1) - t200 * t224 - t202 * t222 - t239;
t151 = 0.1e1 / t153;
t150 = t221 * t193 * t165;
t146 = (-t163 + (t164 * t244 + t163) * t165) * t180;
t145 = -t154 * t269 + (t231 * t270 + t267) * t201;
t144 = t164 * t199 * t202 - t163 * t179 + (-t163 * t262 - t269) * t150;
t143 = t175 * t227 * t282 + (t227 * t273 + (t159 * t261 + (t199 * t230 - t200 * t254) * t180) * t175) * t169;
t141 = t220 * t252 + (t161 * t242 + t256 + (-t241 * t263 + (-0.2e1 * t195 * t200 ^ 2 - t193) * t196 * qJD(2)) * t177) * t165;
t139 = -0.2e1 * t221 * t249 + (-t221 * t234 + ((-t177 * t203 - t162) * t196 + (t240 * t280 + (t179 * t203 + t161) * t197) * t202) * t193) * t165;
t138 = (t144 * t275 - t155 * t181) * t253 + (t144 * t274 + t160 * t155 + (t144 * t250 - t156 * t181) * t142 - (t202 * t255 - t199 * t260 - t139 * t177 - t150 * t161 + (-t150 * t262 - t179) * t147) * t156 * t268 - (-t162 + (-t139 * t201 - t147 * t202) * t199 - t215 * t150) * t248) * t151;
t1 = [t165 * t281 + t226 * t279, t141, 0, t139, t139, 0; t177 * t229 + (-t161 * t155 + (t142 * t177 + t146 * t159) * t156) * t151 + (t146 * t228 + (0.2e1 * t146 * t277 + (t159 * t165 - t159 - (-t147 * t165 * t244 + t252) * t180) * t156 * t163 + (-(t226 * t280 - t147) * t275 + (-(t147 + t246) * t180 + t281 * t177) * t156 * t165) * t164) * t151) * t180, t145 * t180 * t228 + (-(-t141 * t269 + (t147 * t272 - t161 * t164) * t154) * t275 + (-t155 * t261 - (-t154 * t270 + t163 * t238 + t267) * t275) * t259 + (t142 * t250 + t274) * t145) * t151 + (t229 * t261 + ((-t155 * t254 - (t231 * qJD(2) - t147) * t248) * t200 + (t155 * t230 + (t207 * t142 - (-t141 + t256) * t271 - (t231 * t147 - qJD(2)) * t268) * t156) * t199) * t151) * t201, 0, t138, t138, 0; t219 * t199 * t251 + (-t219 * t255 + ((qJD(1) * t174 + 0.2e1 * t179 * t273) * t207 + (-t278 * t160 - t162 * t207 + t179 * t230) * t175) * t199) * t169 (t174 * t200 * t207 + t202 * t245) * t251 + (0.2e1 * t202 * t225 + t217 * t174 + ((t160 * t207 - 0.2e1 * t202 * t233) * t200 + (t205 * t260 + 0.2e1 * t207 * t223) * t192) * t175) * t169, 0, t143, t143, 0;];
JaD_rot  = t1;
