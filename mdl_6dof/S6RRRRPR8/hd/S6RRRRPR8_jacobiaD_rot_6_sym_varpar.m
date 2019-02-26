% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR8_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR8_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:34:43
% EndTime: 2019-02-26 22:34:44
% DurationCPUTime: 1.12s
% Computational Cost: add. (2960->117), mult. (4961->248), div. (518->12), fcn. (5877->11), ass. (0->113)
t209 = qJD(3) + qJD(4);
t297 = qJD(6) - t209;
t219 = sin(qJ(2));
t211 = t219 ^ 2;
t222 = cos(qJ(2));
t214 = 0.1e1 / t222 ^ 2;
t273 = t211 * t214;
t220 = sin(qJ(1));
t212 = t220 ^ 2;
t205 = t212 * t273 + 0.1e1;
t213 = 0.1e1 / t222;
t270 = t213 * t219;
t295 = t219 * t273;
t232 = qJD(2) * (t213 * t295 + t270);
t223 = cos(qJ(1));
t262 = qJD(1) * t223;
t271 = t211 * t220;
t240 = t262 * t271;
t279 = (t212 * t232 + t214 * t240) / t205 ^ 2;
t296 = -0.2e1 * t279;
t247 = 0.1e1 + t273;
t294 = t220 * t247;
t217 = qJ(3) + qJ(4);
t208 = cos(t217);
t264 = t222 * t223;
t207 = sin(t217);
t268 = t220 * t207;
t197 = t208 * t264 + t268;
t221 = cos(qJ(6));
t293 = t297 * t221;
t218 = sin(qJ(6));
t292 = t297 * t218;
t265 = t220 * t222;
t233 = t207 * t265 + t208 * t223;
t259 = qJD(2) * t223;
t248 = t219 * t259;
t167 = t233 * qJD(1) - t197 * t209 + t207 * t248;
t242 = -qJD(1) * t222 + t209;
t243 = t209 * t222 - qJD(1);
t276 = t207 * t223;
t168 = -t243 * t276 + (t242 * t220 - t248) * t208;
t267 = t220 * t208;
t196 = t207 * t264 - t267;
t236 = t196 * t221 - t197 * t218;
t161 = t236 * qJD(6) - t167 * t218 + t168 * t221;
t181 = t196 * t218 + t197 * t221;
t173 = 0.1e1 / t181;
t234 = t207 * t218 + t208 * t221;
t235 = t207 * t221 - t208 * t218;
t174 = 0.1e1 / t181 ^ 2;
t282 = t174 * t236;
t291 = t235 * t173 - t234 * t282;
t266 = t220 * t219;
t204 = atan2(t266, t222);
t200 = cos(t204);
t199 = sin(t204);
t252 = t199 * t266;
t188 = t200 * t222 + t252;
t185 = 0.1e1 / t188;
t186 = 0.1e1 / t188 ^ 2;
t290 = 0.2e1 * t219;
t202 = 0.1e1 / t205;
t289 = t202 - 0.1e1;
t160 = t181 * qJD(6) + t167 * t221 + t168 * t218;
t172 = t236 ^ 2;
t165 = t172 * t174 + 0.1e1;
t175 = t173 * t174;
t285 = t161 * t175;
t288 = (-t160 * t282 - t172 * t285) / t165 ^ 2;
t216 = t223 ^ 2;
t272 = t211 * t216;
t184 = t186 * t272 + 0.1e1;
t260 = qJD(2) * t222;
t249 = t219 * t262;
t261 = qJD(2) * t220;
t171 = ((t220 * t260 + t249) * t213 + t261 * t273) * t202;
t277 = t200 * t219;
t158 = (t171 * t220 - qJD(2)) * t277 + (t249 + (-t171 + t261) * t222) * t199;
t286 = t158 * t185 * t186;
t287 = (-t272 * t286 + (t216 * t219 * t260 - t240) * t186) / t184 ^ 2;
t284 = t171 * t199;
t283 = t171 * t219;
t269 = t219 * t223;
t192 = t234 * t269;
t281 = t174 * t192;
t280 = t186 * t219;
t190 = t202 * t294;
t278 = t190 * t220;
t263 = qJD(1) * t220;
t256 = 0.2e1 * t288;
t255 = -0.2e1 * t286;
t254 = -0.2e1 * t175 * t236;
t253 = t186 * t269;
t251 = t202 * t211 * t213;
t246 = -0.2e1 * t219 * t287;
t245 = t161 * t254;
t244 = t213 * t296;
t241 = t220 * t251;
t239 = t247 * t223;
t195 = -t208 * t265 + t276;
t237 = -t195 * t218 - t221 * t233;
t177 = t195 * t221 - t218 * t233;
t231 = t219 * t261 + t242 * t223;
t191 = t235 * t269;
t182 = 0.1e1 / t184;
t170 = t231 * t208 + t243 * t268;
t169 = t231 * t207 - t243 * t267;
t166 = (-t289 * t219 * t199 + t200 * t241) * t223;
t163 = 0.1e1 / t165;
t162 = t199 * t265 - t277 + (-t199 * t222 + t200 * t266) * t190;
t159 = t294 * t296 + (qJD(1) * t239 + 0.2e1 * t220 * t232) * t202;
t155 = (t173 * t181 + t236 * t282) * t256 + (-t161 * t173 - t236 * t245 + (0.2e1 * t236 * t160 + t181 * t161) * t174) * t163;
t1 = [t244 * t269 + (qJD(2) * t239 - t263 * t270) * t202, t159, 0, 0, 0, 0; (t185 * t246 + (t185 * t260 + (-qJD(1) * t166 - t158) * t280) * t182) * t220 + (t186 * t246 * t166 + (((-t171 * t241 - t289 * t260 + t279 * t290) * t199 + (t244 * t271 + t283 + (-t283 + (t290 + t295) * t261) * t202) * t200) * t253 + (t186 * t260 + t219 * t255) * t166 + (t185 + ((-t212 + t216) * t200 * t251 + t289 * t252) * t186) * t219 * qJD(1)) * t182) * t223, 0.2e1 * (-t162 * t280 + t185 * t222) * t223 * t287 + ((t185 * t263 + (qJD(2) * t162 + t158) * t223 * t186) * t222 + (t185 * t259 + (t159 * t200 * t220 - t199 * t261 - t278 * t284 + t284 + (qJD(2) * t199 + t200 * t262) * t190) * t253 + (-t186 * t263 + t223 * t255) * t162 + ((-t159 + t262) * t199 + ((-0.1e1 + t278) * qJD(2) + (-t190 + t220) * t171) * t200) * t186 * t264) * t219) * t182, 0, 0, 0, 0; (t173 * t237 - t177 * t282) * t256 + ((t177 * qJD(6) - t169 * t221 + t170 * t218) * t173 + t177 * t245 + (t237 * t161 + (t237 * qJD(6) + t169 * t218 + t170 * t221) * t236 - t177 * t160) * t174) * t163 (-t173 * t191 + t236 * t281) * t256 + (t160 * t281 + (-t174 * t191 - t192 * t254) * t161 + t291 * t222 * t259 + (-t291 * t263 + ((-t293 * t173 + t292 * t282) * t208 + (-t292 * t173 - t293 * t282) * t207) * t223) * t219) * t163, t155, t155, 0, -0.2e1 * t288 - 0.2e1 * (t160 * t163 * t174 - (-t163 * t285 - t174 * t288) * t236) * t236;];
JaD_rot  = t1;
