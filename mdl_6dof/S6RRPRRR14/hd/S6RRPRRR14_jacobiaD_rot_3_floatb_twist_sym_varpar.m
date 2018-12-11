% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPRRR14
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR14_jacobiaD_rot_3_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_rot_3_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_rot_3_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiaD_rot_3_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:21
% EndTime: 2018-12-10 18:38:22
% DurationCPUTime: 1.05s
% Computational Cost: add. (6585->112), mult. (8580->230), div. (424->12), fcn. (8986->21), ass. (0->117)
t284 = pkin(6) + qJ(2);
t276 = sin(t284);
t303 = -qJD(2) / 0.2e1;
t229 = t276 * t303;
t285 = pkin(6) - qJ(2);
t277 = sin(t285);
t268 = qJD(2) * t277;
t216 = t268 / 0.2e1 + t229;
t250 = cos(qJ(2));
t251 = cos(qJ(1));
t237 = cos(t284) / 0.2e1;
t241 = cos(t285);
t225 = t241 / 0.2e1 + t237;
t248 = sin(qJ(2));
t249 = sin(qJ(1));
t265 = t249 * t225 + t251 * t248;
t287 = qJD(2) * t249;
t196 = qJD(1) * t265 - t251 * t216 + t250 * t287;
t246 = sin(pkin(7));
t301 = sin(pkin(6));
t302 = cos(pkin(7));
t269 = t302 * t301;
t264 = qJD(1) * t269;
t184 = -t196 * t246 - t249 * t264;
t209 = -t251 * t225 + t249 * t248;
t259 = -t209 * t246 + t251 * t269;
t199 = t259 ^ 2;
t270 = t276 / 0.2e1;
t208 = -(t270 + t277 / 0.2e1) * t246 + cos(pkin(6)) * t302;
t206 = 0.1e1 / t208 ^ 2;
t191 = t199 * t206 + 0.1e1;
t205 = 0.1e1 / t208;
t207 = t205 * t206;
t218 = qJD(2) * t237 + t241 * t303;
t288 = t218 * t246;
t298 = (t184 * t206 * t259 + t199 * t207 * t288) / t191 ^ 2;
t306 = -0.2e1 * t298;
t223 = t270 - t277 / 0.2e1;
t266 = t251 * t223 + t249 * t250;
t224 = t237 - t241 / 0.2e1;
t290 = t259 * t224;
t305 = t246 * (-t205 * t266 + t206 * t290);
t217 = t225 * qJD(2);
t286 = qJD(2) * t251;
t194 = qJD(1) * t266 + t249 * t217 + t248 * t286;
t192 = atan2(t259, t208);
t185 = sin(t192);
t186 = cos(t192);
t169 = t185 * t259 + t186 * t208;
t166 = 0.1e1 / t169;
t214 = t249 * t223 - t251 * t250;
t242 = pkin(7) + pkin(14);
t233 = sin(t242) / 0.2e1;
t243 = pkin(7) - pkin(14);
t238 = sin(t243);
t220 = t233 - t238 / 0.2e1;
t234 = cos(t243) / 0.2e1;
t239 = cos(t242);
t222 = t234 - t239 / 0.2e1;
t247 = cos(pkin(14));
t279 = t249 * t301;
t182 = -t214 * t247 - t220 * t265 + t222 * t279;
t176 = 0.1e1 / t182;
t167 = 0.1e1 / t169 ^ 2;
t177 = 0.1e1 / t182 ^ 2;
t203 = -t246 * t265 - t249 * t269;
t200 = t203 ^ 2;
t164 = t200 * t167 + 0.1e1;
t193 = qJD(1) * t209 - t249 * t216 - t250 * t286;
t183 = t193 * t246 - t251 * t264;
t294 = t183 * t167;
t189 = 0.1e1 / t191;
t280 = t206 * t288;
t261 = t184 * t205 + t259 * t280;
t160 = t261 * t189;
t267 = -t185 * t208 + t186 * t259;
t156 = t160 * t267 + t185 * t184 - t186 * t288;
t299 = t156 * t166 * t167;
t300 = (-t200 * t299 + t203 * t294) / t164 ^ 2;
t297 = t167 * t203;
t275 = qJD(1) * t301;
t271 = t251 * t275;
t174 = t193 * t220 - t194 * t247 + t222 * t271;
t296 = t174 * t176 * t177;
t219 = t233 + t238 / 0.2e1;
t221 = t234 + t239 / 0.2e1;
t245 = sin(pkin(14));
t181 = -t214 * t245 - t219 * t279 + t221 * t265;
t295 = t177 * t181;
t293 = t185 * t203;
t292 = t186 * t203;
t291 = t259 * t205;
t289 = t218 * t246 ^ 2;
t283 = -0.2e1 * t300;
t282 = -0.2e1 * t299;
t175 = t181 ^ 2;
t172 = t175 * t177 + 0.1e1;
t173 = -t193 * t221 - t194 * t245 - t219 * t271;
t281 = 0.2e1 * (t173 * t295 - t175 * t296) / t172 ^ 2;
t278 = t251 * t301;
t274 = t205 * t306;
t273 = 0.2e1 * t181 * t296;
t272 = t249 * t275;
t260 = t185 + (t186 * t291 - t185) * t189;
t258 = qJD(1) * t214 - t251 * t217 + t248 * t287;
t215 = t229 - t268 / 0.2e1;
t188 = t214 * t220 - t247 * t265;
t187 = -t214 * t221 - t245 * t265;
t180 = t209 * t220 + t222 * t278 - t247 * t266;
t179 = -t209 * t221 - t219 * t278 - t245 * t266;
t170 = 0.1e1 / t172;
t162 = 0.1e1 / t164;
t161 = t189 * t305;
t159 = t260 * t203;
t157 = (-t185 * t266 - t186 * t224) * t246 + t267 * t161;
t155 = t305 * t306 + (0.2e1 * t207 * t289 * t290 + t258 * t205 * t246 + (-t266 * t289 + (t184 * t224 + t215 * t259) * t246) * t206) * t189;
t1 = [t203 * t274 + (t183 * t205 + t203 * t280) * t189, t155, 0, 0, 0, 0; t259 * t166 * t283 + (t184 * t166 + (-t156 * t259 + t159 * t183) * t167) * t162 + ((t159 * t282 + t260 * t294) * t162 + (t159 * t283 + ((-t160 * t189 * t291 + 0.2e1 * t298) * t293 + (t259 * t274 + t160 + (-t160 + t261) * t189) * t292) * t162) * t167) * t203, 0.2e1 * (t166 * t214 * t246 - t157 * t297) * t300 + ((t267 * t155 + (-t160 * t169 + t184 * t186) * t161) * t297 + (t203 * t282 + t294) * t157 + (-t194 * t166 + (t214 * t156 + (-t160 * t266 - t215) * t292 + (t160 * t224 + t161 * t218 + t258) * t293) * t167) * t246) * t162, 0, 0, 0, 0; (-t176 * t179 + t180 * t295) * t281 + ((-t196 * t221 + t219 * t272 + t245 * t258) * t176 + t180 * t273 + (-t179 * t174 - (t196 * t220 - t222 * t272 + t247 * t258) * t181 - t180 * t173) * t177) * t170 (-t176 * t187 + t188 * t295) * t281 + ((t193 * t245 - t194 * t221) * t176 + t188 * t273 + (-t187 * t174 - (t193 * t247 + t194 * t220) * t181 - t188 * t173) * t177) * t170, 0, 0, 0, 0;];
JaD_rot  = t1;
