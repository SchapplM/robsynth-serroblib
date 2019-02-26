% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR5
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRR5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:06:28
% EndTime: 2019-02-26 20:06:29
% DurationCPUTime: 0.95s
% Computational Cost: add. (4054->111), mult. (9671->227), div. (577->12), fcn. (12382->13), ass. (0->107)
t245 = sin(pkin(11));
t247 = cos(pkin(11));
t250 = sin(qJ(2));
t248 = cos(pkin(6));
t252 = cos(qJ(2));
t279 = t248 * t252;
t231 = -t245 * t250 + t247 * t279;
t224 = t231 * qJD(2);
t280 = t248 * t250;
t232 = t245 * t252 + t247 * t280;
t249 = sin(qJ(3));
t246 = sin(pkin(6));
t283 = t246 * t249;
t270 = t247 * t283;
t251 = cos(qJ(3));
t276 = qJD(3) * t251;
t202 = -qJD(3) * t270 + t224 * t249 + t232 * t276;
t282 = t246 * t251;
t216 = t232 * t249 + t247 * t282;
t214 = t216 ^ 2;
t235 = -t248 * t251 + t250 * t283;
t229 = 0.1e1 / t235 ^ 2;
t212 = t214 * t229 + 0.1e1;
t210 = 0.1e1 / t212;
t236 = t248 * t249 + t250 * t282;
t277 = qJD(2) * t252;
t269 = t246 * t277;
t221 = t236 * qJD(3) + t249 * t269;
t228 = 0.1e1 / t235;
t287 = t216 * t229;
t180 = (-t202 * t228 + t221 * t287) * t210;
t213 = atan2(-t216, t235);
t208 = sin(t213);
t209 = cos(t213);
t264 = -t208 * t235 - t209 * t216;
t176 = t264 * t180 - t208 * t202 + t209 * t221;
t192 = -t208 * t216 + t209 * t235;
t189 = 0.1e1 / t192;
t190 = 0.1e1 / t192 ^ 2;
t300 = t176 * t189 * t190;
t271 = t245 * t280;
t234 = t247 * t252 - t271;
t261 = -t234 * t249 + t245 * t282;
t299 = -0.2e1 * t261 * t300;
t281 = t246 * t252;
t260 = -t228 * t231 + t281 * t287;
t298 = t249 * t260;
t286 = t221 * t228 * t229;
t297 = -0.2e1 * (t202 * t287 - t214 * t286) / t212 ^ 2;
t220 = t234 * t251 + t245 * t283;
t233 = t245 * t279 + t247 * t250;
t243 = pkin(12) + qJ(5) + qJ(6);
t241 = sin(t243);
t242 = cos(t243);
t201 = t220 * t242 + t233 * t241;
t197 = 0.1e1 / t201;
t198 = 0.1e1 / t201 ^ 2;
t227 = -qJD(2) * t271 + t247 * t277;
t244 = qJD(5) + qJD(6);
t266 = t220 * t244 - t227;
t226 = t233 * qJD(2);
t205 = t261 * qJD(3) - t226 * t251;
t267 = t233 * t244 + t205;
t187 = t267 * t241 + t266 * t242;
t200 = t220 * t241 - t233 * t242;
t196 = t200 ^ 2;
t195 = t196 * t198 + 0.1e1;
t292 = t198 * t200;
t188 = -t266 * t241 + t267 * t242;
t295 = t188 * t197 * t198;
t296 = (t187 * t292 - t196 * t295) / t195 ^ 2;
t294 = t190 * t261;
t293 = t197 * t241;
t291 = t200 * t242;
t204 = t220 * qJD(3) - t226 * t249;
t290 = t204 * t190;
t289 = t208 * t261;
t288 = t209 * t261;
t285 = t233 * t249;
t284 = t233 * t251;
t278 = qJD(2) * t250;
t215 = t261 ^ 2;
t186 = t190 * t215 + 0.1e1;
t275 = 0.2e1 * (-t215 * t300 - t261 * t290) / t186 ^ 2;
t274 = -0.2e1 * t296;
t272 = t200 * t295;
t268 = -0.2e1 * t216 * t286;
t265 = t244 * t284 - t226;
t263 = t198 * t291 - t293;
t218 = t232 * t251 - t270;
t262 = -t218 * t228 + t236 * t287;
t259 = qJD(3) * t285 - t227 * t251 + t234 * t244;
t225 = t232 * qJD(2);
t222 = -t235 * qJD(3) + t251 * t269;
t207 = t234 * t241 - t242 * t284;
t206 = -t234 * t242 - t241 * t284;
t203 = -t216 * qJD(3) + t224 * t251;
t193 = 0.1e1 / t195;
t183 = 0.1e1 / t186;
t182 = t210 * t298;
t181 = t262 * t210;
t178 = (-t208 * t231 + t209 * t281) * t249 + t264 * t182;
t177 = t264 * t181 - t208 * t218 + t209 * t236;
t175 = t262 * t297 + (t236 * t268 - t203 * t228 + (t202 * t236 + t216 * t222 + t218 * t221) * t229) * t210;
t173 = t297 * t298 + (t260 * t276 + (t268 * t281 + t225 * t228 + (t221 * t231 + (t202 * t252 - t216 * t278) * t246) * t229) * t249) * t210;
t172 = t274 + 0.2e1 * (t187 * t198 * t193 + (-t193 * t295 - t198 * t296) * t200) * t200;
t1 = [0, t173, t175, 0, 0, 0; 0 (-t178 * t294 + t189 * t285) * t275 + ((-t227 * t249 - t233 * t276) * t189 + (-t290 + t299) * t178 + (t285 * t176 + (-t173 * t216 - t182 * t202 + (-t249 * t278 + t252 * t276) * t246 + (-t182 * t235 - t231 * t249) * t180) * t288 + (-t231 * t276 - t173 * t235 - t182 * t221 + t225 * t249 + (t182 * t216 - t249 * t281) * t180) * t289) * t190) * t183 (-t177 * t294 - t189 * t220) * t275 + (t177 * t299 + t205 * t189 + (-t220 * t176 - t177 * t204 + (-t175 * t216 - t181 * t202 + t222 + (-t181 * t235 - t218) * t180) * t288 + (-t175 * t235 - t181 * t221 - t203 + (t181 * t216 - t236) * t180) * t289) * t190) * t183, 0, 0, 0; 0, 0.2e1 * (-t197 * t206 + t207 * t292) * t296 + (0.2e1 * t207 * t272 - t265 * t197 * t242 + t259 * t293 + (-t265 * t200 * t241 - t207 * t187 - t206 * t188 - t259 * t291) * t198) * t193, -t263 * t261 * t274 + (t263 * t204 - ((-t197 * t244 - 0.2e1 * t272) * t242 + (t187 * t242 + (-t200 * t244 + t188) * t241) * t198) * t261) * t193, 0, t172, t172;];
JaD_rot  = t1;
