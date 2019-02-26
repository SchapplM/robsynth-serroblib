% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRPR4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:12:19
% EndTime: 2019-02-26 20:12:20
% DurationCPUTime: 0.97s
% Computational Cost: add. (4054->111), mult. (9671->227), div. (577->12), fcn. (12382->13), ass. (0->107)
t246 = sin(pkin(11));
t248 = cos(pkin(11));
t251 = sin(qJ(2));
t249 = cos(pkin(6));
t253 = cos(qJ(2));
t280 = t249 * t253;
t232 = -t246 * t251 + t248 * t280;
t225 = t232 * qJD(2);
t281 = t249 * t251;
t233 = t246 * t253 + t248 * t281;
t250 = sin(qJ(3));
t247 = sin(pkin(6));
t284 = t247 * t250;
t271 = t248 * t284;
t252 = cos(qJ(3));
t277 = qJD(3) * t252;
t203 = -qJD(3) * t271 + t225 * t250 + t233 * t277;
t283 = t247 * t252;
t217 = t233 * t250 + t248 * t283;
t215 = t217 ^ 2;
t236 = -t249 * t252 + t251 * t284;
t230 = 0.1e1 / t236 ^ 2;
t213 = t215 * t230 + 0.1e1;
t211 = 0.1e1 / t213;
t237 = t249 * t250 + t251 * t283;
t278 = qJD(2) * t253;
t270 = t247 * t278;
t222 = t237 * qJD(3) + t250 * t270;
t229 = 0.1e1 / t236;
t288 = t217 * t230;
t181 = (-t203 * t229 + t222 * t288) * t211;
t214 = atan2(-t217, t236);
t209 = sin(t214);
t210 = cos(t214);
t265 = -t209 * t236 - t210 * t217;
t177 = t265 * t181 - t209 * t203 + t210 * t222;
t193 = -t209 * t217 + t210 * t236;
t190 = 0.1e1 / t193;
t191 = 0.1e1 / t193 ^ 2;
t301 = t177 * t190 * t191;
t272 = t246 * t281;
t235 = t248 * t253 - t272;
t262 = -t235 * t250 + t246 * t283;
t300 = -0.2e1 * t262 * t301;
t282 = t247 * t253;
t261 = -t229 * t232 + t282 * t288;
t299 = t250 * t261;
t287 = t222 * t229 * t230;
t298 = -0.2e1 * (t203 * t288 - t215 * t287) / t213 ^ 2;
t221 = t235 * t252 + t246 * t284;
t234 = t246 * t280 + t248 * t251;
t244 = qJ(4) + pkin(12) + qJ(6);
t242 = sin(t244);
t243 = cos(t244);
t202 = t221 * t243 + t234 * t242;
t198 = 0.1e1 / t202;
t199 = 0.1e1 / t202 ^ 2;
t228 = -qJD(2) * t272 + t248 * t278;
t245 = qJD(4) + qJD(6);
t267 = t221 * t245 - t228;
t227 = t234 * qJD(2);
t206 = t262 * qJD(3) - t227 * t252;
t268 = t234 * t245 + t206;
t188 = t268 * t242 + t267 * t243;
t201 = t221 * t242 - t234 * t243;
t197 = t201 ^ 2;
t196 = t197 * t199 + 0.1e1;
t293 = t199 * t201;
t189 = -t267 * t242 + t268 * t243;
t296 = t189 * t198 * t199;
t297 = (t188 * t293 - t197 * t296) / t196 ^ 2;
t295 = t191 * t262;
t294 = t198 * t242;
t292 = t201 * t243;
t205 = t221 * qJD(3) - t227 * t250;
t291 = t205 * t191;
t290 = t209 * t262;
t289 = t210 * t262;
t286 = t234 * t250;
t285 = t234 * t252;
t279 = qJD(2) * t251;
t216 = t262 ^ 2;
t187 = t216 * t191 + 0.1e1;
t276 = 0.2e1 * (-t216 * t301 - t262 * t291) / t187 ^ 2;
t275 = -0.2e1 * t297;
t273 = t201 * t296;
t269 = -0.2e1 * t217 * t287;
t266 = t245 * t285 - t227;
t264 = t199 * t292 - t294;
t219 = t233 * t252 - t271;
t263 = -t219 * t229 + t237 * t288;
t260 = qJD(3) * t286 - t228 * t252 + t235 * t245;
t226 = t233 * qJD(2);
t223 = -t236 * qJD(3) + t252 * t270;
t208 = t235 * t242 - t243 * t285;
t207 = -t235 * t243 - t242 * t285;
t204 = -t217 * qJD(3) + t225 * t252;
t194 = 0.1e1 / t196;
t184 = 0.1e1 / t187;
t183 = t211 * t299;
t182 = t263 * t211;
t179 = (-t209 * t232 + t210 * t282) * t250 + t265 * t183;
t178 = t265 * t182 - t209 * t219 + t210 * t237;
t176 = t263 * t298 + (t237 * t269 - t204 * t229 + (t203 * t237 + t217 * t223 + t219 * t222) * t230) * t211;
t174 = t298 * t299 + (t261 * t277 + (t269 * t282 + t226 * t229 + (t222 * t232 + (t203 * t253 - t217 * t279) * t247) * t230) * t250) * t211;
t173 = t275 + 0.2e1 * (t188 * t199 * t194 + (-t194 * t296 - t199 * t297) * t201) * t201;
t1 = [0, t174, t176, 0, 0, 0; 0 (-t179 * t295 + t190 * t286) * t276 + ((-t228 * t250 - t234 * t277) * t190 + (-t291 + t300) * t179 + (t286 * t177 + (-t174 * t217 - t183 * t203 + (-t250 * t279 + t253 * t277) * t247 + (-t183 * t236 - t232 * t250) * t181) * t289 + (-t232 * t277 - t174 * t236 - t183 * t222 + t226 * t250 + (t183 * t217 - t250 * t282) * t181) * t290) * t191) * t184 (-t178 * t295 - t190 * t221) * t276 + (t178 * t300 + t206 * t190 + (-t221 * t177 - t178 * t205 + (-t176 * t217 - t182 * t203 + t223 + (-t182 * t236 - t219) * t181) * t289 + (-t176 * t236 - t182 * t222 - t204 + (t182 * t217 - t237) * t181) * t290) * t191) * t184, 0, 0, 0; 0, 0.2e1 * (-t198 * t207 + t208 * t293) * t297 + (0.2e1 * t208 * t273 - t266 * t198 * t243 + t260 * t294 + (-t266 * t201 * t242 - t208 * t188 - t207 * t189 - t260 * t292) * t199) * t194, -t264 * t262 * t275 + (t264 * t205 - ((-t198 * t245 - 0.2e1 * t273) * t243 + (t188 * t243 + (-t201 * t245 + t189) * t242) * t199) * t262) * t194, t173, 0, t173;];
JaD_rot  = t1;
