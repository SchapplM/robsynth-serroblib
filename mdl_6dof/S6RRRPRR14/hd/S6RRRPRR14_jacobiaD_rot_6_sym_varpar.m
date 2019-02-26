% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR14
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR14_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR14_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:40
% EndTime: 2019-02-26 22:23:42
% DurationCPUTime: 1.39s
% Computational Cost: add. (5342->152), mult. (14155->300), div. (744->12), fcn. (17850->13), ass. (0->133)
t288 = sin(qJ(1));
t285 = cos(pkin(6));
t306 = qJD(2) * t285 + qJD(1);
t287 = sin(qJ(2));
t334 = t288 * t287;
t318 = t285 * t334;
t330 = qJD(2) * t287;
t290 = cos(qJ(2));
t291 = cos(qJ(1));
t332 = t290 * t291;
t284 = sin(pkin(6));
t337 = t284 * t291;
t362 = -qJD(1) * t318 - qJD(3) * t337 - t288 * t330 + t306 * t332;
t333 = t288 * t290;
t335 = t287 * t291;
t268 = t285 * t335 + t333;
t289 = cos(qJ(3));
t286 = sin(qJ(3));
t336 = t286 * t291;
t254 = t268 * t289 - t284 * t336;
t339 = t284 * t289;
t266 = t285 * t286 + t287 * t339;
t241 = atan2(-t254, t266);
t236 = sin(t241);
t237 = cos(t241);
t219 = -t236 * t254 + t237 * t266;
t217 = 0.1e1 / t219 ^ 2;
t270 = -t318 + t332;
t340 = t284 * t286;
t259 = t270 * t289 + t288 * t340;
t250 = t259 ^ 2;
t215 = t217 * t250 + 0.1e1;
t269 = t285 * t333 + t335;
t245 = -t268 * qJD(1) - t269 * qJD(2);
t327 = qJD(3) * t289;
t328 = qJD(3) * t286;
t224 = t245 * t289 - t270 * t328 + (qJD(1) * t336 + t288 * t327) * t284;
t352 = t217 * t259;
t249 = t254 ^ 2;
t263 = 0.1e1 / t266 ^ 2;
t240 = t249 * t263 + 0.1e1;
t238 = 0.1e1 / t240;
t307 = t362 * t289;
t331 = qJD(1) * t284;
t317 = t288 * t331;
t226 = -t268 * t328 + t286 * t317 + t307;
t265 = t285 * t289 - t287 * t340;
t329 = qJD(2) * t290;
t315 = t284 * t329;
t252 = t265 * qJD(3) + t289 * t315;
t262 = 0.1e1 / t266;
t343 = t254 * t263;
t303 = -t226 * t262 + t252 * t343;
t207 = t303 * t238;
t304 = -t236 * t266 - t237 * t254;
t202 = t304 * t207 - t236 * t226 + t237 * t252;
t216 = 0.1e1 / t219;
t218 = t216 * t217;
t356 = t202 * t218;
t326 = 0.2e1 * (t224 * t352 - t250 * t356) / t215 ^ 2;
t361 = t252 * t263;
t319 = t285 * t332;
t267 = t319 - t334;
t338 = t284 * t290;
t300 = -t262 * t267 + t338 * t343;
t360 = t289 * t300;
t258 = t270 * t286 - t288 * t339;
t283 = qJ(5) + qJ(6);
t280 = sin(t283);
t281 = cos(t283);
t235 = t258 * t280 + t269 * t281;
t229 = 0.1e1 / t235;
t230 = 0.1e1 / t235 ^ 2;
t359 = -0.2e1 * t254;
t358 = 0.2e1 * t259;
t345 = t262 * t361;
t355 = (t226 * t343 - t249 * t345) / t240 ^ 2;
t244 = -qJD(1) * t319 - t291 * t329 + t306 * t334;
t282 = qJD(5) + qJD(6);
t309 = t258 * t282 - t244;
t316 = t289 * t331;
t223 = t259 * qJD(3) + t245 * t286 - t291 * t316;
t311 = t269 * t282 - t223;
t211 = -t311 * t280 + t309 * t281;
t354 = t211 * t229 * t230;
t353 = t217 * t224;
t210 = t309 * t280 + t311 * t281;
t234 = -t258 * t281 + t269 * t280;
t228 = t234 ^ 2;
t222 = t228 * t230 + 0.1e1;
t349 = t230 * t234;
t351 = 0.1e1 / t222 ^ 2 * (t210 * t349 - t228 * t354);
t350 = t229 * t281;
t348 = t234 * t280;
t347 = t236 * t259;
t346 = t237 * t259;
t344 = t254 * t262;
t342 = t269 * t286;
t341 = t269 * t289;
t325 = -0.2e1 * t355;
t324 = t218 * t358;
t323 = 0.2e1 * t351;
t322 = t262 * t355;
t321 = t217 * t347;
t320 = t217 * t346;
t313 = 0.2e1 * t234 * t354;
t312 = t345 * t359;
t225 = t268 * t327 + t362 * t286 - t288 * t316;
t310 = t267 * t282 + t225;
t246 = -t269 * qJD(1) - t268 * qJD(2);
t253 = t268 * t286 + t289 * t337;
t308 = -t253 * t282 + t246;
t305 = -t282 * t342 + t245;
t302 = t230 * t348 + t350;
t301 = t253 * t262 + t265 * t343;
t299 = -t236 + (t237 * t344 + t236) * t238;
t298 = -t244 * t286 + t269 * t327 + t270 * t282;
t251 = -t266 * qJD(3) - t286 * t315;
t243 = t270 * t281 - t280 * t342;
t242 = t270 * t280 + t281 * t342;
t233 = -t253 * t280 + t267 * t281;
t232 = t253 * t281 + t267 * t280;
t220 = 0.1e1 / t222;
t213 = 0.1e1 / t215;
t212 = t238 * t360;
t209 = t301 * t238;
t206 = t299 * t259;
t204 = (-t236 * t267 + t237 * t338) * t289 + t304 * t212;
t203 = t304 * t209 + t236 * t253 + t237 * t265;
t201 = t301 * t325 + (t265 * t312 + t225 * t262 + (t226 * t265 + t251 * t254 - t252 * t253) * t263) * t238;
t199 = t325 * t360 + (-t300 * t328 + (t312 * t338 - t246 * t262 + (t252 * t267 + (t226 * t290 - t254 * t330) * t284) * t263) * t289) * t238;
t198 = -0.2e1 * t351 + 0.2e1 * (t210 * t230 * t220 + (-t220 * t354 - t230 * t351) * t234) * t234;
t1 = [t322 * t358 + (-t224 * t262 + t259 * t361) * t238, t199, t201, 0, 0, 0; t254 * t216 * t326 + (((qJD(3) * t268 - t317) * t286 - t307) * t216 + (t202 * t254 - t206 * t224) * t217) * t213 + (t206 * t217 * t326 + (0.2e1 * t206 * t356 - (-t207 * t238 * t344 + t325) * t321 - (t322 * t359 - t207 + (t207 - t303) * t238) * t320 - t299 * t353) * t213) * t259 (t204 * t352 + t216 * t341) * t326 + (-t204 * t353 + (t244 * t289 + t269 * t328) * t216 + (t204 * t324 + t217 * t341) * t202 - (-t199 * t254 - t212 * t226 + (-t289 * t330 - t290 * t328) * t284 + (-t212 * t266 - t267 * t289) * t207) * t320 - (t267 * t328 - t199 * t266 - t212 * t252 - t246 * t289 + (t212 * t254 - t289 * t338) * t207) * t321) * t213 (t203 * t352 + t216 * t258) * t326 + (t203 * t202 * t324 - t223 * t216 + (t258 * t202 - t203 * t224 - (-t201 * t254 - t209 * t226 + t251 + (-t209 * t266 + t253) * t207) * t346 - (-t201 * t266 - t209 * t252 + t225 + (t209 * t254 - t265) * t207) * t347) * t217) * t213, 0, 0, 0; (-t229 * t232 + t233 * t349) * t323 + ((t308 * t280 + t310 * t281) * t229 + t233 * t313 + (-t232 * t211 - (-t310 * t280 + t308 * t281) * t234 - t233 * t210) * t230) * t220 (-t229 * t242 + t243 * t349) * t323 + (t243 * t313 + t305 * t229 * t280 + t298 * t350 + (-t305 * t234 * t281 - t243 * t210 - t242 * t211 + t298 * t348) * t230) * t220, t302 * t259 * t323 + (-t302 * t224 + ((t229 * t282 + t313) * t280 + (-t210 * t280 + (-t234 * t282 + t211) * t281) * t230) * t259) * t220, 0, t198, t198;];
JaD_rot  = t1;
