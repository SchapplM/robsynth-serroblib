% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR11_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR11_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:36:21
% EndTime: 2019-02-26 22:36:23
% DurationCPUTime: 1.54s
% Computational Cost: add. (5847->150), mult. (14155->298), div. (744->12), fcn. (17850->13), ass. (0->132)
t286 = cos(pkin(6));
t288 = sin(qJ(2));
t362 = sin(qJ(1));
t325 = t362 * t288;
t311 = t286 * t325;
t320 = qJD(2) * t362;
t290 = cos(qJ(2));
t291 = cos(qJ(1));
t340 = t291 * t290;
t285 = sin(pkin(6));
t342 = t285 * t291;
t367 = -qJD(1) * t311 - t288 * t320 + (qJD(2) * t286 + qJD(1)) * t340 - qJD(3) * t342;
t287 = sin(qJ(3));
t289 = cos(qJ(3));
t324 = t362 * t290;
t341 = t291 * t288;
t303 = -t286 * t341 - t324;
t254 = -t287 * t303 + t289 * t342;
t344 = t285 * t288;
t265 = -t286 * t289 + t287 * t344;
t245 = atan2(-t254, t265);
t240 = sin(t245);
t241 = cos(t245);
t221 = -t240 * t254 + t241 * t265;
t219 = 0.1e1 / t221 ^ 2;
t270 = -t311 + t340;
t326 = t285 * t362;
t302 = -t270 * t287 + t289 * t326;
t251 = t302 ^ 2;
t217 = t219 * t251 + 0.1e1;
t301 = -t286 * t324 - t341;
t247 = t303 * qJD(1) + t301 * qJD(2);
t260 = t270 * t289 + t287 * t326;
t323 = qJD(1) * t342;
t225 = t260 * qJD(3) + t247 * t287 - t289 * t323;
t356 = t219 * t302;
t250 = t254 ^ 2;
t263 = 0.1e1 / t265 ^ 2;
t244 = t250 * t263 + 0.1e1;
t242 = 0.1e1 / t244;
t319 = t362 * qJD(1);
t310 = t285 * t319;
t337 = qJD(3) * t289;
t227 = t367 * t287 - t289 * t310 - t303 * t337;
t266 = t286 * t287 + t289 * t344;
t338 = qJD(2) * t290;
t322 = t285 * t338;
t252 = t266 * qJD(3) + t287 * t322;
t262 = 0.1e1 / t265;
t347 = t254 * t263;
t307 = -t227 * t262 + t252 * t347;
t209 = t307 * t242;
t308 = -t240 * t265 - t241 * t254;
t204 = t308 * t209 - t240 * t227 + t241 * t252;
t218 = 0.1e1 / t221;
t220 = t218 * t219;
t360 = t204 * t220;
t336 = 0.2e1 * (-t225 * t356 - t251 * t360) / t217 ^ 2;
t366 = t252 * t263;
t327 = t286 * t340;
t267 = -t325 + t327;
t343 = t285 * t290;
t304 = -t262 * t267 + t343 * t347;
t365 = t287 * t304;
t228 = (qJD(3) * t303 + t310) * t287 + t367 * t289;
t283 = qJ(4) + pkin(12) + qJ(6);
t281 = sin(t283);
t282 = cos(t283);
t237 = t260 * t282 - t281 * t301;
t231 = 0.1e1 / t237;
t232 = 0.1e1 / t237 ^ 2;
t364 = -0.2e1 * t254;
t363 = -0.2e1 * t302;
t349 = t262 * t366;
t359 = (t227 * t347 - t250 * t349) / t244 ^ 2;
t246 = -qJD(1) * t327 - t291 * t338 + (t286 * t320 + t319) * t288;
t284 = qJD(4) + qJD(6);
t314 = t260 * t284 + t246;
t226 = t302 * qJD(3) + t247 * t289 + t287 * t323;
t316 = -t284 * t301 + t226;
t213 = -t314 * t281 + t316 * t282;
t358 = t213 * t231 * t232;
t357 = t219 * t225;
t212 = t316 * t281 + t314 * t282;
t236 = t260 * t281 + t282 * t301;
t230 = t236 ^ 2;
t224 = t230 * t232 + 0.1e1;
t353 = t232 * t236;
t355 = 0.1e1 / t224 ^ 2 * (t212 * t353 - t230 * t358);
t354 = t231 * t281;
t352 = t236 * t282;
t351 = t240 * t302;
t350 = t241 * t302;
t348 = t254 * t262;
t346 = t301 * t287;
t345 = t301 * t289;
t339 = qJD(2) * t288;
t335 = -0.2e1 * t359;
t334 = t220 * t363;
t333 = -0.2e1 * t355;
t332 = 0.2e1 * t355;
t331 = t262 * t359;
t330 = t236 * t358;
t329 = t219 * t351;
t328 = t219 * t350;
t318 = 0.2e1 * t330;
t317 = t349 * t364;
t315 = t267 * t284 - t228;
t248 = t301 * qJD(1) + t303 * qJD(2);
t256 = -t287 * t342 - t289 * t303;
t313 = -t256 * t284 - t248;
t309 = -t284 * t345 + t247;
t306 = t232 * t352 - t354;
t305 = -t256 * t262 + t266 * t347;
t299 = -t240 + (t241 * t348 + t240) * t242;
t298 = -qJD(3) * t346 + t246 * t289 + t270 * t284;
t253 = -t265 * qJD(3) + t289 * t322;
t239 = t270 * t281 + t282 * t345;
t238 = -t270 * t282 + t281 * t345;
t235 = -t256 * t282 + t267 * t281;
t234 = -t256 * t281 - t267 * t282;
t222 = 0.1e1 / t224;
t215 = 0.1e1 / t217;
t214 = t242 * t365;
t211 = t305 * t242;
t208 = t299 * t302;
t206 = (-t240 * t267 + t241 * t343) * t287 + t308 * t214;
t205 = t308 * t211 - t240 * t256 + t241 * t266;
t203 = t305 * t335 + (t266 * t317 - t228 * t262 + (t227 * t266 + t252 * t256 + t253 * t254) * t263) * t242;
t201 = t335 * t365 + (t304 * t337 + (t317 * t343 - t248 * t262 + (t252 * t267 + (t227 * t290 - t254 * t339) * t285) * t263) * t287) * t242;
t200 = t333 + 0.2e1 * (t212 * t232 * t222 + (-t222 * t358 - t232 * t355) * t236) * t236;
t1 = [t331 * t363 + (-t225 * t262 - t302 * t366) * t242, t201, t203, 0, 0, 0; t254 * t218 * t336 + (-t227 * t218 + (t204 * t254 + t208 * t225) * t219) * t215 - (-t208 * t219 * t336 + (-0.2e1 * t208 * t360 + (-t209 * t242 * t348 + t335) * t329 + (t331 * t364 - t209 + (t209 - t307) * t242) * t328 - t299 * t357) * t215) * t302 (-t206 * t356 - t218 * t346) * t336 + (-t206 * t357 + (t246 * t287 + t301 * t337) * t218 + (t206 * t334 - t219 * t346) * t204 + (-t201 * t254 - t214 * t227 + (-t287 * t339 + t290 * t337) * t285 + (-t214 * t265 - t267 * t287) * t209) * t328 + (-t267 * t337 - t201 * t265 - t214 * t252 - t248 * t287 + (t214 * t254 - t287 * t343) * t209) * t329) * t215 (-t205 * t356 - t218 * t260) * t336 + (t205 * t204 * t334 + t226 * t218 + (-t260 * t204 - t205 * t225 + (-t203 * t254 - t211 * t227 + t253 + (-t211 * t265 - t256) * t209) * t350 + (-t203 * t265 - t211 * t252 - t228 + (t211 * t254 - t266) * t209) * t351) * t219) * t215, 0, 0, 0; (-t231 * t234 + t235 * t353) * t332 + ((t315 * t281 + t313 * t282) * t231 + t235 * t318 + (-t234 * t213 - (-t313 * t281 + t315 * t282) * t236 - t235 * t212) * t232) * t222 (-t231 * t238 + t239 * t353) * t332 + (t239 * t318 - t309 * t231 * t282 + t298 * t354 + (-t309 * t236 * t281 - t239 * t212 - t238 * t213 - t298 * t352) * t232) * t222, -t306 * t302 * t333 + (t306 * t225 - ((-t231 * t284 - 0.2e1 * t330) * t282 + (t212 * t282 + (-t236 * t284 + t213) * t281) * t232) * t302) * t222, t200, 0, t200;];
JaD_rot  = t1;
