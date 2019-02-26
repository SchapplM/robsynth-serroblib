% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR8_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR8_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:19:58
% EndTime: 2019-02-26 22:19:59
% DurationCPUTime: 1.50s
% Computational Cost: add. (9346->151), mult. (14155->297), div. (744->12), fcn. (17850->13), ass. (0->135)
t298 = cos(pkin(6));
t299 = sin(qJ(2));
t374 = sin(qJ(1));
t335 = t374 * t299;
t321 = t298 * t335;
t330 = qJD(2) * t374;
t300 = cos(qJ(2));
t301 = cos(qJ(1));
t350 = t301 * t300;
t297 = sin(pkin(6));
t352 = t297 * t301;
t379 = -qJD(1) * t321 - t299 * t330 + (qJD(2) * t298 + qJD(1)) * t350 - qJD(3) * t352;
t295 = qJ(3) + pkin(12);
t290 = sin(t295);
t291 = cos(t295);
t334 = t374 * t300;
t351 = t301 * t299;
t313 = -t298 * t351 - t334;
t263 = -t290 * t313 + t291 * t352;
t354 = t297 * t299;
t274 = t290 * t354 - t291 * t298;
t252 = atan2(-t263, t274);
t247 = sin(t252);
t248 = cos(t252);
t230 = -t247 * t263 + t248 * t274;
t228 = 0.1e1 / t230 ^ 2;
t279 = -t321 + t350;
t336 = t297 * t374;
t312 = -t279 * t290 + t291 * t336;
t262 = t312 ^ 2;
t224 = t228 * t262 + 0.1e1;
t311 = -t298 * t334 - t351;
t256 = t313 * qJD(1) + t311 * qJD(2);
t269 = t279 * t291 + t290 * t336;
t333 = qJD(1) * t352;
t234 = t269 * qJD(3) + t256 * t290 - t291 * t333;
t367 = t234 * t228;
t261 = t263 ^ 2;
t271 = 0.1e1 / t274 ^ 2;
t251 = t261 * t271 + 0.1e1;
t249 = 0.1e1 / t251;
t329 = t374 * qJD(1);
t320 = t297 * t329;
t347 = qJD(3) * t291;
t236 = t379 * t290 - t291 * t320 - t313 * t347;
t275 = t290 * t298 + t291 * t354;
t348 = qJD(2) * t300;
t332 = t297 * t348;
t259 = t275 * qJD(3) + t290 * t332;
t270 = 0.1e1 / t274;
t359 = t263 * t271;
t317 = -t236 * t270 + t259 * t359;
t218 = t317 * t249;
t318 = -t247 * t274 - t248 * t263;
t213 = t318 * t218 - t236 * t247 + t248 * t259;
t227 = 0.1e1 / t230;
t229 = t227 * t228;
t372 = t213 * t229;
t346 = 0.2e1 * (-t262 * t372 - t312 * t367) / t224 ^ 2;
t378 = t259 * t271;
t337 = t298 * t350;
t276 = -t335 + t337;
t353 = t297 * t300;
t314 = -t270 * t276 + t353 * t359;
t377 = t290 * t314;
t237 = (qJD(3) * t313 + t320) * t290 + t379 * t291;
t296 = qJ(5) + qJ(6);
t293 = cos(t296);
t292 = sin(t296);
t357 = t311 * t292;
t246 = t269 * t293 - t357;
t240 = 0.1e1 / t246;
t241 = 0.1e1 / t246 ^ 2;
t376 = -0.2e1 * t263;
t375 = -0.2e1 * t312;
t255 = -qJD(1) * t337 - t301 * t348 + (t298 * t330 + t329) * t299;
t294 = qJD(5) + qJD(6);
t324 = t269 * t294 + t255;
t235 = t312 * qJD(3) + t256 * t291 + t290 * t333;
t355 = t311 * t294;
t326 = t235 - t355;
t225 = t326 * t292 + t324 * t293;
t356 = t311 * t293;
t245 = t269 * t292 + t356;
t239 = t245 ^ 2;
t233 = t239 * t241 + 0.1e1;
t365 = t241 * t245;
t226 = -t324 * t292 + t326 * t293;
t369 = t226 * t240 * t241;
t371 = (t225 * t365 - t239 * t369) / t233 ^ 2;
t361 = t270 * t378;
t370 = (t236 * t359 - t261 * t361) / t251 ^ 2;
t368 = t228 * t312;
t366 = t240 * t292;
t364 = t245 * t293;
t363 = t247 * t312;
t362 = t248 * t312;
t360 = t263 * t270;
t358 = t311 * t290;
t349 = qJD(2) * t299;
t345 = -0.2e1 * t371;
t344 = 0.2e1 * t371;
t343 = -0.2e1 * t370;
t342 = t229 * t375;
t341 = t270 * t370;
t340 = t245 * t369;
t339 = t228 * t363;
t338 = t228 * t362;
t328 = 0.2e1 * t340;
t327 = t361 * t376;
t325 = t276 * t294 - t237;
t257 = t311 * qJD(1) + t313 * qJD(2);
t265 = -t290 * t352 - t291 * t313;
t323 = -t265 * t294 - t257;
t319 = -t291 * t355 + t256;
t316 = t241 * t364 - t366;
t315 = -t265 * t270 + t275 * t359;
t309 = -t247 + (t248 * t360 + t247) * t249;
t308 = -qJD(3) * t358 + t255 * t291 + t279 * t294;
t260 = -t274 * qJD(3) + t291 * t332;
t254 = t279 * t292 + t291 * t356;
t253 = -t279 * t293 + t291 * t357;
t244 = -t265 * t293 + t276 * t292;
t243 = -t265 * t292 - t276 * t293;
t231 = 0.1e1 / t233;
t222 = 0.1e1 / t224;
t221 = t249 * t377;
t219 = t315 * t249;
t217 = t309 * t312;
t215 = (-t247 * t276 + t248 * t353) * t290 + t318 * t221;
t214 = t318 * t219 - t247 * t265 + t248 * t275;
t212 = t315 * t343 + (t275 * t327 - t237 * t270 + (t236 * t275 + t259 * t265 + t260 * t263) * t271) * t249;
t210 = t343 * t377 + (t314 * t347 + (t327 * t353 - t257 * t270 + (t259 * t276 + (t236 * t300 - t263 * t349) * t297) * t271) * t290) * t249;
t209 = t345 + 0.2e1 * (t225 * t231 * t241 + (-t231 * t369 - t241 * t371) * t245) * t245;
t1 = [t341 * t375 + (-t234 * t270 - t312 * t378) * t249, t210, t212, 0, 0, 0; t263 * t227 * t346 + (-t236 * t227 + (t213 * t263 + t217 * t234) * t228) * t222 - (-t217 * t228 * t346 + (-0.2e1 * t217 * t372 + (-t218 * t249 * t360 + t343) * t339 + (t341 * t376 - t218 + (t218 - t317) * t249) * t338 - t309 * t367) * t222) * t312 (-t215 * t368 - t227 * t358) * t346 + (-t215 * t367 + (t255 * t290 + t311 * t347) * t227 + (t215 * t342 - t228 * t358) * t213 + (-t210 * t263 - t221 * t236 + (-t290 * t349 + t300 * t347) * t297 + (-t221 * t274 - t276 * t290) * t218) * t338 + (-t276 * t347 - t210 * t274 - t221 * t259 - t257 * t290 + (t221 * t263 - t290 * t353) * t218) * t339) * t222 (-t214 * t368 - t227 * t269) * t346 + (t214 * t213 * t342 + t235 * t227 + (-t269 * t213 - t214 * t234 + (-t212 * t263 - t219 * t236 + t260 + (-t219 * t274 - t265) * t218) * t362 + (-t212 * t274 - t219 * t259 - t237 + (t219 * t263 - t275) * t218) * t363) * t228) * t222, 0, 0, 0; (-t240 * t243 + t244 * t365) * t344 + ((t325 * t292 + t323 * t293) * t240 + t244 * t328 + (-t243 * t226 - (-t323 * t292 + t325 * t293) * t245 - t244 * t225) * t241) * t231 (-t240 * t253 + t254 * t365) * t344 + (t254 * t328 - t319 * t240 * t293 + t308 * t366 + (-t319 * t245 * t292 - t254 * t225 - t253 * t226 - t308 * t364) * t241) * t231, -t316 * t312 * t345 + (t316 * t234 - ((-t240 * t294 - 0.2e1 * t340) * t293 + (t225 * t293 + (-t245 * t294 + t226) * t292) * t241) * t312) * t231, 0, t209, t209;];
JaD_rot  = t1;
