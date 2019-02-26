% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR9_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR9_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobiaD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:53:24
% EndTime: 2019-02-26 20:53:26
% DurationCPUTime: 1.60s
% Computational Cost: add. (8530->127), mult. (25589->258), div. (448->12), fcn. (32987->17), ass. (0->119)
t327 = sin(pkin(7));
t325 = sin(pkin(13));
t329 = cos(pkin(13));
t334 = sin(qJ(3));
t337 = cos(qJ(3));
t358 = t325 * t337 + t329 * t334;
t351 = qJD(3) * t358;
t300 = t327 * t351;
t331 = cos(pkin(7));
t302 = t331 * t351;
t316 = t325 * t334 - t329 * t337;
t314 = t316 * qJD(3);
t326 = sin(pkin(12));
t328 = sin(pkin(6));
t330 = cos(pkin(12));
t332 = cos(pkin(6));
t284 = t332 * t300 + (t302 * t330 - t314 * t326) * t328;
t304 = t316 * t331;
t352 = t316 * t327;
t290 = (t304 * t330 + t326 * t358) * t328 + t332 * t352;
t287 = 0.1e1 / t290 ^ 2;
t380 = t284 * t287;
t371 = qJD(3) * t337;
t372 = qJD(3) * t334;
t391 = t325 * t372 - t329 * t371;
t338 = cos(qJ(1));
t375 = t332 * t338;
t335 = sin(qJ(1));
t378 = t326 * t335;
t310 = -t330 * t375 + t378;
t374 = t335 * t330;
t311 = t326 * t375 + t374;
t350 = t328 * t352;
t353 = t304 * t310 - t311 * t358 + t338 * t350;
t265 = atan2(t353, t290);
t260 = sin(t265);
t261 = cos(t265);
t249 = t260 * t353 + t261 * t290;
t246 = 0.1e1 / t249;
t354 = t326 * t338 + t332 * t374;
t377 = t328 * t335;
t295 = t327 * t354 + t331 * t377;
t333 = sin(qJ(5));
t336 = cos(qJ(5));
t303 = t358 * t327;
t305 = t358 * t331;
t366 = t332 * t378;
t313 = t330 * t338 - t366;
t347 = t303 * t377 - t305 * t354 - t313 * t316;
t273 = t295 * t333 + t336 * t347;
t267 = 0.1e1 / t273;
t286 = 0.1e1 / t290;
t247 = 0.1e1 / t249 ^ 2;
t268 = 0.1e1 / t273 ^ 2;
t281 = t304 * t354 - t313 * t358 - t335 * t350;
t275 = t281 ^ 2;
t245 = t247 * t275 + 0.1e1;
t306 = t310 * qJD(1);
t307 = t311 * qJD(1);
t346 = qJD(1) * t350;
t255 = -t300 * t377 + t302 * t354 - t306 * t304 + t307 * t358 + t313 * t314 - t338 * t346;
t386 = t247 * t281;
t274 = t353 ^ 2;
t264 = t274 * t287 + 0.1e1;
t262 = 0.1e1 / t264;
t308 = t354 * qJD(1);
t373 = qJD(1) * t338;
t309 = -qJD(1) * t366 + t330 * t373;
t376 = t328 * t338;
t258 = t300 * t376 + t310 * t302 + t308 * t304 - t309 * t358 + t311 * t314 - t335 * t346;
t357 = t258 * t286 - t353 * t380;
t240 = t357 * t262;
t359 = -t260 * t290 + t261 * t353;
t236 = t240 * t359 + t260 * t258 + t261 * t284;
t389 = t236 * t246 * t247;
t390 = (t255 * t386 - t275 * t389) / t245 ^ 2;
t379 = t286 * t380;
t388 = (t258 * t287 * t353 - t274 * t379) / t264 ^ 2;
t243 = 0.1e1 / t245;
t387 = t243 * t247;
t299 = t391 * t327;
t301 = t391 * t331;
t315 = -t325 * t371 - t329 * t372;
t256 = t354 * t301 + t306 * t305 + t307 * t316 + t313 * t315 + (-t299 * t335 + t303 * t373) * t328;
t365 = qJD(1) * t328 * t331;
t291 = -t306 * t327 + t338 * t365;
t272 = -t295 * t336 + t333 * t347;
t370 = qJD(5) * t272;
t251 = t256 * t336 + t291 * t333 - t370;
t385 = t251 * t267 * t268;
t250 = qJD(5) * t273 + t256 * t333 - t291 * t336;
t266 = t272 ^ 2;
t254 = t266 * t268 + 0.1e1;
t383 = t268 * t272;
t384 = 0.1e1 / t254 ^ 2 * (t250 * t383 - t266 * t385);
t382 = t353 * t286;
t289 = t332 * t303 + (t305 * t330 - t316 * t326) * t328;
t381 = t353 * t289;
t369 = -0.2e1 * t390;
t368 = 0.2e1 * t388;
t367 = 0.2e1 * t384;
t362 = -0.2e1 * t281 * t389;
t361 = -0.2e1 * t286 * t388;
t360 = 0.2e1 * t272 * t385;
t294 = -t310 * t327 + t331 * t376;
t348 = t303 * t376 + t305 * t310 + t311 * t316;
t271 = t294 * t333 + t336 * t348;
t270 = -t294 * t336 + t333 * t348;
t356 = -t267 * t333 + t336 * t383;
t355 = -t286 * t348 + t287 * t381;
t349 = t260 + (t261 * t382 - t260) * t262;
t345 = -t310 * t301 + t308 * t305 + t309 * t316 - t311 * t315 + (-qJD(1) * t303 * t335 - t299 * t338) * t328;
t292 = -t308 * t327 - t335 * t365;
t285 = -t332 * t299 + (-t301 * t330 + t315 * t326) * t328;
t252 = 0.1e1 / t254;
t241 = t355 * t262;
t237 = -t241 * t359 + t260 * t348 + t261 * t289;
t235 = t355 * t368 + (0.2e1 * t379 * t381 + t345 * t286 + (-t258 * t289 - t284 * t348 - t285 * t353) * t287) * t262;
t1 = [t281 * t361 + (t255 * t286 - t281 * t380) * t262, 0, t235, 0, 0, 0; (-t236 * t387 + t246 * t369) * t353 + (t258 * t246 + (t349 * t255 + ((-t240 * t262 * t382 + t368) * t260 + (t353 * t361 + t240 + (-t240 + t357) * t262) * t261) * t281) * t386) * t243 + (t243 * t362 + t255 * t387 + t369 * t386) * t349 * t281, 0, 0.2e1 * (-t237 * t386 - t246 * t347) * t390 + (t256 * t246 + t237 * t362 + (-t347 * t236 + t237 * t255 + ((t235 * t353 - t241 * t258 + t285 + (t241 * t290 + t348) * t240) * t261 + (-t235 * t290 + t241 * t284 + t345 + (t241 * t353 - t289) * t240) * t260) * t281) * t247) * t243, 0, 0, 0; (-t267 * t270 + t271 * t383) * t367 + ((qJD(5) * t271 - t292 * t336 + t333 * t345) * t267 + t271 * t360 + (-t270 * t251 - (-qJD(5) * t270 + t292 * t333 + t336 * t345) * t272 - t271 * t250) * t268) * t252, 0, t356 * t281 * t367 + (-t356 * t255 + ((qJD(5) * t267 + t360) * t336 + (-t250 * t336 + (-t251 + t370) * t333) * t268) * t281) * t252, 0, -0.2e1 * t384 + 0.2e1 * (t250 * t268 * t252 + (-t252 * t385 - t268 * t384) * t272) * t272, 0;];
JaD_rot  = t1;
