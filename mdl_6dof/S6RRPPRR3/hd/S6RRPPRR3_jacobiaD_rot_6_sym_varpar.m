% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:29:34
% EndTime: 2019-02-26 21:29:36
% DurationCPUTime: 1.82s
% Computational Cost: add. (12327->154), mult. (24081->307), div. (726->12), fcn. (31139->15), ass. (0->132)
t326 = sin(pkin(11));
t327 = cos(pkin(11));
t330 = sin(qJ(2));
t332 = cos(qJ(2));
t313 = t330 * t326 - t332 * t327;
t328 = cos(pkin(6));
t310 = t313 * t328;
t305 = qJD(2) * t310;
t354 = t332 * t326 + t330 * t327;
t312 = t354 * qJD(2);
t333 = cos(qJ(1));
t399 = sin(pkin(6));
t363 = t333 * t399;
t311 = t354 * t328;
t400 = sin(qJ(1));
t366 = t400 * t311;
t405 = -t400 * t312 - qJD(1) * t366 + (-qJD(1) * t313 - t305) * t333 - qJD(5) * t363;
t325 = pkin(12) + qJ(5);
t323 = sin(t325);
t324 = cos(t325);
t347 = -t333 * t311 + t313 * t400;
t283 = -t323 * t347 + t324 * t363;
t364 = t332 * t399;
t365 = t330 * t399;
t341 = -t326 * t364 - t327 * t365;
t299 = -t323 * t341 - t328 * t324;
t271 = atan2(-t283, t299);
t266 = sin(t271);
t267 = cos(t271);
t249 = -t266 * t283 + t267 * t299;
t247 = 0.1e1 / t249 ^ 2;
t348 = -t333 * t313 - t366;
t359 = t400 * t399;
t344 = -t323 * t348 + t324 * t359;
t282 = t344 ^ 2;
t245 = t282 * t247 + 0.1e1;
t289 = t323 * t359 + t324 * t348;
t340 = qJD(1) * t347 + t305 * t400 - t333 * t312;
t357 = qJD(1) * t363;
t253 = qJD(5) * t289 + t323 * t340 - t324 * t357;
t392 = t253 * t247;
t281 = t283 ^ 2;
t297 = 0.1e1 / t299 ^ 2;
t270 = t281 * t297 + 0.1e1;
t268 = 0.1e1 / t270;
t353 = qJD(1) * t359;
t377 = qJD(5) * t324;
t255 = t323 * t405 - t324 * t353 - t347 * t377;
t300 = t328 * t323 - t324 * t341;
t308 = -t326 * t365 + t327 * t364;
t304 = t308 * qJD(2);
t279 = qJD(5) * t300 + t304 * t323;
t296 = 0.1e1 / t299;
t386 = t283 * t297;
t352 = -t255 * t296 + t279 * t386;
t237 = t352 * t268;
t355 = -t266 * t299 - t267 * t283;
t232 = t237 * t355 - t266 * t255 + t267 * t279;
t246 = 0.1e1 / t249;
t248 = t246 * t247;
t397 = t232 * t248;
t375 = 0.2e1 * (-t282 * t397 - t344 * t392) / t245 ^ 2;
t404 = t279 * t297;
t291 = -t333 * t310 - t354 * t400;
t349 = -t291 * t296 + t308 * t386;
t403 = t323 * t349;
t256 = t323 * (qJD(5) * t347 + t353) + t405 * t324;
t331 = cos(qJ(6));
t294 = t310 * t400 - t333 * t354;
t329 = sin(qJ(6));
t384 = t294 * t329;
t265 = t289 * t331 - t384;
t259 = 0.1e1 / t265;
t260 = 0.1e1 / t265 ^ 2;
t402 = -0.2e1 * t283;
t401 = -0.2e1 * t344;
t254 = qJD(5) * t344 + t323 * t357 + t324 * t340;
t306 = t328 * t312;
t346 = t313 * qJD(2);
t275 = qJD(1) * t291 - t306 * t400 - t333 * t346;
t241 = qJD(6) * t265 + t254 * t329 - t275 * t331;
t383 = t294 * t331;
t264 = t289 * t329 + t383;
t258 = t264 ^ 2;
t252 = t258 * t260 + 0.1e1;
t391 = t260 * t264;
t376 = qJD(6) * t264;
t242 = t254 * t331 + t275 * t329 - t376;
t394 = t242 * t259 * t260;
t396 = (t241 * t391 - t258 * t394) / t252 ^ 2;
t388 = t296 * t404;
t395 = (t255 * t386 - t281 * t388) / t270 ^ 2;
t393 = t247 * t344;
t390 = t266 * t344;
t389 = t267 * t344;
t387 = t283 * t296;
t385 = t294 * t323;
t382 = t329 * t259;
t380 = t331 * t264;
t374 = -0.2e1 * t396;
t373 = 0.2e1 * t396;
t372 = -0.2e1 * t395;
t371 = t248 * t401;
t370 = t296 * t395;
t369 = t247 * t390;
t368 = t247 * t389;
t367 = t264 * t394;
t362 = 0.2e1 * t367;
t361 = t388 * t402;
t285 = -t323 * t363 - t324 * t347;
t356 = qJD(6) * t294 * t324 - t340;
t263 = -t285 * t331 + t291 * t329;
t262 = -t285 * t329 - t291 * t331;
t351 = t380 * t260 - t382;
t350 = -t285 * t296 + t300 * t386;
t345 = -t266 + (t267 * t387 + t266) * t268;
t342 = -qJD(5) * t385 + qJD(6) * t348 - t275 * t324;
t303 = t341 * qJD(2);
t280 = -qJD(5) * t299 + t304 * t324;
t277 = qJD(1) * t294 - t333 * t306 + t346 * t400;
t273 = t324 * t383 + t329 * t348;
t272 = t324 * t384 - t331 * t348;
t250 = 0.1e1 / t252;
t243 = 0.1e1 / t245;
t240 = t268 * t403;
t239 = t350 * t268;
t236 = t345 * t344;
t234 = (-t266 * t291 + t267 * t308) * t323 + t355 * t240;
t233 = t239 * t355 - t266 * t285 + t267 * t300;
t230 = t350 * t372 + (t300 * t361 - t256 * t296 + (t255 * t300 + t279 * t285 + t280 * t283) * t297) * t268;
t229 = t372 * t403 + (t349 * t377 + (t308 * t361 - t277 * t296 + (t255 * t308 + t279 * t291 + t283 * t303) * t297) * t323) * t268;
t1 = [t370 * t401 + (-t253 * t296 - t344 * t404) * t268, t229, 0, 0, t230, 0; t283 * t246 * t375 + (-t255 * t246 + (t232 * t283 + t236 * t253) * t247) * t243 - (-t236 * t247 * t375 + (-0.2e1 * t236 * t397 + (-t237 * t268 * t387 + t372) * t369 + (t370 * t402 - t237 + (t237 - t352) * t268) * t368 - t345 * t392) * t243) * t344 (-t234 * t393 - t246 * t385) * t375 + (-t234 * t392 + (-t275 * t323 + t294 * t377) * t246 + (t234 * t371 - t247 * t385) * t232 + (t308 * t377 - t229 * t283 - t240 * t255 + t303 * t323 + (-t240 * t299 - t291 * t323) * t237) * t368 + (-t291 * t377 - t229 * t299 - t240 * t279 - t277 * t323 + (t240 * t283 - t308 * t323) * t237) * t369) * t243, 0, 0 (-t233 * t393 - t246 * t289) * t375 + (t233 * t232 * t371 + t254 * t246 + (-t289 * t232 - t233 * t253 + (-t230 * t283 - t239 * t255 + t280 + (-t239 * t299 - t285) * t237) * t389 + (-t230 * t299 - t239 * t279 - t256 + (t239 * t283 - t300) * t237) * t390) * t247) * t243, 0; (-t259 * t262 + t263 * t391) * t373 + ((t263 * qJD(6) - t256 * t329 - t277 * t331) * t259 + t263 * t362 + (-t262 * t242 - (-t262 * qJD(6) - t256 * t331 + t277 * t329) * t264 - t263 * t241) * t260) * t250 (-t259 * t272 + t273 * t391) * t373 + (t273 * t362 + t356 * t259 * t331 + t342 * t382 + (t264 * t329 * t356 - t273 * t241 - t272 * t242 - t342 * t380) * t260) * t250, 0, 0, -t351 * t344 * t374 + (t351 * t253 - ((-qJD(6) * t259 - 0.2e1 * t367) * t331 + (t241 * t331 + (t242 - t376) * t329) * t260) * t344) * t250, t374 + 0.2e1 * (t241 * t260 * t250 + (-t250 * t394 - t260 * t396) * t264) * t264;];
JaD_rot  = t1;
