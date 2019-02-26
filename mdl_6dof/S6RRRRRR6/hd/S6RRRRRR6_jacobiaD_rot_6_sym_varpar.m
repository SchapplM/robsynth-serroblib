% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRR6_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR6_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:50:05
% EndTime: 2019-02-26 22:50:06
% DurationCPUTime: 1.66s
% Computational Cost: add. (13010->155), mult. (18761->302), div. (1007->12), fcn. (23678->13), ass. (0->141)
t329 = qJ(3) + qJ(4);
t323 = sin(t329);
t331 = cos(pkin(6));
t333 = cos(qJ(2));
t407 = sin(qJ(1));
t366 = t407 * t333;
t332 = sin(qJ(2));
t334 = cos(qJ(1));
t384 = t334 * t332;
t345 = -t331 * t384 - t366;
t325 = cos(t329);
t330 = sin(pkin(6));
t385 = t330 * t334;
t371 = t325 * t385;
t295 = -t323 * t345 + t371;
t387 = t330 * t332;
t306 = t323 * t387 - t325 * t331;
t284 = atan2(-t295, t306);
t279 = sin(t284);
t280 = cos(t284);
t262 = -t279 * t295 + t280 * t306;
t260 = 0.1e1 / t262 ^ 2;
t367 = t407 * t332;
t355 = t331 * t367;
t383 = t334 * t333;
t311 = -t355 + t383;
t368 = t330 * t407;
t300 = t311 * t323 - t325 * t368;
t294 = t300 ^ 2;
t258 = t260 * t294 + 0.1e1;
t344 = -t331 * t366 - t384;
t290 = t345 * qJD(1) + t344 * qJD(2);
t327 = qJD(3) + qJD(4);
t351 = t327 * t368 + t290;
t365 = qJD(1) * t385;
t388 = t325 * t327;
t266 = t311 * t388 + t351 * t323 - t325 * t365;
t400 = t260 * t300;
t293 = t295 ^ 2;
t304 = 0.1e1 / t306 ^ 2;
t283 = t293 * t304 + 0.1e1;
t281 = 0.1e1 / t283;
t364 = qJD(2) * t407;
t292 = -qJD(1) * t355 - t332 * t364 + (qJD(2) * t331 + qJD(1)) * t383;
t317 = t323 * t385;
t363 = t407 * qJD(1);
t354 = t330 * t363;
t268 = t292 * t323 - t327 * t317 - t325 * t354 - t345 * t388;
t381 = qJD(2) * t333;
t347 = t327 * t331 + t330 * t381;
t370 = t327 * t387;
t287 = t347 * t323 + t325 * t370;
t303 = 0.1e1 / t306;
t392 = t295 * t304;
t350 = -t268 * t303 + t287 * t392;
t250 = t350 * t281;
t352 = -t279 * t306 - t280 * t295;
t245 = t352 * t250 - t279 * t268 + t280 * t287;
t259 = 0.1e1 / t262;
t261 = t259 * t260;
t405 = t245 * t261;
t380 = 0.2e1 * (t266 * t400 - t294 * t405) / t258 ^ 2;
t411 = t287 * t304;
t369 = t331 * t383;
t308 = -t367 + t369;
t386 = t330 * t333;
t346 = -t303 * t308 + t386 * t392;
t410 = t323 * t346;
t269 = (t327 * t345 + t354) * t323 + t292 * t325 - t327 * t371;
t301 = t311 * t325 + t323 * t368;
t328 = qJ(5) + qJ(6);
t322 = sin(t328);
t324 = cos(t328);
t278 = t301 * t324 - t322 * t344;
t272 = 0.1e1 / t278;
t273 = 0.1e1 / t278 ^ 2;
t409 = -0.2e1 * t295;
t408 = 0.2e1 * t300;
t289 = -qJD(1) * t369 - t334 * t381 + (t331 * t364 + t363) * t332;
t326 = qJD(5) + qJD(6);
t358 = t301 * t326 + t289;
t267 = t351 * t325 + (-t311 * t327 + t365) * t323;
t389 = t344 * t326;
t360 = t267 - t389;
t254 = t360 * t322 + t358 * t324;
t277 = t301 * t322 + t324 * t344;
t271 = t277 ^ 2;
t265 = t271 * t273 + 0.1e1;
t398 = t273 * t277;
t255 = -t358 * t322 + t360 * t324;
t402 = t255 * t272 * t273;
t404 = (t254 * t398 - t271 * t402) / t265 ^ 2;
t394 = t303 * t411;
t403 = (t268 * t392 - t293 * t394) / t283 ^ 2;
t401 = t260 * t266;
t399 = t272 * t322;
t397 = t277 * t324;
t396 = t279 * t300;
t395 = t280 * t300;
t393 = t295 * t303;
t391 = t344 * t323;
t390 = t344 * t325;
t382 = qJD(2) * t332;
t379 = -0.2e1 * t404;
t378 = 0.2e1 * t404;
t377 = -0.2e1 * t403;
t376 = t261 * t408;
t375 = t303 * t403;
t374 = t277 * t402;
t373 = t260 * t396;
t372 = t260 * t395;
t362 = 0.2e1 * t374;
t361 = t394 * t409;
t359 = t308 * t326 - t269;
t291 = t344 * qJD(1) + t345 * qJD(2);
t297 = -t325 * t345 - t317;
t357 = -t297 * t326 - t291;
t353 = -t325 * t389 + t290;
t349 = t273 * t397 - t399;
t307 = t323 * t331 + t325 * t387;
t348 = -t297 * t303 + t307 * t392;
t342 = -t279 + (t280 * t393 + t279) * t281;
t341 = t289 * t325 + t311 * t326 - t327 * t391;
t288 = -t323 * t370 + t347 * t325;
t286 = t311 * t322 + t324 * t390;
t285 = -t311 * t324 + t322 * t390;
t276 = -t297 * t324 + t308 * t322;
t275 = -t297 * t322 - t308 * t324;
t263 = 0.1e1 / t265;
t256 = 0.1e1 / t258;
t253 = t281 * t410;
t252 = t348 * t281;
t249 = t342 * t300;
t247 = (-t279 * t308 + t280 * t386) * t323 + t352 * t253;
t246 = t352 * t252 - t279 * t297 + t280 * t307;
t243 = t348 * t377 + (t307 * t361 - t269 * t303 + (t268 * t307 + t287 * t297 + t288 * t295) * t304) * t281;
t242 = t377 * t410 + (t346 * t388 + (t361 * t386 - t291 * t303 + (t287 * t308 + (t268 * t333 - t295 * t382) * t330) * t304) * t323) * t281;
t241 = t379 + 0.2e1 * (t254 * t263 * t273 + (-t263 * t402 - t273 * t404) * t277) * t277;
t240 = t349 * t300 * t379 + (t349 * t266 + ((-t272 * t326 - 0.2e1 * t374) * t324 + (t254 * t324 + (-t277 * t326 + t255) * t322) * t273) * t300) * t263;
t239 = (t246 * t400 - t259 * t301) * t380 + (t246 * t245 * t376 + t267 * t259 + (-t301 * t245 - t246 * t266 - (-t243 * t295 - t252 * t268 + t288 + (-t252 * t306 - t297) * t250) * t395 - (-t243 * t306 - t252 * t287 - t269 + (t252 * t295 - t307) * t250) * t396) * t260) * t256;
t1 = [t375 * t408 + (-t266 * t303 + t300 * t411) * t281, t242, t243, t243, 0, 0; t295 * t259 * t380 + (-t268 * t259 + (t245 * t295 - t249 * t266) * t260) * t256 + (t249 * t260 * t380 + (0.2e1 * t249 * t405 - (-t250 * t281 * t393 + t377) * t373 - (t375 * t409 - t250 + (t250 - t350) * t281) * t372 - t342 * t401) * t256) * t300 (t247 * t400 - t259 * t391) * t380 + (-t247 * t401 + (t289 * t323 + t344 * t388) * t259 + (t247 * t376 - t260 * t391) * t245 - (-t242 * t295 - t253 * t268 + (-t323 * t382 + t333 * t388) * t330 + (-t253 * t306 - t308 * t323) * t250) * t372 - (-t308 * t388 - t242 * t306 - t253 * t287 - t291 * t323 + (t253 * t295 - t323 * t386) * t250) * t373) * t256, t239, t239, 0, 0; (-t272 * t275 + t276 * t398) * t378 + ((t359 * t322 + t357 * t324) * t272 + t276 * t362 + (-t275 * t255 - (-t357 * t322 + t359 * t324) * t277 - t276 * t254) * t273) * t263 (-t272 * t285 + t286 * t398) * t378 + (t286 * t362 - t353 * t272 * t324 + t341 * t399 + (-t353 * t277 * t322 - t286 * t254 - t285 * t255 - t341 * t397) * t273) * t263, t240, t240, t241, t241;];
JaD_rot  = t1;
