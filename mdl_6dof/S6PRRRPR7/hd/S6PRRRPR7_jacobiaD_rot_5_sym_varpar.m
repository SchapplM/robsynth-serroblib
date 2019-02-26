% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRPR7_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR7_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_jacobiaD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:05
% EndTime: 2019-02-26 20:14:07
% DurationCPUTime: 2.58s
% Computational Cost: add. (13649->192), mult. (41070->372), div. (792->12), fcn. (52744->17), ass. (0->155)
t356 = sin(qJ(3));
t359 = cos(qJ(3));
t351 = sin(pkin(7));
t354 = cos(pkin(7));
t357 = sin(qJ(2));
t360 = cos(qJ(2));
t430 = cos(pkin(12));
t431 = cos(pkin(6));
t387 = t431 * t430;
t429 = sin(pkin(12));
t372 = t429 * t357 - t360 * t387;
t352 = sin(pkin(6));
t396 = t352 * t430;
t367 = -t351 * t396 - t354 * t372;
t371 = -t357 * t387 - t360 * t429;
t322 = t356 * t371 + t359 * t367;
t339 = t372 * qJD(2);
t340 = t371 * qJD(2);
t410 = t354 * t356;
t301 = qJD(3) * t322 - t339 * t359 + t340 * t410;
t323 = t356 * t367 - t359 * t371;
t355 = sin(qJ(4));
t358 = cos(qJ(4));
t368 = t351 * t372 - t354 * t396;
t311 = t323 * t358 + t355 * t368;
t412 = t351 * t358;
t276 = qJD(4) * t311 + t301 * t355 + t340 * t412;
t309 = t323 * t355 - t358 * t368;
t307 = t309 ^ 2;
t406 = t357 * t359;
t407 = t356 * t360;
t377 = t354 * t407 + t406;
t397 = t351 * t431;
t335 = t352 * t377 + t356 * t397;
t411 = t351 * t360;
t343 = -t352 * t411 + t354 * t431;
t326 = t335 * t355 - t343 * t358;
t320 = 0.1e1 / t326 ^ 2;
t292 = t307 * t320 + 0.1e1;
t290 = 0.1e1 / t292;
t405 = t359 * t360;
t408 = t356 * t357;
t376 = -t354 * t408 + t405;
t379 = t354 * t405 - t408;
t388 = qJD(3) * t397;
t318 = t359 * t388 + (qJD(2) * t376 + qJD(3) * t379) * t352;
t327 = t335 * t358 + t343 * t355;
t413 = t351 * t357;
t398 = t352 * t413;
t390 = qJD(2) * t398;
t294 = qJD(4) * t327 + t318 * t355 - t358 * t390;
t319 = 0.1e1 / t326;
t419 = t309 * t320;
t259 = (-t276 * t319 + t294 * t419) * t290;
t293 = atan2(-t309, t326);
t288 = sin(t293);
t289 = cos(t293);
t385 = -t288 * t326 - t289 * t309;
t254 = t259 * t385 - t276 * t288 + t289 * t294;
t270 = -t288 * t309 + t289 * t326;
t267 = 0.1e1 / t270;
t268 = 0.1e1 / t270 ^ 2;
t434 = t254 * t267 * t268;
t386 = t431 * t429;
t370 = t357 * t386 - t360 * t430;
t369 = t357 * t430 + t360 * t386;
t395 = t352 * t429;
t389 = t351 * t395;
t373 = -t354 * t369 + t389;
t325 = t356 * t373 - t359 * t370;
t374 = t351 * t369 + t354 * t395;
t312 = t325 * t355 - t358 * t374;
t391 = 0.2e1 * t312 * t434;
t334 = t352 * t379 + t359 * t397;
t381 = -t319 * t322 + t334 * t419;
t433 = t355 * t381;
t420 = t294 * t319 * t320;
t432 = -0.2e1 * (t276 * t419 - t307 * t420) / t292 ^ 2;
t313 = t325 * t358 + t355 * t374;
t409 = t354 * t359;
t415 = t370 * t356;
t324 = -t359 * t389 + t369 * t409 - t415;
t350 = sin(pkin(13));
t353 = cos(pkin(13));
t287 = t313 * t353 + t324 * t350;
t283 = 0.1e1 / t287;
t284 = 0.1e1 / t287 ^ 2;
t428 = t268 * t312;
t341 = t369 * qJD(2);
t342 = t370 * qJD(2);
t303 = t342 * t410 - t341 * t359 + (t359 * t373 + t415) * qJD(3);
t414 = t351 * t355;
t279 = -qJD(4) * t312 + t303 * t358 - t342 * t414;
t302 = qJD(3) * t325 - t341 * t356 - t342 * t409;
t272 = t279 * t353 + t302 * t350;
t427 = t272 * t283 * t284;
t278 = qJD(4) * t313 + t303 * t355 + t342 * t412;
t426 = t278 * t268;
t425 = t283 * t350;
t286 = t313 * t350 - t324 * t353;
t424 = t284 * t286;
t423 = t286 * t353;
t422 = t288 * t312;
t421 = t289 * t312;
t418 = t324 * t355;
t417 = t324 * t358;
t404 = qJD(4) * t355;
t403 = qJD(4) * t358;
t308 = t312 ^ 2;
t266 = t268 * t308 + 0.1e1;
t402 = 0.2e1 * (-t308 * t434 + t312 * t426) / t266 ^ 2;
t271 = t279 * t350 - t302 * t353;
t282 = t286 ^ 2;
t275 = t282 * t284 + 0.1e1;
t401 = 0.2e1 * (t271 * t424 - t282 * t427) / t275 ^ 2;
t399 = t286 * t427;
t393 = 0.2e1 * t399;
t392 = -0.2e1 * t309 * t420;
t383 = -t311 * t319 + t327 * t419;
t328 = -t359 * t372 + t371 * t410;
t314 = t328 * t355 + t371 * t412;
t338 = t376 * t352;
t331 = t338 * t355 - t358 * t398;
t382 = -t314 * t319 + t331 * t419;
t330 = -t359 * t369 + t370 * t410;
t380 = -t330 * t355 - t370 * t412;
t316 = t330 * t358 - t370 * t414;
t329 = -t356 * t369 - t370 * t409;
t378 = -t354 * t406 - t407;
t375 = -t302 * t358 + t324 * t404;
t317 = -t356 * t388 + (qJD(2) * t378 - qJD(3) * t377) * t352;
t306 = -qJD(3) * t329 + t341 * t410 + t342 * t359;
t305 = qJD(3) * t330 - t341 * t409 + t342 * t356;
t304 = t338 * t403 + ((qJD(3) * t378 + qJD(4) * t413) * t355 + (-t355 * t377 - t358 * t411) * qJD(2)) * t352;
t300 = -t323 * qJD(3) + t339 * t356 + t340 * t409;
t299 = t316 * t353 + t329 * t350;
t298 = t316 * t350 - t329 * t353;
t297 = t325 * t350 - t353 * t417;
t296 = -t325 * t353 - t350 * t417;
t295 = -qJD(4) * t326 + t318 * t358 + t355 * t390;
t281 = qJD(4) * t380 + t306 * t358 - t341 * t414;
t280 = (t339 * t410 + t340 * t359 + (t356 * t372 + t371 * t409) * qJD(3)) * t355 + t328 * t403 + t339 * t412 - t371 * t351 * t404;
t277 = -qJD(4) * t309 + t301 * t358 - t340 * t414;
t273 = 0.1e1 / t275;
t264 = 0.1e1 / t266;
t263 = t290 * t433;
t262 = t382 * t290;
t261 = t383 * t290;
t257 = (-t288 * t322 + t289 * t334) * t355 + t385 * t263;
t256 = t262 * t385 - t288 * t314 + t289 * t331;
t255 = t261 * t385 - t288 * t311 + t289 * t327;
t253 = t382 * t432 + (t331 * t392 - t280 * t319 + (t276 * t331 + t294 * t314 + t304 * t309) * t320) * t290;
t251 = t383 * t432 + (t327 * t392 - t277 * t319 + (t276 * t327 + t294 * t311 + t295 * t309) * t320) * t290;
t250 = t432 * t433 + (t381 * t403 + (t334 * t392 - t300 * t319 + (t276 * t334 + t294 * t322 + t309 * t317) * t320) * t355) * t290;
t1 = [0, t253, t250, t251, 0, 0; 0 (t256 * t428 + t267 * t380) * t402 + ((qJD(4) * t316 + t306 * t355 + t341 * t412) * t267 + t256 * t391 + (t380 * t254 - t256 * t278 - (-t253 * t309 - t262 * t276 + t304 + (-t262 * t326 - t314) * t259) * t421 - (-t253 * t326 - t262 * t294 - t280 + (t262 * t309 - t331) * t259) * t422) * t268) * t264 (t257 * t428 + t267 * t418) * t402 + ((-t302 * t355 - t324 * t403) * t267 + (-t426 + t391) * t257 + (t418 * t254 - (t334 * t403 - t250 * t309 - t263 * t276 + t317 * t355 + (-t263 * t326 - t322 * t355) * t259) * t421 - (-t322 * t403 - t250 * t326 - t263 * t294 - t300 * t355 + (t263 * t309 - t334 * t355) * t259) * t422) * t268) * t264 (t255 * t428 - t267 * t313) * t402 + (t255 * t391 + t279 * t267 + (-t313 * t254 - t255 * t278 - (-t251 * t309 - t261 * t276 + t295 + (-t261 * t326 - t311) * t259) * t421 - (-t251 * t326 - t261 * t294 - t277 + (t261 * t309 - t327) * t259) * t422) * t268) * t264, 0, 0; 0 (-t283 * t298 + t299 * t424) * t401 + ((t281 * t350 - t305 * t353) * t283 + t299 * t393 + (-t298 * t272 - (t281 * t353 + t305 * t350) * t286 - t299 * t271) * t284) * t273 (-t283 * t296 + t297 * t424) * t401 + ((-t303 * t353 + t350 * t375) * t283 + t297 * t393 + (-t296 * t272 - (t303 * t350 + t353 * t375) * t286 - t297 * t271) * t284) * t273 (-t284 * t423 + t425) * t312 * t401 + (-0.2e1 * t312 * t353 * t399 - t278 * t425 + (t278 * t423 + (t271 * t353 + t272 * t350) * t312) * t284) * t273, 0, 0;];
JaD_rot  = t1;
