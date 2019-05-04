% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR10V2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRR10V2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10V2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_jacobiaD_rot_6_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:56:32
% EndTime: 2019-04-11 14:56:35
% DurationCPUTime: 2.63s
% Computational Cost: add. (17109->206), mult. (24753->390), div. (1249->12), fcn. (30668->13), ass. (0->164)
t351 = qJ(2) + qJ(3);
t348 = sin(t351);
t353 = sin(qJ(5));
t357 = cos(qJ(5));
t350 = qJD(2) + qJD(3);
t358 = cos(qJ(4));
t389 = qJD(5) * t358 - t350;
t354 = sin(qJ(4));
t419 = qJD(4) * t354;
t460 = (t389 * t353 + t357 * t419) * t348;
t349 = cos(t351);
t355 = sin(qJ(1));
t427 = t350 * t355;
t359 = cos(qJ(1));
t429 = t348 * t359;
t459 = qJD(1) * t429 + t349 * t427;
t421 = t359 * t354;
t423 = t355 * t358;
t333 = t349 * t423 - t421;
t430 = t348 * t357;
t315 = t333 * t353 - t355 * t430;
t425 = t353 * t358;
t428 = t349 * t357;
t329 = t348 * t425 + t428;
t308 = atan2(-t315, t329);
t299 = sin(t308);
t300 = cos(t308);
t277 = -t299 * t315 + t300 * t329;
t275 = 0.1e1 / t277 ^ 2;
t422 = t358 * t359;
t424 = t355 * t354;
t335 = t349 * t422 + t424;
t320 = t335 * t353 - t357 * t429;
t314 = t320 ^ 2;
t273 = t275 * t314 + 0.1e1;
t387 = -qJD(1) * t349 + qJD(4);
t388 = -qJD(4) * t349 + qJD(1);
t426 = t350 * t359;
t304 = t388 * t421 + (-t348 * t426 + t387 * t355) * t358;
t321 = t335 * t357 + t353 * t429;
t420 = qJD(1) * t355;
t455 = t348 * t420 - t349 * t426;
t278 = t321 * qJD(5) + t304 * t353 + t455 * t357;
t444 = t275 * t320;
t313 = t315 ^ 2;
t327 = 0.1e1 / t329 ^ 2;
t307 = t313 * t327 + 0.1e1;
t301 = 0.1e1 / t307;
t418 = qJD(4) * t358;
t394 = t359 * t418;
t395 = t355 * t419;
t404 = t348 * t427;
t306 = t335 * qJD(1) - t349 * t395 - t358 * t404 - t394;
t317 = t348 * t353 * t355 + t333 * t357;
t280 = t317 * qJD(5) + t306 * t353 - t357 * t459;
t390 = t350 * t358 - qJD(5);
t382 = t349 * t390;
t367 = -t353 * t419 + t389 * t357;
t456 = t367 * t348;
t289 = t353 * t382 + t456;
t326 = 0.1e1 / t329;
t434 = t315 * t327;
t379 = -t280 * t326 + t289 * t434;
t264 = t379 * t301;
t383 = -t299 * t329 - t300 * t315;
t258 = t383 * t264 - t280 * t299 + t289 * t300;
t274 = 0.1e1 / t277;
t276 = t274 * t275;
t449 = t258 * t276;
t414 = 0.2e1 * (t278 * t444 - t314 * t449) / t273 ^ 2;
t458 = t289 * t327;
t375 = t349 * t424 + t422;
t431 = t348 * t354;
t405 = t315 * t431;
t374 = -t326 * t375 + t327 * t405;
t457 = t353 * t374;
t330 = -t349 * t353 + t358 * t430;
t324 = t330 * t359;
t400 = t350 * t421;
t454 = (t354 * t420 - t394) * t348 + qJD(6) * t324 - t349 * t400;
t417 = qJD(5) * t348;
t281 = (-qJD(5) * t333 + t459) * t353 + (t355 * t417 + t306) * t357;
t334 = t349 * t421 - t423;
t352 = sin(qJ(6));
t356 = cos(qJ(6));
t298 = t321 * t356 + t334 * t352;
t292 = 0.1e1 / t298;
t293 = 0.1e1 / t298 ^ 2;
t452 = -0.2e1 * t315;
t451 = 0.2e1 * t320;
t279 = (t359 * t417 + t304) * t357 + (-qJD(5) * t335 - t455) * t353;
t303 = t375 * qJD(1) + t348 * t400 - t349 * t394 - t395;
t266 = t298 * qJD(6) + t279 * t352 + t303 * t356;
t297 = t321 * t352 - t334 * t356;
t291 = t297 ^ 2;
t285 = t291 * t293 + 0.1e1;
t440 = t293 * t297;
t415 = qJD(6) * t297;
t267 = t279 * t356 - t303 * t352 - t415;
t446 = t267 * t292 * t293;
t448 = (t266 * t440 - t291 * t446) / t285 ^ 2;
t443 = t326 * t458;
t447 = (t280 * t434 - t313 * t443) / t307 ^ 2;
t445 = t275 * t278;
t442 = t292 * t352;
t441 = t292 * t356;
t439 = t297 * t352;
t438 = t297 * t356;
t437 = t299 * t320;
t436 = t300 * t320;
t435 = t315 * t326;
t433 = t334 * t353;
t432 = t334 * t357;
t416 = qJD(5) * t357;
t413 = -0.2e1 * t448;
t412 = 0.2e1 * t448;
t411 = -0.2e1 * t447;
t410 = t276 * t451;
t409 = t326 * t447;
t408 = t297 * t446;
t407 = t275 * t437;
t406 = t275 * t436;
t403 = t348 * t421;
t393 = t258 * t410;
t392 = 0.2e1 * t408;
t391 = t443 * t452;
t384 = qJD(6) * t432 + t304;
t296 = -t317 * t356 - t352 * t375;
t295 = -t317 * t352 + t356 * t375;
t381 = t390 * t353;
t380 = -qJD(6) * t403 + t330 * t420 + (-t390 * t428 + t460) * t359;
t378 = t293 * t438 - t442;
t377 = -t317 * t326 + t330 * t434;
t322 = t329 * t355;
t331 = t349 * t425 - t430;
t376 = t322 * t326 + t331 * t434;
t373 = -t349 * t350 * t354 - t348 * t418;
t370 = -t299 + (t300 * t435 + t299) * t301;
t369 = qJD(1) * t329;
t368 = qJD(5) * t433 + qJD(6) * t335 + t303 * t357;
t323 = t329 * t359;
t312 = -t324 * t356 - t352 * t403;
t311 = -t324 * t352 + t356 * t403;
t310 = t335 * t352 - t356 * t432;
t309 = -t335 * t356 - t352 * t432;
t305 = t388 * t423 + (t387 * t359 + t404) * t354;
t290 = t357 * t382 - t460;
t288 = -t348 * t381 + t367 * t349;
t287 = t289 * t355 + t359 * t369;
t283 = 0.1e1 / t285;
t271 = 0.1e1 / t273;
t270 = t301 * t457;
t269 = t376 * t301;
t268 = t377 * t301;
t263 = t370 * t320;
t261 = (t299 * t375 - t300 * t431) * t353 - t383 * t270;
t260 = t383 * t269 + t299 * t322 + t300 * t331;
t259 = t383 * t268 - t299 * t317 + t300 * t330;
t256 = t376 * t411 + (t331 * t391 + t287 * t326 + (t280 * t331 + t288 * t315 - t289 * t322) * t327) * t301;
t255 = t377 * t411 + (t330 * t391 - t281 * t326 + (t280 * t330 + t289 * t317 + t290 * t315) * t327) * t301;
t254 = 0.2e1 * t447 * t457 + (-t374 * t416 + (0.2e1 * t405 * t443 - t305 * t326 + (-t280 * t431 - t289 * t375 + t373 * t315) * t327) * t353) * t301;
t253 = (-t292 * t311 + t312 * t440) * t412 + (t312 * t392 + t380 * t442 - t454 * t441 + (-t312 * t266 - t311 * t267 - t380 * t438 - t454 * t439) * t293) * t283;
t252 = (t260 * t444 + t274 * t323) * t414 + (t260 * t393 + (t323 * t258 - t260 * t278 - (-t256 * t315 - t269 * t280 + t288 + (-t269 * t329 + t322) * t264) * t436 - (-t256 * t329 - t269 * t289 + t287 + (t269 * t315 - t331) * t264) * t437) * t275 + (t355 * t369 + (-t349 * t381 - t456) * t359) * t274) * t271;
t1 = [t409 * t451 + (-t278 * t326 + t320 * t458) * t301, t256, t256, t254, t255, 0; t315 * t274 * t414 + (-t280 * t274 + (t258 * t315 - t263 * t278) * t275) * t271 + (t263 * t275 * t414 + (0.2e1 * t263 * t449 - (-t264 * t301 * t435 + t411) * t407 - (t409 * t452 - t264 + (t264 - t379) * t301) * t406 - t370 * t445) * t271) * t320, t252, t252 (t261 * t444 + t274 * t433) * t414 + (-t261 * t445 + (t303 * t353 - t334 * t416) * t274 + (t261 * t410 + t275 * t433) * t258 - (t375 * t416 - t254 * t329 + t270 * t289 - t305 * t353 + (-t270 * t315 + t353 * t431) * t264) * t407 - (-t416 * t431 - t254 * t315 - (-t264 * t329 - t280) * t270 + (t264 * t375 + t373) * t353) * t406) * t271 (t259 * t444 - t274 * t321) * t414 + (t259 * t393 + t279 * t274 + (-t321 * t258 - t259 * t278 - (-t255 * t315 - t268 * t280 + t290 + (-t268 * t329 - t317) * t264) * t436 - (-t255 * t329 - t268 * t289 - t281 + (t268 * t315 - t330) * t264) * t437) * t275) * t271, 0; (-t292 * t295 + t296 * t440) * t412 + ((t296 * qJD(6) - t281 * t352 - t305 * t356) * t292 + t296 * t392 + (-t295 * t267 - (-t295 * qJD(6) - t281 * t356 + t305 * t352) * t297 - t296 * t266) * t293) * t283, t253, t253 (-t292 * t309 + t310 * t440) * t412 + (t310 * t392 - t384 * t441 + t368 * t442 + (-t310 * t266 - t309 * t267 - t368 * t438 - t384 * t439) * t293) * t283, t378 * t320 * t413 + (t378 * t278 + ((-qJD(6) * t292 - 0.2e1 * t408) * t356 + (t266 * t356 + (t267 - t415) * t352) * t293) * t320) * t283, t413 + 0.2e1 * (t266 * t293 * t283 + (-t283 * t446 - t293 * t448) * t297) * t297;];
JaD_rot  = t1;
