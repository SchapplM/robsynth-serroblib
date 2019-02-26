% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRP8_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP8_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:43:49
% EndTime: 2019-02-26 22:43:52
% DurationCPUTime: 3.18s
% Computational Cost: add. (20345->181), mult. (36096->349), div. (1229->12), fcn. (45886->13), ass. (0->153)
t355 = qJ(3) + qJ(4);
t352 = sin(t355);
t353 = cos(t355);
t357 = cos(pkin(6));
t360 = sin(qJ(1));
t362 = cos(qJ(2));
t422 = t360 * t362;
t359 = sin(qJ(2));
t363 = cos(qJ(1));
t423 = t359 * t363;
t345 = t357 * t423 + t422;
t354 = qJD(3) + qJD(4);
t356 = sin(pkin(6));
t419 = qJD(1) * t360;
t380 = -t345 * t354 + t356 * t419;
t418 = qJD(2) * t359;
t399 = t360 * t418;
t401 = t359 * t419;
t420 = t362 * t363;
t326 = -t357 * t401 - t399 + (qJD(2) * t357 + qJD(1)) * t420;
t426 = t356 * t363;
t388 = t354 * t426 - t326;
t291 = t380 * t352 - t388 * t353;
t336 = -t345 * t353 + t352 * t426;
t358 = sin(qJ(5));
t361 = cos(qJ(5));
t424 = t359 * t360;
t395 = -t357 * t420 + t424;
t313 = t336 * t358 + t395 * t361;
t381 = t357 * t422 + t423;
t370 = t381 * qJD(1) + t345 * qJD(2);
t286 = t313 * qJD(5) + t291 * t361 + t370 * t358;
t387 = t395 * t358;
t314 = t336 * t361 - t387;
t463 = t314 * qJD(5) - t291 * t358 + t370 * t361;
t305 = t313 ^ 2;
t428 = t356 * t359;
t341 = t352 * t357 + t353 * t428;
t421 = t361 * t362;
t331 = t341 * t358 + t356 * t421;
t328 = 0.1e1 / t331 ^ 2;
t299 = t305 * t328 + 0.1e1;
t295 = 0.1e1 / t299;
t378 = qJD(2) * t356 * t362 + t354 * t357;
t402 = t354 * t428;
t323 = -t352 * t402 + t378 * t353;
t425 = t358 * t362;
t332 = t341 * t361 - t356 * t425;
t400 = t356 * t418;
t302 = t332 * qJD(5) + t323 * t358 - t361 * t400;
t327 = 0.1e1 / t331;
t431 = t313 * t328;
t385 = -t302 * t431 + t327 * t463;
t270 = t385 * t295;
t300 = atan2(t313, t331);
t293 = sin(t300);
t294 = cos(t300);
t386 = -t293 * t331 + t294 * t313;
t265 = t386 * t270 + t293 * t463 + t294 * t302;
t282 = t293 * t313 + t294 * t331;
t280 = 0.1e1 / t282 ^ 2;
t462 = t265 * t280;
t346 = -t357 * t424 + t420;
t379 = qJD(1) * t426 - t346 * t354;
t325 = -t345 * qJD(1) - t381 * qJD(2);
t427 = t356 * t360;
t389 = t354 * t427 + t325;
t289 = t379 * t352 + t389 * t353;
t338 = t346 * t353 + t352 * t427;
t316 = t338 * t361 + t381 * t358;
t324 = t357 * t399 + t401 + (-qJD(1) * t357 - qJD(2)) * t420;
t283 = t316 * qJD(5) + t289 * t358 + t324 * t361;
t315 = t338 * t358 - t381 * t361;
t447 = 0.2e1 * t315;
t279 = 0.1e1 / t282;
t457 = t279 * t462;
t397 = t447 * t457;
t461 = -t280 * t283 + t397;
t288 = -t389 * t352 + t379 * t353;
t308 = 0.1e1 / t316 ^ 2;
t337 = -t346 * t352 + t353 * t427;
t330 = t337 ^ 2;
t434 = t308 * t330;
t301 = 0.1e1 + t434;
t284 = -t315 * qJD(5) + t289 * t361 - t324 * t358;
t307 = 0.1e1 / t316;
t440 = t284 * t307 * t308;
t406 = t330 * t440;
t433 = t308 * t337;
t445 = (t288 * t433 - t406) / t301 ^ 2;
t458 = 0.2e1 * t445;
t456 = -0.2e1 * t315;
t455 = t302 * t328;
t335 = t345 * t352 + t353 * t426;
t340 = -t352 * t428 + t353 * t357;
t382 = t327 * t335 - t340 * t431;
t454 = t358 * t382;
t297 = 0.1e1 / t301;
t436 = t297 * t308;
t452 = t284 * t436 + t307 * t458;
t306 = t315 ^ 2;
t278 = t280 * t306 + 0.1e1;
t276 = 0.1e1 / t278;
t441 = t280 * t315;
t446 = (t283 * t441 - t306 * t457) / t278 ^ 2;
t451 = -t276 * t462 - 0.2e1 * t279 * t446;
t394 = t433 * t445;
t405 = t337 * t440;
t450 = -t288 * t436 + 0.2e1 * t297 * t405 + 0.2e1 * t394;
t415 = 0.2e1 * t446;
t449 = t276 * t461 + t415 * t441;
t290 = t388 * t352 + t380 * t353;
t448 = 0.2e1 * t313;
t435 = t327 * t455;
t444 = (-t305 * t435 + t431 * t463) / t299 ^ 2;
t443 = t276 * t279;
t439 = t293 * t315;
t438 = t294 * t315;
t437 = t297 * t307;
t432 = t313 * t327;
t430 = t337 * t358;
t429 = t352 * t354;
t417 = qJD(5) * t353;
t416 = qJD(5) * t361;
t414 = -0.2e1 * t444;
t410 = t327 * t444;
t408 = t276 * t441;
t403 = t297 * t433;
t396 = t435 * t448;
t384 = t314 * t327 - t332 * t431;
t317 = -t345 * t361 - t353 * t387;
t339 = (t353 * t425 - t359 * t361) * t356;
t383 = -t317 * t327 - t339 * t431;
t376 = t353 * t381;
t375 = t354 * t381;
t374 = -t293 + (-t294 * t432 + t293) * t295;
t373 = -qJD(5) * t376 - t325;
t372 = t346 * qJD(5) + t324 * t353 + t352 * t375;
t322 = -t378 * t352 - t353 * t402;
t304 = ((-qJD(2) + t417) * t421 + (-t362 * t429 + (-qJD(2) * t353 + qJD(5)) * t359) * t358) * t356;
t303 = -t331 * qJD(5) + t323 * t361 + t358 * t400;
t287 = (-t395 * t417 - t326) * t361 + (t345 * qJD(5) - t353 * t370 + t395 * t429) * t358;
t275 = t295 * t454;
t274 = t383 * t295;
t273 = t384 * t295;
t268 = (t293 * t335 + t294 * t340) * t358 + t386 * t275;
t266 = t386 * t273 + t293 * t314 + t294 * t332;
t264 = t383 * t414 + (t339 * t396 - t287 * t327 + (t302 * t317 - t304 * t313 - t339 * t463) * t328) * t295;
t263 = (t307 * t338 + t361 * t434) * t458 + (0.2e1 * t361 * t406 - t289 * t307 + (qJD(5) * t330 * t358 - 0.2e1 * t288 * t337 * t361 + t284 * t338) * t308) * t297;
t261 = t384 * t414 + (t332 * t396 - t286 * t327 + (-t302 * t314 - t303 * t313 - t332 * t463) * t328) * t295;
t260 = t414 * t454 + (t382 * t416 + (t340 * t396 - t290 * t327 + (-t302 * t335 - t313 * t322 - t340 * t463) * t328) * t358) * t295;
t259 = (t268 * t441 - t279 * t430) * t415 + ((t288 * t358 + t337 * t416) * t279 + t461 * t268 + (-t430 * t265 - (t340 * t416 + t260 * t313 + t275 * t463 + t322 * t358 + (-t275 * t331 + t335 * t358) * t270) * t438 - (t335 * t416 - t260 * t331 - t275 * t302 - t290 * t358 + (-t275 * t313 - t340 * t358) * t270) * t439) * t280) * t276;
t1 = [t410 * t447 + (-t283 * t327 + t315 * t455) * t295, t264, t260, t260, t261, 0; t463 * t443 - (t374 * t283 + ((t270 * t295 * t432 + t414) * t293 + (t410 * t448 - t270 + (t270 - t385) * t295) * t294) * t315) * t408 + t451 * t313 + t449 * t374 * t315 (t372 * t358 + t373 * t361) * t443 - ((t264 * t313 + t274 * t463 + t304 + (-t274 * t331 - t317) * t270) * t294 + (-t264 * t331 - t274 * t302 - t287 + (-t274 * t313 - t339) * t270) * t293) * t408 + t451 * (-t346 * t361 - t358 * t376) + t449 * (t386 * t274 - t293 * t317 + t294 * t339) t259, t259 (t266 * t441 - t279 * t316) * t415 + (t266 * t397 + t284 * t279 + (-t316 * t265 - t266 * t283 - (t261 * t313 + t273 * t463 + t303 + (-t273 * t331 + t314) * t270) * t438 - (-t261 * t331 - t273 * t302 - t286 + (-t273 * t313 - t332) * t270) * t439) * t280) * t276, 0; t286 * t403 - t290 * t437 + t450 * t314 - t452 * t335 -(-t373 * t358 + t372 * t361) * t403 + (-t324 * t352 + t353 * t375) * t437 - t452 * t352 * t381 + t450 * (t346 * t358 - t361 * t376) t263, t263, t394 * t456 + (t405 * t456 + (t283 * t337 + t288 * t315) * t308) * t297, 0;];
JaD_rot  = t1;
