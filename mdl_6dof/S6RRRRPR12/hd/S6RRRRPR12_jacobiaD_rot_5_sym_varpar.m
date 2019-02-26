% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR12_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR12_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_jacobiaD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:36:59
% EndTime: 2019-02-26 22:37:02
% DurationCPUTime: 2.08s
% Computational Cost: add. (8441->175), mult. (25165->321), div. (705->12), fcn. (31370->15), ass. (0->161)
t349 = cos(pkin(7));
t350 = sin(qJ(3));
t446 = cos(pkin(6));
t447 = sin(qJ(2));
t398 = t446 * t447;
t448 = sin(qJ(1));
t384 = t448 * t398;
t351 = cos(qJ(2));
t352 = cos(qJ(1));
t427 = t352 * t351;
t367 = t384 - t427;
t407 = t351 * t446;
t368 = t352 * t447 + t448 * t407;
t348 = sin(pkin(6));
t445 = sin(pkin(7));
t408 = t348 * t445;
t392 = t448 * t408;
t449 = cos(qJ(3));
t309 = -t367 * t449 + (-t368 * t349 + t392) * t350;
t430 = t348 * t349;
t402 = t448 * t430;
t328 = t368 * t445 + t402;
t347 = qJ(4) + pkin(13);
t345 = sin(t347);
t346 = cos(t347);
t289 = t309 * t345 - t328 * t346;
t457 = 0.2e1 * t289;
t401 = t449 * t447;
t428 = t350 * t351;
t375 = t349 * t401 + t428;
t377 = t349 * t428 + t401;
t394 = t446 * t445;
t391 = t350 * t394;
t299 = qJD(3) * t391 + (t375 * qJD(2) + t377 * qJD(3)) * t348;
t411 = t447 * t350;
t412 = t449 * t351;
t376 = t349 * t412 - t411;
t383 = t449 * t394;
t325 = -t376 * t348 - t383;
t323 = 0.1e1 / t325 ^ 2;
t456 = t299 * t323;
t344 = t448 * t447;
t390 = -t352 * t407 + t344;
t333 = t448 * t351 + t352 * t398;
t400 = t352 * t408;
t385 = t449 * t400;
t369 = -t333 * t350 - t385;
t373 = t390 * t449;
t371 = t349 * t373;
t303 = t371 - t369;
t301 = t303 ^ 2;
t295 = t301 * t323 + 0.1e1;
t293 = 0.1e1 / t295;
t381 = t390 * t350;
t414 = t333 * t449;
t363 = -t349 * t381 + t414;
t319 = t368 * qJD(1) + t333 * qJD(2);
t320 = -qJD(1) * t384 - qJD(2) * t344 + (qJD(2) * t446 + qJD(1)) * t427;
t340 = t350 * t400;
t379 = t449 * t392;
t413 = t349 * t449;
t372 = -qJD(1) * t379 - qJD(3) * t340 + t319 * t413 + t320 * t350;
t278 = t363 * qJD(3) + t372;
t322 = 0.1e1 / t325;
t432 = t303 * t323;
t389 = -t278 * t322 + t299 * t432;
t260 = t389 * t293;
t296 = atan2(-t303, t325);
t291 = sin(t296);
t292 = cos(t296);
t393 = -t291 * t325 - t292 * t303;
t255 = t393 * t260 - t291 * t278 + t292 * t299;
t272 = -t291 * t303 + t292 * t325;
t269 = 0.1e1 / t272;
t270 = 0.1e1 / t272 ^ 2;
t366 = t368 * t449;
t454 = -t349 * t366 + t350 * t367 + t379;
t302 = t454 ^ 2;
t266 = t302 * t270 + 0.1e1;
t264 = 0.1e1 / t266;
t439 = t264 * t270;
t318 = t333 * qJD(1) + t368 * qJD(2);
t361 = qJD(1) * t390 + t367 * qJD(2);
t359 = t361 * t449;
t276 = -qJD(1) * t385 + t309 * qJD(3) - t318 * t350 - t349 * t359;
t437 = t270 * t454;
t443 = t255 * t269 * t270;
t444 = (-t276 * t437 - t302 * t443) / t266 ^ 2;
t455 = -t255 * t439 - 0.2e1 * t269 * t444;
t450 = -0.2e1 * t454;
t405 = t443 * t450;
t425 = 0.2e1 * t444;
t453 = t264 * t405 - t276 * t439 - t425 * t437;
t452 = -(qJD(1) * t392 - t333 * qJD(3) - t319 * t349) * t350 + qJD(3) * t385 - t320 * t449;
t290 = t309 * t346 + t328 * t345;
t284 = 0.1e1 / t290;
t285 = 0.1e1 / t290 ^ 2;
t451 = -0.2e1 * t303;
t360 = t361 * t350;
t277 = qJD(1) * t340 + qJD(3) * t454 - t318 * t449 + t349 * t360;
t415 = t352 * t430;
t310 = qJD(1) * t415 - t361 * t445;
t267 = t290 * qJD(4) + t277 * t345 - t310 * t346;
t283 = t289 ^ 2;
t275 = t283 * t285 + 0.1e1;
t435 = t285 * t289;
t426 = qJD(4) * t289;
t268 = t277 * t346 + t310 * t345 - t426;
t438 = t268 * t284 * t285;
t442 = (t267 * t435 - t283 * t438) / t275 ^ 2;
t434 = t322 * t456;
t441 = (t278 * t432 - t301 * t434) / t295 ^ 2;
t440 = t264 * t269;
t273 = 0.1e1 / t275;
t436 = t273 * t285;
t433 = t303 * t322;
t429 = t349 * t350;
t424 = -0.2e1 * t442;
t423 = -0.2e1 * t441;
t421 = t285 * t442;
t420 = t322 * t441;
t418 = t264 * t437;
t417 = t267 * t436;
t416 = t289 * t438;
t410 = t345 * t445;
t409 = t346 * t445;
t404 = 0.2e1 * t416;
t403 = t434 * t451;
t382 = t349 * t390;
t364 = t350 * t382 - t414;
t307 = t340 + t364;
t327 = -t390 * t445 + t415;
t288 = t307 * t346 + t327 * t345;
t287 = t307 * t345 - t327 * t346;
t388 = -t345 * t284 + t346 * t435;
t305 = -t340 + t363;
t326 = t377 * t348 + t391;
t387 = -t305 * t322 + t326 * t432;
t315 = t333 * t413 - t381;
t332 = t375 * t348;
t386 = -t315 * t322 + t332 * t432;
t317 = t367 * t429 - t366;
t298 = t317 * t346 - t367 * t410;
t380 = -t317 * t345 - t367 * t409;
t378 = -t291 + (t292 * t433 + t291) * t293;
t374 = -t349 * t411 + t412;
t370 = t449 * t382;
t316 = -t368 * t350 - t367 * t413;
t312 = (t376 * qJD(2) + t374 * qJD(3)) * t348;
t311 = -qJD(1) * t402 - t319 * t445;
t300 = qJD(3) * t383 + (t374 * qJD(2) + t376 * qJD(3)) * t348;
t282 = t320 * t413 - t319 * t350 + (-t333 * t429 - t373) * qJD(3);
t281 = -t316 * qJD(3) + t318 * t429 + t359;
t280 = qJD(3) * t370 + t452;
t279 = -qJD(3) * t371 - t452;
t263 = t386 * t293;
t262 = t387 * t293;
t256 = t393 * t262 - t291 * t305 + t292 * t326;
t254 = t386 * t423 + (t332 * t403 - t282 * t322 + (t278 * t332 + t299 * t315 + t303 * t312) * t323) * t293;
t253 = t387 * t423 + (t326 * t403 - t279 * t322 + (t278 * t326 + t299 * t305 + t300 * t303) * t323) * t293;
t1 = [t420 * t450 + (-t276 * t322 - t454 * t456) * t293, t254, t253, 0, 0, 0; (t364 * qJD(3) - t372) * t440 + (t378 * t276 - ((-t260 * t293 * t433 + t423) * t291 + (t420 * t451 - t260 + (t260 - t389) * t293) * t292) * t454) * t418 + t455 * (-t370 + t369) - t453 * t378 * t454 (t317 * qJD(3) - t318 * t413 + t360) * t440 + ((-t254 * t303 - t263 * t278 + t312 + (-t263 * t325 - t315) * t260) * t292 + (-t254 * t325 - t263 * t299 - t282 + (t263 * t303 - t332) * t260) * t291) * t418 + t455 * t316 + t453 * (t393 * t263 - t291 * t315 + t292 * t332) (-t256 * t437 - t269 * t309) * t425 + (t256 * t405 + t277 * t269 + (-t309 * t255 - t256 * t276 - (-(-t253 * t303 - t262 * t278 + t300 + (-t262 * t325 - t305) * t260) * t292 - (-t253 * t325 - t262 * t299 - t279 + (t262 * t303 - t326) * t260) * t291) * t454) * t270) * t264, 0, 0, 0; 0.2e1 * (-t284 * t287 + t288 * t435) * t442 + ((t288 * qJD(4) + t280 * t345 - t311 * t346) * t284 + t288 * t404 + (-t287 * t268 - (-t287 * qJD(4) + t280 * t346 + t311 * t345) * t289 - t288 * t267) * t285) * t273 (t421 * t457 - t417) * t298 - (-t268 * t436 + t284 * t424) * t380 + ((t298 * qJD(4) + t281 * t345 + t318 * t409) * t284 - (t380 * qJD(4) + t281 * t346 - t318 * t410) * t435 + t298 * t404) * t273, -t388 * t454 * t424 + (t388 * t276 - ((-qJD(4) * t284 - 0.2e1 * t416) * t346 + (t267 * t346 + (t268 - t426) * t345) * t285) * t454) * t273, t424 + (t417 + (-t273 * t438 - t421) * t289) * t457, 0, 0;];
JaD_rot  = t1;
