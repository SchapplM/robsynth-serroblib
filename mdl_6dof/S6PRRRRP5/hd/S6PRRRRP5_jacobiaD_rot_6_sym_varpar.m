% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRRP5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:17:29
% EndTime: 2019-02-26 20:17:32
% DurationCPUTime: 2.78s
% Computational Cost: add. (14691->202), mult. (44088->385), div. (816->12), fcn. (56595->17), ass. (0->161)
t400 = sin(qJ(3));
t404 = cos(qJ(3));
t395 = sin(pkin(7));
t397 = cos(pkin(7));
t401 = sin(qJ(2));
t405 = cos(qJ(2));
t481 = cos(pkin(12));
t482 = cos(pkin(6));
t434 = t482 * t481;
t480 = sin(pkin(12));
t420 = t480 * t401 - t405 * t434;
t396 = sin(pkin(6));
t444 = t396 * t481;
t412 = -t395 * t444 - t420 * t397;
t419 = -t401 * t434 - t480 * t405;
t368 = t400 * t419 + t412 * t404;
t385 = t420 * qJD(2);
t386 = t419 * qJD(2);
t462 = t397 * t400;
t346 = t368 * qJD(3) - t385 * t404 + t386 * t462;
t369 = t412 * t400 - t404 * t419;
t399 = sin(qJ(4));
t403 = cos(qJ(4));
t413 = t420 * t395 - t397 * t444;
t356 = t369 * t403 + t413 * t399;
t464 = t395 * t403;
t321 = t356 * qJD(4) + t346 * t399 + t386 * t464;
t354 = t369 * t399 - t413 * t403;
t352 = t354 ^ 2;
t457 = t401 * t404;
t458 = t400 * t405;
t423 = t397 * t458 + t457;
t445 = t395 * t482;
t381 = t423 * t396 + t400 * t445;
t463 = t395 * t405;
t389 = -t396 * t463 + t482 * t397;
t372 = t381 * t399 - t389 * t403;
t366 = 0.1e1 / t372 ^ 2;
t337 = t352 * t366 + 0.1e1;
t335 = 0.1e1 / t337;
t455 = t404 * t405;
t459 = t400 * t401;
t422 = -t397 * t459 + t455;
t425 = t397 * t455 - t459;
t436 = qJD(3) * t445;
t363 = t404 * t436 + (t422 * qJD(2) + t425 * qJD(3)) * t396;
t373 = t381 * t403 + t389 * t399;
t465 = t395 * t401;
t446 = t396 * t465;
t438 = qJD(2) * t446;
t339 = t373 * qJD(4) + t363 * t399 - t403 * t438;
t365 = 0.1e1 / t372;
t471 = t354 * t366;
t304 = (-t321 * t365 + t339 * t471) * t335;
t338 = atan2(-t354, t372);
t333 = sin(t338);
t334 = cos(t338);
t432 = -t333 * t372 - t334 * t354;
t299 = t432 * t304 - t321 * t333 + t334 * t339;
t317 = -t333 * t354 + t334 * t372;
t314 = 0.1e1 / t317;
t315 = 0.1e1 / t317 ^ 2;
t485 = t299 * t314 * t315;
t433 = t482 * t480;
t417 = t481 * t401 + t405 * t433;
t443 = t396 * t480;
t437 = t395 * t443;
t414 = -t417 * t397 + t437;
t418 = t401 * t433 - t481 * t405;
t371 = t414 * t400 - t404 * t418;
t415 = t417 * t395 + t397 * t443;
t357 = t371 * t399 - t415 * t403;
t439 = 0.2e1 * t357 * t485;
t380 = t425 * t396 + t404 * t445;
t427 = -t365 * t368 + t380 * t471;
t484 = t399 * t427;
t472 = t339 * t365 * t366;
t483 = -0.2e1 * (t321 * t471 - t352 * t472) / t337 ^ 2;
t358 = t371 * t403 + t415 * t399;
t416 = t417 * t404;
t467 = t418 * t400;
t370 = t397 * t416 - t404 * t437 - t467;
t398 = sin(qJ(5));
t402 = cos(qJ(5));
t332 = t358 * t402 + t370 * t398;
t328 = 0.1e1 / t332;
t329 = 0.1e1 / t332 ^ 2;
t387 = t417 * qJD(2);
t388 = t418 * qJD(2);
t348 = t388 * t462 - t387 * t404 + (t414 * t404 + t467) * qJD(3);
t466 = t395 * t399;
t324 = -t357 * qJD(4) + t348 * t403 - t388 * t466;
t461 = t397 * t404;
t347 = t371 * qJD(3) - t387 * t400 - t388 * t461;
t312 = t332 * qJD(5) + t324 * t398 - t347 * t402;
t331 = t358 * t398 - t370 * t402;
t327 = t331 ^ 2;
t320 = t327 * t329 + 0.1e1;
t475 = t329 * t331;
t452 = qJD(5) * t331;
t313 = t324 * t402 + t347 * t398 - t452;
t478 = t313 * t328 * t329;
t479 = (t312 * t475 - t327 * t478) / t320 ^ 2;
t477 = t315 * t357;
t323 = t358 * qJD(4) + t348 * t399 + t388 * t464;
t476 = t323 * t315;
t474 = t333 * t357;
t473 = t334 * t357;
t470 = t370 * t399;
t469 = t370 * t403;
t460 = t398 * t328;
t456 = t402 * t331;
t454 = qJD(4) * t399;
t453 = qJD(4) * t403;
t353 = t357 ^ 2;
t311 = t315 * t353 + 0.1e1;
t451 = 0.2e1 * (-t353 * t485 + t357 * t476) / t311 ^ 2;
t450 = -0.2e1 * t479;
t449 = 0.2e1 * t479;
t447 = t331 * t478;
t441 = 0.2e1 * t447;
t440 = -0.2e1 * t354 * t472;
t435 = qJD(5) * t469 + t348;
t376 = t418 * t462 - t416;
t361 = t376 * t403 - t418 * t466;
t375 = -t417 * t400 - t418 * t461;
t344 = t361 * t402 + t375 * t398;
t343 = t361 * t398 - t375 * t402;
t430 = t329 * t456 - t460;
t429 = -t356 * t365 + t373 * t471;
t374 = -t420 * t404 + t419 * t462;
t359 = t374 * t399 + t419 * t464;
t384 = t422 * t396;
t377 = t384 * t399 - t403 * t446;
t428 = -t359 * t365 + t377 * t471;
t426 = -t376 * t399 - t418 * t464;
t424 = -t397 * t457 - t458;
t421 = qJD(5) * t371 - t347 * t403 + t370 * t454;
t362 = -t400 * t436 + (t424 * qJD(2) - t423 * qJD(3)) * t396;
t351 = -t375 * qJD(3) + t387 * t462 + t388 * t404;
t350 = t376 * qJD(3) - t387 * t461 + t388 * t400;
t349 = t384 * t453 + ((t424 * qJD(3) + qJD(4) * t465) * t399 + (-t423 * t399 - t403 * t463) * qJD(2)) * t396;
t345 = -qJD(3) * t369 + t385 * t400 + t386 * t461;
t342 = t371 * t398 - t402 * t469;
t341 = -t371 * t402 - t398 * t469;
t340 = -t372 * qJD(4) + t363 * t403 + t399 * t438;
t326 = t426 * qJD(4) + t351 * t403 - t387 * t466;
t325 = (t385 * t462 + t386 * t404 + (t420 * t400 + t419 * t461) * qJD(3)) * t399 + t374 * t453 + t385 * t464 - t419 * t395 * t454;
t322 = -t354 * qJD(4) + t346 * t403 - t386 * t466;
t318 = 0.1e1 / t320;
t309 = 0.1e1 / t311;
t308 = t335 * t484;
t307 = t428 * t335;
t306 = t429 * t335;
t302 = (-t333 * t368 + t334 * t380) * t399 + t432 * t308;
t301 = t432 * t307 - t333 * t359 + t334 * t377;
t300 = t432 * t306 - t333 * t356 + t334 * t373;
t298 = t428 * t483 + (t377 * t440 - t325 * t365 + (t321 * t377 + t339 * t359 + t349 * t354) * t366) * t335;
t296 = t429 * t483 + (t373 * t440 - t322 * t365 + (t321 * t373 + t339 * t356 + t340 * t354) * t366) * t335;
t295 = t483 * t484 + (t427 * t453 + (t380 * t440 - t345 * t365 + (t321 * t380 + t339 * t368 + t354 * t362) * t366) * t399) * t335;
t1 = [0, t298, t295, t296, 0, 0; 0 (t301 * t477 + t314 * t426) * t451 + ((t361 * qJD(4) + t351 * t399 + t387 * t464) * t314 + t301 * t439 + (t426 * t299 - t301 * t323 - (-t298 * t354 - t307 * t321 + t349 + (-t307 * t372 - t359) * t304) * t473 - (-t298 * t372 - t307 * t339 - t325 + (t307 * t354 - t377) * t304) * t474) * t315) * t309 (t302 * t477 + t314 * t470) * t451 + ((-t347 * t399 - t370 * t453) * t314 + (-t476 + t439) * t302 + (t470 * t299 - (t380 * t453 - t295 * t354 - t308 * t321 + t362 * t399 + (-t308 * t372 - t368 * t399) * t304) * t473 - (-t368 * t453 - t295 * t372 - t308 * t339 - t345 * t399 + (t308 * t354 - t380 * t399) * t304) * t474) * t315) * t309 (t300 * t477 - t314 * t358) * t451 + (t300 * t439 + t324 * t314 + (-t358 * t299 - t300 * t323 - (-t296 * t354 - t306 * t321 + t340 + (-t306 * t372 - t356) * t304) * t473 - (-t296 * t372 - t306 * t339 - t322 + (t306 * t354 - t373) * t304) * t474) * t315) * t309, 0, 0; 0 (-t328 * t343 + t344 * t475) * t449 + ((t344 * qJD(5) + t326 * t398 - t350 * t402) * t328 + t344 * t441 + (-t343 * t313 - (-t343 * qJD(5) + t326 * t402 + t350 * t398) * t331 - t344 * t312) * t329) * t318 (-t328 * t341 + t342 * t475) * t449 + (t342 * t441 - t435 * t328 * t402 + t421 * t460 + (-t435 * t331 * t398 - t342 * t312 - t341 * t313 - t421 * t456) * t329) * t318, t430 * t357 * t450 + (t430 * t323 + ((-qJD(5) * t328 - 0.2e1 * t447) * t402 + (t312 * t402 + (t313 - t452) * t398) * t329) * t357) * t318, t450 + 0.2e1 * (t312 * t329 * t318 + (-t318 * t478 - t329 * t479) * t331) * t331, 0;];
JaD_rot  = t1;
