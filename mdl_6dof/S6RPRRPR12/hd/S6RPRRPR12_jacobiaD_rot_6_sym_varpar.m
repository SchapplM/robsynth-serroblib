% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPR12_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR12_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:07:17
% EndTime: 2019-02-26 21:07:20
% DurationCPUTime: 3.17s
% Computational Cost: add. (13786->169), mult. (41237->336), div. (726->12), fcn. (53008->17), ass. (0->149)
t402 = sin(pkin(12));
t407 = cos(pkin(6));
t415 = cos(qJ(1));
t405 = cos(pkin(12));
t411 = sin(qJ(1));
t463 = t411 * t405;
t432 = t402 * t415 + t407 * t463;
t390 = t432 * qJD(1);
t473 = t402 * t411;
t395 = t405 * t415 - t407 * t473;
t391 = t395 * qJD(1);
t403 = sin(pkin(7));
t406 = cos(pkin(7));
t410 = sin(qJ(3));
t414 = cos(qJ(3));
t404 = sin(pkin(6));
t470 = t404 * t411;
t447 = qJD(1) * t470;
t465 = t407 * t415;
t394 = t402 * t465 + t463;
t433 = t405 * t465 - t473;
t430 = t433 * t406;
t469 = t404 * t415;
t452 = t403 * t469;
t435 = -t430 + t452;
t496 = t394 * t410 + t414 * t435;
t348 = t496 * qJD(3) - (-t390 * t406 + t403 * t447) * t410 - t391 * t414;
t475 = t394 * t414;
t506 = t410 * t435;
t372 = -t475 + t506;
t384 = t403 * t433 + t406 * t469;
t409 = sin(qJ(4));
t413 = cos(qJ(4));
t360 = t372 * t413 + t384 * t409;
t379 = t390 * t403 + t406 * t447;
t328 = qJD(4) * t360 + t348 * t409 + t379 * t413;
t423 = t372 * t409;
t359 = -t384 * t413 + t423;
t498 = -t403 * t470 + t432 * t406;
t373 = t395 * t410 + t498 * t414;
t505 = t379 * t409;
t374 = t395 * t414 - t410 * t498;
t386 = t403 * t432 + t406 * t470;
t361 = t374 * t409 - t386 * t413;
t408 = sin(qJ(6));
t412 = cos(qJ(6));
t440 = t361 * t412 - t373 * t408;
t501 = qJD(6) * t440;
t467 = t406 * t410;
t471 = t403 * t407;
t383 = t404 * (t402 * t414 + t405 * t467) + t410 * t471;
t472 = t403 * t404;
t392 = -t405 * t472 + t406 * t407;
t368 = -t383 * t409 + t392 * t413;
t466 = t406 * t414;
t382 = t414 * t471 + (-t402 * t410 + t405 * t466) * t404;
t375 = t382 * qJD(3);
t352 = qJD(4) * t368 + t375 * t413;
t369 = t383 * t413 + t392 * t409;
t366 = 0.1e1 / t369 ^ 2;
t500 = t352 * t366;
t365 = 0.1e1 / t369;
t481 = t360 * t366;
t436 = t365 * t496 - t382 * t481;
t499 = t413 * t436;
t342 = atan2(t360, t369);
t337 = sin(t342);
t338 = cos(t342);
t320 = t337 * t360 + t338 * t369;
t317 = 0.1e1 / t320;
t479 = t373 * t412;
t336 = t361 * t408 + t479;
t330 = 0.1e1 / t336;
t318 = 0.1e1 / t320 ^ 2;
t331 = 0.1e1 / t336 ^ 2;
t495 = 0.2e1 * t360;
t362 = t374 * t413 + t386 * t409;
t494 = 0.2e1 * t362;
t355 = t362 ^ 2;
t314 = t318 * t355 + 0.1e1;
t389 = t394 * qJD(1);
t344 = qJD(1) * t506 - t373 * qJD(3) - t389 * t414;
t377 = t384 * qJD(1);
t460 = qJD(4) * t409;
t325 = t377 * t409 - t374 * t460 + (qJD(4) * t386 + t344) * t413;
t487 = t325 * t318;
t354 = t360 ^ 2;
t341 = t354 * t366 + 0.1e1;
t339 = 0.1e1 / t341;
t327 = t505 + qJD(4) * t423 + (-qJD(4) * t384 - t348) * t413;
t439 = -t327 * t365 - t352 * t481;
t308 = t439 * t339;
t442 = -t337 * t369 + t338 * t360;
t303 = t308 * t442 - t337 * t327 + t338 * t352;
t319 = t317 * t318;
t492 = t303 * t319;
t493 = (-t355 * t492 + t362 * t487) / t314 ^ 2;
t324 = qJD(4) * t362 + t344 * t409 - t377 * t413;
t461 = qJD(1) * t414;
t444 = t461 * t472;
t343 = qJD(3) * t374 - t389 * t410 - t415 * t444 + t430 * t461;
t315 = qJD(6) * t336 - t324 * t412 + t343 * t408;
t329 = t440 ^ 2;
t323 = t329 * t331 + 0.1e1;
t486 = t331 * t440;
t316 = t324 * t408 + t343 * t412 + t501;
t489 = t316 * t330 * t331;
t491 = (-t315 * t486 - t329 * t489) / t323 ^ 2;
t483 = t365 * t500;
t490 = (-t327 * t481 - t354 * t483) / t341 ^ 2;
t488 = t318 * t362;
t485 = t337 * t362;
t484 = t338 * t362;
t482 = t360 * t365;
t480 = t373 * t409;
t478 = t373 * t413;
t464 = t408 * t440;
t462 = t412 * t330;
t459 = 0.2e1 * t493;
t458 = 0.2e1 * t491;
t457 = -0.2e1 * t490;
t456 = t319 * t494;
t455 = t365 * t490;
t454 = t318 * t485;
t453 = t318 * t484;
t446 = -0.2e1 * t440 * t489;
t445 = t483 * t495;
t443 = -qJD(6) * t480 + t344;
t441 = t359 * t412 + t408 * t496;
t334 = t359 * t408 - t412 * t496;
t438 = -t331 * t464 + t462;
t437 = -t359 * t365 - t368 * t481;
t428 = -t337 + (-t338 * t482 + t337) * t339;
t427 = qJD(4) * t478 + qJD(6) * t374 + t343 * t409;
t422 = -t390 * t466 - t391 * t410 + t411 * t444 + (t410 * t452 - t433 * t467 - t475) * qJD(3);
t376 = t383 * qJD(3);
t351 = -qJD(4) * t369 - t375 * t409;
t350 = t374 * t412 - t408 * t480;
t349 = t374 * t408 + t409 * t479;
t321 = 0.1e1 / t323;
t312 = 0.1e1 / t314;
t311 = t339 * t499;
t309 = t437 * t339;
t307 = t428 * t362;
t305 = (t337 * t496 + t338 * t382) * t413 + t442 * t311;
t304 = t309 * t442 - t337 * t359 + t338 * t368;
t301 = t437 * t457 + (t368 * t445 - t328 * t365 + (t327 * t368 - t351 * t360 + t352 * t359) * t366) * t339;
t300 = t457 * t499 + (-t436 * t460 + (t382 * t445 - t422 * t365 + (t327 * t382 - t352 * t496 + t360 * t376) * t366) * t413) * t339;
t1 = [t455 * t494 + (-t325 * t365 + t362 * t500) * t339, 0, t300, t301, 0, 0; -0.2e1 * t360 * t317 * t493 + ((-qJD(4) * t359 + t348 * t413 - t505) * t317 + (-t303 * t360 - t307 * t325) * t318) * t312 + (t307 * t318 * t459 + (0.2e1 * t307 * t492 - (t308 * t339 * t482 + t457) * t454 - (t455 * t495 - t308 + (t308 - t439) * t339) * t453 - t428 * t487) * t312) * t362, 0 (t305 * t488 + t317 * t478) * t459 + (-t305 * t487 + (-t343 * t413 + t373 * t460) * t317 + (t305 * t456 + t318 * t478) * t303 - (-t382 * t460 + t300 * t360 - t311 * t327 - t376 * t413 + (-t311 * t369 + t413 * t496) * t308) * t453 - (-t496 * t460 - t300 * t369 - t311 * t352 - t422 * t413 + (-t311 * t360 - t382 * t413) * t308) * t454) * t312 (t304 * t488 + t317 * t361) * t459 + (t304 * t303 * t456 - t324 * t317 + (t361 * t303 - t304 * t325 - (t301 * t360 - t309 * t327 + t351 + (-t309 * t369 - t359) * t308) * t484 - (-t301 * t369 - t309 * t352 - t328 + (-t309 * t360 - t368) * t308) * t485) * t318) * t312, 0, 0; (t330 * t441 - t334 * t486) * t458 + ((qJD(6) * t334 - t328 * t412 + t408 * t422) * t330 + t334 * t446 + (t441 * t316 + (qJD(6) * t441 + t328 * t408 + t412 * t422) * t440 - t334 * t315) * t331) * t321, 0 (-t330 * t349 - t350 * t486) * t458 + (t350 * t446 + t443 * t330 * t408 + t427 * t462 + (t412 * t440 * t443 - t350 * t315 - t349 * t316 - t427 * t464) * t331) * t321, t438 * t362 * t458 + (-t438 * t325 + ((qJD(6) * t330 + t446) * t408 + (-t315 * t408 + (t316 + t501) * t412) * t331) * t362) * t321, 0, -0.2e1 * t491 - 0.2e1 * (t315 * t331 * t321 - (-t321 * t489 - t331 * t491) * t440) * t440;];
JaD_rot  = t1;
