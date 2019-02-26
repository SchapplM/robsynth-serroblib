% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% JaD_rot [3x7]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S7RRRRRRR1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobiaD_rot_6_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_jacobiaD_rot_6_sym_varpar: qJD has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobiaD_rot_6_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:54:22
% EndTime: 2019-02-26 22:54:27
% DurationCPUTime: 4.96s
% Computational Cost: add. (15691->282), mult. (45225->517), div. (1240->12), fcn. (56295->15), ass. (0->213)
t448 = sin(qJ(3));
t454 = cos(qJ(3));
t455 = cos(qJ(2));
t504 = qJD(1) * t455 + qJD(3);
t449 = sin(qJ(2));
t450 = sin(qJ(1));
t559 = t449 * t450;
t519 = qJD(2) * t559;
t456 = cos(qJ(1));
t551 = t454 * t456;
t542 = qJD(3) * t455;
t549 = qJD(1) * t450;
t596 = -t450 * t542 - t549;
t398 = t596 * t448 - t454 * t519 + t504 * t551;
t552 = t454 * t455;
t560 = t448 * t456;
t427 = t450 * t552 + t560;
t447 = sin(qJ(4));
t453 = cos(qJ(4));
t540 = qJD(4) * t449;
t547 = qJD(2) * t455;
t548 = qJD(1) * t456;
t595 = t449 * t548 + t450 * t547;
t357 = (-qJD(4) * t427 + t595) * t447 + (t450 * t540 + t398) * t453;
t412 = t427 * t453 + t447 * t559;
t446 = sin(qJ(5));
t452 = cos(qJ(5));
t555 = t450 * t448;
t506 = -t455 * t555 + t551;
t385 = t412 * t446 - t506 * t452;
t550 = t455 * t456;
t523 = t448 * t550;
t545 = qJD(3) * t448;
t473 = -qJD(1) * t523 + t448 * t519 + t596 * t454 - t456 * t545;
t343 = t385 * qJD(5) - t357 * t452 - t473 * t446;
t495 = t506 * t446;
t387 = t412 * t452 + t495;
t341 = t387 * qJD(5) + t357 * t446 - t473 * t452;
t383 = t385 ^ 2;
t557 = t449 * t454;
t426 = -t447 * t455 + t453 * t557;
t561 = t448 * t452;
t522 = t449 * t561;
t409 = t426 * t446 + t522;
t407 = 0.1e1 / t409 ^ 2;
t375 = t383 * t407 + 0.1e1;
t373 = 0.1e1 / t375;
t543 = qJD(3) * t454;
t515 = t449 * t543;
t469 = qJD(5) * t426 + t448 * t547 + t515;
t502 = qJD(2) * t454 - qJD(4);
t553 = t453 * t455;
t544 = qJD(3) * t453;
t562 = t448 * t449;
t503 = t454 * qJD(4) - qJD(2);
t587 = t503 * t447;
t599 = qJD(5) * t562 + (t448 * t544 + t587) * t449;
t497 = t502 * t553 - t599;
t352 = t497 * t446 + t469 * t452;
t406 = 0.1e1 / t409;
t567 = t385 * t407;
t489 = -t341 * t406 + t352 * t567;
t324 = t489 * t373;
t376 = atan2(-t385, t409);
t368 = sin(t376);
t369 = cos(t376);
t494 = -t368 * t409 - t369 * t385;
t317 = t494 * t324 - t341 * t368 + t352 * t369;
t338 = -t368 * t385 + t369 * t409;
t336 = 0.1e1 / t338 ^ 2;
t604 = t317 * t336;
t429 = t454 * t550 - t555;
t556 = t449 * t456;
t417 = t429 * t453 + t447 * t556;
t481 = t450 * t454 + t523;
t476 = t481 * t452;
t390 = t417 * t446 + t476;
t581 = 0.2e1 * t390;
t335 = 0.1e1 / t338;
t597 = t335 * t604;
t510 = t581 * t597;
t546 = qJD(2) * t456;
t471 = t449 * t546 + t504 * t450;
t505 = qJD(1) + t542;
t397 = t471 * t454 + t505 * t560;
t480 = t449 * t549 - t455 * t546;
t355 = (t456 * t540 - t397) * t453 + (-qJD(4) * t429 - t480) * t447;
t477 = t481 * t446;
t391 = t417 * t452 - t477;
t396 = t471 * t448 - t505 * t551;
t339 = t391 * qJD(5) + t355 * t446 - t396 * t452;
t574 = t339 * t336;
t603 = -t574 + t510;
t501 = qJD(5) + t544;
t600 = (t501 * t448 + t587) * t446 - t452 * t543;
t416 = t429 * t447 - t453 * t556;
t445 = sin(qJ(6));
t451 = cos(qJ(6));
t366 = t391 * t445 - t416 * t451;
t598 = 0.2e1 * t366;
t592 = t352 * t407;
t482 = t447 * t557 + t553;
t558 = t449 * t453;
t507 = -t427 * t447 + t450 * t558;
t484 = -t406 * t507 - t482 * t567;
t591 = t446 * t484;
t384 = t390 ^ 2;
t334 = t336 * t384 + 0.1e1;
t332 = 0.1e1 / t334;
t580 = (-t384 * t597 + t390 * t574) / t334 ^ 2;
t585 = -t332 * t604 - 0.2e1 * t335 * t580;
t536 = 0.2e1 * t580;
t575 = t336 * t390;
t584 = t332 * t603 + t536 * t575;
t421 = t426 * t456;
t583 = -qJD(5) * t421 + t480 * t448 - t456 * t515;
t367 = t391 * t451 + t416 * t445;
t361 = 0.1e1 / t367;
t362 = 0.1e1 / t367 ^ 2;
t582 = -0.2e1 * t385;
t340 = -t390 * qJD(5) + t355 * t452 + t396 * t446;
t354 = t417 * qJD(4) - t397 * t447 + t480 * t453;
t326 = t367 * qJD(6) + t340 * t445 - t354 * t451;
t360 = t366 ^ 2;
t346 = t360 * t362 + 0.1e1;
t571 = t362 * t366;
t537 = qJD(6) * t366;
t327 = t340 * t451 + t354 * t445 - t537;
t577 = t327 * t361 * t362;
t579 = (t326 * t571 - t360 * t577) / t346 ^ 2;
t572 = t406 * t592;
t578 = (t341 * t567 - t383 * t572) / t375 ^ 2;
t576 = t332 * t335;
t344 = 0.1e1 / t346;
t573 = t344 * t362;
t570 = t368 * t390;
t569 = t369 * t390;
t568 = t385 * t406;
t566 = t416 * t446;
t565 = t416 * t452;
t564 = t445 * t361;
t563 = t446 * t453;
t554 = t451 * t366;
t541 = qJD(4) * t447;
t539 = qJD(5) * t452;
t538 = qJD(5) * t453;
t535 = -0.2e1 * t579;
t534 = 0.2e1 * t579;
t533 = -0.2e1 * t578;
t530 = t362 * t579;
t529 = t406 * t578;
t528 = t326 * t573;
t526 = t332 * t575;
t525 = t366 * t577;
t524 = t446 * t562;
t509 = t572 * t582;
t508 = 0.2e1 * t525;
t496 = qJD(6) * t565 + t355;
t365 = -t387 * t451 + t445 * t507;
t364 = -t387 * t445 - t451 * t507;
t404 = -t421 * t452 + t456 * t524;
t420 = t482 * t456;
t382 = t404 * t451 - t420 * t445;
t381 = t404 * t445 + t420 * t451;
t491 = t502 * t455;
t490 = t426 * t549 + (-t453 * t491 + t599) * t456;
t488 = t362 * t554 - t564;
t410 = t426 * t452 - t524;
t487 = -t387 * t406 + t410 * t567;
t399 = t427 * t452 + t453 * t495;
t483 = -t448 * t563 + t452 * t454;
t418 = t483 * t449;
t486 = -t399 * t406 + t418 * t567;
t419 = t426 * t450;
t402 = -t419 * t446 - t450 * t522;
t428 = t447 * t449 + t453 * t552;
t415 = t428 * t446 + t455 * t561;
t485 = -t402 * t406 + t415 * t567;
t478 = t447 * t481;
t475 = -t368 + (t369 * t568 + t368) * t373;
t474 = qJD(4) * t481;
t472 = qJD(5) * t566 + qJD(6) * t417 - t354 * t452;
t464 = -t429 * qJD(5) + t396 * t453 + t447 * t474;
t467 = -t481 * t538 - t397;
t468 = -qJD(6) * t478 - t467 * t446 + t464 * t452;
t466 = -qJD(2) * t561 - t502 * t563;
t465 = -t412 * qJD(4) - t398 * t447 + t595 * t453;
t401 = -t429 * t446 - t453 * t476;
t463 = t401 * qJD(6) - t396 * t447 + t453 * t474;
t394 = -t503 * t558 + (t449 * t545 - t491) * t447;
t403 = -t421 * t446 - t456 * t522;
t380 = t401 * t451 - t445 * t478;
t378 = t417 * t445 - t451 * t565;
t377 = -t417 * t451 - t445 * t565;
t371 = t394 * t456 + t482 * t549;
t370 = t483 * t547 + ((-qJD(3) - t538) * t561 + (t448 * t541 - t501 * t454) * t446) * t449;
t353 = -t469 * t446 + t497 * t452;
t351 = t428 * t539 + t466 * t449 - t455 * t600;
t350 = -t419 * t539 - t409 * t548 + (t449 * t600 + t466 * t455) * t450;
t349 = -t583 * t446 + t490 * t452;
t348 = (t506 * t538 + t398) * t452 + (-t427 * qJD(5) + t473 * t453 - t506 * t541) * t446;
t331 = t373 * t591;
t330 = t486 * t373;
t329 = t485 * t373;
t328 = t487 * t373;
t322 = (-t368 * t507 - t369 * t482) * t446 + t494 * t331;
t320 = t494 * t329 - t368 * t402 + t369 * t415;
t318 = t494 * t328 - t368 * t387 + t369 * t410;
t316 = t486 * t533 + (t418 * t509 - t348 * t406 + (t341 * t418 + t352 * t399 + t370 * t385) * t407) * t373;
t314 = t485 * t533 + (t415 * t509 - t350 * t406 + (t341 * t415 + t351 * t385 + t352 * t402) * t407) * t373;
t313 = t487 * t533 + (t410 * t509 + t343 * t406 + (t341 * t410 + t352 * t387 + t353 * t385) * t407) * t373;
t312 = t533 * t591 + (t484 * t539 + (-t482 * t509 - t465 * t406 + (-t341 * t482 + t352 * t507 + t385 * t394) * t407) * t446) * t373;
t1 = [t529 * t581 + (-t339 * t406 + t390 * t592) * t373, t314, t316, t312, t313, 0, 0; -t341 * t576 - (t475 * t339 + ((-t324 * t373 * t568 + t533) * t368 + (t529 * t582 - t324 + (t324 - t489) * t373) * t369) * t390) * t526 - t585 * t385 + t584 * t475 * t390 (t320 * t575 - t335 * t403) * t536 + (t320 * t510 + (-t403 * t317 - t320 * t339 - (-t314 * t385 - t329 * t341 + t351 + (-t329 * t409 - t402) * t324) * t569 - (-t314 * t409 - t329 * t352 - t350 + (t329 * t385 - t415) * t324) * t570) * t336 + (t490 * t446 + t583 * t452) * t335) * t332 (t464 * t446 + t467 * t452) * t576 - ((-t316 * t385 - t330 * t341 + t370 + (-t330 * t409 - t399) * t324) * t369 + (-t316 * t409 - t330 * t352 - t348 + (t330 * t385 - t418) * t324) * t368) * t526 + t585 * (t429 * t452 - t453 * t477) + t584 * (t494 * t330 - t368 * t399 + t369 * t418) (t322 * t575 + t335 * t566) * t536 + ((-t354 * t446 - t416 * t539) * t335 + t603 * t322 + (t566 * t317 - (-t482 * t539 - t312 * t385 - t331 * t341 + t394 * t446 + (-t331 * t409 - t446 * t507) * t324) * t569 - (-t507 * t539 - t312 * t409 - t331 * t352 - t465 * t446 + (t331 * t385 + t446 * t482) * t324) * t570) * t336) * t332 (t318 * t575 - t335 * t391) * t536 + (t318 * t510 + t340 * t335 + (-t391 * t317 - t318 * t339 - (-t313 * t385 - t328 * t341 + t353 + (-t328 * t409 - t387) * t324) * t569 - (-t313 * t409 - t328 * t352 + t343 + (t328 * t385 - t410) * t324) * t570) * t336) * t332, 0, 0; (-t361 * t364 + t365 * t571) * t534 + ((t365 * qJD(6) + t343 * t445 - t451 * t465) * t361 + t365 * t508 + (-t364 * t327 - (-t364 * qJD(6) + t343 * t451 + t445 * t465) * t366 - t365 * t326) * t362) * t344 (-t361 * t381 + t382 * t571) * t534 + ((t382 * qJD(6) + t349 * t445 - t371 * t451) * t361 + t382 * t508 + (-t381 * t327 - (-t381 * qJD(6) + t349 * t451 + t371 * t445) * t366 - t382 * t326) * t362) * t344 (t530 * t598 - t528) * t380 + (-t327 * t573 + t361 * t535) * (t401 * t445 + t451 * t478) + ((t468 * t445 + t463 * t451) * t361 - (-t463 * t445 + t468 * t451) * t571 + t380 * t508) * t344 (-t361 * t377 + t378 * t571) * t534 + (t378 * t508 - t496 * t361 * t451 + t472 * t564 + (-t496 * t366 * t445 - t378 * t326 - t377 * t327 - t472 * t554) * t362) * t344, t488 * t390 * t535 + (t488 * t339 + ((-qJD(6) * t361 - 0.2e1 * t525) * t451 + (t326 * t451 + (t327 - t537) * t445) * t362) * t390) * t344, t535 + (t528 + (-t344 * t577 - t530) * t366) * t598, 0;];
JaD_rot  = t1;
