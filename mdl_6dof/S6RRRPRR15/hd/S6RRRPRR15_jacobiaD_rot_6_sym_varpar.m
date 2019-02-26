% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR15_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR15_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:24:16
% EndTime: 2019-02-26 22:24:21
% DurationCPUTime: 5.06s
% Computational Cost: add. (19869->249), mult. (59163->466), div. (983->12), fcn. (74815->17), ass. (0->203)
t462 = sin(qJ(3));
t572 = cos(pkin(6));
t573 = sin(qJ(2));
t510 = t572 * t573;
t465 = cos(qJ(2));
t574 = sin(qJ(1));
t527 = t574 * t465;
t576 = cos(qJ(1));
t485 = t510 * t576 + t527;
t575 = cos(qJ(3));
t478 = t485 * t575;
t459 = cos(pkin(7));
t456 = t574 * t573;
t511 = t572 * t576;
t494 = -t465 * t511 + t456;
t487 = t494 * t459;
t457 = sin(pkin(7));
t458 = sin(pkin(6));
t531 = t458 * t576;
t518 = t457 * t531;
t428 = -t478 + (t487 + t518) * t462;
t484 = t527 * t572 + t573 * t576;
t440 = qJD(1) * t484 + qJD(2) * t485;
t495 = t574 * t510;
t524 = t576 * qJD(1);
t441 = -qJD(1) * t495 - qJD(2) * t456 + (qJD(2) * t511 + t524) * t465;
t530 = t458 * t574;
t516 = t457 * t530;
t452 = t575 * t516;
t529 = t459 * t575;
t379 = qJD(1) * t452 + qJD(3) * t428 - t440 * t529 - t441 * t462;
t513 = qJD(1) * t530;
t432 = t440 * t457 + t459 * t513;
t461 = sin(qJ(5));
t464 = cos(qJ(5));
t591 = -t379 * t461 + t432 * t464;
t590 = -t379 * t464 - t432 * t461;
t445 = t457 * t494 - t459 * t531;
t589 = t445 * t461;
t588 = t445 * t464;
t446 = t457 * t484 + t459 * t530;
t479 = t484 * t575;
t483 = -t465 * t576 + t495;
t508 = t459 * t479 - t462 * t483 - t452;
t412 = t446 * t464 + t461 * t508;
t460 = sin(qJ(6));
t429 = (-t459 * t484 + t516) * t462 - t483 * t575;
t463 = cos(qJ(6));
t556 = t429 * t463;
t387 = t412 * t460 - t556;
t587 = 0.2e1 * t387;
t481 = t485 * t462;
t496 = t575 * t518;
t584 = t481 + t496;
t486 = t494 * t575;
t473 = t459 * t486 + t584;
t408 = -t464 * t473 + t589;
t526 = t573 * t462;
t528 = t575 * t465;
t491 = t459 * t528 - t526;
t523 = t457 * t572;
t503 = t575 * t523;
t442 = -t458 * t491 - t503;
t554 = t457 * t458;
t534 = t465 * t554;
t449 = t459 * t572 - t534;
t505 = t442 * t464 - t449 * t461;
t394 = atan2(-t408, -t505);
t389 = sin(t394);
t390 = cos(t394);
t364 = -t389 * t408 - t390 * t505;
t362 = 0.1e1 / t364 ^ 2;
t411 = t446 * t461 - t464 * t508;
t405 = t411 ^ 2;
t360 = t362 * t405 + 0.1e1;
t438 = qJD(1) * t494 + qJD(2) * t483;
t514 = t458 * t524;
t430 = -t438 * t457 + t459 * t514;
t439 = qJD(1) * t485 + qJD(2) * t484;
t474 = -qJD(1) * t496 + qJD(3) * t429 - t438 * t529 - t439 * t462;
t368 = qJD(5) * t412 + t430 * t461 - t464 * t474;
t565 = t362 * t411;
t404 = t408 ^ 2;
t422 = 0.1e1 / t505 ^ 2;
t393 = t404 * t422 + 0.1e1;
t391 = 0.1e1 / t393;
t410 = t461 * t473 + t588;
t371 = qJD(5) * t410 - t590;
t515 = t573 * t575;
t550 = t462 * t465;
t489 = t459 * t515 + t550;
t492 = t459 * t550 + t515;
t512 = t462 * t523;
t416 = qJD(3) * t512 + (qJD(2) * t489 + qJD(3) * t492) * t458;
t425 = t442 * t461 + t449 * t464;
t519 = t573 * t554;
t504 = qJD(2) * t519;
t395 = qJD(5) * t425 - t416 * t464 + t461 * t504;
t421 = 0.1e1 / t505;
t558 = t408 * t422;
t501 = t371 * t421 + t395 * t558;
t351 = t501 * t391;
t507 = t389 * t505 - t390 * t408;
t345 = t351 * t507 - t371 * t389 + t390 * t395;
t361 = 0.1e1 / t364;
t363 = t361 * t362;
t570 = t345 * t363;
t546 = 0.2e1 * (t368 * t565 - t405 * t570) / t360 ^ 2;
t582 = t395 * t422;
t443 = t458 * t492 + t512;
t497 = -t421 * t428 + t443 * t558;
t581 = t464 * t497;
t557 = t429 * t461;
t579 = qJD(6) * t557 + t474;
t578 = -(-t440 * t459 + t457 * t513) * t462 - t441 * t575;
t388 = t412 * t463 + t429 * t460;
t382 = 0.1e1 / t388;
t383 = 0.1e1 / t388 ^ 2;
t577 = 0.2e1 * t411;
t369 = -qJD(5) * t411 + t430 * t464 + t461 * t474;
t377 = -t439 * t575 + (t438 * t459 + t457 * t514) * t462 - t508 * qJD(3);
t354 = qJD(6) * t388 + t369 * t460 - t377 * t463;
t381 = t387 ^ 2;
t367 = t381 * t383 + 0.1e1;
t563 = t383 * t387;
t547 = qJD(6) * t387;
t355 = t369 * t463 + t377 * t460 - t547;
t567 = t355 * t382 * t383;
t569 = (t354 * t563 - t381 * t567) / t367 ^ 2;
t560 = t421 * t582;
t568 = (t371 * t558 + t404 * t560) / t393 ^ 2;
t566 = t362 * t368;
t365 = 0.1e1 / t367;
t564 = t365 * t383;
t562 = t389 * t411;
t561 = t390 * t411;
t559 = t408 * t421;
t555 = t429 * t464;
t553 = t457 * t461;
t552 = t457 * t464;
t551 = t459 * t462;
t549 = qJD(5) * t461;
t548 = qJD(5) * t464;
t545 = -0.2e1 * t569;
t544 = 0.2e1 * t569;
t543 = -0.2e1 * t568;
t542 = t363 * t577;
t541 = t383 * t569;
t540 = t421 * t568;
t539 = t354 * t564;
t538 = t387 * t567;
t537 = t362 * t562;
t536 = t362 * t561;
t535 = t408 * t560;
t522 = t345 * t542;
t521 = 0.2e1 * t538;
t520 = 0.2e1 * t535;
t427 = -t487 * t575 - t584;
t407 = t427 * t461 - t588;
t386 = t407 * t463 + t428 * t460;
t385 = t407 * t460 - t428 * t463;
t435 = -t462 * t484 - t483 * t529;
t415 = t435 * t461 - t483 * t552;
t436 = t483 * t551 - t479;
t402 = t415 * t463 + t436 * t460;
t401 = t415 * t460 - t436 * t463;
t506 = t427 * t464 + t589;
t500 = -t382 * t460 + t463 * t563;
t499 = t410 * t421 + t425 * t558;
t434 = t459 * t478 - t462 * t494;
t482 = t457 * t485;
t413 = t434 * t464 - t461 * t482;
t448 = t489 * t458;
t437 = -t448 * t464 + t461 * t519;
t498 = -t413 * t421 + t437 * t558;
t414 = -t435 * t464 - t483 * t553;
t493 = -t389 + (-t390 * t559 + t389) * t391;
t490 = -t459 * t526 + t528;
t480 = -qJD(6) * t508 + t377 * t461 + t429 * t548;
t417 = qJD(3) * t503 + (qJD(2) * t490 + qJD(3) * t491) * t458;
t403 = qJD(2) * t461 * t534 + t519 * t548 - (qJD(2) * t491 + qJD(3) * t490) * t458 * t464 + t448 * t549;
t400 = -t460 * t508 + t461 * t556;
t398 = -qJD(3) * t435 + t438 * t575 + t439 * t551;
t397 = qJD(3) * t436 + t438 * t462 - t439 * t529;
t396 = qJD(5) * t505 + t416 * t461 + t464 * t504;
t380 = -qJD(3) * t427 + t578;
t378 = -qJD(3) * t473 - t578;
t374 = -t441 * t553 - t482 * t548 + (t441 * t529 - t440 * t462 + (-t459 * t481 - t486) * qJD(3)) * t464 - t434 * t549;
t373 = -qJD(5) * t414 + t397 * t461 - t439 * t552;
t372 = -qJD(5) * t408 + t591;
t370 = qJD(5) * t506 - t591;
t358 = 0.1e1 / t360;
t357 = t391 * t581;
t356 = t498 * t391;
t353 = t499 * t391;
t350 = t493 * t411;
t348 = (-t389 * t428 - t390 * t443) * t464 - t507 * t357;
t347 = t356 * t507 + t389 * t413 + t390 * t437;
t346 = t353 * t507 - t389 * t410 + t390 * t425;
t344 = t498 * t543 + (t437 * t520 - t374 * t421 + (t371 * t437 - t395 * t413 + t403 * t408) * t422) * t391;
t342 = t499 * t543 + (t425 * t520 + t372 * t421 + (t371 * t425 + t395 * t410 + t396 * t408) * t422) * t391;
t341 = 0.2e1 * t568 * t581 + (t497 * t549 + (-0.2e1 * t443 * t535 - t378 * t421 + (-t371 * t443 + t395 * t428 - t408 * t417) * t422) * t464) * t391;
t1 = [-t540 * t577 + (t368 * t421 + t411 * t582) * t391, t344, t341, 0, t342, 0; t506 * t361 * t546 + ((qJD(5) * t407 + t590) * t361 + (t345 * t506 - t350 * t368) * t362) * t358 + (t350 * t362 * t546 + (0.2e1 * t350 * t570 - (t351 * t391 * t559 + t543) * t537 - (0.2e1 * t408 * t540 - t351 + (t351 - t501) * t391) * t536 - t493 * t566) * t358) * t411 (t347 * t565 - t361 * t414) * t546 + ((qJD(5) * t415 - t397 * t464 - t439 * t553) * t361 + t347 * t522 + (-t414 * t345 - t347 * t368 - (-t344 * t408 - t356 * t371 + t403 + (t356 * t505 + t413) * t351) * t561 - (t344 * t505 - t356 * t395 + t374 + (t356 * t408 - t437) * t351) * t562) * t362) * t358 (t348 * t565 + t361 * t555) * t546 + (-t348 * t566 + (-t377 * t464 + t429 * t549) * t361 + (t348 * t542 + t362 * t555) * t345 - (t443 * t549 - t341 * t408 + t357 * t371 - t417 * t464 + (-t357 * t505 - t428 * t464) * t351) * t536 - (t428 * t549 + t341 * t505 + t357 * t395 + t378 * t464 + (-t357 * t408 + t443 * t464) * t351) * t537) * t358, 0 (t346 * t565 - t361 * t412) * t546 + (t346 * t522 + t369 * t361 + (-t412 * t345 - t346 * t368 - (-t342 * t408 - t353 * t371 + t396 + (t353 * t505 - t410) * t351) * t561 - (t342 * t505 - t353 * t395 - t372 + (t353 * t408 - t425) * t351) * t562) * t362) * t358, 0; (-t382 * t385 + t386 * t563) * t544 + ((qJD(6) * t386 + t370 * t460 - t380 * t463) * t382 + t386 * t521 + (-t385 * t355 - (-qJD(6) * t385 + t370 * t463 + t380 * t460) * t387 - t386 * t354) * t383) * t365 (-t382 * t401 + t402 * t563) * t544 + ((qJD(6) * t402 + t373 * t460 - t398 * t463) * t382 + t402 * t521 + (-t401 * t355 - (-qJD(6) * t401 + t373 * t463 + t398 * t460) * t387 - t402 * t354) * t383) * t365 (t541 * t587 - t539) * t400 + (-t355 * t564 + t382 * t545) * (t460 * t557 + t463 * t508) + ((t480 * t460 + t463 * t579) * t382 - (-t460 * t579 + t480 * t463) * t563 + t400 * t521) * t365, 0, t500 * t411 * t545 + (t500 * t368 + ((-qJD(6) * t382 - 0.2e1 * t538) * t463 + (t354 * t463 + (t355 - t547) * t460) * t383) * t411) * t365, t545 + (t539 + (-t365 * t567 - t541) * t387) * t587;];
JaD_rot  = t1;
