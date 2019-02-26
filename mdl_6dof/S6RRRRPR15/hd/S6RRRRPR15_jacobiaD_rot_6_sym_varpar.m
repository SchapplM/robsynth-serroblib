% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR15_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR15_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:48
% EndTime: 2019-02-26 22:38:53
% DurationCPUTime: 4.68s
% Computational Cost: add. (19869->247), mult. (59163->462), div. (983->12), fcn. (74815->17), ass. (0->188)
t461 = cos(pkin(6));
t469 = cos(qJ(2));
t564 = sin(qJ(1));
t515 = t564 * t469;
t465 = sin(qJ(2));
t470 = cos(qJ(1));
t529 = t470 * t465;
t448 = t461 * t529 + t515;
t486 = t461 * t515 + t529;
t434 = qJD(1) * t486 + qJD(2) * t448;
t516 = t564 * t465;
t508 = t461 * t516;
t528 = t470 * t469;
t435 = -qJD(1) * t508 - qJD(2) * t516 + (qJD(2) * t461 + qJD(1)) * t528;
t460 = cos(pkin(7));
t464 = sin(qJ(3));
t468 = cos(qJ(3));
t458 = sin(pkin(7));
t459 = sin(pkin(6));
t517 = t459 * t564;
t505 = qJD(1) * t517;
t499 = t458 * t505;
t447 = -t461 * t528 + t516;
t538 = t459 * t470;
t492 = t447 * t460 + t458 * t538;
t567 = t448 * t464 + t468 * t492;
t373 = qJD(3) * t567 - (-t434 * t460 + t499) * t464 - t435 * t468;
t545 = t448 * t468;
t421 = t464 * t492 - t545;
t440 = -t447 * t458 + t460 * t538;
t463 = sin(qJ(4));
t467 = cos(qJ(4));
t404 = t421 * t467 + t440 * t463;
t426 = t434 * t458 + t460 * t505;
t365 = qJD(4) * t404 + t373 * t463 + t426 * t467;
t479 = t421 * t463;
t403 = -t440 * t467 + t479;
t575 = t426 * t463;
t509 = t458 * t517;
t478 = -t460 * t486 + t509;
t487 = t508 - t528;
t423 = t464 * t478 - t468 * t487;
t442 = t458 * t486 + t460 * t517;
t405 = t423 * t463 - t442 * t467;
t483 = t486 * t468;
t544 = t487 * t464;
t422 = t460 * t483 - t468 * t509 - t544;
t462 = sin(qJ(6));
t466 = cos(qJ(6));
t501 = t405 * t466 - t422 * t462;
t571 = t501 * qJD(6);
t530 = t468 * t469;
t534 = t464 * t465;
t488 = -t460 * t534 + t530;
t491 = t460 * t530 - t534;
t518 = t461 * t458 * t468;
t411 = qJD(3) * t518 + (qJD(2) * t488 + qJD(3) * t491) * t459;
t532 = t465 * t468;
t533 = t464 * t469;
t489 = t460 * t533 + t532;
t542 = t458 * t464;
t439 = t459 * t489 + t461 * t542;
t539 = t458 * t469;
t446 = -t459 * t539 + t460 * t461;
t417 = -t439 * t463 + t446 * t467;
t541 = t458 * t465;
t519 = t459 * t541;
t506 = qJD(2) * t519;
t389 = qJD(4) * t417 + t411 * t467 + t463 * t506;
t418 = t439 * t467 + t446 * t463;
t415 = 0.1e1 / t418 ^ 2;
t570 = t389 * t415;
t414 = 0.1e1 / t418;
t438 = t459 * t491 + t518;
t551 = t404 * t415;
t494 = t414 * t567 - t438 * t551;
t569 = t467 * t494;
t387 = atan2(t404, t418);
t382 = sin(t387);
t383 = cos(t387);
t357 = t382 * t404 + t383 * t418;
t354 = 0.1e1 / t357;
t549 = t422 * t466;
t381 = t405 * t462 + t549;
t375 = 0.1e1 / t381;
t355 = 0.1e1 / t357 ^ 2;
t376 = 0.1e1 / t381 ^ 2;
t566 = 0.2e1 * t404;
t406 = t423 * t467 + t442 * t463;
t565 = 0.2e1 * t406;
t399 = t406 ^ 2;
t353 = t355 * t399 + 0.1e1;
t432 = qJD(1) * t447 + qJD(2) * t487;
t433 = qJD(1) * t448 + qJD(2) * t486;
t514 = qJD(1) * t538;
t507 = t458 * t514;
t369 = -t433 * t468 + (t432 * t460 + t507) * t464 + (t468 * t478 + t544) * qJD(3);
t424 = -t432 * t458 + t460 * t514;
t527 = qJD(4) * t463;
t362 = t424 * t463 - t423 * t527 + (qJD(4) * t442 + t369) * t467;
t557 = t362 * t355;
t398 = t404 ^ 2;
t386 = t398 * t415 + 0.1e1;
t384 = 0.1e1 / t386;
t364 = t575 + qJD(4) * t479 + (-qJD(4) * t440 - t373) * t467;
t498 = -t364 * t414 - t389 * t551;
t344 = t498 * t384;
t503 = -t382 * t418 + t383 * t404;
t338 = t344 * t503 - t364 * t382 + t383 * t389;
t356 = t354 * t355;
t562 = t338 * t356;
t563 = (-t399 * t562 + t406 * t557) / t353 ^ 2;
t361 = qJD(4) * t406 + t369 * t463 - t424 * t467;
t536 = t460 * t468;
t368 = qJD(3) * t423 - t432 * t536 - t433 * t464 - t468 * t507;
t347 = qJD(6) * t381 - t361 * t466 + t368 * t462;
t374 = t501 ^ 2;
t360 = t374 * t376 + 0.1e1;
t556 = t376 * t501;
t348 = t361 * t462 + t368 * t466 + t571;
t559 = t348 * t375 * t376;
t561 = (-t347 * t556 - t374 * t559) / t360 ^ 2;
t553 = t414 * t570;
t560 = (-t364 * t551 - t398 * t553) / t386 ^ 2;
t558 = t355 * t406;
t555 = t382 * t406;
t554 = t383 * t406;
t552 = t404 * t414;
t550 = t422 * t463;
t548 = t422 * t467;
t543 = t458 * t463;
t540 = t458 * t467;
t537 = t460 * t464;
t535 = t462 * t501;
t531 = t466 * t375;
t526 = 0.2e1 * t563;
t525 = 0.2e1 * t561;
t524 = -0.2e1 * t560;
t523 = t356 * t565;
t522 = t414 * t560;
t521 = t355 * t555;
t520 = t355 * t554;
t513 = qJD(3) * t542;
t512 = t338 * t523;
t511 = -0.2e1 * t501 * t559;
t510 = t553 * t566;
t504 = -qJD(6) * t550 + t369;
t502 = t403 * t466 + t462 * t567;
t379 = t403 * t462 - t466 * t567;
t429 = -t464 * t486 - t487 * t536;
t430 = t487 * t537 - t483;
t493 = -t430 * t463 - t487 * t540;
t500 = -t429 * t462 - t466 * t493;
t395 = t429 * t466 - t462 * t493;
t497 = -t376 * t535 + t531;
t496 = -t403 * t414 - t417 * t551;
t428 = -t447 * t468 - t448 * t537;
t407 = t428 * t467 + t448 * t543;
t445 = t488 * t459;
t431 = t445 * t467 + t463 * t519;
t495 = -t407 * t414 - t431 * t551;
t409 = t430 * t467 - t487 * t543;
t490 = -t460 * t532 - t533;
t485 = -t382 + (-t383 * t552 + t382) * t384;
t482 = qJD(4) * t548 + qJD(6) * t423 + t368 * t463;
t477 = -t434 * t536 - t435 * t464 + t513 * t538 + t468 * t499 + (t447 * t537 - t545) * qJD(3);
t410 = -t461 * t513 + (qJD(2) * t490 - qJD(3) * t489) * t459;
t396 = -t445 * t527 + ((qJD(3) * t490 + qJD(4) * t541) * t467 + (t463 * t539 - t467 * t489) * qJD(2)) * t459;
t393 = t423 * t466 - t462 * t550;
t392 = t423 * t462 + t463 * t549;
t391 = -qJD(3) * t429 + t432 * t468 + t433 * t537;
t390 = qJD(3) * t430 + t432 * t464 - t433 * t536;
t388 = -qJD(4) * t418 - t411 * t463 + t467 * t506;
t367 = (-t435 * t537 - t434 * t468 + (t447 * t464 - t448 * t536) * qJD(3)) * t467 + t435 * t543 + (-t428 * t463 + t448 * t540) * qJD(4);
t366 = qJD(4) * t409 + t391 * t463 + t433 * t540;
t358 = 0.1e1 / t360;
t351 = 0.1e1 / t353;
t350 = t384 * t569;
t349 = t495 * t384;
t346 = t496 * t384;
t343 = t485 * t406;
t341 = (t382 * t567 + t383 * t438) * t467 + t503 * t350;
t340 = t349 * t503 - t382 * t407 + t383 * t431;
t339 = t346 * t503 - t382 * t403 + t383 * t417;
t337 = t495 * t524 + (t431 * t510 - t367 * t414 + (t364 * t431 + t389 * t407 - t396 * t404) * t415) * t384;
t335 = t496 * t524 + (t417 * t510 - t365 * t414 + (t364 * t417 - t388 * t404 + t389 * t403) * t415) * t384;
t334 = t524 * t569 + (-t494 * t527 + (t438 * t510 - t477 * t414 + (t364 * t438 - t389 * t567 - t404 * t410) * t415) * t467) * t384;
t1 = [t522 * t565 + (-t362 * t414 + t406 * t570) * t384, t337, t334, t335, 0, 0; -0.2e1 * t404 * t354 * t563 + ((-qJD(4) * t403 + t373 * t467 - t575) * t354 + (-t338 * t404 - t343 * t362) * t355) * t351 + (t343 * t355 * t526 + (0.2e1 * t343 * t562 - (t344 * t384 * t552 + t524) * t521 - (t522 * t566 - t344 + (t344 - t498) * t384) * t520 - t485 * t557) * t351) * t406 (t340 * t558 - t354 * t409) * t526 + ((qJD(4) * t493 + t391 * t467 - t433 * t543) * t354 + t340 * t512 + (-t409 * t338 - t340 * t362 - (t337 * t404 - t349 * t364 + t396 + (-t349 * t418 - t407) * t344) * t554 - (-t337 * t418 - t349 * t389 - t367 + (-t349 * t404 - t431) * t344) * t555) * t355) * t351 (t341 * t558 + t354 * t548) * t526 + (-t341 * t557 + (-t368 * t467 + t422 * t527) * t354 + (t341 * t523 + t355 * t548) * t338 - (-t438 * t527 + t334 * t404 - t350 * t364 + t410 * t467 + (-t350 * t418 + t467 * t567) * t344) * t520 - (-t567 * t527 - t334 * t418 - t350 * t389 - t477 * t467 + (-t350 * t404 - t438 * t467) * t344) * t521) * t351 (t339 * t558 + t354 * t405) * t526 + (t339 * t512 - t361 * t354 + (t405 * t338 - t339 * t362 - (t335 * t404 - t346 * t364 + t388 + (-t346 * t418 - t403) * t344) * t554 - (-t335 * t418 - t346 * t389 - t365 + (-t346 * t404 - t417) * t344) * t555) * t355) * t351, 0, 0; (t375 * t502 - t379 * t556) * t525 + ((qJD(6) * t379 - t365 * t466 + t462 * t477) * t375 + t379 * t511 + (t502 * t348 + (qJD(6) * t502 + t365 * t462 + t466 * t477) * t501 - t379 * t347) * t376) * t358 (t375 * t500 - t395 * t556) * t525 + ((qJD(6) * t395 - t366 * t466 + t390 * t462) * t375 + t395 * t511 + (t500 * t348 + (qJD(6) * t500 + t366 * t462 + t390 * t466) * t501 - t395 * t347) * t376) * t358 (-t375 * t392 - t393 * t556) * t525 + (t393 * t511 + t504 * t375 * t462 + t482 * t531 + (t466 * t501 * t504 - t393 * t347 - t392 * t348 - t482 * t535) * t376) * t358, t497 * t406 * t525 + (-t497 * t362 + ((qJD(6) * t375 + t511) * t462 + (-t347 * t462 + (t348 + t571) * t466) * t376) * t406) * t358, 0, -0.2e1 * t561 - 0.2e1 * (t347 * t376 * t358 - (-t358 * t559 - t376 * t561) * t501) * t501;];
JaD_rot  = t1;
