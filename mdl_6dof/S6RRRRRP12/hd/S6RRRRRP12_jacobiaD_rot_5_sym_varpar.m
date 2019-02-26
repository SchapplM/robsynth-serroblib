% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRP12_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP12_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:46:24
% EndTime: 2019-02-26 22:46:30
% DurationCPUTime: 5.78s
% Computational Cost: add. (19869->241), mult. (59163->461), div. (983->12), fcn. (74815->17), ass. (0->187)
t566 = cos(pkin(6));
t567 = sin(qJ(2));
t515 = t566 * t567;
t462 = cos(qJ(2));
t568 = sin(qJ(1));
t530 = t568 * t462;
t569 = cos(qJ(1));
t445 = t515 * t569 + t530;
t481 = t530 * t566 + t567 * t569;
t435 = qJD(1) * t481 + qJD(2) * t445;
t454 = t568 * t567;
t496 = t568 * t515;
t516 = t566 * t569;
t436 = -qJD(1) * t496 - qJD(2) * t454 + (qJD(1) * t569 + qJD(2) * t516) * t462;
t458 = sin(qJ(3));
t461 = cos(qJ(3));
t455 = sin(pkin(6));
t564 = sin(pkin(7));
t528 = t455 * t564;
t504 = t568 * t528;
t494 = qJD(1) * t504;
t565 = cos(pkin(7));
t507 = t569 * t528;
t492 = -t462 * t516 + t454;
t576 = t492 * t565;
t473 = t576 + t507;
t572 = t445 * t458 + t461 * t473;
t376 = t572 * qJD(3) - (-t435 * t565 + t494) * t458 - t436 * t461;
t543 = t445 * t461;
t425 = t458 * t473 - t543;
t457 = sin(qJ(4));
t460 = cos(qJ(4));
t529 = t455 * t565;
t508 = t569 * t529;
t474 = t492 * t564 - t508;
t405 = t425 * t457 + t460 * t474;
t505 = t568 * t529;
t479 = qJD(1) * t505 + t435 * t564;
t367 = qJD(4) * t405 - t376 * t460 + t457 * t479;
t406 = t425 * t460 - t457 * t474;
t585 = qJD(4) * t406 + t376 * t457 + t460 * t479;
t582 = -t481 * t565 + t504;
t524 = t461 * t565;
t542 = qJD(3) * t458;
t579 = (t458 * t576 - t543) * qJD(3) - t435 * t524 - t436 * t458 + t461 * t494 + t507 * t542;
t514 = t565 * t567;
t487 = -t458 * t514 + t461 * t462;
t523 = t462 * t565;
t488 = -t458 * t567 + t461 * t523;
t511 = t566 * t564;
t503 = qJD(3) * t511;
t413 = t461 * t503 + (qJD(2) * t487 + qJD(3) * t488) * t455;
t489 = t458 * t523 + t461 * t567;
t440 = t455 * t489 + t458 * t511;
t522 = t462 * t564;
t444 = -t455 * t522 + t565 * t566;
t419 = t440 * t460 + t444 * t457;
t513 = t564 * t567;
t509 = t455 * t513;
t493 = qJD(2) * t509;
t391 = qJD(4) * t419 + t413 * t457 - t460 * t493;
t418 = t440 * t457 - t444 * t460;
t416 = 0.1e1 / t418 ^ 2;
t578 = t391 * t416;
t415 = 0.1e1 / t418;
t439 = t455 * t488 + t461 * t511;
t548 = t405 * t416;
t498 = t415 * t572 - t439 * t548;
t577 = t457 * t498;
t480 = -t462 * t569 + t496;
t471 = t492 * qJD(1) + t480 * qJD(2);
t575 = qJD(1) * t507 + t471 * t565;
t574 = t582 * t461;
t390 = atan2(t405, t418);
t385 = sin(t390);
t386 = cos(t390);
t360 = t385 * t405 + t386 * t418;
t357 = 0.1e1 / t360;
t427 = t582 * t458 - t480 * t461;
t472 = t481 * t564 + t505;
t408 = t427 * t460 + t457 * t472;
t426 = -t458 * t480 - t574;
t456 = sin(qJ(5));
t459 = cos(qJ(5));
t384 = t408 * t459 + t426 * t456;
t378 = 0.1e1 / t384;
t358 = 0.1e1 / t360 ^ 2;
t379 = 0.1e1 / t384 ^ 2;
t571 = 0.2e1 * t405;
t407 = t427 * t457 - t460 * t472;
t570 = 0.2e1 * t407;
t401 = t407 ^ 2;
t356 = t358 * t401 + 0.1e1;
t434 = qJD(1) * t445 + qJD(2) * t481;
t372 = t574 * qJD(3) - t434 * t461 + t575 * t458 + t480 * t542;
t469 = qJD(1) * t508 - t471 * t564;
t364 = qJD(4) * t408 + t372 * t457 - t460 * t469;
t557 = t358 * t407;
t400 = t405 ^ 2;
t389 = t400 * t416 + 0.1e1;
t387 = 0.1e1 / t389;
t502 = -t391 * t548 + t415 * t585;
t347 = t502 * t387;
t510 = -t385 * t418 + t386 * t405;
t341 = t347 * t510 + t385 * t585 + t386 * t391;
t359 = t357 * t358;
t562 = t341 * t359;
t563 = (t364 * t557 - t401 * t562) / t356 ^ 2;
t550 = t415 * t578;
t560 = (-t400 * t550 + t548 * t585) / t389 ^ 2;
t365 = -qJD(4) * t407 + t372 * t460 + t457 * t469;
t371 = t427 * qJD(3) - t434 * t458 - t575 * t461;
t383 = t408 * t456 - t426 * t459;
t540 = qJD(5) * t383;
t351 = t365 * t459 + t371 * t456 - t540;
t559 = t351 * t378 * t379;
t558 = t358 * t364;
t350 = qJD(5) * t384 + t365 * t456 - t371 * t459;
t377 = t383 ^ 2;
t363 = t377 * t379 + 0.1e1;
t554 = t379 * t383;
t556 = 0.1e1 / t363 ^ 2 * (t350 * t554 - t377 * t559);
t555 = t378 * t456;
t553 = t383 * t459;
t552 = t385 * t407;
t551 = t386 * t407;
t549 = t405 * t415;
t547 = t426 * t457;
t546 = t426 * t460;
t541 = qJD(4) * t460;
t539 = 0.2e1 * t563;
t538 = -0.2e1 * t560;
t537 = t359 * t570;
t536 = -0.2e1 * t556;
t535 = 0.2e1 * t556;
t534 = t415 * t560;
t533 = t383 * t559;
t532 = t358 * t552;
t531 = t358 * t551;
t527 = t457 * t564;
t526 = t458 * t565;
t525 = t460 * t564;
t521 = -0.2e1 * t357 * t563;
t520 = t358 * t539;
t519 = t341 * t537;
t518 = 0.2e1 * t533;
t517 = t550 * t571;
t512 = qJD(5) * t546 + t372;
t382 = t406 * t459 - t456 * t572;
t381 = t406 * t456 + t459 * t572;
t432 = -t461 * t481 + t480 * t526;
t411 = t432 * t460 - t480 * t527;
t431 = -t458 * t481 - t480 * t524;
t398 = t411 * t459 + t431 * t456;
t397 = t411 * t456 - t431 * t459;
t501 = t379 * t553 - t555;
t500 = t406 * t415 - t419 * t548;
t430 = -t445 * t526 - t461 * t492;
t409 = t430 * t457 - t445 * t525;
t443 = t487 * t455;
t433 = t443 * t457 - t460 * t509;
t499 = -t409 * t415 - t433 * t548;
t491 = -t432 * t457 - t480 * t525;
t490 = -t385 + (-t386 * t549 + t385) * t387;
t486 = -t458 * t462 - t461 * t514;
t485 = qJD(4) * t547 + qJD(5) * t427 - t371 * t460;
t412 = -t458 * t503 + (qJD(2) * t486 - qJD(3) * t489) * t455;
t399 = t443 * t541 + ((qJD(3) * t486 + qJD(4) * t513) * t457 + (-t457 * t489 - t460 * t522) * qJD(2)) * t455;
t396 = t427 * t456 - t459 * t546;
t395 = -t427 * t459 - t456 * t546;
t394 = -qJD(3) * t431 + t434 * t526 + t461 * t471;
t393 = qJD(3) * t432 - t434 * t524 + t458 * t471;
t392 = -qJD(4) * t418 + t413 * t460 + t457 * t493;
t370 = (-t436 * t526 - t435 * t461 + (-t445 * t524 + t458 * t492) * qJD(3)) * t457 + t430 * t541 - t436 * t525 + t445 * qJD(4) * t527;
t369 = qJD(4) * t491 + t394 * t460 - t434 * t527;
t361 = 0.1e1 / t363;
t354 = 0.1e1 / t356;
t353 = t387 * t577;
t352 = t499 * t387;
t349 = t500 * t387;
t346 = t490 * t407;
t344 = (t385 * t572 + t386 * t439) * t457 + t510 * t353;
t342 = t349 * t510 + t385 * t406 + t386 * t419;
t340 = t499 * t538 + (t433 * t517 - t370 * t415 + (t391 * t409 - t399 * t405 - t433 * t585) * t416) * t387;
t338 = t500 * t538 + (t419 * t517 - t367 * t415 + (-t391 * t406 - t392 * t405 - t419 * t585) * t416) * t387;
t337 = t538 * t577 + (t498 * t541 + (t439 * t517 - t579 * t415 + (-t391 * t572 - t405 * t412 - t439 * t585) * t416) * t457) * t387;
t1 = [t534 * t570 + (-t364 * t415 + t407 * t578) * t387, t340, t337, t338, 0, 0; t405 * t521 + (t585 * t357 + (-t341 * t405 - t346 * t364) * t358) * t354 + (t346 * t520 + (0.2e1 * t346 * t562 - (t347 * t387 * t549 + t538) * t532 - (t534 * t571 - t347 + (t347 - t502) * t387) * t531 - t490 * t558) * t354) * t407, -t491 * t521 + ((qJD(4) * t411 + t394 * t457 + t434 * t525) * t357 + t491 * t358 * t341 - ((t340 * t405 + t352 * t585 + t399 + (-t352 * t418 - t409) * t347) * t386 + (-t340 * t418 - t352 * t391 - t370 + (-t352 * t405 - t433) * t347) * t385) * t557) * t354 + (t407 * t520 + (-t558 + t519) * t354) * (t352 * t510 - t385 * t409 + t386 * t433) (t344 * t557 + t357 * t547) * t539 + (-t344 * t558 + (-t371 * t457 - t426 * t541) * t357 + (t344 * t537 + t358 * t547) * t341 - (t439 * t541 + t337 * t405 + t353 * t585 + t412 * t457 + (-t353 * t418 + t457 * t572) * t347) * t531 - (t572 * t541 - t337 * t418 - t353 * t391 - t579 * t457 + (-t353 * t405 - t439 * t457) * t347) * t532) * t354 (t342 * t557 - t357 * t408) * t539 + (t342 * t519 + t365 * t357 + (-t408 * t341 - t342 * t364 - (t338 * t405 + t349 * t585 + t392 + (-t349 * t418 + t406) * t347) * t551 - (-t338 * t418 - t349 * t391 - t367 + (-t349 * t405 - t419) * t347) * t552) * t358) * t354, 0, 0; (-t378 * t381 + t382 * t554) * t535 + ((qJD(5) * t382 - t367 * t456 - t459 * t579) * t378 + t382 * t518 + (-t381 * t351 - (-qJD(5) * t381 - t367 * t459 + t456 * t579) * t383 - t382 * t350) * t379) * t361 (-t378 * t397 + t398 * t554) * t535 + ((qJD(5) * t398 + t369 * t456 - t393 * t459) * t378 + t398 * t518 + (-t397 * t351 - (-qJD(5) * t397 + t369 * t459 + t393 * t456) * t383 - t398 * t350) * t379) * t361 (-t378 * t395 + t396 * t554) * t535 + (t396 * t518 - t512 * t378 * t459 + t485 * t555 + (-t383 * t456 * t512 - t396 * t350 - t395 * t351 - t485 * t553) * t379) * t361, t501 * t407 * t536 + (t501 * t364 + ((-qJD(5) * t378 - 0.2e1 * t533) * t459 + (t350 * t459 + (t351 - t540) * t456) * t379) * t407) * t361, t536 + 0.2e1 * (t350 * t379 * t361 + (-t361 * t559 - t379 * t556) * t383) * t383, 0;];
JaD_rot  = t1;
