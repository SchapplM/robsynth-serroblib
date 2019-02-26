% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP11
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
% Datum: 2019-02-26 22:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRP11_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP11_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:45:40
% EndTime: 2019-02-26 22:45:45
% DurationCPUTime: 5.69s
% Computational Cost: add. (19869->241), mult. (59163->461), div. (983->12), fcn. (74815->17), ass. (0->187)
t567 = cos(pkin(6));
t568 = sin(qJ(2));
t516 = t567 * t568;
t463 = cos(qJ(2));
t569 = sin(qJ(1));
t531 = t569 * t463;
t570 = cos(qJ(1));
t446 = t570 * t516 + t531;
t482 = t567 * t531 + t570 * t568;
t436 = qJD(1) * t482 + qJD(2) * t446;
t455 = t569 * t568;
t497 = t569 * t516;
t517 = t567 * t570;
t437 = -qJD(1) * t497 - qJD(2) * t455 + (t570 * qJD(1) + qJD(2) * t517) * t463;
t459 = sin(qJ(3));
t462 = cos(qJ(3));
t456 = sin(pkin(6));
t565 = sin(pkin(7));
t529 = t456 * t565;
t505 = t569 * t529;
t495 = qJD(1) * t505;
t566 = cos(pkin(7));
t508 = t570 * t529;
t493 = -t463 * t517 + t455;
t577 = t493 * t566;
t474 = t577 + t508;
t573 = t446 * t459 + t462 * t474;
t377 = qJD(3) * t573 - (-t436 * t566 + t495) * t459 - t437 * t462;
t546 = t446 * t462;
t426 = t459 * t474 - t546;
t458 = sin(qJ(4));
t461 = cos(qJ(4));
t530 = t456 * t566;
t509 = t570 * t530;
t475 = t493 * t565 - t509;
t406 = t426 * t458 + t461 * t475;
t506 = t569 * t530;
t480 = qJD(1) * t506 + t436 * t565;
t368 = qJD(4) * t406 - t377 * t461 + t458 * t480;
t407 = t426 * t461 - t458 * t475;
t586 = qJD(4) * t407 + t377 * t458 + t461 * t480;
t583 = -t482 * t566 + t505;
t525 = t462 * t566;
t543 = qJD(3) * t459;
t580 = (t459 * t577 - t546) * qJD(3) - t436 * t525 - t437 * t459 + t462 * t495 + t508 * t543;
t515 = t566 * t568;
t488 = -t459 * t515 + t462 * t463;
t524 = t463 * t566;
t489 = -t568 * t459 + t462 * t524;
t512 = t567 * t565;
t504 = qJD(3) * t512;
t414 = t462 * t504 + (qJD(2) * t488 + qJD(3) * t489) * t456;
t490 = t459 * t524 + t568 * t462;
t441 = t456 * t490 + t459 * t512;
t523 = t463 * t565;
t445 = -t456 * t523 + t567 * t566;
t420 = t441 * t461 + t445 * t458;
t514 = t565 * t568;
t510 = t456 * t514;
t494 = qJD(2) * t510;
t392 = qJD(4) * t420 + t414 * t458 - t461 * t494;
t419 = t441 * t458 - t445 * t461;
t417 = 0.1e1 / t419 ^ 2;
t579 = t392 * t417;
t416 = 0.1e1 / t419;
t440 = t456 * t489 + t462 * t512;
t551 = t406 * t417;
t499 = t416 * t573 - t440 * t551;
t578 = t458 * t499;
t481 = -t570 * t463 + t497;
t472 = qJD(1) * t493 + t481 * qJD(2);
t576 = qJD(1) * t508 + t472 * t566;
t575 = t583 * t462;
t391 = atan2(t406, t419);
t386 = sin(t391);
t387 = cos(t391);
t361 = t386 * t406 + t387 * t419;
t358 = 0.1e1 / t361;
t428 = t583 * t459 - t481 * t462;
t473 = t482 * t565 + t506;
t409 = t428 * t461 + t458 * t473;
t427 = -t459 * t481 - t575;
t457 = sin(qJ(5));
t460 = cos(qJ(5));
t385 = t409 * t460 + t427 * t457;
t379 = 0.1e1 / t385;
t359 = 0.1e1 / t361 ^ 2;
t380 = 0.1e1 / t385 ^ 2;
t572 = 0.2e1 * t406;
t408 = t428 * t458 - t461 * t473;
t571 = 0.2e1 * t408;
t402 = t408 ^ 2;
t357 = t359 * t402 + 0.1e1;
t435 = qJD(1) * t446 + qJD(2) * t482;
t373 = qJD(3) * t575 - t435 * t462 + t459 * t576 + t481 * t543;
t470 = qJD(1) * t509 - t472 * t565;
t365 = qJD(4) * t409 + t373 * t458 - t461 * t470;
t557 = t365 * t359;
t401 = t406 ^ 2;
t390 = t401 * t417 + 0.1e1;
t388 = 0.1e1 / t390;
t503 = -t392 * t551 + t416 * t586;
t348 = t503 * t388;
t511 = -t386 * t419 + t387 * t406;
t342 = t348 * t511 + t386 * t586 + t387 * t392;
t360 = t358 * t359;
t563 = t342 * t360;
t564 = (-t402 * t563 + t408 * t557) / t357 ^ 2;
t366 = -qJD(4) * t408 + t373 * t461 + t458 * t470;
t372 = qJD(3) * t428 - t435 * t459 - t462 * t576;
t351 = qJD(5) * t385 + t366 * t457 - t372 * t460;
t384 = t409 * t457 - t427 * t460;
t378 = t384 ^ 2;
t364 = t378 * t380 + 0.1e1;
t556 = t380 * t384;
t541 = qJD(5) * t384;
t352 = t366 * t460 + t372 * t457 - t541;
t559 = t352 * t379 * t380;
t561 = (t351 * t556 - t378 * t559) / t364 ^ 2;
t553 = t416 * t579;
t560 = (-t401 * t553 + t551 * t586) / t390 ^ 2;
t558 = t359 * t408;
t555 = t386 * t408;
t554 = t387 * t408;
t552 = t406 * t416;
t550 = t427 * t458;
t549 = t427 * t461;
t545 = t457 * t379;
t544 = t460 * t384;
t542 = qJD(4) * t461;
t540 = 0.2e1 * t564;
t539 = -0.2e1 * t561;
t538 = 0.2e1 * t561;
t537 = -0.2e1 * t560;
t536 = t360 * t571;
t535 = t416 * t560;
t534 = t359 * t555;
t533 = t359 * t554;
t532 = t384 * t559;
t528 = t458 * t565;
t527 = t459 * t566;
t526 = t461 * t565;
t522 = -0.2e1 * t358 * t564;
t521 = t359 * t540;
t520 = 0.2e1 * t532;
t519 = t553 * t572;
t518 = t342 * t536;
t513 = qJD(5) * t549 + t373;
t383 = t407 * t460 - t457 * t573;
t382 = t407 * t457 + t460 * t573;
t433 = -t482 * t462 + t481 * t527;
t412 = t433 * t461 - t481 * t528;
t432 = -t482 * t459 - t481 * t525;
t399 = t412 * t460 + t432 * t457;
t398 = t412 * t457 - t432 * t460;
t502 = t380 * t544 - t545;
t501 = t407 * t416 - t420 * t551;
t431 = -t446 * t527 - t493 * t462;
t410 = t431 * t458 - t446 * t526;
t444 = t488 * t456;
t434 = t444 * t458 - t461 * t510;
t500 = -t410 * t416 - t434 * t551;
t492 = -t433 * t458 - t481 * t526;
t491 = -t386 + (-t387 * t552 + t386) * t388;
t487 = -t459 * t463 - t462 * t515;
t486 = qJD(4) * t550 + qJD(5) * t428 - t372 * t461;
t413 = -t459 * t504 + (qJD(2) * t487 - qJD(3) * t490) * t456;
t400 = t444 * t542 + ((qJD(3) * t487 + qJD(4) * t514) * t458 + (-t490 * t458 - t461 * t523) * qJD(2)) * t456;
t397 = t428 * t457 - t460 * t549;
t396 = -t428 * t460 - t457 * t549;
t395 = -t432 * qJD(3) + t435 * t527 + t472 * t462;
t394 = t433 * qJD(3) - t435 * t525 + t472 * t459;
t393 = -qJD(4) * t419 + t414 * t461 + t458 * t494;
t371 = (-t437 * t527 - t436 * t462 + (-t446 * t525 + t493 * t459) * qJD(3)) * t458 + t431 * t542 - t437 * t526 + t446 * qJD(4) * t528;
t370 = t492 * qJD(4) + t395 * t461 - t435 * t528;
t362 = 0.1e1 / t364;
t355 = 0.1e1 / t357;
t354 = t388 * t578;
t353 = t500 * t388;
t350 = t501 * t388;
t347 = t491 * t408;
t345 = (t386 * t573 + t387 * t440) * t458 + t511 * t354;
t343 = t350 * t511 + t386 * t407 + t387 * t420;
t341 = t500 * t537 + (t434 * t519 - t371 * t416 + (t392 * t410 - t400 * t406 - t434 * t586) * t417) * t388;
t339 = t501 * t537 + (t420 * t519 - t368 * t416 + (-t392 * t407 - t393 * t406 - t420 * t586) * t417) * t388;
t338 = t537 * t578 + (t499 * t542 + (t440 * t519 - t580 * t416 + (-t392 * t573 - t406 * t413 - t440 * t586) * t417) * t458) * t388;
t1 = [t535 * t571 + (-t365 * t416 + t408 * t579) * t388, t341, t338, t339, 0, 0; t406 * t522 + (t586 * t358 + (-t342 * t406 - t347 * t365) * t359) * t355 + (t347 * t521 + (0.2e1 * t347 * t563 - (t348 * t388 * t552 + t537) * t534 - (t535 * t572 - t348 + (t348 - t503) * t388) * t533 - t491 * t557) * t355) * t408, -t492 * t522 + ((t412 * qJD(4) + t395 * t458 + t435 * t526) * t358 + t492 * t359 * t342 - ((t341 * t406 + t353 * t586 + t400 + (-t353 * t419 - t410) * t348) * t387 + (-t341 * t419 - t353 * t392 - t371 + (-t353 * t406 - t434) * t348) * t386) * t558) * t355 + (t408 * t521 + (-t557 + t518) * t355) * (t353 * t511 - t386 * t410 + t387 * t434) (t345 * t558 + t358 * t550) * t540 + (-t345 * t557 + (-t372 * t458 - t427 * t542) * t358 + (t345 * t536 + t359 * t550) * t342 - (t440 * t542 + t338 * t406 + t354 * t586 + t413 * t458 + (-t354 * t419 + t458 * t573) * t348) * t533 - (t573 * t542 - t338 * t419 - t354 * t392 - t580 * t458 + (-t354 * t406 - t440 * t458) * t348) * t534) * t355 (t343 * t558 - t358 * t409) * t540 + (t343 * t518 + t366 * t358 + (-t409 * t342 - t343 * t365 - (t339 * t406 + t350 * t586 + t393 + (-t350 * t419 + t407) * t348) * t554 - (-t339 * t419 - t350 * t392 - t368 + (-t350 * t406 - t420) * t348) * t555) * t359) * t355, 0, 0; (-t379 * t382 + t383 * t556) * t538 + ((qJD(5) * t383 - t368 * t457 - t460 * t580) * t379 + t383 * t520 + (-t382 * t352 - (-qJD(5) * t382 - t368 * t460 + t457 * t580) * t384 - t383 * t351) * t380) * t362 (-t379 * t398 + t399 * t556) * t538 + ((qJD(5) * t399 + t370 * t457 - t394 * t460) * t379 + t399 * t520 + (-t398 * t352 - (-qJD(5) * t398 + t370 * t460 + t394 * t457) * t384 - t399 * t351) * t380) * t362 (-t379 * t396 + t397 * t556) * t538 + (t397 * t520 - t513 * t379 * t460 + t486 * t545 + (-t384 * t457 * t513 - t397 * t351 - t396 * t352 - t486 * t544) * t380) * t362, t502 * t408 * t539 + (t502 * t365 + ((-qJD(5) * t379 - 0.2e1 * t532) * t460 + (t351 * t460 + (t352 - t541) * t457) * t380) * t408) * t362, t539 + 0.2e1 * (t351 * t380 * t362 + (-t362 * t559 - t380 * t561) * t384) * t384, 0;];
JaD_rot  = t1;
