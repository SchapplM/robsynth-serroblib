% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR14
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
% Datum: 2019-02-26 22:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR14_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR14_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_jacobiaD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:13
% EndTime: 2019-02-26 22:38:19
% DurationCPUTime: 5.89s
% Computational Cost: add. (20400->242), mult. (59163->460), div. (983->12), fcn. (74815->17), ass. (0->187)
t569 = cos(pkin(6));
t570 = sin(qJ(2));
t519 = t569 * t570;
t466 = cos(qJ(2));
t571 = sin(qJ(1));
t533 = t571 * t466;
t572 = cos(qJ(1));
t448 = t572 * t519 + t533;
t485 = t569 * t533 + t572 * t570;
t438 = t485 * qJD(1) + t448 * qJD(2);
t457 = t571 * t570;
t500 = t571 * t519;
t520 = t569 * t572;
t439 = -qJD(1) * t500 - qJD(2) * t457 + (t572 * qJD(1) + qJD(2) * t520) * t466;
t463 = sin(qJ(3));
t465 = cos(qJ(3));
t461 = sin(pkin(6));
t567 = sin(pkin(7));
t531 = t461 * t567;
t508 = t571 * t531;
t498 = qJD(1) * t508;
t568 = cos(pkin(7));
t511 = t572 * t531;
t496 = -t466 * t520 + t457;
t579 = t496 * t568;
t477 = t579 + t511;
t575 = t448 * t463 + t477 * t465;
t387 = t575 * qJD(3) - (-t438 * t568 + t498) * t463 - t439 * t465;
t547 = t448 * t465;
t428 = t477 * t463 - t547;
t462 = sin(qJ(4));
t464 = cos(qJ(4));
t532 = t461 * t568;
t512 = t572 * t532;
t478 = t496 * t567 - t512;
t408 = t428 * t462 + t478 * t464;
t509 = t571 * t532;
t483 = qJD(1) * t509 + t438 * t567;
t370 = t408 * qJD(4) - t387 * t464 + t483 * t462;
t409 = t428 * t464 - t478 * t462;
t588 = t409 * qJD(4) + t387 * t462 + t483 * t464;
t585 = -t485 * t568 + t508;
t528 = t465 * t568;
t545 = qJD(3) * t463;
t582 = (t463 * t579 - t547) * qJD(3) - t438 * t528 - t439 * t463 + t465 * t498 + t511 * t545;
t518 = t568 * t570;
t491 = -t463 * t518 + t465 * t466;
t527 = t466 * t568;
t492 = -t570 * t463 + t465 * t527;
t515 = t569 * t567;
t507 = qJD(3) * t515;
t417 = t465 * t507 + (t491 * qJD(2) + t492 * qJD(3)) * t461;
t493 = t463 * t527 + t570 * t465;
t443 = t493 * t461 + t463 * t515;
t447 = -t466 * t531 + t569 * t568;
t422 = t443 * t464 + t447 * t462;
t517 = t567 * t570;
t513 = t461 * t517;
t497 = qJD(2) * t513;
t394 = t422 * qJD(4) + t417 * t462 - t464 * t497;
t421 = t443 * t462 - t447 * t464;
t419 = 0.1e1 / t421 ^ 2;
t581 = t394 * t419;
t418 = 0.1e1 / t421;
t442 = t492 * t461 + t465 * t515;
t552 = t408 * t419;
t502 = t418 * t575 - t442 * t552;
t580 = t462 * t502;
t484 = -t572 * t466 + t500;
t475 = t496 * qJD(1) + t484 * qJD(2);
t578 = qJD(1) * t511 + t475 * t568;
t577 = t585 * t465;
t393 = atan2(t408, t421);
t388 = sin(t393);
t389 = cos(t393);
t363 = t388 * t408 + t389 * t421;
t360 = 0.1e1 / t363;
t430 = t585 * t463 - t484 * t465;
t476 = t485 * t567 + t509;
t411 = t430 * t464 + t476 * t462;
t429 = -t463 * t484 - t577;
t460 = pkin(13) + qJ(6);
t458 = sin(t460);
t459 = cos(t460);
t381 = t411 * t459 + t429 * t458;
t375 = 0.1e1 / t381;
t361 = 0.1e1 / t363 ^ 2;
t376 = 0.1e1 / t381 ^ 2;
t574 = 0.2e1 * t408;
t410 = t430 * t462 - t476 * t464;
t573 = 0.2e1 * t410;
t404 = t410 ^ 2;
t359 = t404 * t361 + 0.1e1;
t437 = t448 * qJD(1) + t485 * qJD(2);
t383 = t577 * qJD(3) - t437 * t465 + t578 * t463 + t484 * t545;
t473 = qJD(1) * t512 - t475 * t567;
t367 = t411 * qJD(4) + t383 * t462 - t473 * t464;
t559 = t367 * t361;
t403 = t408 ^ 2;
t392 = t403 * t419 + 0.1e1;
t390 = 0.1e1 / t392;
t506 = -t394 * t552 + t418 * t588;
t350 = t506 * t390;
t514 = -t388 * t421 + t389 * t408;
t344 = t514 * t350 + t388 * t588 + t389 * t394;
t362 = t360 * t361;
t565 = t344 * t362;
t566 = (-t404 * t565 + t410 * t559) / t359 ^ 2;
t368 = -t410 * qJD(4) + t383 * t464 + t473 * t462;
t382 = t430 * qJD(3) - t437 * t463 - t578 * t465;
t353 = t381 * qJD(6) + t368 * t458 - t382 * t459;
t380 = t411 * t458 - t429 * t459;
t374 = t380 ^ 2;
t366 = t374 * t376 + 0.1e1;
t558 = t376 * t380;
t543 = qJD(6) * t380;
t354 = t368 * t459 + t382 * t458 - t543;
t561 = t354 * t375 * t376;
t563 = (t353 * t558 - t374 * t561) / t366 ^ 2;
t554 = t418 * t581;
t562 = (-t403 * t554 + t552 * t588) / t392 ^ 2;
t560 = t361 * t410;
t557 = t380 * t459;
t556 = t388 * t410;
t555 = t389 * t410;
t553 = t408 * t418;
t551 = t429 * t462;
t550 = t429 * t464;
t546 = t458 * t375;
t544 = qJD(4) * t464;
t542 = 0.2e1 * t566;
t541 = -0.2e1 * t563;
t540 = 0.2e1 * t563;
t539 = -0.2e1 * t562;
t538 = t362 * t573;
t537 = t418 * t562;
t536 = t361 * t556;
t535 = t361 * t555;
t534 = t380 * t561;
t530 = t462 * t567;
t529 = t463 * t568;
t526 = t567 * t464;
t525 = -0.2e1 * t360 * t566;
t524 = t361 * t542;
t523 = t344 * t538;
t522 = 0.2e1 * t534;
t521 = t554 * t574;
t516 = qJD(6) * t550 + t383;
t379 = t409 * t459 - t458 * t575;
t378 = t409 * t458 + t459 * t575;
t435 = -t485 * t465 + t484 * t529;
t414 = t435 * t464 - t484 * t530;
t434 = -t485 * t463 - t484 * t528;
t401 = t414 * t459 + t434 * t458;
t400 = t414 * t458 - t434 * t459;
t505 = t376 * t557 - t546;
t504 = t409 * t418 - t422 * t552;
t433 = -t448 * t529 - t496 * t465;
t412 = t433 * t462 - t448 * t526;
t446 = t491 * t461;
t436 = t446 * t462 - t464 * t513;
t503 = -t412 * t418 - t436 * t552;
t495 = -t435 * t462 - t484 * t526;
t494 = -t388 + (-t389 * t553 + t388) * t390;
t490 = -t463 * t466 - t465 * t518;
t489 = qJD(4) * t551 + qJD(6) * t430 - t382 * t464;
t416 = -t463 * t507 + (t490 * qJD(2) - t493 * qJD(3)) * t461;
t402 = t446 * t544 + ((t490 * qJD(3) + qJD(4) * t517) * t462 + (-t493 * t462 - t466 * t526) * qJD(2)) * t461;
t399 = t430 * t458 - t459 * t550;
t398 = -t430 * t459 - t458 * t550;
t397 = -t434 * qJD(3) + t437 * t529 + t475 * t465;
t396 = t435 * qJD(3) - t437 * t528 + t475 * t463;
t395 = -t421 * qJD(4) + t417 * t464 + t462 * t497;
t373 = (-t439 * t529 - t438 * t465 + (-t448 * t528 + t496 * t463) * qJD(3)) * t462 + t433 * t544 - t439 * t526 + t448 * qJD(4) * t530;
t372 = t495 * qJD(4) + t397 * t464 - t437 * t530;
t364 = 0.1e1 / t366;
t357 = 0.1e1 / t359;
t356 = t390 * t580;
t355 = t503 * t390;
t352 = t504 * t390;
t349 = t494 * t410;
t347 = (t388 * t575 + t389 * t442) * t462 + t514 * t356;
t345 = t514 * t352 + t388 * t409 + t389 * t422;
t343 = t503 * t539 + (t436 * t521 - t373 * t418 + (t394 * t412 - t402 * t408 - t436 * t588) * t419) * t390;
t341 = t504 * t539 + (t422 * t521 - t370 * t418 + (-t394 * t409 - t395 * t408 - t422 * t588) * t419) * t390;
t340 = t539 * t580 + (t502 * t544 + (t442 * t521 - t582 * t418 + (-t394 * t575 - t408 * t416 - t442 * t588) * t419) * t462) * t390;
t1 = [t537 * t573 + (-t367 * t418 + t410 * t581) * t390, t343, t340, t341, 0, 0; t408 * t525 + (t588 * t360 + (-t408 * t344 - t349 * t367) * t361) * t357 + (t349 * t524 + (0.2e1 * t349 * t565 - (t350 * t390 * t553 + t539) * t536 - (t537 * t574 - t350 + (t350 - t506) * t390) * t535 - t494 * t559) * t357) * t410, -t495 * t525 + ((t414 * qJD(4) + t397 * t462 + t437 * t526) * t360 + t495 * t361 * t344 - ((t343 * t408 + t355 * t588 + t402 + (-t355 * t421 - t412) * t350) * t389 + (-t343 * t421 - t355 * t394 - t373 + (-t355 * t408 - t436) * t350) * t388) * t560) * t357 + (t410 * t524 + (-t559 + t523) * t357) * (t514 * t355 - t388 * t412 + t389 * t436) (t347 * t560 + t360 * t551) * t542 + (-t347 * t559 + (-t382 * t462 - t429 * t544) * t360 + (t347 * t538 + t361 * t551) * t344 - (t442 * t544 + t340 * t408 + t356 * t588 + t416 * t462 + (-t356 * t421 + t462 * t575) * t350) * t535 - (t575 * t544 - t340 * t421 - t356 * t394 - t582 * t462 + (-t356 * t408 - t442 * t462) * t350) * t536) * t357 (t345 * t560 - t360 * t411) * t542 + (t345 * t523 + t368 * t360 + (-t411 * t344 - t345 * t367 - (t341 * t408 + t352 * t588 + t395 + (-t352 * t421 + t409) * t350) * t555 - (-t341 * t421 - t352 * t394 - t370 + (-t352 * t408 - t422) * t350) * t556) * t361) * t357, 0, 0; (-t375 * t378 + t379 * t558) * t540 + ((t379 * qJD(6) - t370 * t458 - t459 * t582) * t375 + t379 * t522 + (-t378 * t354 - (-t378 * qJD(6) - t370 * t459 + t458 * t582) * t380 - t379 * t353) * t376) * t364 (-t375 * t400 + t401 * t558) * t540 + ((t401 * qJD(6) + t372 * t458 - t396 * t459) * t375 + t401 * t522 + (-t400 * t354 - (-t400 * qJD(6) + t372 * t459 + t396 * t458) * t380 - t401 * t353) * t376) * t364 (-t375 * t398 + t399 * t558) * t540 + (t399 * t522 - t516 * t375 * t459 + t489 * t546 + (-t516 * t380 * t458 - t399 * t353 - t398 * t354 - t489 * t557) * t376) * t364, t505 * t410 * t541 + (t505 * t367 + ((-qJD(6) * t375 - 0.2e1 * t534) * t459 + (t353 * t459 + (t354 - t543) * t458) * t376) * t410) * t364, 0, t541 + 0.2e1 * (t353 * t376 * t364 + (-t364 * t561 - t376 * t563) * t380) * t380;];
JaD_rot  = t1;
