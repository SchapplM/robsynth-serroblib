% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR13_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR13_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_jacobiaD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:06
% EndTime: 2019-02-26 22:23:12
% DurationCPUTime: 5.71s
% Computational Cost: add. (25275->242), mult. (59163->461), div. (983->12), fcn. (74815->17), ass. (0->189)
t580 = cos(pkin(6));
t581 = sin(qJ(2));
t528 = t580 * t581;
t475 = cos(qJ(2));
t582 = sin(qJ(1));
t543 = t582 * t475;
t583 = cos(qJ(1));
t457 = t583 * t528 + t543;
t494 = t580 * t543 + t583 * t581;
t447 = t494 * qJD(1) + t457 * qJD(2);
t466 = t582 * t581;
t509 = t582 * t528;
t529 = t580 * t583;
t448 = -qJD(1) * t509 - qJD(2) * t466 + (t583 * qJD(1) + qJD(2) * t529) * t475;
t472 = sin(qJ(3));
t474 = cos(qJ(3));
t470 = sin(pkin(6));
t578 = sin(pkin(7));
t539 = t470 * t578;
t517 = t582 * t539;
t507 = qJD(1) * t517;
t579 = cos(pkin(7));
t520 = t583 * t539;
t505 = -t475 * t529 + t466;
t590 = t505 * t579;
t486 = t590 + t520;
t586 = t457 * t472 + t486 * t474;
t402 = t586 * qJD(3) - (-t447 * t579 + t507) * t472 - t448 * t474;
t556 = t457 * t474;
t437 = t486 * t472 - t556;
t469 = pkin(13) + qJ(5);
t467 = sin(t469);
t468 = cos(t469);
t540 = t470 * t579;
t521 = t583 * t540;
t487 = t505 * t578 - t521;
t417 = t437 * t467 + t487 * t468;
t518 = t582 * t540;
t492 = qJD(1) * t518 + t447 * t578;
t379 = t417 * qJD(5) - t402 * t468 + t492 * t467;
t418 = t437 * t468 - t487 * t467;
t599 = t418 * qJD(5) + t402 * t467 + t492 * t468;
t596 = -t494 * t579 + t517;
t537 = t474 * t579;
t555 = qJD(3) * t472;
t593 = (t472 * t590 - t556) * qJD(3) - t447 * t537 - t448 * t472 + t474 * t507 + t520 * t555;
t536 = t475 * t579;
t502 = t472 * t536 + t581 * t474;
t524 = t580 * t578;
t452 = t502 * t470 + t472 * t524;
t535 = t475 * t578;
t456 = -t470 * t535 + t580 * t579;
t428 = t452 * t468 + t456 * t467;
t527 = t579 * t581;
t500 = -t472 * t527 + t474 * t475;
t501 = -t581 * t472 + t474 * t536;
t516 = qJD(3) * t524;
t430 = t474 * t516 + (t500 * qJD(2) + t501 * qJD(3)) * t470;
t526 = t578 * t581;
t522 = t470 * t526;
t506 = qJD(2) * t522;
t403 = t428 * qJD(5) + t430 * t467 - t468 * t506;
t427 = t452 * t467 - t456 * t468;
t425 = 0.1e1 / t427 ^ 2;
t592 = t403 * t425;
t424 = 0.1e1 / t427;
t451 = t501 * t470 + t474 * t524;
t562 = t417 * t425;
t511 = t424 * t586 - t451 * t562;
t591 = t467 * t511;
t493 = -t583 * t475 + t509;
t484 = t505 * qJD(1) + t493 * qJD(2);
t589 = qJD(1) * t520 + t484 * t579;
t588 = t596 * t474;
t392 = atan2(t417, t427);
t383 = sin(t392);
t384 = cos(t392);
t372 = t383 * t417 + t384 * t427;
t369 = 0.1e1 / t372;
t439 = t596 * t472 - t493 * t474;
t485 = t494 * t578 + t518;
t420 = t439 * t468 + t485 * t467;
t473 = cos(qJ(6));
t438 = -t472 * t493 - t588;
t471 = sin(qJ(6));
t560 = t438 * t471;
t396 = t420 * t473 + t560;
t389 = 0.1e1 / t396;
t370 = 0.1e1 / t372 ^ 2;
t390 = 0.1e1 / t396 ^ 2;
t585 = 0.2e1 * t417;
t419 = t439 * t467 - t485 * t468;
t584 = 0.2e1 * t419;
t413 = t419 ^ 2;
t368 = t413 * t370 + 0.1e1;
t446 = t457 * qJD(1) + t494 * qJD(2);
t398 = t588 * qJD(3) - t446 * t474 + t589 * t472 + t493 * t555;
t482 = qJD(1) * t521 - t484 * t578;
t376 = t420 * qJD(5) + t398 * t467 - t482 * t468;
t571 = t370 * t419;
t412 = t417 ^ 2;
t387 = t412 * t425 + 0.1e1;
t385 = 0.1e1 / t387;
t515 = -t403 * t562 + t424 * t599;
t359 = t515 * t385;
t523 = -t383 * t427 + t384 * t417;
t353 = t523 * t359 + t383 * t599 + t384 * t403;
t371 = t369 * t370;
t576 = t353 * t371;
t577 = (t376 * t571 - t413 * t576) / t368 ^ 2;
t564 = t424 * t592;
t574 = (-t412 * t564 + t562 * t599) / t387 ^ 2;
t377 = -t419 * qJD(5) + t398 * t468 + t482 * t467;
t397 = t439 * qJD(3) - t446 * t472 - t589 * t474;
t559 = t438 * t473;
t395 = t420 * t471 - t559;
t553 = qJD(6) * t395;
t364 = t377 * t473 + t397 * t471 - t553;
t573 = t364 * t389 * t390;
t572 = t370 * t376;
t363 = t396 * qJD(6) + t377 * t471 - t397 * t473;
t388 = t395 ^ 2;
t375 = t388 * t390 + 0.1e1;
t566 = t390 * t395;
t570 = 0.1e1 / t375 ^ 2 * (t363 * t566 - t388 * t573);
t569 = t383 * t419;
t568 = t384 * t419;
t567 = t389 * t471;
t565 = t395 * t473;
t563 = t417 * t424;
t561 = t438 * t467;
t554 = qJD(5) * t468;
t552 = 0.2e1 * t577;
t551 = -0.2e1 * t574;
t550 = t371 * t584;
t549 = -0.2e1 * t570;
t548 = 0.2e1 * t570;
t547 = t424 * t574;
t546 = t395 * t573;
t545 = t370 * t569;
t544 = t370 * t568;
t542 = t467 * t578;
t541 = t468 * t578;
t538 = t472 * t579;
t534 = -0.2e1 * t369 * t577;
t533 = t370 * t552;
t532 = t353 * t550;
t531 = 0.2e1 * t546;
t530 = t564 * t585;
t525 = qJD(6) * t438 * t468 + t398;
t394 = t418 * t473 - t471 * t586;
t393 = t418 * t471 + t473 * t586;
t445 = -t494 * t474 + t493 * t538;
t423 = t445 * t468 - t493 * t542;
t444 = -t494 * t472 - t493 * t537;
t408 = t423 * t473 + t444 * t471;
t407 = t423 * t471 - t444 * t473;
t514 = t390 * t565 - t567;
t513 = t418 * t424 - t428 * t562;
t443 = -t457 * t538 - t505 * t474;
t421 = t443 * t467 - t457 * t541;
t455 = t500 * t470;
t442 = t455 * t467 - t468 * t522;
t512 = -t421 * t424 - t442 * t562;
t504 = -t445 * t467 - t493 * t541;
t503 = -t383 + (-t384 * t563 + t383) * t385;
t499 = -t472 * t475 - t474 * t527;
t498 = qJD(5) * t561 + qJD(6) * t439 - t397 * t468;
t429 = -t472 * t516 + (t499 * qJD(2) - t502 * qJD(3)) * t470;
t411 = t455 * t554 + ((t499 * qJD(3) + qJD(5) * t526) * t467 + (-t502 * t467 - t468 * t535) * qJD(2)) * t470;
t410 = t439 * t471 - t468 * t559;
t409 = -t439 * t473 - t468 * t560;
t406 = -t444 * qJD(3) + t446 * t538 + t484 * t474;
t405 = t445 * qJD(3) - t446 * t537 + t484 * t472;
t404 = -t427 * qJD(5) + t430 * t468 + t467 * t506;
t382 = (-t448 * t538 - t447 * t474 + (-t457 * t537 + t505 * t472) * qJD(3)) * t467 + t443 * t554 - t448 * t541 + t457 * qJD(5) * t542;
t381 = t504 * qJD(5) + t406 * t468 - t446 * t542;
t373 = 0.1e1 / t375;
t366 = 0.1e1 / t368;
t365 = t385 * t591;
t362 = t512 * t385;
t361 = t513 * t385;
t358 = t503 * t419;
t356 = (t383 * t586 + t384 * t451) * t467 + t523 * t365;
t354 = t523 * t361 + t383 * t418 + t384 * t428;
t352 = t512 * t551 + (t442 * t530 - t382 * t424 + (t403 * t421 - t411 * t417 - t442 * t599) * t425) * t385;
t350 = t513 * t551 + (t428 * t530 - t379 * t424 + (-t403 * t418 - t404 * t417 - t428 * t599) * t425) * t385;
t349 = t551 * t591 + (t511 * t554 + (t451 * t530 - t593 * t424 + (-t403 * t586 - t417 * t429 - t451 * t599) * t425) * t467) * t385;
t1 = [t547 * t584 + (-t376 * t424 + t419 * t592) * t385, t352, t349, 0, t350, 0; t417 * t534 + (t599 * t369 + (-t353 * t417 - t358 * t376) * t370) * t366 + (t358 * t533 + (0.2e1 * t358 * t576 - (t359 * t385 * t563 + t551) * t545 - (t547 * t585 - t359 + (t359 - t515) * t385) * t544 - t503 * t572) * t366) * t419, -t504 * t534 + ((t423 * qJD(5) + t406 * t467 + t446 * t541) * t369 + t504 * t370 * t353 - ((t352 * t417 + t362 * t599 + t411 + (-t362 * t427 - t421) * t359) * t384 + (-t352 * t427 - t362 * t403 - t382 + (-t362 * t417 - t442) * t359) * t383) * t571) * t366 + (t419 * t533 + (-t572 + t532) * t366) * (t523 * t362 - t383 * t421 + t384 * t442) (t356 * t571 + t369 * t561) * t552 + (-t356 * t572 + (-t397 * t467 - t438 * t554) * t369 + (t356 * t550 + t370 * t561) * t353 - (t451 * t554 + t349 * t417 + t365 * t599 + t429 * t467 + (-t365 * t427 + t467 * t586) * t359) * t544 - (t586 * t554 - t349 * t427 - t365 * t403 - t593 * t467 + (-t365 * t417 - t451 * t467) * t359) * t545) * t366, 0 (t354 * t571 - t369 * t420) * t552 + (t354 * t532 + t377 * t369 + (-t420 * t353 - t354 * t376 - (t350 * t417 + t361 * t599 + t404 + (-t361 * t427 + t418) * t359) * t568 - (-t350 * t427 - t361 * t403 - t379 + (-t361 * t417 - t428) * t359) * t569) * t370) * t366, 0; (-t389 * t393 + t394 * t566) * t548 + ((t394 * qJD(6) - t379 * t471 - t473 * t593) * t389 + t394 * t531 + (-t393 * t364 - (-t393 * qJD(6) - t379 * t473 + t471 * t593) * t395 - t394 * t363) * t390) * t373 (-t389 * t407 + t408 * t566) * t548 + ((t408 * qJD(6) + t381 * t471 - t405 * t473) * t389 + t408 * t531 + (-t407 * t364 - (-t407 * qJD(6) + t381 * t473 + t405 * t471) * t395 - t408 * t363) * t390) * t373 (-t389 * t409 + t410 * t566) * t548 + (t410 * t531 - t525 * t389 * t473 + t498 * t567 + (-t525 * t395 * t471 - t410 * t363 - t409 * t364 - t498 * t565) * t390) * t373, 0, t514 * t419 * t549 + (t514 * t376 + ((-qJD(6) * t389 - 0.2e1 * t546) * t473 + (t363 * t473 + (t364 - t553) * t471) * t390) * t419) * t373, t549 + 0.2e1 * (t363 * t390 * t373 + (-t373 * t573 - t390 * t570) * t395) * t395;];
JaD_rot  = t1;
