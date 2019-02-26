% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function JaD_rot = S6RRRRPR12_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR12_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_jacobiaD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:37:00
% EndTime: 2019-02-26 22:37:06
% DurationCPUTime: 5.90s
% Computational Cost: add. (25275->242), mult. (59163->461), div. (983->12), fcn. (74815->17), ass. (0->189)
t581 = cos(pkin(6));
t582 = sin(qJ(2));
t529 = t581 * t582;
t476 = cos(qJ(2));
t583 = sin(qJ(1));
t544 = t583 * t476;
t584 = cos(qJ(1));
t458 = t584 * t529 + t544;
t495 = t581 * t544 + t584 * t582;
t448 = qJD(1) * t495 + qJD(2) * t458;
t467 = t583 * t582;
t510 = t583 * t529;
t530 = t581 * t584;
t449 = -qJD(1) * t510 - qJD(2) * t467 + (qJD(1) * t584 + qJD(2) * t530) * t476;
t473 = sin(qJ(3));
t475 = cos(qJ(3));
t471 = sin(pkin(6));
t579 = sin(pkin(7));
t540 = t471 * t579;
t518 = t583 * t540;
t508 = qJD(1) * t518;
t580 = cos(pkin(7));
t521 = t584 * t540;
t506 = -t476 * t530 + t467;
t591 = t506 * t580;
t487 = t591 + t521;
t587 = t458 * t473 + t475 * t487;
t403 = t587 * qJD(3) - (-t448 * t580 + t508) * t473 - t449 * t475;
t557 = t458 * t475;
t438 = t473 * t487 - t557;
t470 = qJ(4) + pkin(13);
t468 = sin(t470);
t469 = cos(t470);
t541 = t471 * t580;
t522 = t584 * t541;
t488 = t506 * t579 - t522;
t418 = t438 * t468 + t488 * t469;
t519 = t583 * t541;
t493 = qJD(1) * t519 + t448 * t579;
t380 = qJD(4) * t418 - t403 * t469 + t468 * t493;
t419 = t438 * t469 - t488 * t468;
t600 = qJD(4) * t419 + t403 * t468 + t469 * t493;
t597 = -t495 * t580 + t518;
t538 = t475 * t580;
t556 = qJD(3) * t473;
t594 = (t473 * t591 - t557) * qJD(3) - t448 * t538 - t449 * t473 + t475 * t508 + t521 * t556;
t537 = t476 * t580;
t503 = t473 * t537 + t582 * t475;
t525 = t581 * t579;
t453 = t471 * t503 + t473 * t525;
t536 = t476 * t579;
t457 = -t471 * t536 + t580 * t581;
t429 = t453 * t469 + t457 * t468;
t528 = t580 * t582;
t501 = -t473 * t528 + t475 * t476;
t502 = -t582 * t473 + t475 * t537;
t517 = qJD(3) * t525;
t431 = t475 * t517 + (qJD(2) * t501 + qJD(3) * t502) * t471;
t527 = t579 * t582;
t523 = t471 * t527;
t507 = qJD(2) * t523;
t404 = qJD(4) * t429 + t431 * t468 - t469 * t507;
t428 = t453 * t468 - t457 * t469;
t426 = 0.1e1 / t428 ^ 2;
t593 = t404 * t426;
t425 = 0.1e1 / t428;
t452 = t471 * t502 + t475 * t525;
t563 = t418 * t426;
t512 = t425 * t587 - t452 * t563;
t592 = t468 * t512;
t494 = -t584 * t476 + t510;
t485 = t506 * qJD(1) + t494 * qJD(2);
t590 = qJD(1) * t521 + t485 * t580;
t589 = t597 * t475;
t393 = atan2(t418, t428);
t384 = sin(t393);
t385 = cos(t393);
t373 = t384 * t418 + t385 * t428;
t370 = 0.1e1 / t373;
t440 = t597 * t473 - t494 * t475;
t486 = t495 * t579 + t519;
t421 = t440 * t469 + t468 * t486;
t474 = cos(qJ(6));
t439 = -t473 * t494 - t589;
t472 = sin(qJ(6));
t561 = t439 * t472;
t397 = t421 * t474 + t561;
t390 = 0.1e1 / t397;
t371 = 0.1e1 / t373 ^ 2;
t391 = 0.1e1 / t397 ^ 2;
t586 = 0.2e1 * t418;
t420 = t440 * t468 - t469 * t486;
t585 = 0.2e1 * t420;
t414 = t420 ^ 2;
t369 = t371 * t414 + 0.1e1;
t447 = qJD(1) * t458 + qJD(2) * t495;
t399 = t589 * qJD(3) - t447 * t475 + t590 * t473 + t494 * t556;
t483 = qJD(1) * t522 - t485 * t579;
t377 = qJD(4) * t421 + t399 * t468 - t469 * t483;
t571 = t371 * t420;
t413 = t418 ^ 2;
t388 = t413 * t426 + 0.1e1;
t386 = 0.1e1 / t388;
t516 = -t404 * t563 + t425 * t600;
t360 = t516 * t386;
t524 = -t384 * t428 + t385 * t418;
t354 = t360 * t524 + t384 * t600 + t385 * t404;
t372 = t370 * t371;
t577 = t354 * t372;
t578 = (t377 * t571 - t414 * t577) / t369 ^ 2;
t378 = -qJD(4) * t420 + t399 * t469 + t468 * t483;
t398 = t440 * qJD(3) - t447 * t473 - t590 * t475;
t364 = qJD(6) * t397 + t378 * t472 - t398 * t474;
t560 = t439 * t474;
t396 = t421 * t472 - t560;
t389 = t396 ^ 2;
t376 = t389 * t391 + 0.1e1;
t567 = t391 * t396;
t554 = qJD(6) * t396;
t365 = t378 * t474 + t398 * t472 - t554;
t573 = t365 * t390 * t391;
t575 = (t364 * t567 - t389 * t573) / t376 ^ 2;
t565 = t425 * t593;
t574 = (-t413 * t565 + t563 * t600) / t388 ^ 2;
t572 = t371 * t377;
t570 = t384 * t420;
t569 = t385 * t420;
t568 = t390 * t472;
t566 = t396 * t474;
t564 = t418 * t425;
t562 = t439 * t468;
t555 = qJD(4) * t469;
t553 = 0.2e1 * t578;
t552 = -0.2e1 * t575;
t551 = 0.2e1 * t575;
t550 = -0.2e1 * t574;
t549 = t372 * t585;
t548 = t425 * t574;
t547 = t396 * t573;
t546 = t371 * t570;
t545 = t371 * t569;
t543 = t468 * t579;
t542 = t469 * t579;
t539 = t473 * t580;
t535 = -0.2e1 * t370 * t578;
t534 = t371 * t553;
t533 = t354 * t549;
t532 = 0.2e1 * t547;
t531 = t565 * t586;
t526 = qJD(6) * t439 * t469 + t399;
t395 = t419 * t474 - t472 * t587;
t394 = t419 * t472 + t474 * t587;
t446 = -t495 * t475 + t494 * t539;
t424 = t446 * t469 - t494 * t543;
t445 = -t495 * t473 - t494 * t538;
t409 = t424 * t474 + t445 * t472;
t408 = t424 * t472 - t445 * t474;
t515 = t391 * t566 - t568;
t514 = t419 * t425 - t429 * t563;
t444 = -t458 * t539 - t475 * t506;
t422 = t444 * t468 - t458 * t542;
t456 = t501 * t471;
t443 = t456 * t468 - t469 * t523;
t513 = -t422 * t425 - t443 * t563;
t505 = -t446 * t468 - t494 * t542;
t504 = -t384 + (-t385 * t564 + t384) * t386;
t500 = -t473 * t476 - t475 * t528;
t499 = qJD(4) * t562 + qJD(6) * t440 - t398 * t469;
t430 = -t473 * t517 + (qJD(2) * t500 - qJD(3) * t503) * t471;
t412 = t456 * t555 + ((qJD(3) * t500 + qJD(4) * t527) * t468 + (-t468 * t503 - t469 * t536) * qJD(2)) * t471;
t411 = t440 * t472 - t469 * t560;
t410 = -t440 * t474 - t469 * t561;
t407 = -qJD(3) * t445 + t447 * t539 + t475 * t485;
t406 = qJD(3) * t446 - t447 * t538 + t473 * t485;
t405 = -qJD(4) * t428 + t431 * t469 + t468 * t507;
t383 = (-t449 * t539 - t448 * t475 + (-t458 * t538 + t473 * t506) * qJD(3)) * t468 + t444 * t555 - t449 * t542 + t458 * qJD(4) * t543;
t382 = qJD(4) * t505 + t407 * t469 - t447 * t543;
t374 = 0.1e1 / t376;
t367 = 0.1e1 / t369;
t366 = t386 * t592;
t363 = t513 * t386;
t362 = t514 * t386;
t359 = t504 * t420;
t357 = (t384 * t587 + t385 * t452) * t468 + t524 * t366;
t355 = t362 * t524 + t384 * t419 + t385 * t429;
t353 = t513 * t550 + (t443 * t531 - t383 * t425 + (t404 * t422 - t412 * t418 - t443 * t600) * t426) * t386;
t351 = t514 * t550 + (t429 * t531 - t380 * t425 + (-t404 * t419 - t405 * t418 - t429 * t600) * t426) * t386;
t350 = t550 * t592 + (t512 * t555 + (t452 * t531 - t594 * t425 + (-t404 * t587 - t418 * t430 - t452 * t600) * t426) * t468) * t386;
t1 = [t548 * t585 + (-t377 * t425 + t420 * t593) * t386, t353, t350, t351, 0, 0; t418 * t535 + (t600 * t370 + (-t354 * t418 - t359 * t377) * t371) * t367 + (t359 * t534 + (0.2e1 * t359 * t577 - (t360 * t386 * t564 + t550) * t546 - (t548 * t586 - t360 + (t360 - t516) * t386) * t545 - t504 * t572) * t367) * t420, -t505 * t535 + ((qJD(4) * t424 + t407 * t468 + t447 * t542) * t370 + t505 * t371 * t354 - ((t353 * t418 + t363 * t600 + t412 + (-t363 * t428 - t422) * t360) * t385 + (-t353 * t428 - t363 * t404 - t383 + (-t363 * t418 - t443) * t360) * t384) * t571) * t367 + (t420 * t534 + (-t572 + t533) * t367) * (t363 * t524 - t384 * t422 + t385 * t443) (t357 * t571 + t370 * t562) * t553 + (-t357 * t572 + (-t398 * t468 - t439 * t555) * t370 + (t357 * t549 + t371 * t562) * t354 - (t452 * t555 + t350 * t418 + t366 * t600 + t430 * t468 + (-t366 * t428 + t468 * t587) * t360) * t545 - (t587 * t555 - t350 * t428 - t366 * t404 - t594 * t468 + (-t366 * t418 - t452 * t468) * t360) * t546) * t367 (t355 * t571 - t370 * t421) * t553 + (t355 * t533 + t378 * t370 + (-t421 * t354 - t355 * t377 - (t351 * t418 + t362 * t600 + t405 + (-t362 * t428 + t419) * t360) * t569 - (-t351 * t428 - t362 * t404 - t380 + (-t362 * t418 - t429) * t360) * t570) * t371) * t367, 0, 0; (-t390 * t394 + t395 * t567) * t551 + ((qJD(6) * t395 - t380 * t472 - t474 * t594) * t390 + t395 * t532 + (-t394 * t365 - (-qJD(6) * t394 - t380 * t474 + t472 * t594) * t396 - t395 * t364) * t391) * t374 (-t390 * t408 + t409 * t567) * t551 + ((qJD(6) * t409 + t382 * t472 - t406 * t474) * t390 + t409 * t532 + (-t408 * t365 - (-qJD(6) * t408 + t382 * t474 + t406 * t472) * t396 - t409 * t364) * t391) * t374 (-t390 * t410 + t411 * t567) * t551 + (t411 * t532 - t526 * t390 * t474 + t499 * t568 + (-t396 * t472 * t526 - t411 * t364 - t410 * t365 - t499 * t566) * t391) * t374, t515 * t420 * t552 + (t515 * t377 + ((-qJD(6) * t390 - 0.2e1 * t547) * t474 + (t364 * t474 + (t365 - t554) * t472) * t391) * t420) * t374, 0, t552 + 0.2e1 * (t364 * t374 * t391 + (-t374 * t573 - t391 * t575) * t396) * t396;];
JaD_rot  = t1;
