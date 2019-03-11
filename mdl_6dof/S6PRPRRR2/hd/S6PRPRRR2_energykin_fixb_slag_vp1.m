% Calculate kinetic energy for
% S6PRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRR2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:26:39
% EndTime: 2019-03-08 20:26:42
% DurationCPUTime: 3.16s
% Computational Cost: add. (3698->374), mult. (8963->585), div. (0->0), fcn. (11426->14), ass. (0->172)
t530 = sin(pkin(11));
t532 = cos(pkin(11));
t536 = sin(qJ(2));
t568 = sin(pkin(12));
t569 = cos(pkin(12));
t574 = cos(qJ(2));
t516 = -t536 * t568 + t574 * t569;
t533 = cos(pkin(6));
t543 = t533 * t516;
t545 = t536 * t569 + t568 * t574;
t493 = -t530 * t545 + t532 * t543;
t509 = t545 * t533;
t494 = t509 * t532 + t516 * t530;
t531 = sin(pkin(6));
t563 = t531 * t532;
t439 = Icges(4,5) * t494 + Icges(4,6) * t493 - Icges(4,3) * t563;
t495 = -t530 * t543 - t532 * t545;
t496 = -t509 * t530 + t516 * t532;
t564 = t530 * t531;
t440 = Icges(4,5) * t496 + Icges(4,6) * t495 + Icges(4,3) * t564;
t507 = t516 * t531;
t508 = t545 * t531;
t476 = Icges(4,5) * t508 + Icges(4,6) * t507 + Icges(4,3) * t533;
t555 = t533 * t574;
t511 = -t530 * t536 + t532 * t555;
t561 = t533 * t536;
t512 = t530 * t574 + t532 * t561;
t480 = Icges(3,5) * t512 + Icges(3,6) * t511 - Icges(3,3) * t563;
t513 = -t530 * t555 - t532 * t536;
t514 = -t530 * t561 + t532 * t574;
t481 = Icges(3,5) * t514 + Icges(3,6) * t513 + Icges(3,3) * t564;
t503 = Icges(3,3) * t533 + (Icges(3,5) * t536 + Icges(3,6) * t574) * t531;
t579 = -(t476 + t503) * t533 + (t480 + t439) * t563 - (t481 + t440) * t564;
t575 = qJD(2) ^ 2;
t573 = cos(qJ(4));
t572 = pkin(2) * t574;
t537 = cos(qJ(5));
t571 = pkin(5) * t537;
t534 = sin(qJ(5));
t567 = t493 * t534;
t566 = t495 * t534;
t565 = t507 * t534;
t535 = sin(qJ(4));
t562 = t531 * t535;
t517 = pkin(2) * t531 * t536 + qJ(3) * t533;
t560 = -pkin(3) * t508 + pkin(8) * t507 - t517;
t559 = qJD(2) * t531;
t521 = t530 * t559;
t471 = -qJD(4) * t495 + t521;
t526 = qJD(2) * t533;
t500 = -qJD(4) * t507 + t526;
t558 = qJD(3) * t532;
t556 = t531 * t573;
t469 = t496 * t535 - t530 * t556;
t428 = qJD(5) * t469 + t471;
t498 = t508 * t535 - t533 * t573;
t462 = qJD(5) * t498 + t500;
t554 = t532 * t559;
t553 = pkin(2) * t561 - qJ(3) * t531;
t551 = (-rSges(4,1) * t508 - rSges(4,2) * t507 - rSges(4,3) * t533 - t517) * t531;
t488 = t530 * t572 + t532 * t553;
t489 = -t530 * t553 + t532 * t572;
t550 = qJD(3) * t533 + t488 * t521 + t489 * t554 + qJD(1);
t472 = -qJD(4) * t493 - t554;
t467 = t494 * t535 + t532 * t556;
t429 = qJD(5) * t467 + t472;
t451 = pkin(3) * t494 - pkin(8) * t493;
t452 = pkin(3) * t496 - pkin(8) * t495;
t546 = t451 * t521 + t452 * t554 + t550;
t468 = t494 * t573 - t532 * t562;
t426 = t468 * pkin(4) + t467 * pkin(9);
t470 = t496 * t573 + t530 * t562;
t427 = t470 * pkin(4) + t469 * pkin(9);
t544 = t471 * t426 - t427 * t472 + t546;
t475 = t489 * t526;
t542 = t452 * t526 + t475 + (qJD(2) * t530 * t560 - t558) * t531;
t520 = qJD(3) * t564;
t541 = t520 + ((-t451 - t488) * t533 + t560 * t563) * qJD(2);
t499 = t508 * t573 + t533 * t535;
t459 = t499 * pkin(4) + t498 * pkin(9);
t540 = t500 * t427 - t459 * t471 + t542;
t539 = -t426 * t500 + t472 * t459 + t541;
t529 = qJ(5) + qJ(6);
t528 = cos(t529);
t527 = sin(t529);
t506 = t533 * rSges(3,3) + (rSges(3,1) * t536 + rSges(3,2) * t574) * t531;
t505 = Icges(3,5) * t533 + (Icges(3,1) * t536 + Icges(3,4) * t574) * t531;
t504 = Icges(3,6) * t533 + (Icges(3,4) * t536 + Icges(3,2) * t574) * t531;
t487 = rSges(3,1) * t514 + rSges(3,2) * t513 + rSges(3,3) * t564;
t486 = rSges(3,1) * t512 + rSges(3,2) * t511 - rSges(3,3) * t563;
t485 = Icges(3,1) * t514 + Icges(3,4) * t513 + Icges(3,5) * t564;
t484 = Icges(3,1) * t512 + Icges(3,4) * t511 - Icges(3,5) * t563;
t483 = Icges(3,4) * t514 + Icges(3,2) * t513 + Icges(3,6) * t564;
t482 = Icges(3,4) * t512 + Icges(3,2) * t511 - Icges(3,6) * t563;
t478 = Icges(4,1) * t508 + Icges(4,4) * t507 + Icges(4,5) * t533;
t477 = Icges(4,4) * t508 + Icges(4,2) * t507 + Icges(4,6) * t533;
t464 = t499 * t537 - t565;
t463 = -t499 * t534 - t507 * t537;
t461 = t499 * t528 - t507 * t527;
t460 = -t499 * t527 - t507 * t528;
t458 = (-t486 * t533 - t506 * t563) * qJD(2);
t457 = (t487 * t533 - t506 * t564) * qJD(2);
t456 = rSges(5,1) * t499 - rSges(5,2) * t498 - rSges(5,3) * t507;
t455 = Icges(5,1) * t499 - Icges(5,4) * t498 - Icges(5,5) * t507;
t454 = Icges(5,4) * t499 - Icges(5,2) * t498 - Icges(5,6) * t507;
t453 = Icges(5,5) * t499 - Icges(5,6) * t498 - Icges(5,3) * t507;
t449 = qJD(6) * t498 + t462;
t446 = rSges(4,1) * t496 + rSges(4,2) * t495 + rSges(4,3) * t564;
t445 = rSges(4,1) * t494 + rSges(4,2) * t493 - rSges(4,3) * t563;
t444 = Icges(4,1) * t496 + Icges(4,4) * t495 + Icges(4,5) * t564;
t443 = Icges(4,1) * t494 + Icges(4,4) * t493 - Icges(4,5) * t563;
t442 = Icges(4,4) * t496 + Icges(4,2) * t495 + Icges(4,6) * t564;
t441 = Icges(4,4) * t494 + Icges(4,2) * t493 - Icges(4,6) * t563;
t438 = t470 * t537 - t566;
t437 = -t470 * t534 - t495 * t537;
t436 = t468 * t537 - t567;
t435 = -t468 * t534 - t493 * t537;
t434 = t470 * t528 - t495 * t527;
t433 = -t470 * t527 - t495 * t528;
t432 = t468 * t528 - t493 * t527;
t431 = -t468 * t527 - t493 * t528;
t430 = qJD(1) + (t486 * t530 + t487 * t532) * t559;
t423 = rSges(6,1) * t464 + rSges(6,2) * t463 + rSges(6,3) * t498;
t422 = Icges(6,1) * t464 + Icges(6,4) * t463 + Icges(6,5) * t498;
t421 = Icges(6,4) * t464 + Icges(6,2) * t463 + Icges(6,6) * t498;
t420 = Icges(6,5) * t464 + Icges(6,6) * t463 + Icges(6,3) * t498;
t419 = rSges(5,1) * t470 - rSges(5,2) * t469 - rSges(5,3) * t495;
t418 = rSges(5,1) * t468 - rSges(5,2) * t467 - rSges(5,3) * t493;
t417 = Icges(5,1) * t470 - Icges(5,4) * t469 - Icges(5,5) * t495;
t416 = Icges(5,1) * t468 - Icges(5,4) * t467 - Icges(5,5) * t493;
t415 = Icges(5,4) * t470 - Icges(5,2) * t469 - Icges(5,6) * t495;
t414 = Icges(5,4) * t468 - Icges(5,2) * t467 - Icges(5,6) * t493;
t413 = Icges(5,5) * t470 - Icges(5,6) * t469 - Icges(5,3) * t495;
t412 = Icges(5,5) * t468 - Icges(5,6) * t467 - Icges(5,3) * t493;
t411 = -pkin(5) * t565 + pkin(10) * t498 + t499 * t571;
t410 = rSges(7,1) * t461 + rSges(7,2) * t460 + rSges(7,3) * t498;
t409 = Icges(7,1) * t461 + Icges(7,4) * t460 + Icges(7,5) * t498;
t408 = Icges(7,4) * t461 + Icges(7,2) * t460 + Icges(7,6) * t498;
t407 = Icges(7,5) * t461 + Icges(7,6) * t460 + Icges(7,3) * t498;
t406 = qJD(6) * t467 + t429;
t405 = qJD(6) * t469 + t428;
t403 = t520 + ((-t445 - t488) * t533 + t532 * t551) * qJD(2);
t402 = -t531 * t558 + t475 + (t446 * t533 + t530 * t551) * qJD(2);
t401 = rSges(6,1) * t438 + rSges(6,2) * t437 + rSges(6,3) * t469;
t400 = rSges(6,1) * t436 + rSges(6,2) * t435 + rSges(6,3) * t467;
t399 = Icges(6,1) * t438 + Icges(6,4) * t437 + Icges(6,5) * t469;
t398 = Icges(6,1) * t436 + Icges(6,4) * t435 + Icges(6,5) * t467;
t397 = Icges(6,4) * t438 + Icges(6,2) * t437 + Icges(6,6) * t469;
t396 = Icges(6,4) * t436 + Icges(6,2) * t435 + Icges(6,6) * t467;
t395 = Icges(6,5) * t438 + Icges(6,6) * t437 + Icges(6,3) * t469;
t394 = Icges(6,5) * t436 + Icges(6,6) * t435 + Icges(6,3) * t467;
t393 = rSges(7,1) * t434 + rSges(7,2) * t433 + rSges(7,3) * t469;
t392 = rSges(7,1) * t432 + rSges(7,2) * t431 + rSges(7,3) * t467;
t391 = Icges(7,1) * t434 + Icges(7,4) * t433 + Icges(7,5) * t469;
t390 = Icges(7,1) * t432 + Icges(7,4) * t431 + Icges(7,5) * t467;
t389 = Icges(7,4) * t434 + Icges(7,2) * t433 + Icges(7,6) * t469;
t388 = Icges(7,4) * t432 + Icges(7,2) * t431 + Icges(7,6) * t467;
t387 = Icges(7,5) * t434 + Icges(7,6) * t433 + Icges(7,3) * t469;
t386 = Icges(7,5) * t432 + Icges(7,6) * t431 + Icges(7,3) * t467;
t385 = -pkin(5) * t566 + pkin(10) * t469 + t470 * t571;
t384 = -pkin(5) * t567 + pkin(10) * t467 + t468 * t571;
t383 = (t445 * t530 + t446 * t532) * t559 + t550;
t382 = -t418 * t500 + t456 * t472 + t541;
t381 = t419 * t500 - t456 * t471 + t542;
t380 = t418 * t471 - t419 * t472 + t546;
t379 = -t400 * t462 + t423 * t429 + t539;
t378 = t401 * t462 - t423 * t428 + t540;
t377 = t400 * t428 - t401 * t429 + t544;
t376 = -t384 * t462 - t392 * t449 + t406 * t410 + t411 * t429 + t539;
t375 = t385 * t462 + t393 * t449 - t405 * t410 - t411 * t428 + t540;
t374 = t384 * t428 - t385 * t429 + t392 * t405 - t393 * t406 + t544;
t1 = m(7) * (t374 ^ 2 + t375 ^ 2 + t376 ^ 2) / 0.2e1 + t405 * ((t469 * t387 + t433 * t389 + t434 * t391) * t405 + (t386 * t469 + t388 * t433 + t390 * t434) * t406 + (t407 * t469 + t408 * t433 + t409 * t434) * t449) / 0.2e1 + t406 * ((t387 * t467 + t389 * t431 + t391 * t432) * t405 + (t467 * t386 + t431 * t388 + t432 * t390) * t406 + (t407 * t467 + t408 * t431 + t409 * t432) * t449) / 0.2e1 + t428 * ((t395 * t469 + t397 * t437 + t399 * t438) * t428 + (t394 * t469 + t396 * t437 + t398 * t438) * t429 + (t420 * t469 + t421 * t437 + t422 * t438) * t462) / 0.2e1 + t429 * ((t395 * t467 + t397 * t435 + t399 * t436) * t428 + (t394 * t467 + t396 * t435 + t398 * t436) * t429 + (t420 * t467 + t421 * t435 + t422 * t436) * t462) / 0.2e1 + t462 * ((t395 * t498 + t397 * t463 + t399 * t464) * t428 + (t394 * t498 + t396 * t463 + t398 * t464) * t429 + (t420 * t498 + t421 * t463 + t422 * t464) * t462) / 0.2e1 + t471 * ((-t413 * t495 - t415 * t469 + t417 * t470) * t471 + (-t412 * t495 - t414 * t469 + t416 * t470) * t472 + (-t453 * t495 - t454 * t469 + t455 * t470) * t500) / 0.2e1 + t472 * ((-t413 * t493 - t415 * t467 + t417 * t468) * t471 + (-t412 * t493 - t414 * t467 + t416 * t468) * t472 + (-t453 * t493 - t454 * t467 + t455 * t468) * t500) / 0.2e1 + t500 * ((-t413 * t507 - t415 * t498 + t417 * t499) * t471 + (-t412 * t507 - t414 * t498 + t416 * t499) * t472 + (-t453 * t507 - t454 * t498 + t455 * t499) * t500) / 0.2e1 + ((t533 * t481 + (t483 * t574 + t485 * t536) * t531) * t521 - (t533 * t480 + (t482 * t574 + t484 * t536) * t531) * t554 + (t533 * t503 + (t504 * t574 + t505 * t536) * t531) * t526) * t526 / 0.2e1 + t449 * ((t387 * t498 + t389 * t460 + t391 * t461) * t405 + (t386 * t498 + t388 * t460 + t390 * t461) * t406 + (t407 * t498 + t408 * t460 + t409 * t461) * t449) / 0.2e1 + m(6) * (t377 ^ 2 + t378 ^ 2 + t379 ^ 2) / 0.2e1 + m(5) * (t380 ^ 2 + t381 ^ 2 + t382 ^ 2) / 0.2e1 + m(4) * (t383 ^ 2 + t402 ^ 2 + t403 ^ 2) / 0.2e1 + m(3) * (t430 ^ 2 + t457 ^ 2 + t458 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 - ((t442 * t493 + t444 * t494 + t483 * t511 + t485 * t512) * t564 + (t477 * t493 + t478 * t494 + t504 * t511 + t505 * t512) * t533 + (-t441 * t493 - t443 * t494 - t482 * t511 - t484 * t512 + t579) * t563) * t575 * t563 / 0.2e1 + (t533 * ((t476 * t533 + t477 * t507 + t478 * t508) * t533 + ((t440 * t533 + t442 * t507 + t444 * t508) * t530 - (t439 * t533 + t441 * t507 + t443 * t508) * t532) * t531) + ((-t441 * t495 - t443 * t496 - t482 * t513 - t484 * t514) * t563 + (t477 * t495 + t478 * t496 + t504 * t513 + t505 * t514) * t533 + (t442 * t495 + t444 * t496 + t483 * t513 + t485 * t514 - t579) * t564) * t564) * t575 / 0.2e1;
T  = t1;
