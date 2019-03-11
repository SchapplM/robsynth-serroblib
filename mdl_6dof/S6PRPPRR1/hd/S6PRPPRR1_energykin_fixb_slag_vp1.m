% Calculate kinetic energy for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
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
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPPRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPPRR1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPPRR1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:13:56
% EndTime: 2019-03-08 19:13:59
% DurationCPUTime: 3.23s
% Computational Cost: add. (3444->375), mult. (7510->571), div. (0->0), fcn. (9440->14), ass. (0->170)
t553 = sin(qJ(2));
t589 = sin(pkin(11));
t590 = cos(pkin(11));
t593 = cos(qJ(2));
t532 = -t553 * t590 - t589 * t593;
t546 = sin(pkin(10));
t549 = cos(pkin(10));
t533 = -t553 * t589 + t593 * t590;
t550 = cos(pkin(6));
t557 = t550 * t533;
t508 = t546 * t532 + t549 * t557;
t526 = t532 * t550;
t509 = -t526 * t549 + t533 * t546;
t547 = sin(pkin(6));
t586 = t547 * t549;
t452 = Icges(4,5) * t509 + Icges(4,6) * t508 - Icges(4,3) * t586;
t510 = t532 * t549 - t546 * t557;
t511 = t526 * t546 + t533 * t549;
t587 = t546 * t547;
t453 = Icges(4,5) * t511 + Icges(4,6) * t510 + Icges(4,3) * t587;
t573 = t550 * t593;
t528 = -t546 * t553 + t549 * t573;
t585 = t550 * t553;
t529 = t546 * t593 + t549 * t585;
t493 = Icges(3,5) * t529 + Icges(3,6) * t528 - Icges(3,3) * t586;
t530 = -t546 * t573 - t549 * t553;
t531 = -t546 * t585 + t549 * t593;
t494 = Icges(3,5) * t531 + Icges(3,6) * t530 + Icges(3,3) * t587;
t600 = (t493 + t452) * t586 - (t494 + t453) * t587;
t469 = pkin(3) * t511 - qJ(4) * t510;
t543 = qJD(2) * t550;
t599 = -qJD(4) * t508 + t469 * t543;
t525 = t532 * t547;
t545 = sin(pkin(12));
t548 = cos(pkin(12));
t514 = t525 * t545 + t548 * t550;
t588 = t545 * t550;
t515 = -t525 * t548 + t588;
t524 = t533 * t547;
t598 = Icges(4,4) * t525 + Icges(5,5) * t515 - Icges(4,6) * t550 + Icges(5,6) * t514 + (-Icges(4,2) - Icges(5,3)) * t524;
t489 = -Icges(4,5) * t525 + Icges(4,6) * t524 + Icges(4,3) * t550;
t520 = Icges(3,3) * t550 + (Icges(3,5) * t553 + Icges(3,6) * t593) * t547;
t597 = t489 + t520;
t594 = qJD(2) ^ 2;
t592 = pkin(2) * t593;
t591 = pkin(4) * t548;
t468 = pkin(3) * t509 - qJ(4) * t508;
t570 = pkin(2) * t585 - qJ(3) * t547;
t501 = t546 * t592 + t549 * t570;
t583 = -t468 - t501;
t534 = pkin(2) * t547 * t553 + qJ(3) * t550;
t582 = pkin(3) * t525 + qJ(4) * t524 - t534;
t538 = qJD(3) * t587;
t581 = -qJD(4) * t510 + t538;
t580 = qJD(2) * t547;
t539 = t546 * t580;
t484 = -qJD(5) * t510 + t539;
t516 = -qJD(5) * t524 + t543;
t579 = qJD(3) * t549;
t577 = pkin(12) + qJ(5);
t576 = t545 * t587;
t575 = t545 * t586;
t574 = -pkin(4) * t588 + pkin(8) * t524 + t525 * t591 + t582;
t572 = t549 * t580;
t568 = cos(t577);
t567 = (rSges(4,1) * t525 - rSges(4,2) * t524 - rSges(4,3) * t550 - t534) * t547;
t502 = -t546 * t570 + t549 * t592;
t566 = qJD(3) * t550 + t501 * t539 + t502 * t572 + qJD(1);
t563 = (-rSges(5,1) * t515 - rSges(5,2) * t514 + rSges(5,3) * t524 + t582) * t547;
t485 = -qJD(5) * t508 - t572;
t488 = t502 * t543;
t561 = -t547 * t579 + t488;
t560 = t547 * t568;
t559 = -qJD(4) * t524 + t468 * t539 + t469 * t572 + t566;
t417 = -pkin(4) * t575 - pkin(8) * t508 + t509 * t591;
t418 = pkin(4) * t576 - pkin(8) * t510 + t511 * t591;
t558 = t417 * t539 + t418 * t572 + t559;
t556 = t418 * t543 + t488 + (qJD(2) * t546 * t574 - t579) * t547 + t599;
t555 = ((-t417 + t583) * t550 + t574 * t586) * qJD(2) + t581;
t554 = cos(qJ(6));
t552 = sin(qJ(6));
t544 = sin(t577);
t523 = t550 * rSges(3,3) + (rSges(3,1) * t553 + rSges(3,2) * t593) * t547;
t522 = Icges(3,5) * t550 + (Icges(3,1) * t553 + Icges(3,4) * t593) * t547;
t521 = Icges(3,6) * t550 + (Icges(3,4) * t553 + Icges(3,2) * t593) * t547;
t513 = -t525 * t568 + t550 * t544;
t512 = -t525 * t544 - t550 * t568;
t500 = rSges(3,1) * t531 + rSges(3,2) * t530 + rSges(3,3) * t587;
t499 = rSges(3,1) * t529 + rSges(3,2) * t528 - rSges(3,3) * t586;
t498 = Icges(3,1) * t531 + Icges(3,4) * t530 + Icges(3,5) * t587;
t497 = Icges(3,1) * t529 + Icges(3,4) * t528 - Icges(3,5) * t586;
t496 = Icges(3,4) * t531 + Icges(3,2) * t530 + Icges(3,6) * t587;
t495 = Icges(3,4) * t529 + Icges(3,2) * t528 - Icges(3,6) * t586;
t491 = -Icges(4,1) * t525 + Icges(4,4) * t524 + Icges(4,5) * t550;
t483 = t511 * t548 + t576;
t482 = -t511 * t545 + t548 * t587;
t481 = t509 * t548 - t575;
t480 = -t509 * t545 - t548 * t586;
t479 = t511 * t568 + t544 * t587;
t478 = t511 * t544 - t546 * t560;
t477 = t509 * t568 - t544 * t586;
t476 = t509 * t544 + t549 * t560;
t475 = t513 * t554 - t524 * t552;
t474 = -t513 * t552 - t524 * t554;
t473 = qJD(6) * t512 + t516;
t472 = pkin(5) * t513 + pkin(9) * t512;
t471 = (-t499 * t550 - t523 * t586) * qJD(2);
t470 = (t500 * t550 - t523 * t587) * qJD(2);
t466 = Icges(5,1) * t515 + Icges(5,4) * t514 - Icges(5,5) * t524;
t465 = Icges(5,4) * t515 + Icges(5,2) * t514 - Icges(5,6) * t524;
t462 = rSges(4,1) * t511 + rSges(4,2) * t510 + rSges(4,3) * t587;
t461 = rSges(4,1) * t509 + rSges(4,2) * t508 - rSges(4,3) * t586;
t458 = rSges(6,1) * t513 - rSges(6,2) * t512 - rSges(6,3) * t524;
t457 = Icges(4,1) * t511 + Icges(4,4) * t510 + Icges(4,5) * t587;
t456 = Icges(4,1) * t509 + Icges(4,4) * t508 - Icges(4,5) * t586;
t455 = Icges(4,4) * t511 + Icges(4,2) * t510 + Icges(4,6) * t587;
t454 = Icges(4,4) * t509 + Icges(4,2) * t508 - Icges(4,6) * t586;
t451 = Icges(6,1) * t513 - Icges(6,4) * t512 - Icges(6,5) * t524;
t450 = Icges(6,4) * t513 - Icges(6,2) * t512 - Icges(6,6) * t524;
t449 = Icges(6,5) * t513 - Icges(6,6) * t512 - Icges(6,3) * t524;
t447 = t479 * t554 - t510 * t552;
t446 = -t479 * t552 - t510 * t554;
t445 = t477 * t554 - t508 * t552;
t444 = -t477 * t552 - t508 * t554;
t443 = qJD(1) + (t499 * t546 + t500 * t549) * t580;
t442 = qJD(6) * t476 + t485;
t441 = qJD(6) * t478 + t484;
t440 = pkin(5) * t479 + pkin(9) * t478;
t439 = pkin(5) * t477 + pkin(9) * t476;
t438 = rSges(5,1) * t483 + rSges(5,2) * t482 - rSges(5,3) * t510;
t437 = rSges(5,1) * t481 + rSges(5,2) * t480 - rSges(5,3) * t508;
t436 = Icges(5,1) * t483 + Icges(5,4) * t482 - Icges(5,5) * t510;
t435 = Icges(5,1) * t481 + Icges(5,4) * t480 - Icges(5,5) * t508;
t434 = Icges(5,4) * t483 + Icges(5,2) * t482 - Icges(5,6) * t510;
t433 = Icges(5,4) * t481 + Icges(5,2) * t480 - Icges(5,6) * t508;
t432 = Icges(5,5) * t483 + Icges(5,6) * t482 - Icges(5,3) * t510;
t431 = Icges(5,5) * t481 + Icges(5,6) * t480 - Icges(5,3) * t508;
t430 = rSges(7,1) * t475 + rSges(7,2) * t474 + rSges(7,3) * t512;
t429 = Icges(7,1) * t475 + Icges(7,4) * t474 + Icges(7,5) * t512;
t428 = Icges(7,4) * t475 + Icges(7,2) * t474 + Icges(7,6) * t512;
t427 = Icges(7,5) * t475 + Icges(7,6) * t474 + Icges(7,3) * t512;
t426 = rSges(6,1) * t479 - rSges(6,2) * t478 - rSges(6,3) * t510;
t425 = rSges(6,1) * t477 - rSges(6,2) * t476 - rSges(6,3) * t508;
t424 = Icges(6,1) * t479 - Icges(6,4) * t478 - Icges(6,5) * t510;
t423 = Icges(6,1) * t477 - Icges(6,4) * t476 - Icges(6,5) * t508;
t422 = Icges(6,4) * t479 - Icges(6,2) * t478 - Icges(6,6) * t510;
t421 = Icges(6,4) * t477 - Icges(6,2) * t476 - Icges(6,6) * t508;
t420 = Icges(6,5) * t479 - Icges(6,6) * t478 - Icges(6,3) * t510;
t419 = Icges(6,5) * t477 - Icges(6,6) * t476 - Icges(6,3) * t508;
t413 = t538 + ((-t461 - t501) * t550 + t549 * t567) * qJD(2);
t412 = (t462 * t550 + t546 * t567) * qJD(2) + t561;
t411 = rSges(7,1) * t447 + rSges(7,2) * t446 + rSges(7,3) * t478;
t410 = rSges(7,1) * t445 + rSges(7,2) * t444 + rSges(7,3) * t476;
t409 = Icges(7,1) * t447 + Icges(7,4) * t446 + Icges(7,5) * t478;
t408 = Icges(7,1) * t445 + Icges(7,4) * t444 + Icges(7,5) * t476;
t407 = Icges(7,4) * t447 + Icges(7,2) * t446 + Icges(7,6) * t478;
t406 = Icges(7,4) * t445 + Icges(7,2) * t444 + Icges(7,6) * t476;
t405 = Icges(7,5) * t447 + Icges(7,6) * t446 + Icges(7,3) * t478;
t404 = Icges(7,5) * t445 + Icges(7,6) * t444 + Icges(7,3) * t476;
t403 = (t461 * t546 + t462 * t549) * t580 + t566;
t402 = ((-t437 + t583) * t550 + t549 * t563) * qJD(2) + t581;
t401 = (t438 * t550 + t546 * t563) * qJD(2) + t561 + t599;
t400 = (t437 * t546 + t438 * t549) * t580 + t559;
t399 = -t425 * t516 + t458 * t485 + t555;
t398 = t426 * t516 - t458 * t484 + t556;
t397 = t425 * t484 - t426 * t485 + t558;
t396 = -t410 * t473 + t430 * t442 - t439 * t516 + t472 * t485 + t555;
t395 = t411 * t473 - t430 * t441 + t440 * t516 - t472 * t484 + t556;
t394 = t410 * t441 - t411 * t442 + t439 * t484 - t440 * t485 + t558;
t1 = t484 * ((-t420 * t510 - t422 * t478 + t424 * t479) * t484 + (-t419 * t510 - t421 * t478 + t423 * t479) * t485 + (-t449 * t510 - t450 * t478 + t451 * t479) * t516) / 0.2e1 + t485 * ((-t420 * t508 - t422 * t476 + t424 * t477) * t484 + (-t419 * t508 - t421 * t476 + t423 * t477) * t485 + (-t449 * t508 - t450 * t476 + t451 * t477) * t516) / 0.2e1 + t516 * ((-t420 * t524 - t422 * t512 + t424 * t513) * t484 + (-t419 * t524 - t421 * t512 + t423 * t513) * t485 + (-t449 * t524 - t450 * t512 + t451 * t513) * t516) / 0.2e1 + t441 * ((t478 * t405 + t446 * t407 + t447 * t409) * t441 + (t404 * t478 + t406 * t446 + t408 * t447) * t442 + (t427 * t478 + t428 * t446 + t429 * t447) * t473) / 0.2e1 + t442 * ((t405 * t476 + t407 * t444 + t409 * t445) * t441 + (t476 * t404 + t444 * t406 + t445 * t408) * t442 + (t427 * t476 + t428 * t444 + t429 * t445) * t473) / 0.2e1 + t473 * ((t405 * t512 + t407 * t474 + t409 * t475) * t441 + (t404 * t512 + t406 * t474 + t408 * t475) * t442 + (t427 * t512 + t428 * t474 + t429 * t475) * t473) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t443 ^ 2 + t470 ^ 2 + t471 ^ 2) / 0.2e1 + m(4) * (t403 ^ 2 + t412 ^ 2 + t413 ^ 2) / 0.2e1 + m(5) * (t400 ^ 2 + t401 ^ 2 + t402 ^ 2) / 0.2e1 + m(6) * (t397 ^ 2 + t398 ^ 2 + t399 ^ 2) / 0.2e1 + m(7) * (t394 ^ 2 + t395 ^ 2 + t396 ^ 2) / 0.2e1 + ((t550 * t494 + (t496 * t593 + t498 * t553) * t547) * t539 - (t550 * t493 + (t495 * t593 + t497 * t553) * t547) * t572 + (t550 * t520 + (t521 * t593 + t522 * t553) * t547) * t543) * t543 / 0.2e1 - (((-t432 * t508 + t434 * t480 + t436 * t481) * t546 - (-t431 * t508 + t433 * t480 + t435 * t481) * t549) * t547 + (t455 * t508 + t457 * t509 + t496 * t528 + t498 * t529) * t587 + (-t454 * t508 - t456 * t509 - t495 * t528 - t497 * t529 + t600) * t586 + (t465 * t480 + t466 * t481 + t491 * t509 - t508 * t598 + t521 * t528 + t522 * t529 - t586 * t597) * t550) * t594 * t586 / 0.2e1 + (((t465 * t514 + t466 * t515 + t489 * t550 - t491 * t525 - t598 * t524) * t550 + ((-t433 * t514 - t435 * t515 - t452 * t550 + t456 * t525 - (-t431 + t454) * t524) * t549 + (t434 * t514 + t436 * t515 + t453 * t550 - t457 * t525 - (t432 - t455) * t524) * t546) * t547) * t550 + (((-t432 * t510 + t434 * t482 + t436 * t483) * t546 - (-t431 * t510 + t433 * t482 + t435 * t483) * t549) * t547 + (-t454 * t510 - t456 * t511 - t495 * t530 - t497 * t531) * t586 + (t455 * t510 + t457 * t511 + t496 * t530 + t498 * t531 - t600) * t587 + (t465 * t482 + t466 * t483 + t491 * t511 - t510 * t598 + t521 * t530 + t522 * t531 + t587 * t597) * t550) * t587) * t594 / 0.2e1;
T  = t1;
