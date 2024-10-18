% Calculate kinetic energy for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% m [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR12_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR12_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR12_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR12_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR12_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:07:13
% EndTime: 2024-09-28 18:07:14
% DurationCPUTime: 0.45s
% Computational Cost: add. (3246->367), mult. (3785->574), div. (0->0), fcn. (4022->26), ass. (0->197)
t579 = sin(pkin(11));
t578 = qJ(2) + qJ(3);
t574 = cos(t578);
t615 = pkin(3) * t574;
t636 = t579 * t615;
t582 = cos(pkin(11));
t635 = t582 * t615;
t633 = qJD(2) ^ 2;
t593 = pkin(7) + pkin(8);
t580 = sin(pkin(6));
t632 = pkin(10) * t580;
t592 = cos(qJ(2));
t631 = t592 * pkin(2);
t591 = cos(qJ(3));
t566 = pkin(3) * t591 + pkin(2);
t630 = -pkin(2) + t566;
t575 = qJ(4) + t578;
t564 = sin(t575);
t629 = t564 * t579;
t581 = sin(pkin(5));
t628 = t579 * t581;
t627 = t580 * t581;
t626 = t580 * t582;
t625 = t582 * t581;
t583 = cos(pkin(6));
t585 = sin(qJ(5));
t624 = t583 * t585;
t589 = cos(qJ(5));
t623 = t583 * t589;
t584 = cos(pkin(5));
t622 = t584 * t579;
t621 = t584 * t582;
t620 = t584 * t585;
t588 = sin(qJ(2));
t619 = t584 * t588;
t618 = t584 * t589;
t617 = t584 * t592;
t577 = pkin(9) + t593;
t587 = sin(qJ(3));
t608 = t592 * t587 * pkin(3);
t512 = -t577 * t581 + (t588 * t566 + t608) * t584;
t609 = pkin(2) * t619;
t544 = -t581 * t593 + t609;
t616 = t512 - t544;
t614 = qJD(2) * t581;
t553 = t579 * t614;
t529 = qJD(3) * t628 + t553;
t570 = qJD(2) * t584;
t548 = qJD(3) * t584 + t570;
t612 = -qJD(2) - qJD(3);
t611 = t579 * t632;
t610 = pkin(10) * t626;
t607 = t585 * t627;
t606 = t589 * t627;
t605 = t583 * t620;
t604 = t583 * t618;
t517 = qJD(4) * t628 + t529;
t535 = qJD(4) * t584 + t548;
t602 = pkin(7) * t581 + t544;
t497 = t579 * t631 + t582 * t602;
t498 = -t579 * t602 + t582 * t631;
t603 = t582 * t498 * t614 + t497 * t553 + qJD(1);
t572 = pkin(5) - t578;
t571 = pkin(5) + t578;
t524 = t581 * t588 * pkin(2) + (-pkin(7) + t593) * t584;
t601 = t498 * t570 - t524 * t553;
t518 = (-qJD(4) + t612) * t625;
t454 = t582 * t616 + t636;
t455 = -t579 * t616 + t635;
t530 = t612 * t625;
t600 = t529 * t454 - t455 * t530 + t603;
t554 = pkin(10) * t583 + t577;
t599 = -t554 * t581 - t512 + t609;
t536 = pkin(4) * t582 + t584 * t611;
t538 = pkin(4) * t622 - t610;
t586 = sin(qJ(4));
t590 = cos(qJ(4));
t598 = pkin(3) * t582 + t536 * t590 - t538 * t586;
t537 = -pkin(4) * t579 + t584 * t610;
t539 = pkin(4) * t621 + t611;
t597 = pkin(3) * t579 - t537 * t590 + t539 * t586;
t596 = (-t497 * t584 - t524 * t625) * qJD(2);
t495 = (t577 - t593) * t584 + (t588 * t630 + t608) * t581;
t595 = t548 * t455 - t495 * t529 + t601;
t594 = -t454 * t548 + t530 * t495 + t596;
t573 = sin(t578);
t565 = cos(t575);
t563 = -qJ(4) + t572;
t562 = qJ(4) + t571;
t561 = cos(t571);
t560 = sin(t572);
t559 = cos(t562);
t558 = sin(t563);
t556 = cos(t572) / 0.2e1;
t555 = sin(t571) / 0.2e1;
t550 = cos(t563) / 0.2e1;
t549 = sin(t562) / 0.2e1;
t546 = -pkin(4) * t586 + t590 * t632;
t545 = pkin(4) * t590 + t586 * t632 + pkin(3);
t543 = t556 - t561 / 0.2e1;
t542 = t556 + t561 / 0.2e1;
t541 = t555 - t560 / 0.2e1;
t540 = t555 + t560 / 0.2e1;
t534 = -t579 * t619 + t582 * t592;
t533 = -t579 * t617 - t582 * t588;
t532 = t579 * t592 + t582 * t619;
t531 = -t579 * t588 + t582 * t617;
t528 = t550 - t559 / 0.2e1;
t527 = t550 + t559 / 0.2e1;
t526 = t549 - t558 / 0.2e1;
t525 = t549 + t558 / 0.2e1;
t523 = -t565 * t627 + t583 * t584;
t522 = rSges(3,3) * t584 + (rSges(3,1) * t588 + rSges(3,2) * t592) * t581;
t521 = Icges(3,5) * t584 + (Icges(3,1) * t588 + Icges(3,4) * t592) * t581;
t520 = Icges(3,6) * t584 + (Icges(3,4) * t588 + Icges(3,2) * t592) * t581;
t519 = Icges(3,3) * t584 + (Icges(3,5) * t588 + Icges(3,6) * t592) * t581;
t516 = -t541 * t579 + t574 * t582;
t515 = -t542 * t579 - t573 * t582;
t514 = t541 * t582 + t574 * t579;
t513 = t542 * t582 - t573 * t579;
t511 = -t526 * t579 + t565 * t582;
t510 = -t527 * t579 - t564 * t582;
t509 = t526 * t582 + t565 * t579;
t508 = t527 * t582 - t629;
t507 = t580 * t620 + (t564 * t589 + t565 * t624) * t581;
t506 = t580 * t618 + (-t564 * t585 + t565 * t623) * t581;
t505 = rSges(4,1) * t543 + rSges(4,2) * t540 + rSges(4,3) * t584;
t504 = qJD(5) * t523 + t535;
t503 = Icges(4,1) * t543 + Icges(4,4) * t540 + Icges(4,5) * t584;
t502 = Icges(4,4) * t543 + Icges(4,2) * t540 + Icges(4,6) * t584;
t501 = Icges(4,5) * t543 + Icges(4,6) * t540 + Icges(4,3) * t584;
t500 = -t583 * t625 + (-t565 * t621 + t629) * t580;
t499 = t564 * t626 + (t565 * t580 * t584 + t581 * t583) * t579;
t496 = rSges(5,1) * t528 + rSges(5,2) * t525 + rSges(5,3) * t584;
t494 = Icges(5,1) * t528 + Icges(5,4) * t525 + Icges(5,5) * t584;
t493 = Icges(5,4) * t528 + Icges(5,2) * t525 + Icges(5,6) * t584;
t492 = Icges(5,5) * t528 + Icges(5,6) * t525 + Icges(5,3) * t584;
t491 = rSges(3,1) * t534 + rSges(3,2) * t533 + rSges(3,3) * t628;
t490 = rSges(3,1) * t532 + rSges(3,2) * t531 - rSges(3,3) * t625;
t488 = Icges(3,1) * t534 + Icges(3,4) * t533 + Icges(3,5) * t628;
t487 = Icges(3,1) * t532 + Icges(3,4) * t531 - Icges(3,5) * t625;
t486 = Icges(3,4) * t534 + Icges(3,2) * t533 + Icges(3,6) * t628;
t485 = Icges(3,4) * t532 + Icges(3,2) * t531 - Icges(3,6) * t625;
t484 = Icges(3,5) * t534 + Icges(3,6) * t533 + Icges(3,3) * t628;
t483 = Icges(3,5) * t532 + Icges(3,6) * t531 - Icges(3,3) * t625;
t482 = pkin(3) * t621 + t537 * t586 + t539 * t590;
t481 = -pkin(3) * t622 - t536 * t586 - t538 * t590;
t477 = (-t579 * t585 + t582 * t604) * t565 + (-t579 * t623 - t582 * t620) * t564 - t582 * t606;
t476 = (t579 * t589 + t582 * t605) * t565 + (-t579 * t624 + t582 * t618) * t564 - t582 * t607;
t475 = (-t579 * t604 - t582 * t585) * t565 + (t579 * t620 - t582 * t623) * t564 + t579 * t606;
t474 = (-t579 * t605 + t582 * t589) * t565 + (-t579 * t618 - t582 * t624) * t564 + t579 * t607;
t473 = qJD(5) * t500 + t518;
t472 = qJD(5) * t499 + t517;
t471 = rSges(4,1) * t516 + rSges(4,2) * t515 + rSges(4,3) * t628;
t470 = rSges(4,1) * t514 + rSges(4,2) * t513 - rSges(4,3) * t625;
t469 = Icges(4,1) * t516 + Icges(4,4) * t515 + Icges(4,5) * t628;
t468 = Icges(4,1) * t514 + Icges(4,4) * t513 - Icges(4,5) * t625;
t467 = Icges(4,4) * t516 + Icges(4,2) * t515 + Icges(4,6) * t628;
t466 = Icges(4,4) * t514 + Icges(4,2) * t513 - Icges(4,6) * t625;
t465 = Icges(4,5) * t516 + Icges(4,6) * t515 + Icges(4,3) * t628;
t464 = Icges(4,5) * t514 + Icges(4,6) * t513 - Icges(4,3) * t625;
t463 = rSges(5,1) * t511 + rSges(5,2) * t510 + rSges(5,3) * t628;
t462 = rSges(5,1) * t509 + rSges(5,2) * t508 - rSges(5,3) * t625;
t461 = Icges(5,1) * t511 + Icges(5,4) * t510 + Icges(5,5) * t628;
t460 = Icges(5,1) * t509 + Icges(5,4) * t508 - Icges(5,5) * t625;
t459 = Icges(5,4) * t511 + Icges(5,2) * t510 + Icges(5,6) * t628;
t458 = Icges(5,4) * t509 + Icges(5,2) * t508 - Icges(5,6) * t625;
t457 = Icges(5,5) * t511 + Icges(5,6) * t510 + Icges(5,3) * t628;
t456 = Icges(5,5) * t509 + Icges(5,6) * t508 - Icges(5,3) * t625;
t453 = (-t490 * t584 - t522 * t625) * qJD(2);
t452 = (t491 * t584 - t522 * t628) * qJD(2);
t449 = rSges(6,1) * t507 + rSges(6,2) * t506 + rSges(6,3) * t523;
t448 = Icges(6,1) * t507 + Icges(6,4) * t506 + Icges(6,5) * t523;
t447 = Icges(6,4) * t507 + Icges(6,2) * t506 + Icges(6,6) * t523;
t446 = Icges(6,5) * t507 + Icges(6,6) * t506 + Icges(6,3) * t523;
t445 = qJD(1) + (t490 * t579 + t491 * t582) * t614;
t444 = (-t577 + t554) * t584 + ((-t546 * t591 + (-pkin(3) + t545) * t587) * t592 + (t545 * t591 + t546 * t587 - t630) * t588) * t581;
t443 = rSges(6,1) * t476 + rSges(6,2) * t477 + rSges(6,3) * t500;
t442 = rSges(6,1) * t474 + rSges(6,2) * t475 + rSges(6,3) * t499;
t441 = Icges(6,1) * t476 + Icges(6,4) * t477 + Icges(6,5) * t500;
t440 = Icges(6,1) * t474 + Icges(6,4) * t475 + Icges(6,5) * t499;
t439 = Icges(6,4) * t476 + Icges(6,2) * t477 + Icges(6,6) * t500;
t438 = Icges(6,4) * t474 + Icges(6,2) * t475 + Icges(6,6) * t499;
t437 = Icges(6,5) * t476 + Icges(6,6) * t477 + Icges(6,3) * t500;
t436 = Icges(6,5) * t474 + Icges(6,6) * t475 + Icges(6,3) * t499;
t435 = -t470 * t548 + t505 * t530 + t596;
t434 = t471 * t548 - t505 * t529 + t601;
t433 = t470 * t529 - t471 * t530 + t603;
t432 = (t481 * t587 + t591 * t598) * t592 + (t481 * t591 - t587 * t598) * t588 - t599 * t579 - t635;
t431 = (t482 * t587 + t591 * t597) * t592 + (t482 * t591 - t587 * t597) * t588 + t599 * t582 - t636;
t430 = -t462 * t535 + t496 * t518 + t594;
t429 = t463 * t535 - t496 * t517 + t595;
t428 = t462 * t517 - t463 * t518 + t600;
t427 = -t431 * t535 - t443 * t504 + t444 * t518 + t449 * t473 + t594;
t426 = t432 * t535 + t442 * t504 - t444 * t517 - t449 * t472 + t595;
t425 = t431 * t517 - t432 * t518 - t442 * t473 + t443 * t472 + t600;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t445 ^ 2 + t452 ^ 2 + t453 ^ 2) / 0.2e1 - t633 * ((-t484 * t625 + t486 * t531 + t488 * t532) * t628 - (-t483 * t625 + t485 * t531 + t487 * t532) * t625 + (-t519 * t625 + t520 * t531 + t521 * t532) * t584) * t625 / 0.2e1 + m(4) * (t433 ^ 2 + t434 ^ 2 + t435 ^ 2) / 0.2e1 + t529 * ((t465 * t628 + t515 * t467 + t516 * t469) * t529 + (t464 * t628 + t466 * t515 + t468 * t516) * t530 + (t501 * t628 + t502 * t515 + t503 * t516) * t548) / 0.2e1 + t530 * ((-t465 * t625 + t467 * t513 + t469 * t514) * t529 + (-t464 * t625 + t513 * t466 + t514 * t468) * t530 + (-t501 * t625 + t502 * t513 + t503 * t514) * t548) / 0.2e1 + t548 * ((t465 * t584 + t467 * t540 + t469 * t543) * t529 + (t464 * t584 + t466 * t540 + t468 * t543) * t530 + (t584 * t501 + t540 * t502 + t543 * t503) * t548) / 0.2e1 + m(5) * (t428 ^ 2 + t429 ^ 2 + t430 ^ 2) / 0.2e1 + t517 * ((t457 * t628 + t510 * t459 + t511 * t461) * t517 + (t456 * t628 + t458 * t510 + t460 * t511) * t518 + (t492 * t628 + t493 * t510 + t494 * t511) * t535) / 0.2e1 + t518 * ((-t457 * t625 + t459 * t508 + t461 * t509) * t517 + (-t456 * t625 + t508 * t458 + t509 * t460) * t518 + (-t492 * t625 + t493 * t508 + t494 * t509) * t535) / 0.2e1 + t535 * ((t457 * t584 + t459 * t525 + t461 * t528) * t517 + (t456 * t584 + t458 * t525 + t460 * t528) * t518 + (t584 * t492 + t525 * t493 + t528 * t494) * t535) / 0.2e1 + m(6) * (t425 ^ 2 + t426 ^ 2 + t427 ^ 2) / 0.2e1 + t472 * ((t499 * t436 + t475 * t438 + t474 * t440) * t472 + (t437 * t499 + t439 * t475 + t441 * t474) * t473 + (t446 * t499 + t447 * t475 + t448 * t474) * t504) / 0.2e1 + t473 * ((t436 * t500 + t438 * t477 + t440 * t476) * t472 + (t500 * t437 + t477 * t439 + t476 * t441) * t473 + (t446 * t500 + t447 * t477 + t448 * t476) * t504) / 0.2e1 + t504 * ((t436 * t523 + t438 * t506 + t440 * t507) * t472 + (t437 * t523 + t439 * t506 + t441 * t507) * t473 + (t523 * t446 + t506 * t447 + t507 * t448) * t504) / 0.2e1 + (((t484 * t628 + t486 * t533 + t488 * t534) * t628 - (t483 * t628 + t485 * t533 + t487 * t534) * t625 + (t519 * t628 + t520 * t533 + t521 * t534) * t584) * t628 + t584 * (t584 ^ 2 * t519 + (((t486 * t592 + t488 * t588) * t579 - (t485 * t592 + t487 * t588) * t582) * t581 + (-t483 * t582 + t484 * t579 + t520 * t592 + t521 * t588) * t584) * t581)) * t633 / 0.2e1;
T = t1;
