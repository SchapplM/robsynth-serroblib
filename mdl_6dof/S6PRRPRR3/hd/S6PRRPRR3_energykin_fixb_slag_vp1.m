% Calculate kinetic energy for
% S6PRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_energykin_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRR3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:02:13
% EndTime: 2019-03-08 22:02:17
% DurationCPUTime: 3.76s
% Computational Cost: add. (5726->392), mult. (15437->602), div. (0->0), fcn. (19987->16), ass. (0->179)
t644 = Icges(4,3) + Icges(5,3);
t592 = sin(pkin(12));
t595 = cos(pkin(12));
t601 = sin(qJ(2));
t597 = cos(pkin(6));
t604 = cos(qJ(2));
t622 = t597 * t604;
t585 = -t592 * t622 - t595 * t601;
t593 = sin(pkin(7));
t596 = cos(pkin(7));
t594 = sin(pkin(6));
t629 = t592 * t594;
t566 = -t585 * t593 + t596 * t629;
t643 = pkin(9) * t566;
t583 = -t592 * t601 + t595 * t622;
t625 = t595 * t594;
t615 = t583 * t593 + t596 * t625;
t642 = pkin(9) * t615;
t600 = sin(qJ(3));
t603 = cos(qJ(3));
t630 = sin(pkin(13));
t631 = cos(pkin(13));
t588 = -t600 * t630 + t603 * t631;
t576 = t588 * t596;
t623 = t597 * t601;
t584 = t592 * t604 + t595 * t623;
t587 = -t600 * t631 - t603 * t630;
t610 = t593 * t588;
t607 = t594 * t610;
t522 = t576 * t583 + t584 * t587 - t595 * t607;
t575 = t587 * t593;
t577 = t587 * t596;
t523 = t575 * t625 - t577 * t583 + t584 * t588;
t614 = t583 * t596 - t593 * t625;
t540 = -t584 * t600 + t603 * t614;
t541 = t584 * t603 + t600 * t614;
t641 = Icges(4,5) * t541 + Icges(5,5) * t523 + Icges(4,6) * t540 + Icges(5,6) * t522 - t644 * t615;
t586 = -t592 * t623 + t595 * t604;
t524 = t585 * t576 + t586 * t587 + t592 * t607;
t525 = -t575 * t629 - t577 * t585 + t586 * t588;
t613 = t585 * t596 + t593 * t629;
t542 = -t586 * t600 + t603 * t613;
t543 = t586 * t603 + t600 * t613;
t640 = Icges(4,5) * t543 + Icges(5,5) * t525 + Icges(4,6) * t542 + Icges(5,6) * t524 + t644 * t566;
t626 = t594 * t604;
t627 = t594 * t601;
t536 = t576 * t626 + t587 * t627 + t597 * t610;
t537 = -t597 * t575 + (-t577 * t604 + t588 * t601) * t594;
t624 = t596 * t604;
t563 = t597 * t593 * t603 + (-t600 * t601 + t603 * t624) * t594;
t628 = t593 * t600;
t564 = t597 * t628 + (t600 * t624 + t601 * t603) * t594;
t582 = -t593 * t626 + t597 * t596;
t639 = Icges(4,5) * t564 + Icges(5,5) * t537 + Icges(4,6) * t563 + Icges(5,6) * t536 + t644 * t582;
t638 = qJD(2) ^ 2;
t634 = cos(qJ(5));
t633 = pkin(3) * t603;
t632 = pkin(9) + qJ(4);
t621 = qJD(2) * t594;
t589 = t592 * t621;
t555 = qJD(3) * t566 + t589;
t591 = qJD(2) * t597;
t569 = qJD(3) * t582 + t591;
t509 = -qJD(5) * t524 + t555;
t519 = -qJD(5) * t536 + t569;
t619 = t595 * t621;
t545 = pkin(2) * t584 - t642;
t546 = pkin(2) * t586 + t643;
t618 = t545 * t589 + t546 * t619 + qJD(1);
t556 = -qJD(3) * t615 - t619;
t580 = pkin(3) * t628 + t596 * t632;
t581 = pkin(3) * t596 * t600 - t593 * t632;
t511 = -t580 * t625 + t581 * t583 + t584 * t633 + t642;
t617 = qJD(4) * t582 + t555 * t511 + t618;
t568 = pkin(2) * t627 + pkin(9) * t582;
t616 = t546 * t591 - t568 * t589;
t510 = -qJD(5) * t522 + t556;
t512 = t580 * t629 + t581 * t585 + t586 * t633 - t643;
t612 = -qJD(4) * t615 + t569 * t512 + t616;
t611 = (-t545 * t597 - t568 * t625) * qJD(2);
t488 = pkin(4) * t523 - pkin(10) * t522;
t489 = pkin(4) * t525 - pkin(10) * t524;
t609 = t555 * t488 + (-t489 - t512) * t556 + t617;
t532 = (-pkin(9) * t596 + t580) * t597 + ((pkin(9) * t593 + t581) * t604 + t633 * t601) * t594;
t608 = qJD(4) * t566 + t556 * t532 + t611;
t508 = pkin(4) * t537 - pkin(10) * t536;
t606 = t569 * t489 + (-t508 - t532) * t555 + t612;
t605 = t556 * t508 + (-t488 - t511) * t569 + t608;
t602 = cos(qJ(6));
t599 = sin(qJ(5));
t598 = sin(qJ(6));
t573 = t597 * rSges(3,3) + (rSges(3,1) * t601 + rSges(3,2) * t604) * t594;
t572 = Icges(3,5) * t597 + (Icges(3,1) * t601 + Icges(3,4) * t604) * t594;
t571 = Icges(3,6) * t597 + (Icges(3,4) * t601 + Icges(3,2) * t604) * t594;
t570 = Icges(3,3) * t597 + (Icges(3,5) * t601 + Icges(3,6) * t604) * t594;
t554 = rSges(3,1) * t586 + rSges(3,2) * t585 + rSges(3,3) * t629;
t553 = rSges(3,1) * t584 + rSges(3,2) * t583 - rSges(3,3) * t625;
t552 = Icges(3,1) * t586 + Icges(3,4) * t585 + Icges(3,5) * t629;
t551 = Icges(3,1) * t584 + Icges(3,4) * t583 - Icges(3,5) * t625;
t550 = Icges(3,4) * t586 + Icges(3,2) * t585 + Icges(3,6) * t629;
t549 = Icges(3,4) * t584 + Icges(3,2) * t583 - Icges(3,6) * t625;
t548 = Icges(3,5) * t586 + Icges(3,6) * t585 + Icges(3,3) * t629;
t547 = Icges(3,5) * t584 + Icges(3,6) * t583 - Icges(3,3) * t625;
t534 = (-t553 * t597 - t573 * t625) * qJD(2);
t533 = (t554 * t597 - t573 * t629) * qJD(2);
t531 = rSges(4,1) * t564 + rSges(4,2) * t563 + rSges(4,3) * t582;
t530 = Icges(4,1) * t564 + Icges(4,4) * t563 + Icges(4,5) * t582;
t529 = Icges(4,4) * t564 + Icges(4,2) * t563 + Icges(4,6) * t582;
t527 = t537 * t634 + t582 * t599;
t526 = t537 * t599 - t582 * t634;
t518 = qJD(1) + (t553 * t592 + t554 * t595) * t621;
t516 = t525 * t634 + t566 * t599;
t515 = t525 * t599 - t566 * t634;
t514 = t523 * t634 - t599 * t615;
t513 = t523 * t599 + t615 * t634;
t507 = rSges(4,1) * t543 + rSges(4,2) * t542 + rSges(4,3) * t566;
t506 = rSges(4,1) * t541 + rSges(4,2) * t540 - rSges(4,3) * t615;
t505 = Icges(4,1) * t543 + Icges(4,4) * t542 + Icges(4,5) * t566;
t504 = Icges(4,1) * t541 + Icges(4,4) * t540 - Icges(4,5) * t615;
t503 = Icges(4,4) * t543 + Icges(4,2) * t542 + Icges(4,6) * t566;
t502 = Icges(4,4) * t541 + Icges(4,2) * t540 - Icges(4,6) * t615;
t498 = rSges(5,1) * t537 + rSges(5,2) * t536 + rSges(5,3) * t582;
t497 = Icges(5,1) * t537 + Icges(5,4) * t536 + Icges(5,5) * t582;
t496 = Icges(5,4) * t537 + Icges(5,2) * t536 + Icges(5,6) * t582;
t494 = t527 * t602 - t536 * t598;
t493 = -t527 * t598 - t536 * t602;
t490 = pkin(5) * t527 + pkin(11) * t526;
t487 = qJD(6) * t526 + t519;
t485 = rSges(5,1) * t525 + rSges(5,2) * t524 + rSges(5,3) * t566;
t484 = rSges(5,1) * t523 + rSges(5,2) * t522 - rSges(5,3) * t615;
t483 = Icges(5,1) * t525 + Icges(5,4) * t524 + Icges(5,5) * t566;
t482 = Icges(5,1) * t523 + Icges(5,4) * t522 - Icges(5,5) * t615;
t481 = Icges(5,4) * t525 + Icges(5,2) * t524 + Icges(5,6) * t566;
t480 = Icges(5,4) * t523 + Icges(5,2) * t522 - Icges(5,6) * t615;
t477 = t516 * t602 - t524 * t598;
t476 = -t516 * t598 - t524 * t602;
t475 = t514 * t602 - t522 * t598;
t474 = -t514 * t598 - t522 * t602;
t472 = rSges(6,1) * t527 - rSges(6,2) * t526 - rSges(6,3) * t536;
t471 = Icges(6,1) * t527 - Icges(6,4) * t526 - Icges(6,5) * t536;
t470 = Icges(6,4) * t527 - Icges(6,2) * t526 - Icges(6,6) * t536;
t469 = Icges(6,5) * t527 - Icges(6,6) * t526 - Icges(6,3) * t536;
t468 = pkin(5) * t516 + pkin(11) * t515;
t467 = pkin(5) * t514 + pkin(11) * t513;
t466 = qJD(6) * t513 + t510;
t465 = qJD(6) * t515 + t509;
t464 = rSges(6,1) * t516 - rSges(6,2) * t515 - rSges(6,3) * t524;
t463 = rSges(6,1) * t514 - rSges(6,2) * t513 - rSges(6,3) * t522;
t462 = Icges(6,1) * t516 - Icges(6,4) * t515 - Icges(6,5) * t524;
t461 = Icges(6,1) * t514 - Icges(6,4) * t513 - Icges(6,5) * t522;
t460 = Icges(6,4) * t516 - Icges(6,2) * t515 - Icges(6,6) * t524;
t459 = Icges(6,4) * t514 - Icges(6,2) * t513 - Icges(6,6) * t522;
t458 = Icges(6,5) * t516 - Icges(6,6) * t515 - Icges(6,3) * t524;
t457 = Icges(6,5) * t514 - Icges(6,6) * t513 - Icges(6,3) * t522;
t456 = -t506 * t569 + t531 * t556 + t611;
t455 = t507 * t569 - t531 * t555 + t616;
t454 = rSges(7,1) * t494 + rSges(7,2) * t493 + rSges(7,3) * t526;
t453 = Icges(7,1) * t494 + Icges(7,4) * t493 + Icges(7,5) * t526;
t452 = Icges(7,4) * t494 + Icges(7,2) * t493 + Icges(7,6) * t526;
t451 = Icges(7,5) * t494 + Icges(7,6) * t493 + Icges(7,3) * t526;
t450 = t506 * t555 - t507 * t556 + t618;
t449 = rSges(7,1) * t477 + rSges(7,2) * t476 + rSges(7,3) * t515;
t448 = rSges(7,1) * t475 + rSges(7,2) * t474 + rSges(7,3) * t513;
t447 = Icges(7,1) * t477 + Icges(7,4) * t476 + Icges(7,5) * t515;
t446 = Icges(7,1) * t475 + Icges(7,4) * t474 + Icges(7,5) * t513;
t445 = Icges(7,4) * t477 + Icges(7,2) * t476 + Icges(7,6) * t515;
t444 = Icges(7,4) * t475 + Icges(7,2) * t474 + Icges(7,6) * t513;
t443 = Icges(7,5) * t477 + Icges(7,6) * t476 + Icges(7,3) * t515;
t442 = Icges(7,5) * t475 + Icges(7,6) * t474 + Icges(7,3) * t513;
t441 = t498 * t556 + (-t484 - t511) * t569 + t608;
t440 = t485 * t569 + (-t498 - t532) * t555 + t612;
t439 = t484 * t555 + (-t485 - t512) * t556 + t617;
t438 = -t463 * t519 + t472 * t510 + t605;
t437 = t464 * t519 - t472 * t509 + t606;
t436 = t463 * t509 - t464 * t510 + t609;
t435 = -t448 * t487 + t454 * t466 - t467 * t519 + t490 * t510 + t605;
t434 = t449 * t487 - t454 * t465 + t468 * t519 - t490 * t509 + t606;
t433 = t448 * t465 - t449 * t466 + t467 * t509 - t468 * t510 + t609;
t1 = -t638 * ((-t548 * t625 + t550 * t583 + t552 * t584) * t629 - (-t547 * t625 + t549 * t583 + t551 * t584) * t625 + (-t570 * t625 + t571 * t583 + t572 * t584) * t597) * t625 / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t518 ^ 2 + t533 ^ 2 + t534 ^ 2) / 0.2e1 + m(4) * (t450 ^ 2 + t455 ^ 2 + t456 ^ 2) / 0.2e1 + m(5) * (t439 ^ 2 + t440 ^ 2 + t441 ^ 2) / 0.2e1 + m(6) * (t436 ^ 2 + t437 ^ 2 + t438 ^ 2) / 0.2e1 + t509 * ((-t524 * t458 - t515 * t460 + t516 * t462) * t509 + (-t457 * t524 - t459 * t515 + t461 * t516) * t510 + (-t469 * t524 - t470 * t515 + t471 * t516) * t519) / 0.2e1 + t510 * ((-t458 * t522 - t460 * t513 + t462 * t514) * t509 + (-t522 * t457 - t513 * t459 + t514 * t461) * t510 + (-t469 * t522 - t470 * t513 + t471 * t514) * t519) / 0.2e1 + t519 * ((-t458 * t536 - t460 * t526 + t462 * t527) * t509 + (-t457 * t536 - t459 * t526 + t461 * t527) * t510 + (-t536 * t469 - t526 * t470 + t527 * t471) * t519) / 0.2e1 + t465 * ((t515 * t443 + t476 * t445 + t477 * t447) * t465 + (t442 * t515 + t444 * t476 + t446 * t477) * t466 + (t451 * t515 + t452 * t476 + t453 * t477) * t487) / 0.2e1 + t466 * ((t443 * t513 + t445 * t474 + t447 * t475) * t465 + (t513 * t442 + t474 * t444 + t475 * t446) * t466 + (t451 * t513 + t452 * t474 + t453 * t475) * t487) / 0.2e1 + t487 * ((t443 * t526 + t445 * t493 + t447 * t494) * t465 + (t442 * t526 + t444 * t493 + t446 * t494) * t466 + (t526 * t451 + t493 * t452 + t494 * t453) * t487) / 0.2e1 + m(7) * (t433 ^ 2 + t434 ^ 2 + t435 ^ 2) / 0.2e1 + ((t496 * t524 + t497 * t525 + t529 * t542 + t530 * t543 + t566 * t639) * t569 + (t480 * t524 + t482 * t525 + t502 * t542 + t504 * t543 + t566 * t641) * t556 + (t524 * t481 + t525 * t483 + t542 * t503 + t543 * t505 + t566 * t640) * t555) * t555 / 0.2e1 + ((t496 * t522 + t497 * t523 + t529 * t540 + t530 * t541 - t615 * t639) * t569 + (t522 * t480 + t523 * t482 + t540 * t502 + t541 * t504 - t615 * t641) * t556 + (t481 * t522 + t483 * t523 + t503 * t540 + t505 * t541 - t615 * t640) * t555) * t556 / 0.2e1 + ((t536 * t496 + t537 * t497 + t563 * t529 + t564 * t530 + t582 * t639) * t569 + (t480 * t536 + t482 * t537 + t502 * t563 + t504 * t564 + t582 * t641) * t556 + (t481 * t536 + t483 * t537 + t503 * t563 + t505 * t564 + t582 * t640) * t555) * t569 / 0.2e1 + (((t548 * t629 + t550 * t585 + t552 * t586) * t629 - (t547 * t629 + t549 * t585 + t551 * t586) * t625 + (t570 * t629 + t571 * t585 + t572 * t586) * t597) * t629 + t597 * (t597 ^ 2 * t570 + (((t550 * t604 + t552 * t601) * t592 - (t549 * t604 + t551 * t601) * t595) * t594 + (-t547 * t595 + t548 * t592 + t571 * t604 + t572 * t601) * t597) * t594)) * t638 / 0.2e1;
T  = t1;
