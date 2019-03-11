% Calculate kinetic energy for
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_energykin_fixb_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR7_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR7_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRR7_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:50:35
% EndTime: 2019-03-08 20:50:38
% DurationCPUTime: 3.40s
% Computational Cost: add. (9194->389), mult. (26023->616), div. (0->0), fcn. (34218->18), ass. (0->176)
t664 = qJD(2) ^ 2;
t663 = cos(qJ(4));
t662 = cos(qJ(5));
t661 = cos(pkin(8));
t660 = sin(pkin(8));
t618 = sin(pkin(13));
t620 = sin(pkin(6));
t659 = t618 * t620;
t619 = sin(pkin(7));
t658 = t619 * t620;
t624 = cos(pkin(6));
t657 = t619 * t624;
t622 = cos(pkin(13));
t656 = t620 * t622;
t623 = cos(pkin(7));
t655 = t620 * t623;
t628 = sin(qJ(2));
t654 = t620 * t628;
t630 = cos(qJ(2));
t653 = t623 * t630;
t652 = t624 * t628;
t651 = t624 * t630;
t613 = -t618 * t651 - t622 * t628;
t602 = -t613 * t619 + t618 * t655;
t614 = -t618 * t652 + t622 * t630;
t585 = pkin(2) * t614 + qJ(3) * t602;
t611 = -t618 * t628 + t622 * t651;
t601 = -t611 * t619 - t622 * t655;
t616 = qJD(2) * t624;
t650 = qJD(3) * t601 + t585 * t616;
t617 = sin(pkin(14));
t621 = cos(pkin(14));
t636 = t613 * t623 + t618 * t658;
t581 = -t614 * t617 + t621 * t636;
t563 = -t581 * t660 + t602 * t661;
t649 = qJD(2) * t620;
t615 = t618 * t649;
t553 = qJD(4) * t563 + t615;
t599 = t621 * t657 + (-t617 * t628 + t621 * t653) * t620;
t610 = t624 * t623 - t630 * t658;
t578 = -t599 * t660 + t610 * t661;
t571 = qJD(4) * t578 + t616;
t582 = t614 * t621 + t617 * t636;
t627 = sin(qJ(4));
t639 = t663 * t660;
t640 = t661 * t663;
t539 = -t581 * t640 + t582 * t627 - t602 * t639;
t517 = qJD(5) * t539 + t553;
t600 = t621 * t654 + (t620 * t653 + t657) * t617;
t558 = -t599 * t640 + t600 * t627 - t610 * t639;
t530 = qJD(5) * t558 + t571;
t647 = t622 * t649;
t603 = pkin(2) * t654 + qJ(3) * t610;
t645 = (-t600 * pkin(3) - pkin(10) * t578 - t603) * t620;
t644 = (-rSges(4,1) * t600 - rSges(4,2) * t599 - rSges(4,3) * t610 - t603) * t620;
t612 = t618 * t630 + t622 * t652;
t584 = pkin(2) * t612 + qJ(3) * t601;
t643 = qJD(3) * t610 + t584 * t615 + t585 * t647 + qJD(1);
t637 = t611 * t623 - t619 * t656;
t579 = -t612 * t617 + t621 * t637;
t562 = -t579 * t660 + t601 * t661;
t554 = qJD(4) * t562 - t647;
t580 = t612 * t621 + t617 * t637;
t537 = -t579 * t640 + t580 * t627 - t601 * t639;
t518 = qJD(5) * t537 + t554;
t542 = t580 * pkin(3) + pkin(10) * t562;
t543 = t582 * pkin(3) + pkin(10) * t563;
t638 = t542 * t615 + t543 * t647 + t643;
t635 = qJD(2) * t618 * t645 + t543 * t616 + t650;
t538 = t580 * t663 + (t579 * t661 + t601 * t660) * t627;
t514 = pkin(4) * t538 + pkin(11) * t537;
t540 = t582 * t663 + (t581 * t661 + t602 * t660) * t627;
t515 = pkin(4) * t540 + pkin(11) * t539;
t634 = t553 * t514 - t515 * t554 + t638;
t598 = qJD(3) * t602;
t633 = t598 + ((-t542 - t584) * t624 + t622 * t645) * qJD(2);
t559 = t600 * t663 + (t599 * t661 + t610 * t660) * t627;
t529 = pkin(4) * t559 + pkin(11) * t558;
t632 = t571 * t515 - t529 * t553 + t635;
t631 = -t514 * t571 + t554 * t529 + t633;
t629 = cos(qJ(6));
t626 = sin(qJ(5));
t625 = sin(qJ(6));
t608 = t624 * rSges(3,3) + (rSges(3,1) * t628 + rSges(3,2) * t630) * t620;
t607 = Icges(3,5) * t624 + (Icges(3,1) * t628 + Icges(3,4) * t630) * t620;
t606 = Icges(3,6) * t624 + (Icges(3,4) * t628 + Icges(3,2) * t630) * t620;
t605 = Icges(3,3) * t624 + (Icges(3,5) * t628 + Icges(3,6) * t630) * t620;
t593 = rSges(3,1) * t614 + rSges(3,2) * t613 + rSges(3,3) * t659;
t592 = rSges(3,1) * t612 + rSges(3,2) * t611 - rSges(3,3) * t656;
t591 = Icges(3,1) * t614 + Icges(3,4) * t613 + Icges(3,5) * t659;
t590 = Icges(3,1) * t612 + Icges(3,4) * t611 - Icges(3,5) * t656;
t589 = Icges(3,4) * t614 + Icges(3,2) * t613 + Icges(3,6) * t659;
t588 = Icges(3,4) * t612 + Icges(3,2) * t611 - Icges(3,6) * t656;
t587 = Icges(3,5) * t614 + Icges(3,6) * t613 + Icges(3,3) * t659;
t586 = Icges(3,5) * t612 + Icges(3,6) * t611 - Icges(3,3) * t656;
t570 = (-t592 * t624 - t608 * t656) * qJD(2);
t569 = (t593 * t624 - t608 * t659) * qJD(2);
t567 = Icges(4,1) * t600 + Icges(4,4) * t599 + Icges(4,5) * t610;
t566 = Icges(4,4) * t600 + Icges(4,2) * t599 + Icges(4,6) * t610;
t565 = Icges(4,5) * t600 + Icges(4,6) * t599 + Icges(4,3) * t610;
t552 = qJD(1) + (t592 * t618 + t593 * t622) * t649;
t551 = rSges(4,1) * t582 + rSges(4,2) * t581 + rSges(4,3) * t602;
t550 = rSges(4,1) * t580 + rSges(4,2) * t579 + rSges(4,3) * t601;
t549 = Icges(4,1) * t582 + Icges(4,4) * t581 + Icges(4,5) * t602;
t548 = Icges(4,1) * t580 + Icges(4,4) * t579 + Icges(4,5) * t601;
t547 = Icges(4,4) * t582 + Icges(4,2) * t581 + Icges(4,6) * t602;
t546 = Icges(4,4) * t580 + Icges(4,2) * t579 + Icges(4,6) * t601;
t545 = Icges(4,5) * t582 + Icges(4,6) * t581 + Icges(4,3) * t602;
t544 = Icges(4,5) * t580 + Icges(4,6) * t579 + Icges(4,3) * t601;
t536 = t559 * t662 + t578 * t626;
t535 = t559 * t626 - t578 * t662;
t528 = rSges(5,1) * t559 - rSges(5,2) * t558 + rSges(5,3) * t578;
t527 = Icges(5,1) * t559 - Icges(5,4) * t558 + Icges(5,5) * t578;
t526 = Icges(5,4) * t559 - Icges(5,2) * t558 + Icges(5,6) * t578;
t525 = Icges(5,5) * t559 - Icges(5,6) * t558 + Icges(5,3) * t578;
t524 = t540 * t662 + t563 * t626;
t523 = t540 * t626 - t563 * t662;
t522 = t538 * t662 + t562 * t626;
t521 = t538 * t626 - t562 * t662;
t520 = t536 * t629 + t558 * t625;
t519 = -t536 * t625 + t558 * t629;
t513 = pkin(5) * t536 + pkin(12) * t535;
t512 = qJD(6) * t535 + t530;
t511 = t598 + ((-t550 - t584) * t624 + t622 * t644) * qJD(2);
t510 = (t551 * t624 + t618 * t644) * qJD(2) + t650;
t508 = rSges(5,1) * t540 - rSges(5,2) * t539 + rSges(5,3) * t563;
t507 = rSges(5,1) * t538 - rSges(5,2) * t537 + rSges(5,3) * t562;
t506 = Icges(5,1) * t540 - Icges(5,4) * t539 + Icges(5,5) * t563;
t505 = Icges(5,1) * t538 - Icges(5,4) * t537 + Icges(5,5) * t562;
t504 = Icges(5,4) * t540 - Icges(5,2) * t539 + Icges(5,6) * t563;
t503 = Icges(5,4) * t538 - Icges(5,2) * t537 + Icges(5,6) * t562;
t502 = Icges(5,5) * t540 - Icges(5,6) * t539 + Icges(5,3) * t563;
t501 = Icges(5,5) * t538 - Icges(5,6) * t537 + Icges(5,3) * t562;
t500 = rSges(6,1) * t536 - rSges(6,2) * t535 + rSges(6,3) * t558;
t499 = Icges(6,1) * t536 - Icges(6,4) * t535 + Icges(6,5) * t558;
t498 = Icges(6,4) * t536 - Icges(6,2) * t535 + Icges(6,6) * t558;
t497 = Icges(6,5) * t536 - Icges(6,6) * t535 + Icges(6,3) * t558;
t496 = t524 * t629 + t539 * t625;
t495 = -t524 * t625 + t539 * t629;
t494 = t522 * t629 + t537 * t625;
t493 = -t522 * t625 + t537 * t629;
t491 = (t550 * t618 + t551 * t622) * t649 + t643;
t490 = pkin(5) * t524 + pkin(12) * t523;
t489 = pkin(5) * t522 + pkin(12) * t521;
t488 = qJD(6) * t521 + t518;
t487 = qJD(6) * t523 + t517;
t486 = rSges(6,1) * t524 - rSges(6,2) * t523 + rSges(6,3) * t539;
t485 = rSges(6,1) * t522 - rSges(6,2) * t521 + rSges(6,3) * t537;
t484 = Icges(6,1) * t524 - Icges(6,4) * t523 + Icges(6,5) * t539;
t483 = Icges(6,1) * t522 - Icges(6,4) * t521 + Icges(6,5) * t537;
t482 = Icges(6,4) * t524 - Icges(6,2) * t523 + Icges(6,6) * t539;
t481 = Icges(6,4) * t522 - Icges(6,2) * t521 + Icges(6,6) * t537;
t480 = Icges(6,5) * t524 - Icges(6,6) * t523 + Icges(6,3) * t539;
t479 = Icges(6,5) * t522 - Icges(6,6) * t521 + Icges(6,3) * t537;
t478 = rSges(7,1) * t520 + rSges(7,2) * t519 + rSges(7,3) * t535;
t477 = Icges(7,1) * t520 + Icges(7,4) * t519 + Icges(7,5) * t535;
t476 = Icges(7,4) * t520 + Icges(7,2) * t519 + Icges(7,6) * t535;
t475 = Icges(7,5) * t520 + Icges(7,6) * t519 + Icges(7,3) * t535;
t474 = rSges(7,1) * t496 + rSges(7,2) * t495 + rSges(7,3) * t523;
t473 = rSges(7,1) * t494 + rSges(7,2) * t493 + rSges(7,3) * t521;
t472 = Icges(7,1) * t496 + Icges(7,4) * t495 + Icges(7,5) * t523;
t471 = Icges(7,1) * t494 + Icges(7,4) * t493 + Icges(7,5) * t521;
t470 = Icges(7,4) * t496 + Icges(7,2) * t495 + Icges(7,6) * t523;
t469 = Icges(7,4) * t494 + Icges(7,2) * t493 + Icges(7,6) * t521;
t468 = Icges(7,5) * t496 + Icges(7,6) * t495 + Icges(7,3) * t523;
t467 = Icges(7,5) * t494 + Icges(7,6) * t493 + Icges(7,3) * t521;
t466 = -t507 * t571 + t528 * t554 + t633;
t465 = t508 * t571 - t528 * t553 + t635;
t464 = t507 * t553 - t508 * t554 + t638;
t463 = -t485 * t530 + t500 * t518 + t631;
t462 = t486 * t530 - t500 * t517 + t632;
t461 = t485 * t517 - t486 * t518 + t634;
t460 = -t473 * t512 + t478 * t488 - t489 * t530 + t513 * t518 + t631;
t459 = t474 * t512 - t478 * t487 + t490 * t530 - t513 * t517 + t632;
t458 = t473 * t487 - t474 * t488 + t489 * t517 - t490 * t518 + t634;
t1 = t517 * ((t480 * t539 - t482 * t523 + t484 * t524) * t517 + (t479 * t539 - t481 * t523 + t483 * t524) * t518 + (t497 * t539 - t498 * t523 + t499 * t524) * t530) / 0.2e1 + t518 * ((t480 * t537 - t482 * t521 + t484 * t522) * t517 + (t479 * t537 - t481 * t521 + t483 * t522) * t518 + (t497 * t537 - t498 * t521 + t499 * t522) * t530) / 0.2e1 + t530 * ((t480 * t558 - t482 * t535 + t484 * t536) * t517 + (t479 * t558 - t481 * t535 + t483 * t536) * t518 + (t497 * t558 - t498 * t535 + t499 * t536) * t530) / 0.2e1 + t571 * ((t502 * t578 - t504 * t558 + t506 * t559) * t553 + (t501 * t578 - t503 * t558 + t505 * t559) * t554 + (t525 * t578 - t526 * t558 + t527 * t559) * t571) / 0.2e1 + t553 * ((t502 * t563 - t504 * t539 + t506 * t540) * t553 + (t501 * t563 - t503 * t539 + t505 * t540) * t554 + (t525 * t563 - t526 * t539 + t527 * t540) * t571) / 0.2e1 + t554 * ((t502 * t562 - t504 * t537 + t506 * t538) * t553 + (t501 * t562 - t503 * t537 + t505 * t538) * t554 + (t525 * t562 - t526 * t537 + t527 * t538) * t571) / 0.2e1 + m(7) * (t458 ^ 2 + t459 ^ 2 + t460 ^ 2) / 0.2e1 + m(6) * (t461 ^ 2 + t462 ^ 2 + t463 ^ 2) / 0.2e1 + m(5) * (t464 ^ 2 + t465 ^ 2 + t466 ^ 2) / 0.2e1 + m(4) * (t491 ^ 2 + t510 ^ 2 + t511 ^ 2) / 0.2e1 + m(3) * (t552 ^ 2 + t569 ^ 2 + t570 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + t512 * ((t468 * t535 + t470 * t519 + t472 * t520) * t487 + (t467 * t535 + t469 * t519 + t471 * t520) * t488 + (t475 * t535 + t476 * t519 + t477 * t520) * t512) / 0.2e1 + t487 * ((t468 * t523 + t470 * t495 + t472 * t496) * t487 + (t467 * t523 + t469 * t495 + t471 * t496) * t488 + (t475 * t523 + t476 * t495 + t477 * t496) * t512) / 0.2e1 + t488 * ((t468 * t521 + t470 * t493 + t472 * t494) * t487 + (t467 * t521 + t469 * t493 + t471 * t494) * t488 + (t475 * t521 + t476 * t493 + t477 * t494) * t512) / 0.2e1 - ((-t587 * t656 + t589 * t611 + t591 * t612) * t659 - (-t586 * t656 + t588 * t611 + t590 * t612) * t656 + ((t545 * t601 + t547 * t579 + t549 * t580) * t618 - (t544 * t601 + t546 * t579 + t548 * t580) * t622) * t620 + (t565 * t601 + t566 * t579 + t567 * t580 - t605 * t656 + t606 * t611 + t607 * t612) * t624) * t664 * t656 / 0.2e1 + (((t565 * t610 + t566 * t599 + t567 * t600 + t624 * t605) * t624 + ((-t586 * t622 + t587 * t618 + t606 * t630 + t607 * t628) * t624 + (t545 * t610 + t547 * t599 + t549 * t600) * t618 - (t544 * t610 + t546 * t599 + t548 * t600) * t622 + ((t589 * t630 + t591 * t628) * t618 - (t588 * t630 + t590 * t628) * t622) * t620) * t620) * t624 + (((t545 * t602 + t547 * t581 + t549 * t582) * t618 - (t544 * t602 + t546 * t581 + t548 * t582) * t622) * t620 + (t587 * t659 + t589 * t613 + t591 * t614) * t659 - (t586 * t659 + t588 * t613 + t590 * t614) * t656 + (t565 * t602 + t566 * t581 + t567 * t582 + t605 * t659 + t606 * t613 + t607 * t614) * t624) * t659) * t664 / 0.2e1;
T  = t1;
