% Calculate kinetic energy for
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 01:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_energykin_fixb_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRR6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRR6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:14:04
% EndTime: 2019-03-09 01:14:06
% DurationCPUTime: 2.75s
% Computational Cost: add. (9450->388), mult. (26599->621), div. (0->0), fcn. (34954->18), ass. (0->178)
t659 = qJD(2) ^ 2;
t658 = cos(qJ(4));
t657 = cos(qJ(5));
t656 = cos(pkin(8));
t655 = sin(pkin(8));
t617 = sin(pkin(14));
t619 = sin(pkin(6));
t654 = t617 * t619;
t618 = sin(pkin(7));
t653 = t618 * t619;
t622 = cos(pkin(6));
t652 = t618 * t622;
t620 = cos(pkin(14));
t651 = t619 * t620;
t621 = cos(pkin(7));
t650 = t619 * t621;
t630 = cos(qJ(2));
t649 = t621 * t630;
t627 = sin(qJ(2));
t648 = t622 * t627;
t647 = t622 * t630;
t613 = -t617 * t647 - t620 * t627;
t601 = -t613 * t618 + t617 * t650;
t646 = qJD(2) * t619;
t615 = t617 * t646;
t591 = qJD(3) * t601 + t615;
t610 = t622 * t621 - t630 * t653;
t616 = qJD(2) * t622;
t603 = qJD(3) * t610 + t616;
t614 = -t617 * t648 + t620 * t630;
t626 = sin(qJ(3));
t629 = cos(qJ(3));
t638 = t613 * t621 + t617 * t653;
t578 = -t614 * t626 + t629 * t638;
t560 = -t578 * t655 + t601 * t656;
t548 = qJD(4) * t560 + t591;
t598 = t629 * t652 + (-t626 * t627 + t629 * t649) * t619;
t575 = -t598 * t655 + t610 * t656;
t566 = qJD(4) * t575 + t603;
t644 = t620 * t646;
t611 = -t617 * t627 + t620 * t647;
t600 = -t611 * t618 - t620 * t650;
t612 = t617 * t630 + t620 * t648;
t581 = pkin(2) * t612 + pkin(10) * t600;
t582 = pkin(2) * t614 + pkin(10) * t601;
t643 = t581 * t615 + t582 * t644 + qJD(1);
t579 = t614 * t629 + t626 * t638;
t625 = sin(qJ(4));
t641 = t658 * t655;
t642 = t656 * t658;
t536 = -t578 * t642 + t579 * t625 - t601 * t641;
t514 = qJD(5) * t536 + t548;
t599 = t626 * t652 + (t626 * t649 + t627 * t629) * t619;
t557 = -t598 * t642 + t599 * t625 - t610 * t641;
t528 = qJD(5) * t557 + t566;
t592 = qJD(3) * t600 - t644;
t602 = t619 * t627 * pkin(2) + pkin(10) * t610;
t640 = t582 * t616 - t602 * t615;
t639 = t611 * t621 - t618 * t651;
t576 = -t612 * t626 + t629 * t639;
t559 = -t576 * t655 + t600 * t656;
t549 = qJD(4) * t559 + t592;
t577 = t612 * t629 + t626 * t639;
t534 = -t576 * t642 + t577 * t625 - t600 * t641;
t515 = qJD(5) * t534 + t549;
t538 = t577 * pkin(3) + pkin(11) * t559;
t539 = t579 * pkin(3) + pkin(11) * t560;
t637 = t591 * t538 - t539 * t592 + t643;
t636 = (-t581 * t622 - t602 * t651) * qJD(2);
t561 = t599 * pkin(3) + pkin(11) * t575;
t635 = t603 * t539 - t561 * t591 + t640;
t535 = t577 * t658 + (t576 * t656 + t600 * t655) * t625;
t512 = pkin(4) * t535 + pkin(12) * t534;
t537 = t579 * t658 + (t578 * t656 + t601 * t655) * t625;
t513 = pkin(4) * t537 + pkin(12) * t536;
t634 = t548 * t512 - t513 * t549 + t637;
t633 = -t538 * t603 + t592 * t561 + t636;
t558 = t599 * t658 + (t598 * t656 + t610 * t655) * t625;
t527 = pkin(4) * t558 + pkin(12) * t557;
t632 = t566 * t513 - t527 * t548 + t635;
t631 = -t512 * t566 + t549 * t527 + t633;
t628 = cos(qJ(6));
t624 = sin(qJ(5));
t623 = sin(qJ(6));
t608 = t622 * rSges(3,3) + (rSges(3,1) * t627 + rSges(3,2) * t630) * t619;
t607 = Icges(3,5) * t622 + (Icges(3,1) * t627 + Icges(3,4) * t630) * t619;
t606 = Icges(3,6) * t622 + (Icges(3,4) * t627 + Icges(3,2) * t630) * t619;
t605 = Icges(3,3) * t622 + (Icges(3,5) * t627 + Icges(3,6) * t630) * t619;
t590 = rSges(3,1) * t614 + rSges(3,2) * t613 + rSges(3,3) * t654;
t589 = rSges(3,1) * t612 + rSges(3,2) * t611 - rSges(3,3) * t651;
t588 = Icges(3,1) * t614 + Icges(3,4) * t613 + Icges(3,5) * t654;
t587 = Icges(3,1) * t612 + Icges(3,4) * t611 - Icges(3,5) * t651;
t586 = Icges(3,4) * t614 + Icges(3,2) * t613 + Icges(3,6) * t654;
t585 = Icges(3,4) * t612 + Icges(3,2) * t611 - Icges(3,6) * t651;
t584 = Icges(3,5) * t614 + Icges(3,6) * t613 + Icges(3,3) * t654;
t583 = Icges(3,5) * t612 + Icges(3,6) * t611 - Icges(3,3) * t651;
t568 = (-t589 * t622 - t608 * t651) * qJD(2);
t567 = (t590 * t622 - t608 * t654) * qJD(2);
t565 = rSges(4,1) * t599 + rSges(4,2) * t598 + rSges(4,3) * t610;
t564 = Icges(4,1) * t599 + Icges(4,4) * t598 + Icges(4,5) * t610;
t563 = Icges(4,4) * t599 + Icges(4,2) * t598 + Icges(4,6) * t610;
t562 = Icges(4,5) * t599 + Icges(4,6) * t598 + Icges(4,3) * t610;
t551 = qJD(1) + (t589 * t617 + t590 * t620) * t646;
t547 = rSges(4,1) * t579 + rSges(4,2) * t578 + rSges(4,3) * t601;
t546 = rSges(4,1) * t577 + rSges(4,2) * t576 + rSges(4,3) * t600;
t545 = Icges(4,1) * t579 + Icges(4,4) * t578 + Icges(4,5) * t601;
t544 = Icges(4,1) * t577 + Icges(4,4) * t576 + Icges(4,5) * t600;
t543 = Icges(4,4) * t579 + Icges(4,2) * t578 + Icges(4,6) * t601;
t542 = Icges(4,4) * t577 + Icges(4,2) * t576 + Icges(4,6) * t600;
t541 = Icges(4,5) * t579 + Icges(4,6) * t578 + Icges(4,3) * t601;
t540 = Icges(4,5) * t577 + Icges(4,6) * t576 + Icges(4,3) * t600;
t533 = t558 * t657 + t575 * t624;
t532 = t558 * t624 - t575 * t657;
t525 = rSges(5,1) * t558 - rSges(5,2) * t557 + rSges(5,3) * t575;
t524 = Icges(5,1) * t558 - Icges(5,4) * t557 + Icges(5,5) * t575;
t523 = Icges(5,4) * t558 - Icges(5,2) * t557 + Icges(5,6) * t575;
t522 = Icges(5,5) * t558 - Icges(5,6) * t557 + Icges(5,3) * t575;
t521 = t537 * t657 + t560 * t624;
t520 = t537 * t624 - t560 * t657;
t519 = t535 * t657 + t559 * t624;
t518 = t535 * t624 - t559 * t657;
t517 = t533 * t628 + t557 * t623;
t516 = -t533 * t623 + t557 * t628;
t511 = pkin(5) * t533 + pkin(13) * t532;
t509 = qJD(6) * t532 + t528;
t508 = -t546 * t603 + t565 * t592 + t636;
t507 = t547 * t603 - t565 * t591 + t640;
t505 = rSges(5,1) * t537 - rSges(5,2) * t536 + rSges(5,3) * t560;
t504 = rSges(5,1) * t535 - rSges(5,2) * t534 + rSges(5,3) * t559;
t503 = Icges(5,1) * t537 - Icges(5,4) * t536 + Icges(5,5) * t560;
t502 = Icges(5,1) * t535 - Icges(5,4) * t534 + Icges(5,5) * t559;
t501 = Icges(5,4) * t537 - Icges(5,2) * t536 + Icges(5,6) * t560;
t500 = Icges(5,4) * t535 - Icges(5,2) * t534 + Icges(5,6) * t559;
t499 = Icges(5,5) * t537 - Icges(5,6) * t536 + Icges(5,3) * t560;
t498 = Icges(5,5) * t535 - Icges(5,6) * t534 + Icges(5,3) * t559;
t497 = rSges(6,1) * t533 - rSges(6,2) * t532 + rSges(6,3) * t557;
t496 = Icges(6,1) * t533 - Icges(6,4) * t532 + Icges(6,5) * t557;
t495 = Icges(6,4) * t533 - Icges(6,2) * t532 + Icges(6,6) * t557;
t494 = Icges(6,5) * t533 - Icges(6,6) * t532 + Icges(6,3) * t557;
t493 = t521 * t628 + t536 * t623;
t492 = -t521 * t623 + t536 * t628;
t491 = t519 * t628 + t534 * t623;
t490 = -t519 * t623 + t534 * t628;
t488 = t546 * t591 - t547 * t592 + t643;
t487 = pkin(5) * t521 + pkin(13) * t520;
t486 = pkin(5) * t519 + pkin(13) * t518;
t485 = qJD(6) * t518 + t515;
t484 = qJD(6) * t520 + t514;
t483 = rSges(6,1) * t521 - rSges(6,2) * t520 + rSges(6,3) * t536;
t482 = rSges(6,1) * t519 - rSges(6,2) * t518 + rSges(6,3) * t534;
t481 = Icges(6,1) * t521 - Icges(6,4) * t520 + Icges(6,5) * t536;
t480 = Icges(6,1) * t519 - Icges(6,4) * t518 + Icges(6,5) * t534;
t479 = Icges(6,4) * t521 - Icges(6,2) * t520 + Icges(6,6) * t536;
t478 = Icges(6,4) * t519 - Icges(6,2) * t518 + Icges(6,6) * t534;
t477 = Icges(6,5) * t521 - Icges(6,6) * t520 + Icges(6,3) * t536;
t476 = Icges(6,5) * t519 - Icges(6,6) * t518 + Icges(6,3) * t534;
t475 = rSges(7,1) * t517 + rSges(7,2) * t516 + rSges(7,3) * t532;
t474 = Icges(7,1) * t517 + Icges(7,4) * t516 + Icges(7,5) * t532;
t473 = Icges(7,4) * t517 + Icges(7,2) * t516 + Icges(7,6) * t532;
t472 = Icges(7,5) * t517 + Icges(7,6) * t516 + Icges(7,3) * t532;
t471 = rSges(7,1) * t493 + rSges(7,2) * t492 + rSges(7,3) * t520;
t470 = rSges(7,1) * t491 + rSges(7,2) * t490 + rSges(7,3) * t518;
t469 = Icges(7,1) * t493 + Icges(7,4) * t492 + Icges(7,5) * t520;
t468 = Icges(7,1) * t491 + Icges(7,4) * t490 + Icges(7,5) * t518;
t467 = Icges(7,4) * t493 + Icges(7,2) * t492 + Icges(7,6) * t520;
t466 = Icges(7,4) * t491 + Icges(7,2) * t490 + Icges(7,6) * t518;
t465 = Icges(7,5) * t493 + Icges(7,6) * t492 + Icges(7,3) * t520;
t464 = Icges(7,5) * t491 + Icges(7,6) * t490 + Icges(7,3) * t518;
t463 = -t504 * t566 + t525 * t549 + t633;
t462 = t505 * t566 - t525 * t548 + t635;
t461 = t504 * t548 - t505 * t549 + t637;
t460 = -t482 * t528 + t497 * t515 + t631;
t459 = t483 * t528 - t497 * t514 + t632;
t458 = t482 * t514 - t483 * t515 + t634;
t457 = -t470 * t509 + t475 * t485 - t486 * t528 + t511 * t515 + t631;
t456 = t471 * t509 - t475 * t484 + t487 * t528 - t511 * t514 + t632;
t455 = t470 * t484 - t471 * t485 + t486 * t514 - t487 * t515 + t634;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + t509 * ((t465 * t532 + t467 * t516 + t469 * t517) * t484 + (t464 * t532 + t466 * t516 + t468 * t517) * t485 + (t472 * t532 + t473 * t516 + t474 * t517) * t509) / 0.2e1 + t484 * ((t520 * t465 + t492 * t467 + t493 * t469) * t484 + (t464 * t520 + t466 * t492 + t468 * t493) * t485 + (t472 * t520 + t473 * t492 + t474 * t493) * t509) / 0.2e1 + t485 * ((t465 * t518 + t467 * t490 + t469 * t491) * t484 + (t518 * t464 + t490 * t466 + t491 * t468) * t485 + (t472 * t518 + t473 * t490 + t474 * t491) * t509) / 0.2e1 + t514 * ((t477 * t536 - t479 * t520 + t481 * t521) * t514 + (t476 * t536 - t478 * t520 + t480 * t521) * t515 + (t494 * t536 - t495 * t520 + t496 * t521) * t528) / 0.2e1 + t515 * ((t477 * t534 - t479 * t518 + t481 * t519) * t514 + (t476 * t534 - t478 * t518 + t480 * t519) * t515 + (t494 * t534 - t495 * t518 + t496 * t519) * t528) / 0.2e1 + t528 * ((t477 * t557 - t479 * t532 + t481 * t533) * t514 + (t476 * t557 - t478 * t532 + t480 * t533) * t515 + (t494 * t557 - t495 * t532 + t496 * t533) * t528) / 0.2e1 + t548 * ((t499 * t560 - t501 * t536 + t503 * t537) * t548 + (t498 * t560 - t500 * t536 + t502 * t537) * t549 + (t522 * t560 - t523 * t536 + t524 * t537) * t566) / 0.2e1 + t566 * ((t499 * t575 - t501 * t557 + t503 * t558) * t548 + (t498 * t575 - t500 * t557 + t502 * t558) * t549 + (t522 * t575 - t523 * t557 + t524 * t558) * t566) / 0.2e1 + t549 * ((t499 * t559 - t501 * t534 + t503 * t535) * t548 + (t498 * t559 - t500 * t534 + t502 * t535) * t549 + (t522 * t559 - t523 * t534 + t524 * t535) * t566) / 0.2e1 + t603 * ((t541 * t610 + t543 * t598 + t545 * t599) * t591 + (t540 * t610 + t542 * t598 + t544 * t599) * t592 + (t562 * t610 + t563 * t598 + t564 * t599) * t603) / 0.2e1 + t592 * ((t541 * t600 + t543 * t576 + t545 * t577) * t591 + (t540 * t600 + t542 * t576 + t544 * t577) * t592 + (t562 * t600 + t563 * t576 + t564 * t577) * t603) / 0.2e1 + t591 * ((t541 * t601 + t543 * t578 + t545 * t579) * t591 + (t540 * t601 + t542 * t578 + t544 * t579) * t592 + (t562 * t601 + t563 * t578 + t564 * t579) * t603) / 0.2e1 + m(7) * (t455 ^ 2 + t456 ^ 2 + t457 ^ 2) / 0.2e1 + m(6) * (t458 ^ 2 + t459 ^ 2 + t460 ^ 2) / 0.2e1 + m(5) * (t461 ^ 2 + t462 ^ 2 + t463 ^ 2) / 0.2e1 + m(4) * (t488 ^ 2 + t507 ^ 2 + t508 ^ 2) / 0.2e1 + m(3) * (t551 ^ 2 + t567 ^ 2 + t568 ^ 2) / 0.2e1 - t659 * ((-t584 * t651 + t586 * t611 + t588 * t612) * t654 - (-t583 * t651 + t585 * t611 + t587 * t612) * t651 + (-t605 * t651 + t606 * t611 + t607 * t612) * t622) * t651 / 0.2e1 + (t622 * (t622 ^ 2 * t605 + (((t586 * t630 + t588 * t627) * t617 - (t585 * t630 + t587 * t627) * t620) * t619 + (-t583 * t620 + t584 * t617 + t606 * t630 + t607 * t627) * t622) * t619) + ((t584 * t654 + t586 * t613 + t588 * t614) * t654 - (t583 * t654 + t585 * t613 + t587 * t614) * t651 + (t605 * t654 + t606 * t613 + t607 * t614) * t622) * t654) * t659 / 0.2e1;
T  = t1;
