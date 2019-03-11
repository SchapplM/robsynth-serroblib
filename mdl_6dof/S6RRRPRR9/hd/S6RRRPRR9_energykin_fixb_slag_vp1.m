% Calculate kinetic energy for
% S6RRRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 19:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR9_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_energykin_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR9_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR9_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR9_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:56:32
% EndTime: 2019-03-09 18:56:36
% DurationCPUTime: 3.91s
% Computational Cost: add. (5795->399), mult. (15489->608), div. (0->0), fcn. (20021->16), ass. (0->184)
t654 = Icges(4,3) + Icges(5,3);
t600 = cos(pkin(6));
t605 = sin(qJ(1));
t608 = cos(qJ(2));
t630 = t605 * t608;
t604 = sin(qJ(2));
t609 = cos(qJ(1));
t631 = t604 * t609;
t585 = -t600 * t630 - t631;
t597 = sin(pkin(7));
t599 = cos(pkin(7));
t598 = sin(pkin(6));
t636 = t598 * t605;
t566 = -t585 * t597 + t599 * t636;
t653 = pkin(10) * t566;
t629 = t608 * t609;
t632 = t604 * t605;
t583 = t600 * t629 - t632;
t634 = t598 * t609;
t622 = t583 * t597 + t599 * t634;
t652 = pkin(10) * t622;
t584 = t600 * t631 + t630;
t586 = -t600 * t632 + t629;
t623 = (Icges(3,5) * t584 + Icges(3,6) * t583 - Icges(3,3) * t634) * t609 - (Icges(3,5) * t586 + Icges(3,6) * t585 + Icges(3,3) * t636) * t605;
t651 = t598 * t623;
t603 = sin(qJ(3));
t607 = cos(qJ(3));
t639 = sin(pkin(13));
t640 = cos(pkin(13));
t589 = -t603 * t639 + t607 * t640;
t576 = t589 * t599;
t588 = -t603 * t640 - t607 * t639;
t616 = t597 * t589;
t614 = t598 * t616;
t526 = t576 * t583 + t584 * t588 - t609 * t614;
t575 = t588 * t597;
t577 = t588 * t599;
t527 = t575 * t634 - t577 * t583 + t584 * t589;
t621 = t583 * t599 - t597 * t634;
t541 = -t584 * t603 + t607 * t621;
t542 = t584 * t607 + t603 * t621;
t650 = Icges(4,5) * t542 + Icges(5,5) * t527 + Icges(4,6) * t541 + Icges(5,6) * t526 - t654 * t622;
t528 = t585 * t576 + t586 * t588 + t605 * t614;
t529 = -t575 * t636 - t577 * t585 + t586 * t589;
t620 = t585 * t599 + t597 * t636;
t543 = -t586 * t603 + t607 * t620;
t544 = t586 * t607 + t603 * t620;
t649 = Icges(4,5) * t544 + Icges(5,5) * t529 + Icges(4,6) * t543 + Icges(5,6) * t528 + t654 * t566;
t635 = t598 * t608;
t637 = t598 * t604;
t536 = t576 * t635 + t588 * t637 + t600 * t616;
t537 = -t575 * t600 + (-t577 * t608 + t589 * t604) * t598;
t633 = t599 * t608;
t563 = t597 * t600 * t607 + (-t603 * t604 + t607 * t633) * t598;
t638 = t597 * t603;
t564 = t600 * t638 + (t603 * t633 + t604 * t607) * t598;
t582 = -t597 * t635 + t599 * t600;
t648 = Icges(4,5) * t564 + Icges(5,5) * t537 + Icges(4,6) * t563 + Icges(5,6) * t536 + t654 * t582;
t643 = cos(qJ(5));
t642 = pkin(3) * t607;
t641 = pkin(10) + qJ(4);
t545 = pkin(2) * t584 - t652;
t546 = pkin(2) * t586 + t653;
t626 = qJD(2) * t598;
t593 = t605 * t626;
t625 = t609 * t626;
t628 = t545 * t593 + t546 * t625;
t555 = qJD(3) * t566 + t593;
t627 = qJD(1) * (pkin(1) * t605 - pkin(9) * t634);
t594 = qJD(2) * t600 + qJD(1);
t509 = -qJD(5) * t528 + t555;
t567 = qJD(3) * t582 + t594;
t580 = pkin(3) * t638 + t599 * t641;
t581 = pkin(3) * t599 * t603 - t597 * t641;
t515 = -t580 * t634 + t581 * t583 + t584 * t642 + t652;
t624 = qJD(4) * t582 + t555 * t515 + t628;
t519 = -qJD(5) * t536 + t567;
t556 = -qJD(3) * t622 - t625;
t510 = -qJD(5) * t526 + t556;
t568 = pkin(2) * t637 + pkin(10) * t582;
t587 = qJD(1) * (pkin(1) * t609 + pkin(9) * t636);
t619 = t594 * t546 - t568 * t593 + t587;
t516 = t580 * t636 + t581 * t585 + t586 * t642 - t653;
t618 = -qJD(4) * t622 + t567 * t516 + t619;
t489 = pkin(4) * t527 - pkin(11) * t526;
t490 = pkin(4) * t529 - pkin(11) * t528;
t617 = t555 * t489 + (-t490 - t516) * t556 + t624;
t615 = -t545 * t594 - t568 * t625 - t627;
t534 = (-pkin(10) * t599 + t580) * t600 + ((pkin(10) * t597 + t581) * t608 + t642 * t604) * t598;
t613 = qJD(4) * t566 + t556 * t534 + t615;
t508 = pkin(4) * t537 - pkin(11) * t536;
t612 = t567 * t490 + (-t508 - t534) * t555 + t618;
t611 = t556 * t508 + (-t489 - t515) * t567 + t613;
t606 = cos(qJ(6));
t602 = sin(qJ(5));
t601 = sin(qJ(6));
t592 = rSges(2,1) * t609 - rSges(2,2) * t605;
t591 = rSges(2,1) * t605 + rSges(2,2) * t609;
t573 = rSges(3,3) * t600 + (rSges(3,1) * t604 + rSges(3,2) * t608) * t598;
t572 = Icges(3,5) * t600 + (Icges(3,1) * t604 + Icges(3,4) * t608) * t598;
t571 = Icges(3,6) * t600 + (Icges(3,4) * t604 + Icges(3,2) * t608) * t598;
t570 = Icges(3,3) * t600 + (Icges(3,5) * t604 + Icges(3,6) * t608) * t598;
t554 = rSges(3,1) * t586 + rSges(3,2) * t585 + rSges(3,3) * t636;
t553 = rSges(3,1) * t584 + rSges(3,2) * t583 - rSges(3,3) * t634;
t552 = Icges(3,1) * t586 + Icges(3,4) * t585 + Icges(3,5) * t636;
t551 = Icges(3,1) * t584 + Icges(3,4) * t583 - Icges(3,5) * t634;
t550 = Icges(3,4) * t586 + Icges(3,2) * t585 + Icges(3,6) * t636;
t549 = Icges(3,4) * t584 + Icges(3,2) * t583 - Icges(3,6) * t634;
t533 = rSges(4,1) * t564 + rSges(4,2) * t563 + rSges(4,3) * t582;
t532 = Icges(4,1) * t564 + Icges(4,4) * t563 + Icges(4,5) * t582;
t531 = Icges(4,4) * t564 + Icges(4,2) * t563 + Icges(4,6) * t582;
t525 = t537 * t643 + t582 * t602;
t524 = t537 * t602 - t582 * t643;
t523 = t554 * t594 - t573 * t593 + t587;
t522 = -t553 * t594 - t573 * t625 - t627;
t518 = (t553 * t605 + t554 * t609) * t626;
t514 = t529 * t643 + t566 * t602;
t513 = t529 * t602 - t566 * t643;
t512 = t527 * t643 - t602 * t622;
t511 = t527 * t602 + t622 * t643;
t507 = rSges(4,1) * t544 + rSges(4,2) * t543 + rSges(4,3) * t566;
t506 = rSges(4,1) * t542 + rSges(4,2) * t541 - rSges(4,3) * t622;
t505 = Icges(4,1) * t544 + Icges(4,4) * t543 + Icges(4,5) * t566;
t504 = Icges(4,1) * t542 + Icges(4,4) * t541 - Icges(4,5) * t622;
t503 = Icges(4,4) * t544 + Icges(4,2) * t543 + Icges(4,6) * t566;
t502 = Icges(4,4) * t542 + Icges(4,2) * t541 - Icges(4,6) * t622;
t498 = rSges(5,1) * t537 + rSges(5,2) * t536 + rSges(5,3) * t582;
t497 = Icges(5,1) * t537 + Icges(5,4) * t536 + Icges(5,5) * t582;
t496 = Icges(5,4) * t537 + Icges(5,2) * t536 + Icges(5,6) * t582;
t493 = t525 * t606 - t536 * t601;
t492 = -t525 * t601 - t536 * t606;
t488 = pkin(5) * t525 + pkin(12) * t524;
t487 = qJD(6) * t524 + t519;
t485 = rSges(5,1) * t529 + rSges(5,2) * t528 + rSges(5,3) * t566;
t484 = rSges(5,1) * t527 + rSges(5,2) * t526 - rSges(5,3) * t622;
t483 = Icges(5,1) * t529 + Icges(5,4) * t528 + Icges(5,5) * t566;
t482 = Icges(5,1) * t527 + Icges(5,4) * t526 - Icges(5,5) * t622;
t481 = Icges(5,4) * t529 + Icges(5,2) * t528 + Icges(5,6) * t566;
t480 = Icges(5,4) * t527 + Icges(5,2) * t526 - Icges(5,6) * t622;
t477 = t514 * t606 - t528 * t601;
t476 = -t514 * t601 - t528 * t606;
t475 = t512 * t606 - t526 * t601;
t474 = -t512 * t601 - t526 * t606;
t472 = pkin(5) * t514 + pkin(12) * t513;
t471 = pkin(5) * t512 + pkin(12) * t511;
t470 = rSges(6,1) * t525 - rSges(6,2) * t524 - rSges(6,3) * t536;
t469 = Icges(6,1) * t525 - Icges(6,4) * t524 - Icges(6,5) * t536;
t468 = Icges(6,4) * t525 - Icges(6,2) * t524 - Icges(6,6) * t536;
t467 = Icges(6,5) * t525 - Icges(6,6) * t524 - Icges(6,3) * t536;
t466 = qJD(6) * t511 + t510;
t465 = qJD(6) * t513 + t509;
t464 = rSges(6,1) * t514 - rSges(6,2) * t513 - rSges(6,3) * t528;
t463 = rSges(6,1) * t512 - rSges(6,2) * t511 - rSges(6,3) * t526;
t462 = Icges(6,1) * t514 - Icges(6,4) * t513 - Icges(6,5) * t528;
t461 = Icges(6,1) * t512 - Icges(6,4) * t511 - Icges(6,5) * t526;
t460 = Icges(6,4) * t514 - Icges(6,2) * t513 - Icges(6,6) * t528;
t459 = Icges(6,4) * t512 - Icges(6,2) * t511 - Icges(6,6) * t526;
t458 = Icges(6,5) * t514 - Icges(6,6) * t513 - Icges(6,3) * t528;
t457 = Icges(6,5) * t512 - Icges(6,6) * t511 - Icges(6,3) * t526;
t456 = t507 * t567 - t533 * t555 + t619;
t455 = -t506 * t567 + t533 * t556 + t615;
t454 = rSges(7,1) * t493 + rSges(7,2) * t492 + rSges(7,3) * t524;
t453 = Icges(7,1) * t493 + Icges(7,4) * t492 + Icges(7,5) * t524;
t452 = Icges(7,4) * t493 + Icges(7,2) * t492 + Icges(7,6) * t524;
t451 = Icges(7,5) * t493 + Icges(7,6) * t492 + Icges(7,3) * t524;
t450 = t506 * t555 - t507 * t556 + t628;
t449 = rSges(7,1) * t477 + rSges(7,2) * t476 + rSges(7,3) * t513;
t448 = rSges(7,1) * t475 + rSges(7,2) * t474 + rSges(7,3) * t511;
t447 = Icges(7,1) * t477 + Icges(7,4) * t476 + Icges(7,5) * t513;
t446 = Icges(7,1) * t475 + Icges(7,4) * t474 + Icges(7,5) * t511;
t445 = Icges(7,4) * t477 + Icges(7,2) * t476 + Icges(7,6) * t513;
t444 = Icges(7,4) * t475 + Icges(7,2) * t474 + Icges(7,6) * t511;
t443 = Icges(7,5) * t477 + Icges(7,6) * t476 + Icges(7,3) * t513;
t442 = Icges(7,5) * t475 + Icges(7,6) * t474 + Icges(7,3) * t511;
t441 = t485 * t567 + (-t498 - t534) * t555 + t618;
t440 = t498 * t556 + (-t484 - t515) * t567 + t613;
t439 = t484 * t555 + (-t485 - t516) * t556 + t624;
t438 = t464 * t519 - t470 * t509 + t612;
t437 = -t463 * t519 + t470 * t510 + t611;
t436 = t463 * t509 - t464 * t510 + t617;
t435 = t449 * t487 - t454 * t465 + t472 * t519 - t488 * t509 + t612;
t434 = -t448 * t487 + t454 * t466 - t471 * t519 + t488 * t510 + t611;
t433 = t448 * t465 - t449 * t466 + t471 * t509 - t472 * t510 + t617;
t1 = -((-t570 * t634 + t571 * t583 + t572 * t584) * t594 + ((t550 * t583 + t552 * t584) * t605 + (-t583 * t549 - t584 * t551 + t651) * t609) * t626) * t625 / 0.2e1 + ((t570 * t636 + t571 * t585 + t572 * t586) * t594 + (-(t549 * t585 + t551 * t586) * t609 + (t585 * t550 + t586 * t552 - t651) * t605) * t626) * t593 / 0.2e1 + m(6) * (t436 ^ 2 + t437 ^ 2 + t438 ^ 2) / 0.2e1 + m(7) * (t433 ^ 2 + t434 ^ 2 + t435 ^ 2) / 0.2e1 + t594 * ((t600 * t570 + (t571 * t608 + t572 * t604) * t598) * t594 + (((t550 * t608 + t552 * t604) * t605 - (t549 * t608 + t551 * t604) * t609) * t598 - t623 * t600) * t626) / 0.2e1 + t509 * ((-t528 * t458 - t513 * t460 + t514 * t462) * t509 + (-t457 * t528 - t459 * t513 + t461 * t514) * t510 + (-t467 * t528 - t468 * t513 + t469 * t514) * t519) / 0.2e1 + t510 * ((-t458 * t526 - t460 * t511 + t462 * t512) * t509 + (-t526 * t457 - t511 * t459 + t512 * t461) * t510 + (-t467 * t526 - t468 * t511 + t469 * t512) * t519) / 0.2e1 + t519 * ((-t458 * t536 - t460 * t524 + t462 * t525) * t509 + (-t457 * t536 - t459 * t524 + t461 * t525) * t510 + (-t536 * t467 - t524 * t468 + t525 * t469) * t519) / 0.2e1 + t465 * ((t513 * t443 + t476 * t445 + t477 * t447) * t465 + (t442 * t513 + t444 * t476 + t446 * t477) * t466 + (t451 * t513 + t452 * t476 + t453 * t477) * t487) / 0.2e1 + t466 * ((t443 * t511 + t445 * t474 + t447 * t475) * t465 + (t511 * t442 + t474 * t444 + t475 * t446) * t466 + (t451 * t511 + t452 * t474 + t453 * t475) * t487) / 0.2e1 + t487 * ((t443 * t524 + t445 * t492 + t447 * t493) * t465 + (t442 * t524 + t444 * t492 + t446 * t493) * t466 + (t524 * t451 + t492 * t452 + t493 * t453) * t487) / 0.2e1 + m(3) * (t518 ^ 2 + t522 ^ 2 + t523 ^ 2) / 0.2e1 + m(4) * (t450 ^ 2 + t455 ^ 2 + t456 ^ 2) / 0.2e1 + m(5) * (t439 ^ 2 + t440 ^ 2 + t441 ^ 2) / 0.2e1 + ((t496 * t528 + t497 * t529 + t531 * t543 + t532 * t544 + t566 * t648) * t567 + (t480 * t528 + t482 * t529 + t502 * t543 + t504 * t544 + t566 * t650) * t556 + (t528 * t481 + t529 * t483 + t543 * t503 + t544 * t505 + t649 * t566) * t555) * t555 / 0.2e1 + ((t496 * t526 + t497 * t527 + t531 * t541 + t532 * t542 - t622 * t648) * t567 + (t526 * t480 + t527 * t482 + t541 * t502 + t542 * t504 - t622 * t650) * t556 + (t481 * t526 + t483 * t527 + t503 * t541 + t505 * t542 - t622 * t649) * t555) * t556 / 0.2e1 + ((t536 * t496 + t537 * t497 + t563 * t531 + t564 * t532 + t648 * t582) * t567 + (t480 * t536 + t482 * t537 + t502 * t563 + t504 * t564 + t582 * t650) * t556 + (t481 * t536 + t483 * t537 + t503 * t563 + t505 * t564 + t582 * t649) * t555) * t567 / 0.2e1 + (Icges(2,3) + m(2) * (t591 ^ 2 + t592 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
