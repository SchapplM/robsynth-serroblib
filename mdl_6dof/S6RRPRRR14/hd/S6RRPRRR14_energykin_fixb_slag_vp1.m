% Calculate kinetic energy for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (leg links until cut joint, platform)
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2019-01-03 10:26
% Revision: 5fdbc45bcf2cc60deefd7ac2d71d743ed41bf7e4 (2018-12-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S6RRPRRR14_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_energykin_fixb_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR14_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR14_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-01-03 10:07:48
% EndTime: 2019-01-03 10:07:52
% DurationCPUTime: 3.45s
% Computational Cost: add. (9263->396), mult. (26075->612), div. (0->0), fcn. (34252->18), ass. (0->179)
t623 = sin(pkin(6));
t626 = cos(pkin(6));
t633 = cos(qJ(2));
t634 = cos(qJ(1));
t657 = t633 * t634;
t630 = sin(qJ(2));
t631 = sin(qJ(1));
t660 = t630 * t631;
t610 = t626 * t657 - t660;
t658 = t631 * t633;
t659 = t630 * t634;
t611 = t626 * t659 + t658;
t612 = -t626 * t658 - t659;
t613 = -t626 * t660 + t657;
t662 = t623 * t634;
t663 = t623 * t631;
t643 = (Icges(3,5) * t611 + Icges(3,6) * t610 - Icges(3,3) * t662) * t634 - (Icges(3,5) * t613 + Icges(3,6) * t612 + Icges(3,3) * t663) * t631;
t672 = t623 * t643;
t669 = cos(qJ(4));
t668 = cos(qJ(5));
t667 = cos(pkin(8));
t666 = sin(pkin(8));
t622 = sin(pkin(7));
t665 = t622 * t626;
t664 = t623 * t630;
t625 = cos(pkin(7));
t661 = t625 * t633;
t621 = sin(pkin(14));
t624 = cos(pkin(14));
t641 = t612 * t625 + t622 * t663;
t581 = -t613 * t621 + t624 * t641;
t601 = -t612 * t622 + t625 * t663;
t565 = -t581 * t666 + t601 * t667;
t655 = qJD(2) * t623;
t618 = t631 * t655;
t551 = qJD(4) * t565 + t618;
t656 = qJD(1) * (pkin(1) * t631 - pkin(10) * t662);
t619 = qJD(2) * t626 + qJD(1);
t584 = pkin(2) * t613 + qJ(3) * t601;
t600 = -t610 * t622 - t625 * t662;
t614 = qJD(1) * (pkin(1) * t634 + pkin(10) * t663);
t654 = qJD(3) * t600 + t619 * t584 + t614;
t583 = pkin(2) * t611 + qJ(3) * t600;
t609 = -t622 * t623 * t633 + t625 * t626;
t652 = t634 * t655;
t653 = qJD(3) * t609 + t583 * t618 + t584 * t652;
t582 = t613 * t624 + t621 * t641;
t629 = sin(qJ(4));
t647 = t669 * t666;
t648 = t667 * t669;
t539 = -t581 * t648 + t582 * t629 - t601 * t647;
t518 = qJD(5) * t539 + t551;
t596 = t624 * t665 + (-t621 * t630 + t624 * t661) * t623;
t576 = -t596 * t666 + t609 * t667;
t570 = qJD(4) * t576 + t619;
t651 = qJD(3) * t601 - t656;
t597 = t624 * t664 + (t623 * t661 + t665) * t621;
t555 = -t596 * t648 + t597 * t629 - t609 * t647;
t529 = qJD(5) * t555 + t570;
t602 = pkin(2) * t664 + qJ(3) * t609;
t646 = (-t597 * pkin(3) - pkin(11) * t576 - t602) * t655;
t645 = (-rSges(4,1) * t597 - rSges(4,2) * t596 - rSges(4,3) * t609 - t602) * t655;
t642 = t610 * t625 - t622 * t662;
t579 = -t611 * t621 + t624 * t642;
t564 = -t579 * t666 + t600 * t667;
t552 = qJD(4) * t564 - t652;
t580 = t611 * t624 + t621 * t642;
t541 = t580 * pkin(3) + pkin(11) * t564;
t542 = t582 * pkin(3) + pkin(11) * t565;
t644 = t541 * t618 + t542 * t652 + t653;
t537 = -t579 * t648 + t580 * t629 - t600 * t647;
t519 = qJD(5) * t537 + t552;
t538 = t580 * t669 + (t579 * t667 + t600 * t666) * t629;
t513 = pkin(4) * t538 + pkin(12) * t537;
t540 = t582 * t669 + (t581 * t667 + t601 * t666) * t629;
t514 = pkin(4) * t540 + pkin(12) * t539;
t640 = t551 * t513 - t514 * t552 + t644;
t639 = t619 * t542 + t631 * t646 + t654;
t556 = t597 * t669 + (t596 * t667 + t609 * t666) * t629;
t528 = pkin(4) * t556 + pkin(12) * t555;
t638 = t570 * t514 - t528 * t551 + t639;
t637 = (-t541 - t583) * t619 + t634 * t646 + t651;
t636 = -t513 * t570 + t552 * t528 + t637;
t632 = cos(qJ(6));
t628 = sin(qJ(5));
t627 = sin(qJ(6));
t617 = rSges(2,1) * t634 - rSges(2,2) * t631;
t616 = rSges(2,1) * t631 + rSges(2,2) * t634;
t607 = rSges(3,3) * t626 + (rSges(3,1) * t630 + rSges(3,2) * t633) * t623;
t606 = Icges(3,5) * t626 + (Icges(3,1) * t630 + Icges(3,4) * t633) * t623;
t605 = Icges(3,6) * t626 + (Icges(3,4) * t630 + Icges(3,2) * t633) * t623;
t604 = Icges(3,3) * t626 + (Icges(3,5) * t630 + Icges(3,6) * t633) * t623;
t592 = rSges(3,1) * t613 + rSges(3,2) * t612 + rSges(3,3) * t663;
t591 = rSges(3,1) * t611 + rSges(3,2) * t610 - rSges(3,3) * t662;
t590 = Icges(3,1) * t613 + Icges(3,4) * t612 + Icges(3,5) * t663;
t589 = Icges(3,1) * t611 + Icges(3,4) * t610 - Icges(3,5) * t662;
t588 = Icges(3,4) * t613 + Icges(3,2) * t612 + Icges(3,6) * t663;
t587 = Icges(3,4) * t611 + Icges(3,2) * t610 - Icges(3,6) * t662;
t568 = Icges(4,1) * t597 + Icges(4,4) * t596 + Icges(4,5) * t609;
t567 = Icges(4,4) * t597 + Icges(4,2) * t596 + Icges(4,6) * t609;
t566 = Icges(4,5) * t597 + Icges(4,6) * t596 + Icges(4,3) * t609;
t562 = t592 * t619 - t607 * t618 + t614;
t561 = -t591 * t619 - t607 * t652 - t656;
t553 = (t591 * t631 + t592 * t634) * t655;
t550 = rSges(4,1) * t582 + rSges(4,2) * t581 + rSges(4,3) * t601;
t549 = rSges(4,1) * t580 + rSges(4,2) * t579 + rSges(4,3) * t600;
t548 = Icges(4,1) * t582 + Icges(4,4) * t581 + Icges(4,5) * t601;
t547 = Icges(4,1) * t580 + Icges(4,4) * t579 + Icges(4,5) * t600;
t546 = Icges(4,4) * t582 + Icges(4,2) * t581 + Icges(4,6) * t601;
t545 = Icges(4,4) * t580 + Icges(4,2) * t579 + Icges(4,6) * t600;
t544 = Icges(4,5) * t582 + Icges(4,6) * t581 + Icges(4,3) * t601;
t543 = Icges(4,5) * t580 + Icges(4,6) * t579 + Icges(4,3) * t600;
t531 = t556 * t668 + t576 * t628;
t530 = t556 * t628 - t576 * t668;
t527 = t540 * t668 + t565 * t628;
t526 = t540 * t628 - t565 * t668;
t525 = t538 * t668 + t564 * t628;
t524 = t538 * t628 - t564 * t668;
t523 = rSges(5,1) * t556 - rSges(5,2) * t555 + rSges(5,3) * t576;
t522 = Icges(5,1) * t556 - Icges(5,4) * t555 + Icges(5,5) * t576;
t521 = Icges(5,4) * t556 - Icges(5,2) * t555 + Icges(5,6) * t576;
t520 = Icges(5,5) * t556 - Icges(5,6) * t555 + Icges(5,3) * t576;
t517 = t531 * t632 + t555 * t627;
t516 = -t531 * t627 + t555 * t632;
t512 = pkin(5) * t531 + pkin(13) * t530;
t511 = qJD(6) * t530 + t529;
t509 = t550 * t619 + t631 * t645 + t654;
t508 = (-t549 - t583) * t619 + t634 * t645 + t651;
t507 = rSges(5,1) * t540 - rSges(5,2) * t539 + rSges(5,3) * t565;
t506 = rSges(5,1) * t538 - rSges(5,2) * t537 + rSges(5,3) * t564;
t505 = Icges(5,1) * t540 - Icges(5,4) * t539 + Icges(5,5) * t565;
t504 = Icges(5,1) * t538 - Icges(5,4) * t537 + Icges(5,5) * t564;
t503 = Icges(5,4) * t540 - Icges(5,2) * t539 + Icges(5,6) * t565;
t502 = Icges(5,4) * t538 - Icges(5,2) * t537 + Icges(5,6) * t564;
t501 = Icges(5,5) * t540 - Icges(5,6) * t539 + Icges(5,3) * t565;
t500 = Icges(5,5) * t538 - Icges(5,6) * t537 + Icges(5,3) * t564;
t499 = t527 * t632 + t539 * t627;
t498 = -t527 * t627 + t539 * t632;
t497 = t525 * t632 + t537 * t627;
t496 = -t525 * t627 + t537 * t632;
t494 = rSges(6,1) * t531 - rSges(6,2) * t530 + rSges(6,3) * t555;
t493 = Icges(6,1) * t531 - Icges(6,4) * t530 + Icges(6,5) * t555;
t492 = Icges(6,4) * t531 - Icges(6,2) * t530 + Icges(6,6) * t555;
t491 = Icges(6,5) * t531 - Icges(6,6) * t530 + Icges(6,3) * t555;
t490 = (t549 * t631 + t550 * t634) * t655 + t653;
t489 = pkin(5) * t527 + pkin(13) * t526;
t488 = pkin(5) * t525 + pkin(13) * t524;
t487 = qJD(6) * t524 + t519;
t486 = qJD(6) * t526 + t518;
t485 = rSges(6,1) * t527 - rSges(6,2) * t526 + rSges(6,3) * t539;
t484 = rSges(6,1) * t525 - rSges(6,2) * t524 + rSges(6,3) * t537;
t483 = Icges(6,1) * t527 - Icges(6,4) * t526 + Icges(6,5) * t539;
t482 = Icges(6,1) * t525 - Icges(6,4) * t524 + Icges(6,5) * t537;
t481 = Icges(6,4) * t527 - Icges(6,2) * t526 + Icges(6,6) * t539;
t480 = Icges(6,4) * t525 - Icges(6,2) * t524 + Icges(6,6) * t537;
t479 = Icges(6,5) * t527 - Icges(6,6) * t526 + Icges(6,3) * t539;
t478 = Icges(6,5) * t525 - Icges(6,6) * t524 + Icges(6,3) * t537;
t477 = rSges(7,1) * t517 + rSges(7,2) * t516 + rSges(7,3) * t530;
t476 = Icges(7,1) * t517 + Icges(7,4) * t516 + Icges(7,5) * t530;
t475 = Icges(7,4) * t517 + Icges(7,2) * t516 + Icges(7,6) * t530;
t474 = Icges(7,5) * t517 + Icges(7,6) * t516 + Icges(7,3) * t530;
t473 = rSges(7,1) * t499 + rSges(7,2) * t498 + rSges(7,3) * t526;
t472 = rSges(7,1) * t497 + rSges(7,2) * t496 + rSges(7,3) * t524;
t471 = Icges(7,1) * t499 + Icges(7,4) * t498 + Icges(7,5) * t526;
t470 = Icges(7,1) * t497 + Icges(7,4) * t496 + Icges(7,5) * t524;
t469 = Icges(7,4) * t499 + Icges(7,2) * t498 + Icges(7,6) * t526;
t468 = Icges(7,4) * t497 + Icges(7,2) * t496 + Icges(7,6) * t524;
t467 = Icges(7,5) * t499 + Icges(7,6) * t498 + Icges(7,3) * t526;
t466 = Icges(7,5) * t497 + Icges(7,6) * t496 + Icges(7,3) * t524;
t465 = t507 * t570 - t523 * t551 + t639;
t464 = -t506 * t570 + t523 * t552 + t637;
t463 = t506 * t551 - t507 * t552 + t644;
t462 = t485 * t529 - t494 * t518 + t638;
t461 = -t484 * t529 + t494 * t519 + t636;
t460 = t484 * t518 - t485 * t519 + t640;
t459 = t473 * t511 - t477 * t486 + t489 * t529 - t512 * t518 + t638;
t458 = -t472 * t511 + t477 * t487 - t488 * t529 + t512 * t519 + t636;
t457 = t472 * t486 - t473 * t487 + t488 * t518 - t489 * t519 + t640;
t1 = m(3) * (t553 ^ 2 + t561 ^ 2 + t562 ^ 2) / 0.2e1 + t511 * ((t467 * t530 + t469 * t516 + t471 * t517) * t486 + (t466 * t530 + t468 * t516 + t470 * t517) * t487 + (t474 * t530 + t475 * t516 + t476 * t517) * t511) / 0.2e1 + t486 * ((t467 * t526 + t469 * t498 + t471 * t499) * t486 + (t466 * t526 + t468 * t498 + t470 * t499) * t487 + (t474 * t526 + t475 * t498 + t476 * t499) * t511) / 0.2e1 + t529 * ((t479 * t555 - t481 * t530 + t483 * t531) * t518 + (t478 * t555 - t480 * t530 + t482 * t531) * t519 + (t491 * t555 - t492 * t530 + t493 * t531) * t529) / 0.2e1 + t487 * ((t467 * t524 + t469 * t496 + t471 * t497) * t486 + (t466 * t524 + t468 * t496 + t470 * t497) * t487 + (t474 * t524 + t475 * t496 + t476 * t497) * t511) / 0.2e1 + t519 * ((t479 * t537 - t481 * t524 + t483 * t525) * t518 + (t478 * t537 - t480 * t524 + t482 * t525) * t519 + (t491 * t537 - t492 * t524 + t493 * t525) * t529) / 0.2e1 + t518 * ((t479 * t539 - t481 * t526 + t483 * t527) * t518 + (t478 * t539 - t480 * t526 + t482 * t527) * t519 + (t491 * t539 - t492 * t526 + t493 * t527) * t529) / 0.2e1 + t552 * ((t501 * t564 - t503 * t537 + t505 * t538) * t551 + (t500 * t564 - t502 * t537 + t504 * t538) * t552 + (t520 * t564 - t521 * t537 + t522 * t538) * t570) / 0.2e1 + t570 * ((t501 * t576 - t503 * t555 + t505 * t556) * t551 + (t500 * t576 - t502 * t555 + t504 * t556) * t552 + (t520 * t576 - t521 * t555 + t522 * t556) * t570) / 0.2e1 + t551 * ((t501 * t565 - t503 * t539 + t505 * t540) * t551 + (t500 * t565 - t502 * t539 + t504 * t540) * t552 + (t520 * t565 - t521 * t539 + t522 * t540) * t570) / 0.2e1 + m(7) * (t457 ^ 2 + t458 ^ 2 + t459 ^ 2) / 0.2e1 + m(6) * (t460 ^ 2 + t461 ^ 2 + t462 ^ 2) / 0.2e1 + m(5) * (t463 ^ 2 + t464 ^ 2 + t465 ^ 2) / 0.2e1 + m(4) * (t490 ^ 2 + t508 ^ 2 + t509 ^ 2) / 0.2e1 + (((t544 * t609 + t546 * t596 + t548 * t597) * t631 - (t543 * t609 + t545 * t596 + t547 * t597) * t634 + ((t588 * t633 + t590 * t630) * t631 - (t587 * t633 + t589 * t630) * t634) * t623 - t643 * t626) * t655 + (t566 * t609 + t567 * t596 + t568 * t597 + t604 * t626 + (t605 * t633 + t606 * t630) * t623) * t619) * t619 / 0.2e1 + (m(2) * (t616 ^ 2 + t617 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t543 * t601 - t545 * t581 - t547 * t582 - t587 * t612 - t589 * t613) * t634 + (t544 * t601 + t546 * t581 + t548 * t582 + t588 * t612 + t590 * t613 - t672) * t631) * t655 + (t566 * t601 + t567 * t581 + t568 * t582 + t604 * t663 + t605 * t612 + t606 * t613) * t619) * t618 / 0.2e1 - (((-t543 * t600 - t545 * t579 - t547 * t580 - t587 * t610 - t589 * t611 + t672) * t634 + (t544 * t600 + t546 * t579 + t548 * t580 + t588 * t610 + t590 * t611) * t631) * t655 + (t566 * t600 + t567 * t579 + t568 * t580 - t604 * t662 + t605 * t610 + t606 * t611) * t619) * t652 / 0.2e1;
T  = t1;
