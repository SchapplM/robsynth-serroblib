% Calculate kinetic energy for
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 06:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_energykin_fixb_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR10_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR10_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:49:41
% EndTime: 2019-03-10 05:49:44
% DurationCPUTime: 3.16s
% Computational Cost: add. (9519->395), mult. (26651->626), div. (0->0), fcn. (34988->18), ass. (0->181)
t622 = sin(pkin(6));
t624 = cos(pkin(6));
t633 = cos(qJ(2));
t634 = cos(qJ(1));
t653 = t633 * t634;
t629 = sin(qJ(2));
t630 = sin(qJ(1));
t656 = t629 * t630;
t610 = t624 * t653 - t656;
t654 = t630 * t633;
t655 = t629 * t634;
t611 = t624 * t655 + t654;
t612 = -t624 * t654 - t655;
t613 = -t624 * t656 + t653;
t658 = t622 * t634;
t659 = t622 * t630;
t646 = (Icges(3,5) * t611 + Icges(3,6) * t610 - Icges(3,3) * t658) * t634 - (Icges(3,5) * t613 + Icges(3,6) * t612 + Icges(3,3) * t659) * t630;
t666 = t622 * t646;
t664 = cos(qJ(4));
t663 = cos(qJ(5));
t662 = cos(pkin(8));
t661 = sin(pkin(8));
t621 = sin(pkin(7));
t660 = t621 * t624;
t623 = cos(pkin(7));
t657 = t623 * t633;
t599 = -t610 * t621 - t623 * t658;
t580 = pkin(2) * t611 + pkin(11) * t599;
t600 = -t612 * t621 + t623 * t659;
t581 = pkin(2) * t613 + pkin(11) * t600;
t650 = qJD(2) * t622;
t618 = t630 * t650;
t649 = t634 * t650;
t652 = t580 * t618 + t581 * t649;
t590 = qJD(3) * t600 + t618;
t651 = qJD(1) * (pkin(1) * t630 - pkin(10) * t658);
t619 = qJD(2) * t624 + qJD(1);
t628 = sin(qJ(3));
t632 = cos(qJ(3));
t643 = t612 * t623 + t621 * t659;
t578 = -t613 * t628 + t632 * t643;
t561 = -t578 * t661 + t600 * t662;
t547 = qJD(4) * t561 + t590;
t609 = -t621 * t622 * t633 + t623 * t624;
t601 = qJD(3) * t609 + t619;
t579 = t613 * t632 + t628 * t643;
t627 = sin(qJ(4));
t647 = t664 * t661;
t648 = t662 * t664;
t535 = -t578 * t648 + t579 * t627 - t600 * t647;
t513 = qJD(5) * t535 + t547;
t597 = t632 * t660 + (-t628 * t629 + t632 * t657) * t622;
t575 = -t597 * t661 + t609 * t662;
t567 = qJD(4) * t575 + t601;
t591 = qJD(3) * t599 - t649;
t598 = t628 * t660 + (t628 * t657 + t629 * t632) * t622;
t552 = -t597 * t648 + t598 * t627 - t609 * t647;
t527 = qJD(5) * t552 + t567;
t644 = t610 * t623 - t621 * t658;
t576 = -t611 * t628 + t632 * t644;
t560 = -t576 * t661 + t599 * t662;
t577 = t611 * t632 + t628 * t644;
t537 = t577 * pkin(3) + pkin(12) * t560;
t538 = t579 * pkin(3) + pkin(12) * t561;
t645 = t590 * t537 - t538 * t591 + t652;
t548 = qJD(4) * t560 + t591;
t602 = pkin(2) * t622 * t629 + pkin(11) * t609;
t614 = qJD(1) * (pkin(1) * t634 + pkin(10) * t659);
t642 = t619 * t581 - t602 * t618 + t614;
t533 = -t576 * t648 + t577 * t627 - t599 * t647;
t514 = qJD(5) * t533 + t548;
t534 = t577 * t664 + (t576 * t662 + t599 * t661) * t627;
t511 = pkin(4) * t534 + pkin(13) * t533;
t536 = t579 * t664 + (t578 * t662 + t600 * t661) * t627;
t512 = pkin(4) * t536 + pkin(13) * t535;
t641 = t547 * t511 - t512 * t548 + t645;
t640 = -t580 * t619 - t602 * t649 - t651;
t562 = t598 * pkin(3) + pkin(12) * t575;
t639 = t601 * t538 - t562 * t590 + t642;
t638 = -t537 * t601 + t591 * t562 + t640;
t553 = t598 * t664 + (t597 * t662 + t609 * t661) * t627;
t525 = pkin(4) * t553 + pkin(13) * t552;
t637 = t567 * t512 - t525 * t547 + t639;
t636 = -t511 * t567 + t548 * t525 + t638;
t631 = cos(qJ(6));
t626 = sin(qJ(5));
t625 = sin(qJ(6));
t617 = rSges(2,1) * t634 - rSges(2,2) * t630;
t616 = rSges(2,1) * t630 + rSges(2,2) * t634;
t607 = rSges(3,3) * t624 + (rSges(3,1) * t629 + rSges(3,2) * t633) * t622;
t606 = Icges(3,5) * t624 + (Icges(3,1) * t629 + Icges(3,4) * t633) * t622;
t605 = Icges(3,6) * t624 + (Icges(3,4) * t629 + Icges(3,2) * t633) * t622;
t604 = Icges(3,3) * t624 + (Icges(3,5) * t629 + Icges(3,6) * t633) * t622;
t589 = rSges(3,1) * t613 + rSges(3,2) * t612 + rSges(3,3) * t659;
t588 = rSges(3,1) * t611 + rSges(3,2) * t610 - rSges(3,3) * t658;
t587 = Icges(3,1) * t613 + Icges(3,4) * t612 + Icges(3,5) * t659;
t586 = Icges(3,1) * t611 + Icges(3,4) * t610 - Icges(3,5) * t658;
t585 = Icges(3,4) * t613 + Icges(3,2) * t612 + Icges(3,6) * t659;
t584 = Icges(3,4) * t611 + Icges(3,2) * t610 - Icges(3,6) * t658;
t566 = rSges(4,1) * t598 + rSges(4,2) * t597 + rSges(4,3) * t609;
t565 = Icges(4,1) * t598 + Icges(4,4) * t597 + Icges(4,5) * t609;
t564 = Icges(4,4) * t598 + Icges(4,2) * t597 + Icges(4,6) * t609;
t563 = Icges(4,5) * t598 + Icges(4,6) * t597 + Icges(4,3) * t609;
t557 = t589 * t619 - t607 * t618 + t614;
t556 = -t588 * t619 - t607 * t649 - t651;
t550 = (t588 * t630 + t589 * t634) * t650;
t546 = rSges(4,1) * t579 + rSges(4,2) * t578 + rSges(4,3) * t600;
t545 = rSges(4,1) * t577 + rSges(4,2) * t576 + rSges(4,3) * t599;
t544 = Icges(4,1) * t579 + Icges(4,4) * t578 + Icges(4,5) * t600;
t543 = Icges(4,1) * t577 + Icges(4,4) * t576 + Icges(4,5) * t599;
t542 = Icges(4,4) * t579 + Icges(4,2) * t578 + Icges(4,6) * t600;
t541 = Icges(4,4) * t577 + Icges(4,2) * t576 + Icges(4,6) * t599;
t540 = Icges(4,5) * t579 + Icges(4,6) * t578 + Icges(4,3) * t600;
t539 = Icges(4,5) * t577 + Icges(4,6) * t576 + Icges(4,3) * t599;
t530 = t553 * t663 + t575 * t626;
t529 = t553 * t626 - t575 * t663;
t524 = t536 * t663 + t561 * t626;
t523 = t536 * t626 - t561 * t663;
t522 = t534 * t663 + t560 * t626;
t521 = t534 * t626 - t560 * t663;
t520 = rSges(5,1) * t553 - rSges(5,2) * t552 + rSges(5,3) * t575;
t519 = Icges(5,1) * t553 - Icges(5,4) * t552 + Icges(5,5) * t575;
t518 = Icges(5,4) * t553 - Icges(5,2) * t552 + Icges(5,6) * t575;
t517 = Icges(5,5) * t553 - Icges(5,6) * t552 + Icges(5,3) * t575;
t516 = t530 * t631 + t552 * t625;
t515 = -t530 * t625 + t552 * t631;
t510 = pkin(5) * t530 + pkin(14) * t529;
t508 = qJD(6) * t529 + t527;
t506 = t546 * t601 - t566 * t590 + t642;
t505 = -t545 * t601 + t566 * t591 + t640;
t504 = rSges(5,1) * t536 - rSges(5,2) * t535 + rSges(5,3) * t561;
t503 = rSges(5,1) * t534 - rSges(5,2) * t533 + rSges(5,3) * t560;
t502 = Icges(5,1) * t536 - Icges(5,4) * t535 + Icges(5,5) * t561;
t501 = Icges(5,1) * t534 - Icges(5,4) * t533 + Icges(5,5) * t560;
t500 = Icges(5,4) * t536 - Icges(5,2) * t535 + Icges(5,6) * t561;
t499 = Icges(5,4) * t534 - Icges(5,2) * t533 + Icges(5,6) * t560;
t498 = Icges(5,5) * t536 - Icges(5,6) * t535 + Icges(5,3) * t561;
t497 = Icges(5,5) * t534 - Icges(5,6) * t533 + Icges(5,3) * t560;
t496 = t524 * t631 + t535 * t625;
t495 = -t524 * t625 + t535 * t631;
t494 = t522 * t631 + t533 * t625;
t493 = -t522 * t625 + t533 * t631;
t492 = rSges(6,1) * t530 - rSges(6,2) * t529 + rSges(6,3) * t552;
t491 = Icges(6,1) * t530 - Icges(6,4) * t529 + Icges(6,5) * t552;
t490 = Icges(6,4) * t530 - Icges(6,2) * t529 + Icges(6,6) * t552;
t489 = Icges(6,5) * t530 - Icges(6,6) * t529 + Icges(6,3) * t552;
t487 = t545 * t590 - t546 * t591 + t652;
t486 = pkin(5) * t524 + pkin(14) * t523;
t485 = pkin(5) * t522 + pkin(14) * t521;
t484 = qJD(6) * t521 + t514;
t483 = qJD(6) * t523 + t513;
t482 = rSges(6,1) * t524 - rSges(6,2) * t523 + rSges(6,3) * t535;
t481 = rSges(6,1) * t522 - rSges(6,2) * t521 + rSges(6,3) * t533;
t480 = Icges(6,1) * t524 - Icges(6,4) * t523 + Icges(6,5) * t535;
t479 = Icges(6,1) * t522 - Icges(6,4) * t521 + Icges(6,5) * t533;
t478 = Icges(6,4) * t524 - Icges(6,2) * t523 + Icges(6,6) * t535;
t477 = Icges(6,4) * t522 - Icges(6,2) * t521 + Icges(6,6) * t533;
t476 = Icges(6,5) * t524 - Icges(6,6) * t523 + Icges(6,3) * t535;
t475 = Icges(6,5) * t522 - Icges(6,6) * t521 + Icges(6,3) * t533;
t474 = rSges(7,1) * t516 + rSges(7,2) * t515 + rSges(7,3) * t529;
t473 = Icges(7,1) * t516 + Icges(7,4) * t515 + Icges(7,5) * t529;
t472 = Icges(7,4) * t516 + Icges(7,2) * t515 + Icges(7,6) * t529;
t471 = Icges(7,5) * t516 + Icges(7,6) * t515 + Icges(7,3) * t529;
t470 = rSges(7,1) * t496 + rSges(7,2) * t495 + rSges(7,3) * t523;
t469 = rSges(7,1) * t494 + rSges(7,2) * t493 + rSges(7,3) * t521;
t468 = Icges(7,1) * t496 + Icges(7,4) * t495 + Icges(7,5) * t523;
t467 = Icges(7,1) * t494 + Icges(7,4) * t493 + Icges(7,5) * t521;
t466 = Icges(7,4) * t496 + Icges(7,2) * t495 + Icges(7,6) * t523;
t465 = Icges(7,4) * t494 + Icges(7,2) * t493 + Icges(7,6) * t521;
t464 = Icges(7,5) * t496 + Icges(7,6) * t495 + Icges(7,3) * t523;
t463 = Icges(7,5) * t494 + Icges(7,6) * t493 + Icges(7,3) * t521;
t462 = t504 * t567 - t520 * t547 + t639;
t461 = -t503 * t567 + t520 * t548 + t638;
t460 = t503 * t547 - t504 * t548 + t645;
t459 = t482 * t527 - t492 * t513 + t637;
t458 = -t481 * t527 + t492 * t514 + t636;
t457 = t481 * t513 - t482 * t514 + t641;
t456 = t470 * t508 - t474 * t483 + t486 * t527 - t510 * t513 + t637;
t455 = -t469 * t508 + t474 * t484 - t485 * t527 + t510 * t514 + t636;
t454 = t469 * t483 - t470 * t484 + t485 * t513 - t486 * t514 + t641;
t1 = ((t604 * t659 + t605 * t612 + t606 * t613) * t619 + (-(t584 * t612 + t586 * t613) * t634 + (t612 * t585 + t613 * t587 - t666) * t630) * t650) * t618 / 0.2e1 - ((-t604 * t658 + t605 * t610 + t606 * t611) * t619 + ((t585 * t610 + t587 * t611) * t630 + (-t610 * t584 - t611 * t586 + t666) * t634) * t650) * t649 / 0.2e1 + m(7) * (t454 ^ 2 + t455 ^ 2 + t456 ^ 2) / 0.2e1 + m(6) * (t457 ^ 2 + t458 ^ 2 + t459 ^ 2) / 0.2e1 + t527 * ((t476 * t552 - t478 * t529 + t480 * t530) * t513 + (t475 * t552 - t477 * t529 + t479 * t530) * t514 + (t552 * t489 - t529 * t490 + t530 * t491) * t527) / 0.2e1 + t513 * ((t535 * t476 - t523 * t478 + t524 * t480) * t513 + (t475 * t535 - t477 * t523 + t479 * t524) * t514 + (t489 * t535 - t490 * t523 + t491 * t524) * t527) / 0.2e1 + t514 * ((t476 * t533 - t478 * t521 + t480 * t522) * t513 + (t533 * t475 - t521 * t477 + t522 * t479) * t514 + (t489 * t533 - t490 * t521 + t491 * t522) * t527) / 0.2e1 + t548 * ((t498 * t560 - t500 * t533 + t502 * t534) * t547 + (t497 * t560 - t499 * t533 + t501 * t534) * t548 + (t517 * t560 - t518 * t533 + t519 * t534) * t567) / 0.2e1 + t567 * ((t498 * t575 - t500 * t552 + t502 * t553) * t547 + (t497 * t575 - t499 * t552 + t501 * t553) * t548 + (t517 * t575 - t552 * t518 + t519 * t553) * t567) / 0.2e1 + t547 * ((t498 * t561 - t500 * t535 + t502 * t536) * t547 + (t497 * t561 - t499 * t535 + t501 * t536) * t548 + (t517 * t561 - t518 * t535 + t519 * t536) * t567) / 0.2e1 + t590 * ((t540 * t600 + t578 * t542 + t579 * t544) * t590 + (t539 * t600 + t541 * t578 + t543 * t579) * t591 + (t563 * t600 + t564 * t578 + t565 * t579) * t601) / 0.2e1 + t591 * ((t540 * t599 + t542 * t576 + t544 * t577) * t590 + (t539 * t599 + t541 * t576 + t543 * t577) * t591 + (t563 * t599 + t564 * t576 + t565 * t577) * t601) / 0.2e1 + t601 * ((t540 * t609 + t542 * t597 + t544 * t598) * t590 + (t539 * t609 + t541 * t597 + t543 * t598) * t591 + (t563 * t609 + t564 * t597 + t565 * t598) * t601) / 0.2e1 + t619 * ((t624 * t604 + (t605 * t633 + t606 * t629) * t622) * t619 + (((t585 * t633 + t587 * t629) * t630 - (t584 * t633 + t586 * t629) * t634) * t622 - t646 * t624) * t650) / 0.2e1 + t483 * ((t523 * t464 + t495 * t466 + t496 * t468) * t483 + (t463 * t523 + t465 * t495 + t467 * t496) * t484 + (t471 * t523 + t472 * t495 + t473 * t496) * t508) / 0.2e1 + t484 * ((t464 * t521 + t466 * t493 + t468 * t494) * t483 + (t521 * t463 + t493 * t465 + t494 * t467) * t484 + (t471 * t521 + t472 * t493 + t473 * t494) * t508) / 0.2e1 + t508 * ((t464 * t529 + t466 * t515 + t468 * t516) * t483 + (t463 * t529 + t465 * t515 + t467 * t516) * t484 + (t529 * t471 + t515 * t472 + t516 * t473) * t508) / 0.2e1 + m(5) * (t460 ^ 2 + t461 ^ 2 + t462 ^ 2) / 0.2e1 + m(4) * (t487 ^ 2 + t505 ^ 2 + t506 ^ 2) / 0.2e1 + m(3) * (t550 ^ 2 + t556 ^ 2 + t557 ^ 2) / 0.2e1 + (m(2) * (t616 ^ 2 + t617 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
