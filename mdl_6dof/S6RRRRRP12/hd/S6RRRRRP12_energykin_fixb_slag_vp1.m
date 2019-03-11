% Calculate kinetic energy for
% S6RRRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP12_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP12_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP12_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP12_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP12_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:59:37
% EndTime: 2019-03-10 02:59:40
% DurationCPUTime: 3.37s
% Computational Cost: add. (5351->336), mult. (14547->522), div. (0->0), fcn. (18749->14), ass. (0->161)
t654 = Icges(6,1) + Icges(7,1);
t653 = -Icges(6,4) + Icges(7,5);
t652 = Icges(7,4) + Icges(6,5);
t651 = Icges(6,2) + Icges(7,3);
t650 = Icges(7,2) + Icges(6,3);
t649 = -Icges(6,6) + Icges(7,6);
t648 = rSges(7,1) + pkin(5);
t647 = rSges(7,3) + qJ(6);
t590 = sin(pkin(6));
t591 = cos(pkin(6));
t597 = cos(qJ(2));
t598 = cos(qJ(1));
t621 = t597 * t598;
t595 = sin(qJ(2));
t596 = sin(qJ(1));
t624 = t595 * t596;
t576 = t591 * t621 - t624;
t622 = t596 * t597;
t623 = t595 * t598;
t577 = t591 * t623 + t622;
t578 = -t591 * t622 - t623;
t579 = -t591 * t624 + t621;
t625 = t590 * t598;
t626 = t590 * t596;
t609 = (Icges(3,5) * t577 + Icges(3,6) * t576 - Icges(3,3) * t625) * t598 - (Icges(3,5) * t579 + Icges(3,6) * t578 + Icges(3,3) * t626) * t596;
t646 = t590 * t609;
t594 = sin(qJ(3));
t628 = sin(pkin(7));
t612 = t590 * t628;
t629 = cos(pkin(7));
t632 = cos(qJ(3));
t540 = t577 * t632 + (t576 * t629 - t598 * t612) * t594;
t613 = t590 * t629;
t563 = -t576 * t628 - t598 * t613;
t593 = sin(qJ(4));
t631 = cos(qJ(4));
t521 = t540 * t631 + t563 * t593;
t610 = t632 * t628;
t608 = t590 * t610;
t611 = t629 * t632;
t539 = -t576 * t611 + t577 * t594 + t598 * t608;
t592 = sin(qJ(5));
t630 = cos(qJ(5));
t491 = t521 * t592 - t539 * t630;
t492 = t521 * t630 + t539 * t592;
t520 = t540 * t593 - t563 * t631;
t645 = t651 * t491 + t653 * t492 + t649 * t520;
t542 = t579 * t632 + (t578 * t629 + t596 * t612) * t594;
t564 = -t578 * t628 + t596 * t613;
t523 = t542 * t631 + t564 * t593;
t541 = -t578 * t611 + t579 * t594 - t596 * t608;
t493 = t523 * t592 - t541 * t630;
t494 = t523 * t630 + t541 * t592;
t522 = t542 * t593 - t564 * t631;
t644 = t651 * t493 + t653 * t494 + t649 * t522;
t643 = t649 * t491 + t652 * t492 + t650 * t520;
t642 = t649 * t493 + t652 * t494 + t650 * t522;
t641 = t653 * t491 + t654 * t492 + t652 * t520;
t640 = t653 * t493 + t654 * t494 + t652 * t522;
t562 = t591 * t628 * t594 + (t594 * t597 * t629 + t595 * t632) * t590;
t575 = t591 * t629 - t597 * t612;
t538 = t562 * t631 + t575 * t593;
t627 = t590 * t595;
t561 = -t590 * t597 * t611 - t591 * t610 + t594 * t627;
t518 = t538 * t592 - t561 * t630;
t519 = t538 * t630 + t561 * t592;
t537 = t562 * t593 - t575 * t631;
t639 = t651 * t518 + t653 * t519 + t649 * t537;
t638 = t649 * t518 + t652 * t519 + t650 * t537;
t637 = t653 * t518 + t654 * t519 + t652 * t537;
t620 = rSges(7,2) * t520 + t647 * t491 + t492 * t648;
t619 = rSges(7,2) * t522 + t647 * t493 + t494 * t648;
t618 = rSges(7,2) * t537 + t647 * t518 + t519 * t648;
t543 = t577 * pkin(2) + pkin(10) * t563;
t544 = t579 * pkin(2) + pkin(10) * t564;
t615 = qJD(2) * t590;
t587 = t596 * t615;
t614 = t598 * t615;
t617 = t543 * t587 + t544 * t614;
t553 = qJD(3) * t564 + t587;
t616 = qJD(1) * (pkin(1) * t596 - pkin(9) * t625);
t588 = qJD(2) * t591 + qJD(1);
t514 = qJD(4) * t541 + t553;
t565 = qJD(3) * t575 + t588;
t529 = qJD(4) * t561 + t565;
t554 = qJD(3) * t563 - t614;
t510 = pkin(3) * t540 + pkin(11) * t539;
t511 = pkin(3) * t542 + pkin(11) * t541;
t607 = t553 * t510 - t511 * t554 + t617;
t515 = qJD(4) * t539 + t554;
t566 = pkin(2) * t627 + pkin(10) * t575;
t580 = qJD(1) * (pkin(1) * t598 + pkin(9) * t626);
t606 = t588 * t544 - t566 * t587 + t580;
t488 = pkin(4) * t521 + pkin(12) * t520;
t489 = pkin(4) * t523 + pkin(12) * t522;
t605 = t514 * t488 - t489 * t515 + t607;
t604 = -t543 * t588 - t566 * t614 - t616;
t528 = pkin(3) * t562 + pkin(11) * t561;
t603 = t565 * t511 - t528 * t553 + t606;
t602 = -t510 * t565 + t554 * t528 + t604;
t509 = pkin(4) * t538 + pkin(12) * t537;
t601 = t529 * t489 - t509 * t514 + t603;
t600 = -t488 * t529 + t515 * t509 + t602;
t585 = rSges(2,1) * t598 - rSges(2,2) * t596;
t584 = rSges(2,1) * t596 + rSges(2,2) * t598;
t572 = rSges(3,3) * t591 + (rSges(3,1) * t595 + rSges(3,2) * t597) * t590;
t571 = Icges(3,5) * t591 + (Icges(3,1) * t595 + Icges(3,4) * t597) * t590;
t570 = Icges(3,6) * t591 + (Icges(3,4) * t595 + Icges(3,2) * t597) * t590;
t569 = Icges(3,3) * t591 + (Icges(3,5) * t595 + Icges(3,6) * t597) * t590;
t552 = rSges(3,1) * t579 + rSges(3,2) * t578 + rSges(3,3) * t626;
t551 = rSges(3,1) * t577 + rSges(3,2) * t576 - rSges(3,3) * t625;
t550 = Icges(3,1) * t579 + Icges(3,4) * t578 + Icges(3,5) * t626;
t549 = Icges(3,1) * t577 + Icges(3,4) * t576 - Icges(3,5) * t625;
t548 = Icges(3,4) * t579 + Icges(3,2) * t578 + Icges(3,6) * t626;
t547 = Icges(3,4) * t577 + Icges(3,2) * t576 - Icges(3,6) * t625;
t527 = rSges(4,1) * t562 - rSges(4,2) * t561 + rSges(4,3) * t575;
t526 = Icges(4,1) * t562 - Icges(4,4) * t561 + Icges(4,5) * t575;
t525 = Icges(4,4) * t562 - Icges(4,2) * t561 + Icges(4,6) * t575;
t524 = Icges(4,5) * t562 - Icges(4,6) * t561 + Icges(4,3) * t575;
t517 = t552 * t588 - t572 * t587 + t580;
t516 = -t551 * t588 - t572 * t614 - t616;
t513 = (t551 * t596 + t552 * t598) * t615;
t508 = qJD(5) * t537 + t529;
t506 = rSges(4,1) * t542 - rSges(4,2) * t541 + rSges(4,3) * t564;
t505 = rSges(4,1) * t540 - rSges(4,2) * t539 + rSges(4,3) * t563;
t504 = Icges(4,1) * t542 - Icges(4,4) * t541 + Icges(4,5) * t564;
t503 = Icges(4,1) * t540 - Icges(4,4) * t539 + Icges(4,5) * t563;
t502 = Icges(4,4) * t542 - Icges(4,2) * t541 + Icges(4,6) * t564;
t501 = Icges(4,4) * t540 - Icges(4,2) * t539 + Icges(4,6) * t563;
t500 = Icges(4,5) * t542 - Icges(4,6) * t541 + Icges(4,3) * t564;
t499 = Icges(4,5) * t540 - Icges(4,6) * t539 + Icges(4,3) * t563;
t498 = rSges(5,1) * t538 - rSges(5,2) * t537 + rSges(5,3) * t561;
t497 = Icges(5,1) * t538 - Icges(5,4) * t537 + Icges(5,5) * t561;
t496 = Icges(5,4) * t538 - Icges(5,2) * t537 + Icges(5,6) * t561;
t495 = Icges(5,5) * t538 - Icges(5,6) * t537 + Icges(5,3) * t561;
t487 = qJD(5) * t520 + t515;
t486 = qJD(5) * t522 + t514;
t483 = rSges(5,1) * t523 - rSges(5,2) * t522 + rSges(5,3) * t541;
t482 = rSges(5,1) * t521 - rSges(5,2) * t520 + rSges(5,3) * t539;
t481 = Icges(5,1) * t523 - Icges(5,4) * t522 + Icges(5,5) * t541;
t480 = Icges(5,1) * t521 - Icges(5,4) * t520 + Icges(5,5) * t539;
t479 = Icges(5,4) * t523 - Icges(5,2) * t522 + Icges(5,6) * t541;
t478 = Icges(5,4) * t521 - Icges(5,2) * t520 + Icges(5,6) * t539;
t477 = Icges(5,5) * t523 - Icges(5,6) * t522 + Icges(5,3) * t541;
t476 = Icges(5,5) * t521 - Icges(5,6) * t520 + Icges(5,3) * t539;
t474 = rSges(6,1) * t519 - rSges(6,2) * t518 + rSges(6,3) * t537;
t463 = rSges(6,1) * t494 - rSges(6,2) * t493 + rSges(6,3) * t522;
t461 = rSges(6,1) * t492 - rSges(6,2) * t491 + rSges(6,3) * t520;
t459 = t506 * t565 - t527 * t553 + t606;
t458 = -t505 * t565 + t527 * t554 + t604;
t445 = t505 * t553 - t506 * t554 + t617;
t444 = t483 * t529 - t498 * t514 + t603;
t443 = -t482 * t529 + t498 * t515 + t602;
t442 = t482 * t514 - t483 * t515 + t607;
t441 = t463 * t508 - t474 * t486 + t601;
t440 = -t461 * t508 + t474 * t487 + t600;
t439 = t461 * t486 - t463 * t487 + t605;
t438 = qJD(6) * t491 - t486 * t618 + t508 * t619 + t601;
t437 = qJD(6) * t493 + t487 * t618 - t508 * t620 + t600;
t436 = qJD(6) * t518 + t486 * t620 - t487 * t619 + t605;
t1 = ((t569 * t626 + t570 * t578 + t571 * t579) * t588 + (-(t547 * t578 + t549 * t579) * t598 + (t578 * t548 + t579 * t550 - t646) * t596) * t615) * t587 / 0.2e1 - ((-t569 * t625 + t570 * t576 + t571 * t577) * t588 + ((t548 * t576 + t550 * t577) * t596 + (-t576 * t547 - t577 * t549 + t646) * t598) * t615) * t614 / 0.2e1 + t515 * ((t477 * t539 - t479 * t520 + t481 * t521) * t514 + (t476 * t539 - t478 * t520 + t480 * t521) * t515 + (t495 * t539 - t496 * t520 + t497 * t521) * t529) / 0.2e1 + t514 * ((t477 * t541 - t479 * t522 + t481 * t523) * t514 + (t476 * t541 - t478 * t522 + t480 * t523) * t515 + (t495 * t541 - t496 * t522 + t497 * t523) * t529) / 0.2e1 + t529 * ((t477 * t561 - t479 * t537 + t481 * t538) * t514 + (t476 * t561 - t478 * t537 + t480 * t538) * t515 + (t495 * t561 - t496 * t537 + t497 * t538) * t529) / 0.2e1 + t565 * ((t500 * t575 - t502 * t561 + t504 * t562) * t553 + (t499 * t575 - t501 * t561 + t503 * t562) * t554 + (t524 * t575 - t525 * t561 + t526 * t562) * t565) / 0.2e1 + t554 * ((t500 * t563 - t502 * t539 + t504 * t540) * t553 + (t499 * t563 - t501 * t539 + t503 * t540) * t554 + (t524 * t563 - t525 * t539 + t526 * t540) * t565) / 0.2e1 + t553 * ((t500 * t564 - t502 * t541 + t504 * t542) * t553 + (t499 * t564 - t501 * t541 + t503 * t542) * t554 + (t524 * t564 - t525 * t541 + t526 * t542) * t565) / 0.2e1 + t588 * ((t591 * t569 + (t570 * t597 + t571 * t595) * t590) * t588 + (((t548 * t597 + t550 * t595) * t596 - (t547 * t597 + t549 * t595) * t598) * t590 - t609 * t591) * t615) / 0.2e1 + m(6) * (t439 ^ 2 + t440 ^ 2 + t441 ^ 2) / 0.2e1 + m(7) * (t436 ^ 2 + t437 ^ 2 + t438 ^ 2) / 0.2e1 + m(5) * (t442 ^ 2 + t443 ^ 2 + t444 ^ 2) / 0.2e1 + m(4) * (t445 ^ 2 + t458 ^ 2 + t459 ^ 2) / 0.2e1 + m(3) * (t513 ^ 2 + t516 ^ 2 + t517 ^ 2) / 0.2e1 + ((t493 * t639 + t494 * t637 + t522 * t638) * t508 + (t493 * t645 + t494 * t641 + t522 * t643) * t487 + (t644 * t493 + t640 * t494 + t642 * t522) * t486) * t486 / 0.2e1 + ((t491 * t639 + t492 * t637 + t520 * t638) * t508 + (t645 * t491 + t641 * t492 + t643 * t520) * t487 + (t491 * t644 + t492 * t640 + t520 * t642) * t486) * t487 / 0.2e1 + ((t639 * t518 + t637 * t519 + t638 * t537) * t508 + (t518 * t645 + t519 * t641 + t537 * t643) * t487 + (t518 * t644 + t519 * t640 + t537 * t642) * t486) * t508 / 0.2e1 + (Icges(2,3) + m(2) * (t584 ^ 2 + t585 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
