% Calculate kinetic energy for
% S6RRRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 23:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR12_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR12_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_energykin_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR12_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR12_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR12_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:29:10
% EndTime: 2019-03-09 23:29:13
% DurationCPUTime: 3.47s
% Computational Cost: add. (5223->387), mult. (12805->589), div. (0->0), fcn. (16360->16), ass. (0->176)
t646 = Icges(5,3) + Icges(6,3);
t590 = sin(pkin(6));
t591 = cos(pkin(6));
t600 = cos(qJ(2));
t601 = cos(qJ(1));
t623 = t600 * t601;
t596 = sin(qJ(2));
t597 = sin(qJ(1));
t626 = t596 * t597;
t574 = t591 * t623 - t626;
t624 = t597 * t600;
t625 = t596 * t601;
t575 = t591 * t625 + t624;
t576 = -t591 * t624 - t625;
t577 = -t591 * t626 + t623;
t627 = t590 * t601;
t628 = t590 * t597;
t612 = (Icges(3,5) * t575 + Icges(3,6) * t574 - Icges(3,3) * t627) * t601 - (Icges(3,5) * t577 + Icges(3,6) * t576 + Icges(3,3) * t628) * t597;
t645 = t590 * t612;
t595 = sin(qJ(3));
t633 = sin(pkin(7));
t616 = t590 * t633;
t634 = cos(pkin(7));
t637 = cos(qJ(3));
t538 = t575 * t637 + (t574 * t634 - t601 * t616) * t595;
t617 = t590 * t634;
t561 = -t574 * t633 - t601 * t617;
t619 = qJ(4) + pkin(13);
t589 = sin(t619);
t615 = cos(t619);
t510 = t538 * t589 - t561 * t615;
t511 = t538 * t615 + t561 * t589;
t594 = sin(qJ(4));
t599 = cos(qJ(4));
t516 = -t538 * t594 + t561 * t599;
t632 = t561 * t594;
t517 = t538 * t599 + t632;
t613 = t637 * t633;
t611 = t590 * t613;
t614 = t634 * t637;
t537 = -t574 * t614 + t575 * t595 + t601 * t611;
t644 = Icges(5,5) * t517 + Icges(6,5) * t511 + Icges(5,6) * t516 - Icges(6,6) * t510 + t646 * t537;
t540 = t577 * t637 + (t576 * t634 + t597 * t616) * t595;
t562 = -t576 * t633 + t597 * t617;
t512 = t540 * t589 - t562 * t615;
t513 = t540 * t615 + t562 * t589;
t518 = -t540 * t594 + t562 * t599;
t631 = t562 * t594;
t519 = t540 * t599 + t631;
t539 = -t576 * t614 + t577 * t595 - t597 * t611;
t643 = Icges(5,5) * t519 + Icges(6,5) * t513 + Icges(5,6) * t518 - Icges(6,6) * t512 + t646 * t539;
t560 = t591 * t633 * t595 + (t595 * t600 * t634 + t596 * t637) * t590;
t573 = t591 * t634 - t600 * t616;
t525 = t560 * t589 - t573 * t615;
t526 = t560 * t615 + t573 * t589;
t535 = -t560 * t594 + t573 * t599;
t630 = t573 * t594;
t536 = t560 * t599 + t630;
t629 = t590 * t596;
t559 = -t590 * t600 * t614 - t591 * t613 + t595 * t629;
t642 = Icges(5,5) * t536 + Icges(6,5) * t526 + Icges(5,6) * t535 - Icges(6,6) * t525 + t646 * t559;
t636 = pkin(4) * t599;
t541 = t575 * pkin(2) + pkin(10) * t561;
t542 = t577 * pkin(2) + pkin(10) * t562;
t620 = qJD(2) * t590;
t585 = t597 * t620;
t618 = t601 * t620;
t622 = t541 * t585 + t542 * t618;
t551 = qJD(3) * t562 + t585;
t621 = qJD(1) * (pkin(1) * t597 - pkin(9) * t627);
t586 = qJD(2) * t591 + qJD(1);
t508 = qJD(4) * t539 + t551;
t563 = qJD(3) * t573 + t586;
t527 = qJD(4) * t559 + t563;
t552 = qJD(3) * t561 - t618;
t502 = pkin(3) * t538 + pkin(11) * t537;
t503 = pkin(3) * t540 + pkin(11) * t539;
t610 = t551 * t502 - t503 * t552 + t622;
t509 = qJD(4) * t537 + t552;
t564 = pkin(2) * t629 + pkin(10) * t573;
t578 = qJD(1) * (pkin(1) * t601 + pkin(9) * t628);
t609 = t586 * t542 - t564 * t585 + t578;
t450 = pkin(4) * t632 + qJ(5) * t537 + t538 * t636;
t608 = qJD(5) * t559 + t508 * t450 + t610;
t607 = -t541 * t586 - t564 * t618 - t621;
t524 = pkin(3) * t560 + pkin(11) * t559;
t606 = t563 * t503 - t524 * t551 + t609;
t451 = pkin(4) * t631 + qJ(5) * t539 + t540 * t636;
t605 = qJD(5) * t537 + t527 * t451 + t606;
t604 = -t502 * t563 + t552 * t524 + t607;
t477 = pkin(4) * t630 + qJ(5) * t559 + t560 * t636;
t603 = qJD(5) * t539 + t509 * t477 + t604;
t598 = cos(qJ(6));
t593 = sin(qJ(6));
t583 = rSges(2,1) * t601 - rSges(2,2) * t597;
t582 = rSges(2,1) * t597 + rSges(2,2) * t601;
t571 = rSges(3,3) * t591 + (rSges(3,1) * t596 + rSges(3,2) * t600) * t590;
t570 = Icges(3,5) * t591 + (Icges(3,1) * t596 + Icges(3,4) * t600) * t590;
t569 = Icges(3,6) * t591 + (Icges(3,4) * t596 + Icges(3,2) * t600) * t590;
t568 = Icges(3,3) * t591 + (Icges(3,5) * t596 + Icges(3,6) * t600) * t590;
t550 = rSges(3,1) * t577 + rSges(3,2) * t576 + rSges(3,3) * t628;
t549 = rSges(3,1) * t575 + rSges(3,2) * t574 - rSges(3,3) * t627;
t548 = Icges(3,1) * t577 + Icges(3,4) * t576 + Icges(3,5) * t628;
t547 = Icges(3,1) * t575 + Icges(3,4) * t574 - Icges(3,5) * t627;
t546 = Icges(3,4) * t577 + Icges(3,2) * t576 + Icges(3,6) * t628;
t545 = Icges(3,4) * t575 + Icges(3,2) * t574 - Icges(3,6) * t627;
t523 = rSges(4,1) * t560 - rSges(4,2) * t559 + rSges(4,3) * t573;
t522 = Icges(4,1) * t560 - Icges(4,4) * t559 + Icges(4,5) * t573;
t521 = Icges(4,4) * t560 - Icges(4,2) * t559 + Icges(4,6) * t573;
t520 = Icges(4,5) * t560 - Icges(4,6) * t559 + Icges(4,3) * t573;
t515 = t550 * t586 - t571 * t585 + t578;
t514 = -t549 * t586 - t571 * t618 - t621;
t507 = t526 * t598 + t559 * t593;
t506 = -t526 * t593 + t559 * t598;
t505 = (t549 * t597 + t550 * t601) * t620;
t501 = qJD(6) * t525 + t527;
t500 = pkin(5) * t526 + pkin(12) * t525;
t498 = rSges(4,1) * t540 - rSges(4,2) * t539 + rSges(4,3) * t562;
t497 = rSges(4,1) * t538 - rSges(4,2) * t537 + rSges(4,3) * t561;
t496 = Icges(4,1) * t540 - Icges(4,4) * t539 + Icges(4,5) * t562;
t495 = Icges(4,1) * t538 - Icges(4,4) * t537 + Icges(4,5) * t561;
t494 = Icges(4,4) * t540 - Icges(4,2) * t539 + Icges(4,6) * t562;
t493 = Icges(4,4) * t538 - Icges(4,2) * t537 + Icges(4,6) * t561;
t492 = Icges(4,5) * t540 - Icges(4,6) * t539 + Icges(4,3) * t562;
t491 = Icges(4,5) * t538 - Icges(4,6) * t537 + Icges(4,3) * t561;
t490 = rSges(5,1) * t536 + rSges(5,2) * t535 + rSges(5,3) * t559;
t489 = Icges(5,1) * t536 + Icges(5,4) * t535 + Icges(5,5) * t559;
t488 = Icges(5,4) * t536 + Icges(5,2) * t535 + Icges(5,6) * t559;
t485 = t513 * t598 + t539 * t593;
t484 = -t513 * t593 + t539 * t598;
t483 = t511 * t598 + t537 * t593;
t482 = -t511 * t593 + t537 * t598;
t481 = rSges(6,1) * t526 - rSges(6,2) * t525 + rSges(6,3) * t559;
t480 = Icges(6,1) * t526 - Icges(6,4) * t525 + Icges(6,5) * t559;
t479 = Icges(6,4) * t526 - Icges(6,2) * t525 + Icges(6,6) * t559;
t476 = qJD(6) * t510 + t509;
t475 = qJD(6) * t512 + t508;
t474 = pkin(5) * t513 + pkin(12) * t512;
t473 = pkin(5) * t511 + pkin(12) * t510;
t472 = rSges(5,1) * t519 + rSges(5,2) * t518 + rSges(5,3) * t539;
t471 = rSges(5,1) * t517 + rSges(5,2) * t516 + rSges(5,3) * t537;
t470 = Icges(5,1) * t519 + Icges(5,4) * t518 + Icges(5,5) * t539;
t469 = Icges(5,1) * t517 + Icges(5,4) * t516 + Icges(5,5) * t537;
t468 = Icges(5,4) * t519 + Icges(5,2) * t518 + Icges(5,6) * t539;
t467 = Icges(5,4) * t517 + Icges(5,2) * t516 + Icges(5,6) * t537;
t464 = rSges(6,1) * t513 - rSges(6,2) * t512 + rSges(6,3) * t539;
t463 = rSges(6,1) * t511 - rSges(6,2) * t510 + rSges(6,3) * t537;
t462 = Icges(6,1) * t513 - Icges(6,4) * t512 + Icges(6,5) * t539;
t461 = Icges(6,1) * t511 - Icges(6,4) * t510 + Icges(6,5) * t537;
t460 = Icges(6,4) * t513 - Icges(6,2) * t512 + Icges(6,6) * t539;
t459 = Icges(6,4) * t511 - Icges(6,2) * t510 + Icges(6,6) * t537;
t455 = rSges(7,1) * t507 + rSges(7,2) * t506 + rSges(7,3) * t525;
t454 = Icges(7,1) * t507 + Icges(7,4) * t506 + Icges(7,5) * t525;
t453 = Icges(7,4) * t507 + Icges(7,2) * t506 + Icges(7,6) * t525;
t452 = Icges(7,5) * t507 + Icges(7,6) * t506 + Icges(7,3) * t525;
t447 = t498 * t563 - t523 * t551 + t609;
t446 = -t497 * t563 + t523 * t552 + t607;
t445 = rSges(7,1) * t485 + rSges(7,2) * t484 + rSges(7,3) * t512;
t444 = rSges(7,1) * t483 + rSges(7,2) * t482 + rSges(7,3) * t510;
t443 = Icges(7,1) * t485 + Icges(7,4) * t484 + Icges(7,5) * t512;
t442 = Icges(7,1) * t483 + Icges(7,4) * t482 + Icges(7,5) * t510;
t441 = Icges(7,4) * t485 + Icges(7,2) * t484 + Icges(7,6) * t512;
t440 = Icges(7,4) * t483 + Icges(7,2) * t482 + Icges(7,6) * t510;
t439 = Icges(7,5) * t485 + Icges(7,6) * t484 + Icges(7,3) * t512;
t438 = Icges(7,5) * t483 + Icges(7,6) * t482 + Icges(7,3) * t510;
t437 = t497 * t551 - t498 * t552 + t622;
t436 = t472 * t527 - t490 * t508 + t606;
t435 = -t471 * t527 + t490 * t509 + t604;
t434 = t471 * t508 - t472 * t509 + t610;
t433 = t464 * t527 + (-t477 - t481) * t508 + t605;
t432 = t481 * t509 + (-t450 - t463) * t527 + t603;
t431 = t463 * t508 + (-t451 - t464) * t509 + t608;
t430 = t445 * t501 - t455 * t475 + t474 * t527 + (-t477 - t500) * t508 + t605;
t429 = -t444 * t501 + t455 * t476 + t500 * t509 + (-t450 - t473) * t527 + t603;
t428 = t444 * t475 - t445 * t476 + t473 * t508 + (-t451 - t474) * t509 + t608;
t1 = m(7) * (t428 ^ 2 + t429 ^ 2 + t430 ^ 2) / 0.2e1 + m(5) * (t434 ^ 2 + t435 ^ 2 + t436 ^ 2) / 0.2e1 + m(6) * (t431 ^ 2 + t432 ^ 2 + t433 ^ 2) / 0.2e1 + t475 * ((t512 * t439 + t484 * t441 + t485 * t443) * t475 + (t438 * t512 + t440 * t484 + t442 * t485) * t476 + (t452 * t512 + t453 * t484 + t454 * t485) * t501) / 0.2e1 + t476 * ((t439 * t510 + t441 * t482 + t443 * t483) * t475 + (t510 * t438 + t482 * t440 + t483 * t442) * t476 + (t452 * t510 + t453 * t482 + t454 * t483) * t501) / 0.2e1 + t501 * ((t439 * t525 + t441 * t506 + t443 * t507) * t475 + (t438 * t525 + t440 * t506 + t442 * t507) * t476 + (t525 * t452 + t506 * t453 + t507 * t454) * t501) / 0.2e1 + t552 * ((t492 * t561 - t494 * t537 + t496 * t538) * t551 + (t561 * t491 - t537 * t493 + t538 * t495) * t552 + (t520 * t561 - t521 * t537 + t522 * t538) * t563) / 0.2e1 + t563 * ((t492 * t573 - t494 * t559 + t496 * t560) * t551 + (t491 * t573 - t493 * t559 + t495 * t560) * t552 + (t573 * t520 - t559 * t521 + t560 * t522) * t563) / 0.2e1 + t551 * ((t562 * t492 - t539 * t494 + t540 * t496) * t551 + (t491 * t562 - t493 * t539 + t495 * t540) * t552 + (t520 * t562 - t521 * t539 + t522 * t540) * t563) / 0.2e1 + t586 * ((t591 * t568 + (t569 * t600 + t570 * t596) * t590) * t586 + (((t546 * t600 + t548 * t596) * t597 - (t545 * t600 + t547 * t596) * t601) * t590 - t612 * t591) * t620) / 0.2e1 + m(4) * (t437 ^ 2 + t446 ^ 2 + t447 ^ 2) / 0.2e1 + m(3) * (t505 ^ 2 + t514 ^ 2 + t515 ^ 2) / 0.2e1 - ((-t568 * t627 + t569 * t574 + t570 * t575) * t586 + ((t546 * t574 + t548 * t575) * t597 + (-t574 * t545 - t575 * t547 + t645) * t601) * t620) * t618 / 0.2e1 + ((t568 * t628 + t569 * t576 + t570 * t577) * t586 + (-(t545 * t576 + t547 * t577) * t601 + (t576 * t546 + t577 * t548 - t645) * t597) * t620) * t585 / 0.2e1 + ((-t479 * t512 + t480 * t513 + t488 * t518 + t489 * t519 + t539 * t642) * t527 + (-t459 * t512 + t461 * t513 + t467 * t518 + t469 * t519 + t539 * t644) * t509 + (-t512 * t460 + t513 * t462 + t518 * t468 + t519 * t470 + t643 * t539) * t508) * t508 / 0.2e1 + ((-t479 * t510 + t480 * t511 + t488 * t516 + t489 * t517 + t537 * t642) * t527 + (-t510 * t459 + t511 * t461 + t516 * t467 + t517 * t469 + t644 * t537) * t509 + (-t460 * t510 + t462 * t511 + t468 * t516 + t470 * t517 + t537 * t643) * t508) * t509 / 0.2e1 + ((-t525 * t479 + t526 * t480 + t535 * t488 + t536 * t489 + t642 * t559) * t527 + (-t459 * t525 + t461 * t526 + t467 * t535 + t469 * t536 + t559 * t644) * t509 + (-t460 * t525 + t462 * t526 + t468 * t535 + t470 * t536 + t559 * t643) * t508) * t527 / 0.2e1 + (m(2) * (t582 ^ 2 + t583 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
