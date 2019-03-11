% Calculate kinetic energy for
% S6RRRPRR13
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
% Datum: 2019-03-09 20:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR13_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR13_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_energykin_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR13_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR13_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR13_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:51:55
% EndTime: 2019-03-09 19:51:58
% DurationCPUTime: 3.47s
% Computational Cost: add. (5106->387), mult. (12499->589), div. (0->0), fcn. (15964->16), ass. (0->176)
t646 = Icges(4,2) + Icges(5,3);
t591 = sin(pkin(6));
t593 = cos(pkin(6));
t600 = cos(qJ(2));
t601 = cos(qJ(1));
t624 = t600 * t601;
t597 = sin(qJ(2));
t598 = sin(qJ(1));
t627 = t597 * t598;
t574 = t593 * t624 - t627;
t625 = t598 * t600;
t626 = t597 * t601;
t575 = t593 * t626 + t625;
t576 = -t593 * t625 - t626;
t577 = -t593 * t627 + t624;
t628 = t591 * t601;
t629 = t591 * t598;
t611 = (Icges(3,5) * t575 + Icges(3,6) * t574 - Icges(3,3) * t628) * t601 - (Icges(3,5) * t577 + Icges(3,6) * t576 + Icges(3,3) * t629) * t598;
t645 = t591 * t611;
t596 = sin(qJ(3));
t634 = sin(pkin(7));
t616 = t591 * t634;
t635 = cos(pkin(7));
t637 = cos(qJ(3));
t538 = t575 * t637 + (t574 * t635 - t601 * t616) * t596;
t617 = t591 * t635;
t561 = -t574 * t634 - t601 * t617;
t590 = sin(pkin(13));
t592 = cos(pkin(13));
t516 = -t538 * t590 + t561 * t592;
t633 = t561 * t590;
t517 = t538 * t592 + t633;
t612 = t637 * t634;
t610 = t591 * t612;
t613 = t635 * t637;
t537 = -t574 * t613 + t575 * t596 + t601 * t610;
t644 = -Icges(4,4) * t538 + Icges(5,5) * t517 - Icges(4,6) * t561 + Icges(5,6) * t516 + t646 * t537;
t540 = t577 * t637 + (t576 * t635 + t598 * t616) * t596;
t562 = -t576 * t634 + t598 * t617;
t518 = -t540 * t590 + t562 * t592;
t632 = t562 * t590;
t519 = t540 * t592 + t632;
t539 = -t576 * t613 + t577 * t596 - t598 * t610;
t643 = -Icges(4,4) * t540 + Icges(5,5) * t519 - Icges(4,6) * t562 + Icges(5,6) * t518 + t646 * t539;
t560 = t593 * t634 * t596 + (t596 * t600 * t635 + t597 * t637) * t591;
t573 = t593 * t635 - t600 * t616;
t531 = -t560 * t590 + t573 * t592;
t631 = t573 * t590;
t532 = t560 * t592 + t631;
t630 = t591 * t597;
t559 = -t591 * t600 * t613 - t593 * t612 + t596 * t630;
t642 = -Icges(4,4) * t560 + Icges(5,5) * t532 - Icges(4,6) * t573 + Icges(5,6) * t531 + t646 * t559;
t636 = pkin(4) * t592;
t541 = t575 * pkin(2) + pkin(10) * t561;
t542 = t577 * pkin(2) + pkin(10) * t562;
t620 = qJD(2) * t591;
t585 = t598 * t620;
t618 = t601 * t620;
t622 = t541 * t585 + t542 * t618;
t551 = qJD(3) * t562 + t585;
t621 = qJD(1) * (pkin(1) * t598 - pkin(9) * t628);
t586 = qJD(2) * t593 + qJD(1);
t619 = pkin(13) + qJ(5);
t508 = qJD(5) * t539 + t551;
t563 = qJD(3) * t573 + t586;
t615 = cos(t619);
t502 = pkin(3) * t538 + qJ(4) * t537;
t614 = qJD(4) * t559 + t551 * t502 + t622;
t527 = qJD(5) * t559 + t563;
t552 = qJD(3) * t561 - t618;
t509 = qJD(5) * t537 + t552;
t564 = pkin(2) * t630 + pkin(10) * t573;
t578 = qJD(1) * (pkin(1) * t601 + pkin(9) * t629);
t609 = t586 * t542 - t564 * t585 + t578;
t503 = pkin(3) * t540 + qJ(4) * t539;
t608 = qJD(4) * t537 + t563 * t503 + t609;
t450 = pkin(4) * t633 + pkin(11) * t537 + t538 * t636;
t451 = pkin(4) * t632 + pkin(11) * t539 + t540 * t636;
t607 = t551 * t450 + (-t451 - t503) * t552 + t614;
t606 = -t541 * t586 - t564 * t618 - t621;
t524 = pkin(3) * t560 + qJ(4) * t559;
t605 = qJD(4) * t539 + t552 * t524 + t606;
t477 = pkin(4) * t631 + pkin(11) * t559 + t560 * t636;
t604 = t563 * t451 + (-t477 - t524) * t551 + t608;
t603 = t552 * t477 + (-t450 - t502) * t563 + t605;
t599 = cos(qJ(6));
t595 = sin(qJ(6));
t589 = sin(t619);
t583 = rSges(2,1) * t601 - rSges(2,2) * t598;
t582 = rSges(2,1) * t598 + rSges(2,2) * t601;
t571 = rSges(3,3) * t593 + (rSges(3,1) * t597 + rSges(3,2) * t600) * t591;
t570 = Icges(3,5) * t593 + (Icges(3,1) * t597 + Icges(3,4) * t600) * t591;
t569 = Icges(3,6) * t593 + (Icges(3,4) * t597 + Icges(3,2) * t600) * t591;
t568 = Icges(3,3) * t593 + (Icges(3,5) * t597 + Icges(3,6) * t600) * t591;
t550 = rSges(3,1) * t577 + rSges(3,2) * t576 + rSges(3,3) * t629;
t549 = rSges(3,1) * t575 + rSges(3,2) * t574 - rSges(3,3) * t628;
t548 = Icges(3,1) * t577 + Icges(3,4) * t576 + Icges(3,5) * t629;
t547 = Icges(3,1) * t575 + Icges(3,4) * t574 - Icges(3,5) * t628;
t546 = Icges(3,4) * t577 + Icges(3,2) * t576 + Icges(3,6) * t629;
t545 = Icges(3,4) * t575 + Icges(3,2) * t574 - Icges(3,6) * t628;
t526 = t560 * t615 + t573 * t589;
t525 = t560 * t589 - t573 * t615;
t523 = rSges(4,1) * t560 - rSges(4,2) * t559 + rSges(4,3) * t573;
t522 = Icges(4,1) * t560 - Icges(4,4) * t559 + Icges(4,5) * t573;
t520 = Icges(4,5) * t560 - Icges(4,6) * t559 + Icges(4,3) * t573;
t515 = t550 * t586 - t571 * t585 + t578;
t514 = -t549 * t586 - t571 * t618 - t621;
t513 = t540 * t615 + t562 * t589;
t512 = t540 * t589 - t562 * t615;
t511 = t538 * t615 + t561 * t589;
t510 = t538 * t589 - t561 * t615;
t507 = t526 * t599 + t559 * t595;
t506 = -t526 * t595 + t559 * t599;
t505 = (t549 * t598 + t550 * t601) * t620;
t501 = qJD(6) * t525 + t527;
t500 = pkin(5) * t526 + pkin(12) * t525;
t498 = rSges(4,1) * t540 - rSges(4,2) * t539 + rSges(4,3) * t562;
t497 = rSges(4,1) * t538 - rSges(4,2) * t537 + rSges(4,3) * t561;
t496 = Icges(4,1) * t540 - Icges(4,4) * t539 + Icges(4,5) * t562;
t495 = Icges(4,1) * t538 - Icges(4,4) * t537 + Icges(4,5) * t561;
t492 = Icges(4,5) * t540 - Icges(4,6) * t539 + Icges(4,3) * t562;
t491 = Icges(4,5) * t538 - Icges(4,6) * t537 + Icges(4,3) * t561;
t490 = rSges(5,1) * t532 + rSges(5,2) * t531 + rSges(5,3) * t559;
t489 = Icges(5,1) * t532 + Icges(5,4) * t531 + Icges(5,5) * t559;
t488 = Icges(5,4) * t532 + Icges(5,2) * t531 + Icges(5,6) * t559;
t485 = t513 * t599 + t539 * t595;
t484 = -t513 * t595 + t539 * t599;
t483 = t511 * t599 + t537 * t595;
t482 = -t511 * t595 + t537 * t599;
t481 = rSges(6,1) * t526 - rSges(6,2) * t525 + rSges(6,3) * t559;
t480 = Icges(6,1) * t526 - Icges(6,4) * t525 + Icges(6,5) * t559;
t479 = Icges(6,4) * t526 - Icges(6,2) * t525 + Icges(6,6) * t559;
t478 = Icges(6,5) * t526 - Icges(6,6) * t525 + Icges(6,3) * t559;
t476 = qJD(6) * t510 + t509;
t475 = qJD(6) * t512 + t508;
t474 = pkin(5) * t513 + pkin(12) * t512;
t473 = pkin(5) * t511 + pkin(12) * t510;
t471 = rSges(5,1) * t519 + rSges(5,2) * t518 + rSges(5,3) * t539;
t470 = rSges(5,1) * t517 + rSges(5,2) * t516 + rSges(5,3) * t537;
t469 = Icges(5,1) * t519 + Icges(5,4) * t518 + Icges(5,5) * t539;
t468 = Icges(5,1) * t517 + Icges(5,4) * t516 + Icges(5,5) * t537;
t467 = Icges(5,4) * t519 + Icges(5,2) * t518 + Icges(5,6) * t539;
t466 = Icges(5,4) * t517 + Icges(5,2) * t516 + Icges(5,6) * t537;
t463 = rSges(6,1) * t513 - rSges(6,2) * t512 + rSges(6,3) * t539;
t462 = rSges(6,1) * t511 - rSges(6,2) * t510 + rSges(6,3) * t537;
t461 = Icges(6,1) * t513 - Icges(6,4) * t512 + Icges(6,5) * t539;
t460 = Icges(6,1) * t511 - Icges(6,4) * t510 + Icges(6,5) * t537;
t459 = Icges(6,4) * t513 - Icges(6,2) * t512 + Icges(6,6) * t539;
t458 = Icges(6,4) * t511 - Icges(6,2) * t510 + Icges(6,6) * t537;
t457 = Icges(6,5) * t513 - Icges(6,6) * t512 + Icges(6,3) * t539;
t456 = Icges(6,5) * t511 - Icges(6,6) * t510 + Icges(6,3) * t537;
t455 = rSges(7,1) * t507 + rSges(7,2) * t506 + rSges(7,3) * t525;
t454 = Icges(7,1) * t507 + Icges(7,4) * t506 + Icges(7,5) * t525;
t453 = Icges(7,4) * t507 + Icges(7,2) * t506 + Icges(7,6) * t525;
t452 = Icges(7,5) * t507 + Icges(7,6) * t506 + Icges(7,3) * t525;
t447 = t498 * t563 - t523 * t551 + t609;
t446 = -t497 * t563 + t523 * t552 + t606;
t445 = rSges(7,1) * t485 + rSges(7,2) * t484 + rSges(7,3) * t512;
t444 = rSges(7,1) * t483 + rSges(7,2) * t482 + rSges(7,3) * t510;
t443 = Icges(7,1) * t485 + Icges(7,4) * t484 + Icges(7,5) * t512;
t442 = Icges(7,1) * t483 + Icges(7,4) * t482 + Icges(7,5) * t510;
t441 = Icges(7,4) * t485 + Icges(7,2) * t484 + Icges(7,6) * t512;
t440 = Icges(7,4) * t483 + Icges(7,2) * t482 + Icges(7,6) * t510;
t439 = Icges(7,5) * t485 + Icges(7,6) * t484 + Icges(7,3) * t512;
t438 = Icges(7,5) * t483 + Icges(7,6) * t482 + Icges(7,3) * t510;
t437 = t497 * t551 - t498 * t552 + t622;
t436 = t471 * t563 + (-t490 - t524) * t551 + t608;
t435 = t490 * t552 + (-t470 - t502) * t563 + t605;
t434 = t470 * t551 + (-t471 - t503) * t552 + t614;
t433 = t463 * t527 - t481 * t508 + t604;
t432 = -t462 * t527 + t481 * t509 + t603;
t431 = t462 * t508 - t463 * t509 + t607;
t430 = t445 * t501 - t455 * t475 + t474 * t527 - t500 * t508 + t604;
t429 = -t444 * t501 + t455 * t476 - t473 * t527 + t500 * t509 + t603;
t428 = t444 * t475 - t445 * t476 + t473 * t508 - t474 * t509 + t607;
t1 = m(4) * (t437 ^ 2 + t446 ^ 2 + t447 ^ 2) / 0.2e1 + t476 * ((t439 * t510 + t441 * t482 + t443 * t483) * t475 + (t510 * t438 + t482 * t440 + t483 * t442) * t476 + (t452 * t510 + t453 * t482 + t454 * t483) * t501) / 0.2e1 + t501 * ((t439 * t525 + t441 * t506 + t443 * t507) * t475 + (t438 * t525 + t440 * t506 + t442 * t507) * t476 + (t525 * t452 + t506 * t453 + t507 * t454) * t501) / 0.2e1 + t475 * ((t512 * t439 + t484 * t441 + t485 * t443) * t475 + (t438 * t512 + t440 * t484 + t442 * t485) * t476 + (t452 * t512 + t453 * t484 + t454 * t485) * t501) / 0.2e1 + t527 * ((t457 * t559 - t459 * t525 + t461 * t526) * t508 + (t456 * t559 - t458 * t525 + t460 * t526) * t509 + (t559 * t478 - t525 * t479 + t526 * t480) * t527) / 0.2e1 + t509 * ((t457 * t537 - t459 * t510 + t461 * t511) * t508 + (t537 * t456 - t510 * t458 + t511 * t460) * t509 + (t478 * t537 - t479 * t510 + t480 * t511) * t527) / 0.2e1 + t508 * ((t539 * t457 - t512 * t459 + t513 * t461) * t508 + (t456 * t539 - t458 * t512 + t460 * t513) * t509 + (t478 * t539 - t479 * t512 + t480 * t513) * t527) / 0.2e1 + t586 * ((t593 * t568 + (t569 * t600 + t570 * t597) * t591) * t586 + (((t546 * t600 + t548 * t597) * t598 - (t545 * t600 + t547 * t597) * t601) * t591 - t611 * t593) * t620) / 0.2e1 + m(3) * (t505 ^ 2 + t514 ^ 2 + t515 ^ 2) / 0.2e1 + m(7) * (t428 ^ 2 + t429 ^ 2 + t430 ^ 2) / 0.2e1 + m(6) * (t431 ^ 2 + t432 ^ 2 + t433 ^ 2) / 0.2e1 + m(5) * (t434 ^ 2 + t435 ^ 2 + t436 ^ 2) / 0.2e1 - ((-t568 * t628 + t569 * t574 + t570 * t575) * t586 + ((t546 * t574 + t548 * t575) * t598 + (-t574 * t545 - t575 * t547 + t645) * t601) * t620) * t618 / 0.2e1 + ((t568 * t629 + t569 * t576 + t570 * t577) * t586 + (-(t545 * t576 + t547 * t577) * t601 + (t576 * t546 + t577 * t548 - t645) * t598) * t620) * t585 / 0.2e1 + ((t488 * t518 + t489 * t519 + t520 * t562 + t522 * t540 + t539 * t642) * t563 + (t466 * t518 + t468 * t519 + t491 * t562 + t495 * t540 + t539 * t644) * t552 + (t518 * t467 + t519 * t469 + t562 * t492 + t540 * t496 + t643 * t539) * t551) * t551 / 0.2e1 + ((t488 * t516 + t489 * t517 + t520 * t561 + t522 * t538 + t537 * t642) * t563 + (t516 * t466 + t517 * t468 + t561 * t491 + t538 * t495 + t644 * t537) * t552 + (t467 * t516 + t469 * t517 + t492 * t561 + t496 * t538 + t537 * t643) * t551) * t552 / 0.2e1 + ((t531 * t488 + t532 * t489 + t573 * t520 + t560 * t522 + t642 * t559) * t563 + (t466 * t531 + t468 * t532 + t491 * t573 + t495 * t560 + t559 * t644) * t552 + (t467 * t531 + t469 * t532 + t492 * t573 + t496 * t560 + t559 * t643) * t551) * t563 / 0.2e1 + (Icges(2,3) + m(2) * (t582 ^ 2 + t583 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
