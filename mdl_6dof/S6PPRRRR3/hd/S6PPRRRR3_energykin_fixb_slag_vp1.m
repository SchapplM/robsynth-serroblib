% Calculate kinetic energy for
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPRRRR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_energykin_fixb_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRR3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRRRR3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:07:16
% EndTime: 2019-03-08 19:07:18
% DurationCPUTime: 2.07s
% Computational Cost: add. (9098->324), mult. (25735->516), div. (0->0), fcn. (34043->18), ass. (0->154)
t630 = cos(qJ(4));
t629 = cos(qJ(5));
t628 = cos(pkin(8));
t627 = sin(pkin(8));
t594 = sin(pkin(13));
t600 = cos(pkin(6));
t626 = t594 * t600;
t595 = sin(pkin(7));
t596 = sin(pkin(6));
t625 = t595 * t596;
t624 = t595 * t600;
t599 = cos(pkin(7));
t623 = t596 * t599;
t597 = cos(pkin(14));
t622 = t597 * t599;
t598 = cos(pkin(13));
t621 = t598 * t600;
t593 = sin(pkin(14));
t586 = t593 * t621 + t594 * t597;
t604 = sin(qJ(3));
t606 = cos(qJ(3));
t585 = -t593 * t594 + t597 * t621;
t614 = t585 * t599 - t598 * t625;
t568 = -t586 * t604 + t606 * t614;
t580 = -t585 * t595 - t598 * t623;
t556 = -t568 * t627 + t580 * t628;
t576 = qJD(3) * t580;
t546 = qJD(4) * t556 + t576;
t588 = -t593 * t626 + t597 * t598;
t587 = -t593 * t598 - t597 * t626;
t613 = t587 * t599 + t594 * t625;
t570 = -t588 * t604 + t606 * t613;
t581 = -t587 * t595 + t594 * t623;
t557 = -t570 * t627 + t581 * t628;
t577 = qJD(3) * t581;
t547 = qJD(4) * t557 + t577;
t578 = t606 * t624 + (-t593 * t604 + t606 * t622) * t596;
t584 = -t597 * t625 + t599 * t600;
t572 = -t578 * t627 + t584 * t628;
t583 = qJD(3) * t584;
t563 = qJD(4) * t572 + t583;
t620 = qJD(2) * t596;
t590 = qJD(2) * t600 + qJD(1);
t569 = t586 * t606 + t604 * t614;
t603 = sin(qJ(4));
t616 = t630 * t627;
t617 = t628 * t630;
t530 = -t568 * t617 + t569 * t603 - t580 * t616;
t512 = qJD(5) * t530 + t546;
t571 = t588 * t606 + t604 * t613;
t532 = -t570 * t617 + t571 * t603 - t581 * t616;
t513 = qJD(5) * t532 + t547;
t579 = t604 * t624 + (t593 * t606 + t604 * t622) * t596;
t554 = -t578 * t617 + t579 * t603 - t584 * t616;
t526 = qJD(5) * t554 + t563;
t618 = t598 * t620;
t536 = t569 * pkin(3) + pkin(10) * t556;
t558 = t579 * pkin(3) + pkin(10) * t572;
t589 = t594 * t620;
t615 = -t536 * t583 + t558 * t576 + t589;
t537 = t571 * pkin(3) + pkin(10) * t557;
t612 = t536 * t577 - t537 * t576 + t590;
t611 = t537 * t583 - t558 * t577 - t618;
t531 = t569 * t630 + (t568 * t628 + t580 * t627) * t603;
t507 = pkin(4) * t531 + pkin(11) * t530;
t555 = t579 * t630 + (t578 * t628 + t584 * t627) * t603;
t524 = pkin(4) * t555 + pkin(11) * t554;
t610 = -t507 * t563 + t546 * t524 + t615;
t533 = t571 * t630 + (t570 * t628 + t581 * t627) * t603;
t508 = pkin(4) * t533 + pkin(11) * t532;
t609 = t547 * t507 - t508 * t546 + t612;
t608 = t563 * t508 - t524 * t547 + t611;
t605 = cos(qJ(6));
t602 = sin(qJ(5));
t601 = sin(qJ(6));
t562 = rSges(4,1) * t579 + rSges(4,2) * t578 + rSges(4,3) * t584;
t561 = Icges(4,1) * t579 + Icges(4,4) * t578 + Icges(4,5) * t584;
t560 = Icges(4,4) * t579 + Icges(4,2) * t578 + Icges(4,6) * t584;
t559 = Icges(4,5) * t579 + Icges(4,6) * t578 + Icges(4,3) * t584;
t545 = rSges(4,1) * t571 + rSges(4,2) * t570 + rSges(4,3) * t581;
t544 = rSges(4,1) * t569 + rSges(4,2) * t568 + rSges(4,3) * t580;
t543 = Icges(4,1) * t571 + Icges(4,4) * t570 + Icges(4,5) * t581;
t542 = Icges(4,1) * t569 + Icges(4,4) * t568 + Icges(4,5) * t580;
t541 = Icges(4,4) * t571 + Icges(4,2) * t570 + Icges(4,6) * t581;
t540 = Icges(4,4) * t569 + Icges(4,2) * t568 + Icges(4,6) * t580;
t539 = Icges(4,5) * t571 + Icges(4,6) * t570 + Icges(4,3) * t581;
t538 = Icges(4,5) * t569 + Icges(4,6) * t568 + Icges(4,3) * t580;
t535 = t555 * t629 + t572 * t602;
t534 = t555 * t602 - t572 * t629;
t523 = rSges(5,1) * t555 - rSges(5,2) * t554 + rSges(5,3) * t572;
t522 = Icges(5,1) * t555 - Icges(5,4) * t554 + Icges(5,5) * t572;
t521 = Icges(5,4) * t555 - Icges(5,2) * t554 + Icges(5,6) * t572;
t520 = Icges(5,5) * t555 - Icges(5,6) * t554 + Icges(5,3) * t572;
t519 = t533 * t629 + t557 * t602;
t518 = t533 * t602 - t557 * t629;
t517 = t531 * t629 + t556 * t602;
t516 = t531 * t602 - t556 * t629;
t515 = t535 * t605 + t554 * t601;
t514 = -t535 * t601 + t554 * t605;
t511 = -t618 + (t545 * t584 - t562 * t581) * qJD(3);
t510 = t589 + (-t544 * t584 + t562 * t580) * qJD(3);
t509 = pkin(5) * t535 + pkin(12) * t534;
t505 = qJD(6) * t534 + t526;
t504 = (t544 * t581 - t545 * t580) * qJD(3) + t590;
t502 = rSges(5,1) * t533 - rSges(5,2) * t532 + rSges(5,3) * t557;
t501 = rSges(5,1) * t531 - rSges(5,2) * t530 + rSges(5,3) * t556;
t500 = Icges(5,1) * t533 - Icges(5,4) * t532 + Icges(5,5) * t557;
t499 = Icges(5,1) * t531 - Icges(5,4) * t530 + Icges(5,5) * t556;
t498 = Icges(5,4) * t533 - Icges(5,2) * t532 + Icges(5,6) * t557;
t497 = Icges(5,4) * t531 - Icges(5,2) * t530 + Icges(5,6) * t556;
t496 = Icges(5,5) * t533 - Icges(5,6) * t532 + Icges(5,3) * t557;
t495 = Icges(5,5) * t531 - Icges(5,6) * t530 + Icges(5,3) * t556;
t494 = rSges(6,1) * t535 - rSges(6,2) * t534 + rSges(6,3) * t554;
t493 = Icges(6,1) * t535 - Icges(6,4) * t534 + Icges(6,5) * t554;
t492 = Icges(6,4) * t535 - Icges(6,2) * t534 + Icges(6,6) * t554;
t491 = Icges(6,5) * t535 - Icges(6,6) * t534 + Icges(6,3) * t554;
t490 = t519 * t605 + t532 * t601;
t489 = -t519 * t601 + t532 * t605;
t488 = t517 * t605 + t530 * t601;
t487 = -t517 * t601 + t530 * t605;
t485 = pkin(5) * t519 + pkin(12) * t518;
t484 = pkin(5) * t517 + pkin(12) * t516;
t483 = qJD(6) * t518 + t513;
t482 = qJD(6) * t516 + t512;
t481 = rSges(6,1) * t519 - rSges(6,2) * t518 + rSges(6,3) * t532;
t480 = rSges(6,1) * t517 - rSges(6,2) * t516 + rSges(6,3) * t530;
t479 = Icges(6,1) * t519 - Icges(6,4) * t518 + Icges(6,5) * t532;
t478 = Icges(6,1) * t517 - Icges(6,4) * t516 + Icges(6,5) * t530;
t477 = Icges(6,4) * t519 - Icges(6,2) * t518 + Icges(6,6) * t532;
t476 = Icges(6,4) * t517 - Icges(6,2) * t516 + Icges(6,6) * t530;
t475 = Icges(6,5) * t519 - Icges(6,6) * t518 + Icges(6,3) * t532;
t474 = Icges(6,5) * t517 - Icges(6,6) * t516 + Icges(6,3) * t530;
t473 = rSges(7,1) * t515 + rSges(7,2) * t514 + rSges(7,3) * t534;
t472 = Icges(7,1) * t515 + Icges(7,4) * t514 + Icges(7,5) * t534;
t471 = Icges(7,4) * t515 + Icges(7,2) * t514 + Icges(7,6) * t534;
t470 = Icges(7,5) * t515 + Icges(7,6) * t514 + Icges(7,3) * t534;
t469 = rSges(7,1) * t490 + rSges(7,2) * t489 + rSges(7,3) * t518;
t468 = rSges(7,1) * t488 + rSges(7,2) * t487 + rSges(7,3) * t516;
t467 = Icges(7,1) * t490 + Icges(7,4) * t489 + Icges(7,5) * t518;
t466 = Icges(7,1) * t488 + Icges(7,4) * t487 + Icges(7,5) * t516;
t465 = Icges(7,4) * t490 + Icges(7,2) * t489 + Icges(7,6) * t518;
t464 = Icges(7,4) * t488 + Icges(7,2) * t487 + Icges(7,6) * t516;
t463 = Icges(7,5) * t490 + Icges(7,6) * t489 + Icges(7,3) * t518;
t462 = Icges(7,5) * t488 + Icges(7,6) * t487 + Icges(7,3) * t516;
t461 = t502 * t563 - t523 * t547 + t611;
t460 = -t501 * t563 + t523 * t546 + t615;
t459 = t501 * t547 - t502 * t546 + t612;
t458 = t481 * t526 - t494 * t513 + t608;
t457 = -t480 * t526 + t494 * t512 + t610;
t456 = t480 * t513 - t481 * t512 + t609;
t455 = t469 * t505 - t473 * t483 + t485 * t526 - t509 * t513 + t608;
t454 = -t468 * t505 + t473 * t482 - t484 * t526 + t509 * t512 + t610;
t453 = t468 * t483 - t469 * t482 + t484 * t513 - t485 * t512 + t609;
t1 = m(7) * (t453 ^ 2 + t454 ^ 2 + t455 ^ 2) / 0.2e1 + m(6) * (t456 ^ 2 + t457 ^ 2 + t458 ^ 2) / 0.2e1 + m(5) * (t459 ^ 2 + t460 ^ 2 + t461 ^ 2) / 0.2e1 + m(4) * (t504 ^ 2 + t510 ^ 2 + t511 ^ 2) / 0.2e1 + m(3) * (t590 ^ 2 + (t594 ^ 2 + t598 ^ 2) * qJD(2) ^ 2 * t596 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + t505 * ((t463 * t534 + t465 * t514 + t467 * t515) * t483 + (t462 * t534 + t464 * t514 + t466 * t515) * t482 + (t470 * t534 + t471 * t514 + t472 * t515) * t505) / 0.2e1 + t483 * ((t463 * t518 + t465 * t489 + t467 * t490) * t483 + (t462 * t518 + t464 * t489 + t466 * t490) * t482 + (t470 * t518 + t471 * t489 + t472 * t490) * t505) / 0.2e1 + t482 * ((t463 * t516 + t465 * t487 + t467 * t488) * t483 + (t462 * t516 + t464 * t487 + t466 * t488) * t482 + (t470 * t516 + t471 * t487 + t472 * t488) * t505) / 0.2e1 + t513 * ((t475 * t532 - t477 * t518 + t479 * t519) * t513 + (t474 * t532 - t476 * t518 + t478 * t519) * t512 + (t491 * t532 - t492 * t518 + t493 * t519) * t526) / 0.2e1 + t512 * ((t475 * t530 - t477 * t516 + t479 * t517) * t513 + (t474 * t530 - t476 * t516 + t478 * t517) * t512 + (t491 * t530 - t492 * t516 + t493 * t517) * t526) / 0.2e1 + t526 * ((t475 * t554 - t477 * t534 + t479 * t535) * t513 + (t474 * t554 - t476 * t534 + t478 * t535) * t512 + (t491 * t554 - t492 * t534 + t493 * t535) * t526) / 0.2e1 + t546 * ((t496 * t556 - t498 * t530 + t500 * t531) * t547 + (t495 * t556 - t497 * t530 + t499 * t531) * t546 + (t520 * t556 - t521 * t530 + t522 * t531) * t563) / 0.2e1 + t563 * ((t496 * t572 - t498 * t554 + t500 * t555) * t547 + (t495 * t572 - t497 * t554 + t499 * t555) * t546 + (t520 * t572 - t521 * t554 + t522 * t555) * t563) / 0.2e1 + t547 * ((t496 * t557 - t498 * t532 + t500 * t533) * t547 + (t495 * t557 - t497 * t532 + t499 * t533) * t546 + (t520 * t557 - t521 * t532 + t522 * t533) * t563) / 0.2e1 + (t580 * ((t539 * t580 + t541 * t568 + t543 * t569) * t581 + (t538 * t580 + t540 * t568 + t542 * t569) * t580 + (t559 * t580 + t560 * t568 + t561 * t569) * t584) + t584 * ((t584 * t539 + t578 * t541 + t579 * t543) * t581 + (t584 * t538 + t578 * t540 + t579 * t542) * t580 + (t584 * t559 + t578 * t560 + t579 * t561) * t584) + t581 * ((t539 * t581 + t541 * t570 + t543 * t571) * t581 + (t538 * t581 + t540 * t570 + t542 * t571) * t580 + (t559 * t581 + t560 * t570 + t561 * t571) * t584)) * qJD(3) ^ 2 / 0.2e1;
T  = t1;
