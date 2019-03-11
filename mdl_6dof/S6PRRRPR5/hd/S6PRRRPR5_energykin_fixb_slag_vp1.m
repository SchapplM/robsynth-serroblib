% Calculate kinetic energy for
% S6PRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_energykin_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPR5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:22:59
% EndTime: 2019-03-08 23:23:02
% DurationCPUTime: 3.14s
% Computational Cost: add. (5154->380), mult. (12753->583), div. (0->0), fcn. (16326->16), ass. (0->171)
t635 = Icges(5,3) + Icges(6,3);
t584 = sin(pkin(12));
t586 = cos(pkin(12));
t592 = sin(qJ(2));
t587 = cos(pkin(6));
t595 = cos(qJ(2));
t615 = t587 * t595;
t573 = -t584 * t592 + t586 * t615;
t616 = t587 * t592;
t574 = t584 * t595 + t586 * t616;
t591 = sin(qJ(3));
t585 = sin(pkin(6));
t623 = sin(pkin(7));
t608 = t585 * t623;
t624 = cos(pkin(7));
t627 = cos(qJ(3));
t534 = t574 * t627 + (t573 * t624 - t586 * t608) * t591;
t609 = t585 * t624;
t560 = -t573 * t623 - t586 * t609;
t612 = qJ(4) + pkin(13);
t583 = sin(t612);
t607 = cos(t612);
t509 = t534 * t583 - t560 * t607;
t510 = t534 * t607 + t560 * t583;
t590 = sin(qJ(4));
t594 = cos(qJ(4));
t513 = -t534 * t590 + t560 * t594;
t622 = t560 * t590;
t514 = t534 * t594 + t622;
t605 = t627 * t623;
t604 = t585 * t605;
t606 = t624 * t627;
t533 = -t573 * t606 + t574 * t591 + t586 * t604;
t634 = Icges(5,5) * t514 + Icges(6,5) * t510 + Icges(5,6) * t513 - Icges(6,6) * t509 + t635 * t533;
t575 = -t584 * t615 - t586 * t592;
t576 = -t584 * t616 + t586 * t595;
t536 = t576 * t627 + (t575 * t624 + t584 * t608) * t591;
t561 = -t575 * t623 + t584 * t609;
t511 = t536 * t583 - t561 * t607;
t512 = t536 * t607 + t561 * t583;
t515 = -t536 * t590 + t561 * t594;
t621 = t561 * t590;
t516 = t536 * t594 + t621;
t535 = -t575 * t606 + t576 * t591 - t584 * t604;
t633 = Icges(5,5) * t516 + Icges(6,5) * t512 + Icges(5,6) * t515 - Icges(6,6) * t511 + t635 * t535;
t559 = t587 * t623 * t591 + (t591 * t595 * t624 + t592 * t627) * t585;
t572 = t587 * t624 - t595 * t608;
t524 = t559 * t583 - t572 * t607;
t525 = t559 * t607 + t572 * t583;
t538 = -t559 * t590 + t572 * t594;
t620 = t572 * t590;
t539 = t559 * t594 + t620;
t618 = t585 * t592;
t558 = -t585 * t595 * t606 - t587 * t605 + t591 * t618;
t632 = Icges(5,5) * t539 + Icges(6,5) * t525 + Icges(5,6) * t538 - Icges(6,6) * t524 + t635 * t558;
t631 = qJD(2) ^ 2;
t626 = pkin(4) * t594;
t619 = t584 * t585;
t617 = t586 * t585;
t614 = qJD(2) * t585;
t580 = t584 * t614;
t550 = qJD(3) * t561 + t580;
t582 = qJD(2) * t587;
t563 = qJD(3) * t572 + t582;
t505 = qJD(4) * t535 + t550;
t528 = qJD(4) * t558 + t563;
t611 = t586 * t614;
t540 = t574 * pkin(2) + pkin(9) * t560;
t541 = t576 * pkin(2) + pkin(9) * t561;
t610 = t540 * t580 + t541 * t611 + qJD(1);
t551 = qJD(3) * t560 - t611;
t562 = pkin(2) * t618 + pkin(9) * t572;
t603 = t541 * t582 - t562 * t580;
t506 = qJD(4) * t533 + t551;
t501 = pkin(3) * t534 + pkin(10) * t533;
t502 = pkin(3) * t536 + pkin(10) * t535;
t602 = t550 * t501 - t502 * t551 + t610;
t601 = (-t540 * t587 - t562 * t617) * qJD(2);
t449 = pkin(4) * t622 + qJ(5) * t533 + t534 * t626;
t600 = qJD(5) * t558 + t505 * t449 + t602;
t523 = pkin(3) * t559 + pkin(10) * t558;
t599 = t563 * t502 - t523 * t550 + t603;
t450 = pkin(4) * t621 + qJ(5) * t535 + t536 * t626;
t598 = qJD(5) * t533 + t528 * t450 + t599;
t597 = -t501 * t563 + t551 * t523 + t601;
t480 = pkin(4) * t620 + qJ(5) * t558 + t559 * t626;
t596 = qJD(5) * t535 + t506 * t480 + t597;
t593 = cos(qJ(6));
t589 = sin(qJ(6));
t570 = t587 * rSges(3,3) + (rSges(3,1) * t592 + rSges(3,2) * t595) * t585;
t569 = Icges(3,5) * t587 + (Icges(3,1) * t592 + Icges(3,4) * t595) * t585;
t568 = Icges(3,6) * t587 + (Icges(3,4) * t592 + Icges(3,2) * t595) * t585;
t567 = Icges(3,3) * t587 + (Icges(3,5) * t592 + Icges(3,6) * t595) * t585;
t549 = rSges(3,1) * t576 + rSges(3,2) * t575 + rSges(3,3) * t619;
t548 = rSges(3,1) * t574 + rSges(3,2) * t573 - rSges(3,3) * t617;
t547 = Icges(3,1) * t576 + Icges(3,4) * t575 + Icges(3,5) * t619;
t546 = Icges(3,1) * t574 + Icges(3,4) * t573 - Icges(3,5) * t617;
t545 = Icges(3,4) * t576 + Icges(3,2) * t575 + Icges(3,6) * t619;
t544 = Icges(3,4) * t574 + Icges(3,2) * t573 - Icges(3,6) * t617;
t543 = Icges(3,5) * t576 + Icges(3,6) * t575 + Icges(3,3) * t619;
t542 = Icges(3,5) * t574 + Icges(3,6) * t573 - Icges(3,3) * t617;
t522 = (-t548 * t587 - t570 * t617) * qJD(2);
t521 = (t549 * t587 - t570 * t619) * qJD(2);
t520 = rSges(4,1) * t559 - rSges(4,2) * t558 + rSges(4,3) * t572;
t519 = Icges(4,1) * t559 - Icges(4,4) * t558 + Icges(4,5) * t572;
t518 = Icges(4,4) * t559 - Icges(4,2) * t558 + Icges(4,6) * t572;
t517 = Icges(4,5) * t559 - Icges(4,6) * t558 + Icges(4,3) * t572;
t508 = t525 * t593 + t558 * t589;
t507 = -t525 * t589 + t558 * t593;
t504 = qJD(1) + (t548 * t584 + t549 * t586) * t614;
t500 = qJD(6) * t524 + t528;
t499 = pkin(5) * t525 + pkin(11) * t524;
t497 = rSges(5,1) * t539 + rSges(5,2) * t538 + rSges(5,3) * t558;
t496 = rSges(4,1) * t536 - rSges(4,2) * t535 + rSges(4,3) * t561;
t495 = rSges(4,1) * t534 - rSges(4,2) * t533 + rSges(4,3) * t560;
t494 = Icges(5,1) * t539 + Icges(5,4) * t538 + Icges(5,5) * t558;
t493 = Icges(5,4) * t539 + Icges(5,2) * t538 + Icges(5,6) * t558;
t491 = Icges(4,1) * t536 - Icges(4,4) * t535 + Icges(4,5) * t561;
t490 = Icges(4,1) * t534 - Icges(4,4) * t533 + Icges(4,5) * t560;
t489 = Icges(4,4) * t536 - Icges(4,2) * t535 + Icges(4,6) * t561;
t488 = Icges(4,4) * t534 - Icges(4,2) * t533 + Icges(4,6) * t560;
t487 = Icges(4,5) * t536 - Icges(4,6) * t535 + Icges(4,3) * t561;
t486 = Icges(4,5) * t534 - Icges(4,6) * t533 + Icges(4,3) * t560;
t485 = rSges(6,1) * t525 - rSges(6,2) * t524 + rSges(6,3) * t558;
t483 = Icges(6,1) * t525 - Icges(6,4) * t524 + Icges(6,5) * t558;
t482 = Icges(6,4) * t525 - Icges(6,2) * t524 + Icges(6,6) * t558;
t479 = t512 * t593 + t535 * t589;
t478 = -t512 * t589 + t535 * t593;
t477 = t510 * t593 + t533 * t589;
t476 = -t510 * t589 + t533 * t593;
t475 = qJD(6) * t509 + t506;
t474 = qJD(6) * t511 + t505;
t473 = pkin(5) * t512 + pkin(11) * t511;
t472 = pkin(5) * t510 + pkin(11) * t509;
t471 = rSges(5,1) * t516 + rSges(5,2) * t515 + rSges(5,3) * t535;
t470 = rSges(5,1) * t514 + rSges(5,2) * t513 + rSges(5,3) * t533;
t469 = Icges(5,1) * t516 + Icges(5,4) * t515 + Icges(5,5) * t535;
t468 = Icges(5,1) * t514 + Icges(5,4) * t513 + Icges(5,5) * t533;
t467 = Icges(5,4) * t516 + Icges(5,2) * t515 + Icges(5,6) * t535;
t466 = Icges(5,4) * t514 + Icges(5,2) * t513 + Icges(5,6) * t533;
t463 = rSges(6,1) * t512 - rSges(6,2) * t511 + rSges(6,3) * t535;
t462 = rSges(6,1) * t510 - rSges(6,2) * t509 + rSges(6,3) * t533;
t461 = Icges(6,1) * t512 - Icges(6,4) * t511 + Icges(6,5) * t535;
t460 = Icges(6,1) * t510 - Icges(6,4) * t509 + Icges(6,5) * t533;
t459 = Icges(6,4) * t512 - Icges(6,2) * t511 + Icges(6,6) * t535;
t458 = Icges(6,4) * t510 - Icges(6,2) * t509 + Icges(6,6) * t533;
t454 = rSges(7,1) * t508 + rSges(7,2) * t507 + rSges(7,3) * t524;
t453 = Icges(7,1) * t508 + Icges(7,4) * t507 + Icges(7,5) * t524;
t452 = Icges(7,4) * t508 + Icges(7,2) * t507 + Icges(7,6) * t524;
t451 = Icges(7,5) * t508 + Icges(7,6) * t507 + Icges(7,3) * t524;
t446 = -t495 * t563 + t520 * t551 + t601;
t445 = t496 * t563 - t520 * t550 + t603;
t444 = rSges(7,1) * t479 + rSges(7,2) * t478 + rSges(7,3) * t511;
t443 = rSges(7,1) * t477 + rSges(7,2) * t476 + rSges(7,3) * t509;
t442 = Icges(7,1) * t479 + Icges(7,4) * t478 + Icges(7,5) * t511;
t441 = Icges(7,1) * t477 + Icges(7,4) * t476 + Icges(7,5) * t509;
t440 = Icges(7,4) * t479 + Icges(7,2) * t478 + Icges(7,6) * t511;
t439 = Icges(7,4) * t477 + Icges(7,2) * t476 + Icges(7,6) * t509;
t438 = Icges(7,5) * t479 + Icges(7,6) * t478 + Icges(7,3) * t511;
t437 = Icges(7,5) * t477 + Icges(7,6) * t476 + Icges(7,3) * t509;
t436 = t495 * t550 - t496 * t551 + t610;
t435 = -t470 * t528 + t497 * t506 + t597;
t434 = t471 * t528 - t497 * t505 + t599;
t433 = t470 * t505 - t471 * t506 + t602;
t432 = t485 * t506 + (-t449 - t462) * t528 + t596;
t431 = t463 * t528 + (-t480 - t485) * t505 + t598;
t430 = t462 * t505 + (-t450 - t463) * t506 + t600;
t429 = -t443 * t500 + t454 * t475 + t499 * t506 + (-t449 - t472) * t528 + t596;
t428 = t444 * t500 - t454 * t474 + t473 * t528 + (-t480 - t499) * t505 + t598;
t427 = t443 * t474 - t444 * t475 + t472 * t505 + (-t450 - t473) * t506 + t600;
t1 = t563 * ((t487 * t572 - t489 * t558 + t491 * t559) * t550 + (t486 * t572 - t488 * t558 + t490 * t559) * t551 + (t572 * t517 - t558 * t518 + t559 * t519) * t563) / 0.2e1 + t550 * ((t561 * t487 - t535 * t489 + t536 * t491) * t550 + (t486 * t561 - t488 * t535 + t490 * t536) * t551 + (t517 * t561 - t518 * t535 + t519 * t536) * t563) / 0.2e1 + t551 * ((t487 * t560 - t533 * t489 + t534 * t491) * t550 + (t560 * t486 - t533 * t488 + t534 * t490) * t551 + (t517 * t560 - t518 * t533 + t519 * t534) * t563) / 0.2e1 + m(7) * (t427 ^ 2 + t428 ^ 2 + t429 ^ 2) / 0.2e1 + m(6) * (t430 ^ 2 + t431 ^ 2 + t432 ^ 2) / 0.2e1 + m(5) * (t433 ^ 2 + t434 ^ 2 + t435 ^ 2) / 0.2e1 + m(4) * (t436 ^ 2 + t445 ^ 2 + t446 ^ 2) / 0.2e1 + m(3) * (t504 ^ 2 + t521 ^ 2 + t522 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + t500 * ((t438 * t524 + t440 * t507 + t442 * t508) * t474 + (t437 * t524 + t439 * t507 + t441 * t508) * t475 + (t524 * t451 + t507 * t452 + t508 * t453) * t500) / 0.2e1 + t475 * ((t438 * t509 + t440 * t476 + t442 * t477) * t474 + (t509 * t437 + t476 * t439 + t477 * t441) * t475 + (t451 * t509 + t452 * t476 + t453 * t477) * t500) / 0.2e1 + t474 * ((t511 * t438 + t478 * t440 + t479 * t442) * t474 + (t437 * t511 + t439 * t478 + t441 * t479) * t475 + (t451 * t511 + t452 * t478 + t453 * t479) * t500) / 0.2e1 - t631 * ((-t543 * t617 + t545 * t573 + t547 * t574) * t619 - (-t542 * t617 + t544 * t573 + t546 * t574) * t617 + (-t567 * t617 + t568 * t573 + t569 * t574) * t587) * t617 / 0.2e1 + ((-t482 * t511 + t483 * t512 + t493 * t515 + t494 * t516 + t535 * t632) * t528 + (-t458 * t511 + t460 * t512 + t466 * t515 + t468 * t516 + t535 * t634) * t506 + (-t511 * t459 + t512 * t461 + t515 * t467 + t516 * t469 + t633 * t535) * t505) * t505 / 0.2e1 + ((-t482 * t509 + t483 * t510 + t493 * t513 + t494 * t514 + t533 * t632) * t528 + (-t509 * t458 + t510 * t460 + t513 * t466 + t514 * t468 + t634 * t533) * t506 + (-t459 * t509 + t461 * t510 + t467 * t513 + t469 * t514 + t533 * t633) * t505) * t506 / 0.2e1 + ((-t524 * t482 + t525 * t483 + t538 * t493 + t539 * t494 + t632 * t558) * t528 + (-t458 * t524 + t460 * t525 + t466 * t538 + t468 * t539 + t558 * t634) * t506 + (-t459 * t524 + t461 * t525 + t467 * t538 + t469 * t539 + t558 * t633) * t505) * t528 / 0.2e1 + (t587 * (t587 ^ 2 * t567 + (((t545 * t595 + t547 * t592) * t584 - (t544 * t595 + t546 * t592) * t586) * t585 + (-t542 * t586 + t543 * t584 + t568 * t595 + t569 * t592) * t587) * t585) + ((t543 * t619 + t545 * t575 + t547 * t576) * t619 - (t542 * t619 + t544 * t575 + t546 * t576) * t617 + (t567 * t619 + t568 * t575 + t569 * t576) * t587) * t619) * t631 / 0.2e1;
T  = t1;
