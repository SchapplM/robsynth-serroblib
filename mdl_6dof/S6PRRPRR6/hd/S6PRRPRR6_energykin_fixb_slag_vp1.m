% Calculate kinetic energy for
% S6PRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_energykin_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRR6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:22:50
% EndTime: 2019-03-08 22:22:53
% DurationCPUTime: 3.14s
% Computational Cost: add. (5037->380), mult. (12447->583), div. (0->0), fcn. (15930->16), ass. (0->171)
t635 = Icges(4,2) + Icges(5,3);
t585 = sin(pkin(12));
t588 = cos(pkin(12));
t593 = sin(qJ(2));
t589 = cos(pkin(6));
t595 = cos(qJ(2));
t616 = t589 * t595;
t573 = -t585 * t593 + t588 * t616;
t617 = t589 * t593;
t574 = t585 * t595 + t588 * t617;
t592 = sin(qJ(3));
t586 = sin(pkin(6));
t624 = sin(pkin(7));
t608 = t586 * t624;
t625 = cos(pkin(7));
t627 = cos(qJ(3));
t536 = t574 * t627 + (t573 * t625 - t588 * t608) * t592;
t609 = t586 * t625;
t560 = -t573 * t624 - t588 * t609;
t584 = sin(pkin(13));
t587 = cos(pkin(13));
t513 = -t536 * t584 + t560 * t587;
t623 = t560 * t584;
t514 = t536 * t587 + t623;
t605 = t627 * t624;
t603 = t586 * t605;
t606 = t625 * t627;
t535 = -t573 * t606 + t574 * t592 + t588 * t603;
t634 = -Icges(4,4) * t536 + Icges(5,5) * t514 - Icges(4,6) * t560 + Icges(5,6) * t513 + t635 * t535;
t575 = -t585 * t616 - t588 * t593;
t576 = -t585 * t617 + t588 * t595;
t538 = t576 * t627 + (t575 * t625 + t585 * t608) * t592;
t561 = -t575 * t624 + t585 * t609;
t515 = -t538 * t584 + t561 * t587;
t622 = t561 * t584;
t516 = t538 * t587 + t622;
t537 = -t575 * t606 + t576 * t592 - t585 * t603;
t633 = -Icges(4,4) * t538 + Icges(5,5) * t516 - Icges(4,6) * t561 + Icges(5,6) * t515 + t635 * t537;
t559 = t589 * t624 * t592 + (t592 * t595 * t625 + t593 * t627) * t586;
t572 = t589 * t625 - t595 * t608;
t533 = -t559 * t584 + t572 * t587;
t621 = t572 * t584;
t534 = t559 * t587 + t621;
t619 = t586 * t593;
t558 = -t586 * t595 * t606 - t589 * t605 + t592 * t619;
t632 = -Icges(4,4) * t559 + Icges(5,5) * t534 - Icges(4,6) * t572 + Icges(5,6) * t533 + t635 * t558;
t631 = qJD(2) ^ 2;
t626 = pkin(4) * t587;
t620 = t585 * t586;
t618 = t588 * t586;
t614 = qJD(2) * t586;
t580 = t585 * t614;
t550 = qJD(3) * t561 + t580;
t582 = qJD(2) * t589;
t563 = qJD(3) * t572 + t582;
t612 = pkin(13) + qJ(5);
t505 = qJD(5) * t537 + t550;
t528 = qJD(5) * t558 + t563;
t611 = t588 * t614;
t540 = t574 * pkin(2) + pkin(9) * t560;
t541 = t576 * pkin(2) + pkin(9) * t561;
t610 = t540 * t580 + t541 * t611 + qJD(1);
t607 = cos(t612);
t551 = qJD(3) * t560 - t611;
t501 = pkin(3) * t536 + qJ(4) * t535;
t604 = qJD(4) * t558 + t550 * t501 + t610;
t562 = pkin(2) * t619 + pkin(9) * t572;
t602 = t541 * t582 - t562 * t580;
t506 = qJD(5) * t535 + t551;
t502 = pkin(3) * t538 + qJ(4) * t537;
t601 = qJD(4) * t535 + t563 * t502 + t602;
t600 = (-t540 * t589 - t562 * t618) * qJD(2);
t449 = pkin(4) * t623 + pkin(10) * t535 + t536 * t626;
t450 = pkin(4) * t622 + pkin(10) * t537 + t538 * t626;
t599 = t550 * t449 + (-t450 - t502) * t551 + t604;
t523 = pkin(3) * t559 + qJ(4) * t558;
t598 = qJD(4) * t537 + t551 * t523 + t600;
t480 = pkin(4) * t621 + pkin(10) * t558 + t559 * t626;
t597 = t563 * t450 + (-t480 - t523) * t550 + t601;
t596 = t551 * t480 + (-t449 - t501) * t563 + t598;
t594 = cos(qJ(6));
t591 = sin(qJ(6));
t583 = sin(t612);
t570 = t589 * rSges(3,3) + (rSges(3,1) * t593 + rSges(3,2) * t595) * t586;
t569 = Icges(3,5) * t589 + (Icges(3,1) * t593 + Icges(3,4) * t595) * t586;
t568 = Icges(3,6) * t589 + (Icges(3,4) * t593 + Icges(3,2) * t595) * t586;
t567 = Icges(3,3) * t589 + (Icges(3,5) * t593 + Icges(3,6) * t595) * t586;
t549 = rSges(3,1) * t576 + rSges(3,2) * t575 + rSges(3,3) * t620;
t548 = rSges(3,1) * t574 + rSges(3,2) * t573 - rSges(3,3) * t618;
t547 = Icges(3,1) * t576 + Icges(3,4) * t575 + Icges(3,5) * t620;
t546 = Icges(3,1) * t574 + Icges(3,4) * t573 - Icges(3,5) * t618;
t545 = Icges(3,4) * t576 + Icges(3,2) * t575 + Icges(3,6) * t620;
t544 = Icges(3,4) * t574 + Icges(3,2) * t573 - Icges(3,6) * t618;
t543 = Icges(3,5) * t576 + Icges(3,6) * t575 + Icges(3,3) * t620;
t542 = Icges(3,5) * t574 + Icges(3,6) * t573 - Icges(3,3) * t618;
t525 = t559 * t607 + t572 * t583;
t524 = t559 * t583 - t572 * t607;
t522 = (-t548 * t589 - t570 * t618) * qJD(2);
t521 = (t549 * t589 - t570 * t620) * qJD(2);
t520 = rSges(4,1) * t559 - rSges(4,2) * t558 + rSges(4,3) * t572;
t519 = Icges(4,1) * t559 - Icges(4,4) * t558 + Icges(4,5) * t572;
t517 = Icges(4,5) * t559 - Icges(4,6) * t558 + Icges(4,3) * t572;
t512 = t538 * t607 + t561 * t583;
t511 = t538 * t583 - t561 * t607;
t510 = t536 * t607 + t560 * t583;
t509 = t536 * t583 - t560 * t607;
t508 = t525 * t594 + t558 * t591;
t507 = -t525 * t591 + t558 * t594;
t504 = qJD(1) + (t548 * t585 + t549 * t588) * t614;
t500 = qJD(6) * t524 + t528;
t499 = pkin(5) * t525 + pkin(11) * t524;
t497 = rSges(4,1) * t538 - rSges(4,2) * t537 + rSges(4,3) * t561;
t496 = rSges(4,1) * t536 - rSges(4,2) * t535 + rSges(4,3) * t560;
t495 = Icges(4,1) * t538 - Icges(4,4) * t537 + Icges(4,5) * t561;
t494 = Icges(4,1) * t536 - Icges(4,4) * t535 + Icges(4,5) * t560;
t491 = Icges(4,5) * t538 - Icges(4,6) * t537 + Icges(4,3) * t561;
t490 = Icges(4,5) * t536 - Icges(4,6) * t535 + Icges(4,3) * t560;
t489 = rSges(5,1) * t534 + rSges(5,2) * t533 + rSges(5,3) * t558;
t488 = Icges(5,1) * t534 + Icges(5,4) * t533 + Icges(5,5) * t558;
t487 = Icges(5,4) * t534 + Icges(5,2) * t533 + Icges(5,6) * t558;
t485 = rSges(6,1) * t525 - rSges(6,2) * t524 + rSges(6,3) * t558;
t484 = Icges(6,1) * t525 - Icges(6,4) * t524 + Icges(6,5) * t558;
t483 = Icges(6,4) * t525 - Icges(6,2) * t524 + Icges(6,6) * t558;
t482 = Icges(6,5) * t525 - Icges(6,6) * t524 + Icges(6,3) * t558;
t479 = t512 * t594 + t537 * t591;
t478 = -t512 * t591 + t537 * t594;
t477 = t510 * t594 + t535 * t591;
t476 = -t510 * t591 + t535 * t594;
t475 = qJD(6) * t509 + t506;
t474 = qJD(6) * t511 + t505;
t473 = pkin(5) * t512 + pkin(11) * t511;
t472 = pkin(5) * t510 + pkin(11) * t509;
t470 = rSges(5,1) * t516 + rSges(5,2) * t515 + rSges(5,3) * t537;
t469 = rSges(5,1) * t514 + rSges(5,2) * t513 + rSges(5,3) * t535;
t468 = Icges(5,1) * t516 + Icges(5,4) * t515 + Icges(5,5) * t537;
t467 = Icges(5,1) * t514 + Icges(5,4) * t513 + Icges(5,5) * t535;
t466 = Icges(5,4) * t516 + Icges(5,2) * t515 + Icges(5,6) * t537;
t465 = Icges(5,4) * t514 + Icges(5,2) * t513 + Icges(5,6) * t535;
t462 = rSges(6,1) * t512 - rSges(6,2) * t511 + rSges(6,3) * t537;
t461 = rSges(6,1) * t510 - rSges(6,2) * t509 + rSges(6,3) * t535;
t460 = Icges(6,1) * t512 - Icges(6,4) * t511 + Icges(6,5) * t537;
t459 = Icges(6,1) * t510 - Icges(6,4) * t509 + Icges(6,5) * t535;
t458 = Icges(6,4) * t512 - Icges(6,2) * t511 + Icges(6,6) * t537;
t457 = Icges(6,4) * t510 - Icges(6,2) * t509 + Icges(6,6) * t535;
t456 = Icges(6,5) * t512 - Icges(6,6) * t511 + Icges(6,3) * t537;
t455 = Icges(6,5) * t510 - Icges(6,6) * t509 + Icges(6,3) * t535;
t454 = rSges(7,1) * t508 + rSges(7,2) * t507 + rSges(7,3) * t524;
t453 = Icges(7,1) * t508 + Icges(7,4) * t507 + Icges(7,5) * t524;
t452 = Icges(7,4) * t508 + Icges(7,2) * t507 + Icges(7,6) * t524;
t451 = Icges(7,5) * t508 + Icges(7,6) * t507 + Icges(7,3) * t524;
t446 = -t496 * t563 + t520 * t551 + t600;
t445 = t497 * t563 - t520 * t550 + t602;
t444 = rSges(7,1) * t479 + rSges(7,2) * t478 + rSges(7,3) * t511;
t443 = rSges(7,1) * t477 + rSges(7,2) * t476 + rSges(7,3) * t509;
t442 = Icges(7,1) * t479 + Icges(7,4) * t478 + Icges(7,5) * t511;
t441 = Icges(7,1) * t477 + Icges(7,4) * t476 + Icges(7,5) * t509;
t440 = Icges(7,4) * t479 + Icges(7,2) * t478 + Icges(7,6) * t511;
t439 = Icges(7,4) * t477 + Icges(7,2) * t476 + Icges(7,6) * t509;
t438 = Icges(7,5) * t479 + Icges(7,6) * t478 + Icges(7,3) * t511;
t437 = Icges(7,5) * t477 + Icges(7,6) * t476 + Icges(7,3) * t509;
t436 = t496 * t550 - t497 * t551 + t610;
t435 = t489 * t551 + (-t469 - t501) * t563 + t598;
t434 = t470 * t563 + (-t489 - t523) * t550 + t601;
t433 = t469 * t550 + (-t470 - t502) * t551 + t604;
t432 = -t461 * t528 + t485 * t506 + t596;
t431 = t462 * t528 - t485 * t505 + t597;
t430 = t461 * t505 - t462 * t506 + t599;
t429 = -t443 * t500 + t454 * t475 - t472 * t528 + t499 * t506 + t596;
t428 = t444 * t500 - t454 * t474 + t473 * t528 - t499 * t505 + t597;
t427 = t443 * t474 - t444 * t475 + t472 * t505 - t473 * t506 + t599;
t1 = t500 * ((t438 * t524 + t440 * t507 + t442 * t508) * t474 + (t437 * t524 + t439 * t507 + t441 * t508) * t475 + (t524 * t451 + t507 * t452 + t508 * t453) * t500) / 0.2e1 + t506 * ((t456 * t535 - t458 * t509 + t460 * t510) * t505 + (t535 * t455 - t509 * t457 + t510 * t459) * t506 + (t482 * t535 - t483 * t509 + t484 * t510) * t528) / 0.2e1 + t528 * ((t456 * t558 - t458 * t524 + t460 * t525) * t505 + (t455 * t558 - t457 * t524 + t459 * t525) * t506 + (t558 * t482 - t524 * t483 + t525 * t484) * t528) / 0.2e1 + t505 * ((t537 * t456 - t511 * t458 + t512 * t460) * t505 + (t455 * t537 - t457 * t511 + t459 * t512) * t506 + (t482 * t537 - t483 * t511 + t484 * t512) * t528) / 0.2e1 + m(4) * (t436 ^ 2 + t445 ^ 2 + t446 ^ 2) / 0.2e1 + m(3) * (t504 ^ 2 + t521 ^ 2 + t522 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(7) * (t427 ^ 2 + t428 ^ 2 + t429 ^ 2) / 0.2e1 + m(6) * (t430 ^ 2 + t431 ^ 2 + t432 ^ 2) / 0.2e1 + m(5) * (t433 ^ 2 + t434 ^ 2 + t435 ^ 2) / 0.2e1 + t474 * ((t511 * t438 + t478 * t440 + t479 * t442) * t474 + (t437 * t511 + t439 * t478 + t441 * t479) * t475 + (t451 * t511 + t452 * t478 + t453 * t479) * t500) / 0.2e1 + t475 * ((t438 * t509 + t440 * t476 + t442 * t477) * t474 + (t509 * t437 + t476 * t439 + t477 * t441) * t475 + (t451 * t509 + t452 * t476 + t453 * t477) * t500) / 0.2e1 - t631 * ((-t543 * t618 + t545 * t573 + t547 * t574) * t620 - (-t542 * t618 + t544 * t573 + t546 * t574) * t618 + (-t567 * t618 + t568 * t573 + t569 * t574) * t589) * t618 / 0.2e1 + ((t487 * t515 + t488 * t516 + t517 * t561 + t519 * t538 + t537 * t632) * t563 + (t465 * t515 + t467 * t516 + t490 * t561 + t494 * t538 + t537 * t634) * t551 + (t515 * t466 + t516 * t468 + t561 * t491 + t538 * t495 + t537 * t633) * t550) * t550 / 0.2e1 + ((t487 * t513 + t488 * t514 + t517 * t560 + t519 * t536 + t535 * t632) * t563 + (t513 * t465 + t514 * t467 + t560 * t490 + t536 * t494 + t535 * t634) * t551 + (t466 * t513 + t468 * t514 + t491 * t560 + t495 * t536 + t535 * t633) * t550) * t551 / 0.2e1 + ((t533 * t487 + t534 * t488 + t572 * t517 + t559 * t519 + t632 * t558) * t563 + (t465 * t533 + t467 * t534 + t490 * t572 + t494 * t559 + t558 * t634) * t551 + (t466 * t533 + t468 * t534 + t491 * t572 + t495 * t559 + t558 * t633) * t550) * t563 / 0.2e1 + (t589 * (t589 ^ 2 * t567 + (((t545 * t595 + t547 * t593) * t585 - (t544 * t595 + t546 * t593) * t588) * t586 + (-t542 * t588 + t543 * t585 + t568 * t595 + t569 * t593) * t589) * t586) + ((t543 * t620 + t545 * t575 + t547 * t576) * t620 - (t542 * t620 + t544 * t575 + t546 * t576) * t618 + (t567 * t620 + t568 * t575 + t569 * t576) * t589) * t620) * t631 / 0.2e1;
T  = t1;
