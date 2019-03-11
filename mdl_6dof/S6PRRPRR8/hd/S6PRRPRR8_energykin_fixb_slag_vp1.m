% Calculate kinetic energy for
% S6PRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRR8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR8_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR8_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRR8_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:35:59
% EndTime: 2019-03-08 22:36:02
% DurationCPUTime: 2.96s
% Computational Cost: add. (4083->339), mult. (11203->524), div. (0->0), fcn. (14272->14), ass. (0->159)
t623 = Icges(4,1) + Icges(5,2);
t622 = Icges(5,1) + Icges(4,3);
t621 = -Icges(4,4) - Icges(5,6);
t620 = -Icges(5,4) + Icges(4,5);
t619 = Icges(5,5) - Icges(4,6);
t618 = Icges(4,2) + Icges(5,3);
t567 = sin(pkin(12));
t569 = cos(pkin(12));
t574 = sin(qJ(2));
t570 = cos(pkin(6));
t576 = cos(qJ(2));
t596 = t570 * t576;
t557 = -t567 * t574 + t569 * t596;
t597 = t570 * t574;
t558 = t567 * t576 + t569 * t597;
t573 = sin(qJ(3));
t568 = sin(pkin(6));
t601 = sin(pkin(7));
t604 = cos(qJ(3));
t586 = t604 * t601;
t584 = t568 * t586;
t602 = cos(pkin(7));
t587 = t602 * t604;
t518 = -t557 * t587 + t558 * t573 + t569 * t584;
t588 = t573 * t601;
t589 = t573 * t602;
t599 = t568 * t569;
t519 = t557 * t589 + t558 * t604 - t588 * t599;
t591 = t568 * t602;
t545 = -t557 * t601 - t569 * t591;
t617 = t518 * t618 + t519 * t621 + t545 * t619;
t559 = -t567 * t596 - t569 * t574;
t560 = -t567 * t597 + t569 * t576;
t520 = -t559 * t587 + t560 * t573 - t567 * t584;
t590 = t568 * t601;
t521 = t560 * t604 + (t559 * t602 + t567 * t590) * t573;
t546 = -t559 * t601 + t567 * t591;
t616 = t520 * t618 + t521 * t621 + t546 * t619;
t615 = t518 * t619 + t519 * t620 + t545 * t622;
t614 = t520 * t619 + t521 * t620 + t546 * t622;
t613 = t621 * t518 + t519 * t623 + t620 * t545;
t612 = t621 * t520 + t521 * t623 + t620 * t546;
t598 = t568 * t574;
t543 = -t568 * t576 * t587 - t570 * t586 + t573 * t598;
t544 = t570 * t588 + (t574 * t604 + t576 * t589) * t568;
t556 = t570 * t602 - t576 * t590;
t611 = t543 * t618 + t544 * t621 + t556 * t619;
t610 = t543 * t619 + t544 * t620 + t556 * t622;
t609 = t621 * t543 + t544 * t623 + t620 * t556;
t608 = qJD(2) ^ 2;
t603 = cos(qJ(5));
t600 = t567 * t568;
t595 = qJD(2) * t568;
t565 = t567 * t595;
t536 = qJD(3) * t546 + t565;
t566 = qJD(2) * t570;
t548 = qJD(3) * t556 + t566;
t487 = qJD(5) * t521 + t536;
t513 = qJD(5) * t544 + t548;
t593 = t569 * t595;
t526 = t558 * pkin(2) + pkin(9) * t545;
t527 = t560 * pkin(2) + pkin(9) * t546;
t592 = t526 * t565 + t527 * t593 + qJD(1);
t537 = qJD(3) * t545 - t593;
t480 = pkin(3) * t519 + qJ(4) * t518;
t585 = qJD(4) * t543 + t536 * t480 + t592;
t547 = pkin(2) * t598 + pkin(9) * t556;
t583 = t527 * t566 - t547 * t565;
t488 = qJD(5) * t519 + t537;
t481 = pkin(3) * t521 + qJ(4) * t520;
t582 = qJD(4) * t518 + t548 * t481 + t583;
t581 = (-t526 * t570 - t547 * t599) * qJD(2);
t496 = pkin(4) * t545 + pkin(10) * t519;
t497 = pkin(4) * t546 + pkin(10) * t521;
t580 = t536 * t496 + (-t481 - t497) * t537 + t585;
t508 = pkin(3) * t544 + qJ(4) * t543;
t579 = qJD(4) * t520 + t537 * t508 + t581;
t525 = pkin(4) * t556 + pkin(10) * t544;
t578 = t548 * t497 + (-t508 - t525) * t536 + t582;
t577 = t537 * t525 + (-t480 - t496) * t548 + t579;
t575 = cos(qJ(6));
t572 = sin(qJ(5));
t571 = sin(qJ(6));
t554 = t570 * rSges(3,3) + (rSges(3,1) * t574 + rSges(3,2) * t576) * t568;
t553 = Icges(3,5) * t570 + (Icges(3,1) * t574 + Icges(3,4) * t576) * t568;
t552 = Icges(3,6) * t570 + (Icges(3,4) * t574 + Icges(3,2) * t576) * t568;
t551 = Icges(3,3) * t570 + (Icges(3,5) * t574 + Icges(3,6) * t576) * t568;
t535 = rSges(3,1) * t560 + rSges(3,2) * t559 + rSges(3,3) * t600;
t534 = rSges(3,1) * t558 + rSges(3,2) * t557 - rSges(3,3) * t599;
t533 = Icges(3,1) * t560 + Icges(3,4) * t559 + Icges(3,5) * t600;
t532 = Icges(3,1) * t558 + Icges(3,4) * t557 - Icges(3,5) * t599;
t531 = Icges(3,4) * t560 + Icges(3,2) * t559 + Icges(3,6) * t600;
t530 = Icges(3,4) * t558 + Icges(3,2) * t557 - Icges(3,6) * t599;
t529 = Icges(3,5) * t560 + Icges(3,6) * t559 + Icges(3,3) * t600;
t528 = Icges(3,5) * t558 + Icges(3,6) * t557 - Icges(3,3) * t599;
t524 = t543 * t572 + t556 * t603;
t523 = -t543 * t603 + t556 * t572;
t507 = (-t534 * t570 - t554 * t599) * qJD(2);
t506 = (t535 * t570 - t554 * t600) * qJD(2);
t505 = rSges(5,1) * t556 - rSges(5,2) * t544 + rSges(5,3) * t543;
t504 = rSges(4,1) * t544 - rSges(4,2) * t543 + rSges(4,3) * t556;
t495 = t520 * t572 + t546 * t603;
t494 = -t520 * t603 + t546 * t572;
t493 = t518 * t572 + t545 * t603;
t492 = -t518 * t603 + t545 * t572;
t491 = t524 * t575 + t544 * t571;
t490 = -t524 * t571 + t544 * t575;
t486 = qJD(1) + (t534 * t567 + t535 * t569) * t595;
t484 = pkin(5) * t524 + pkin(11) * t523;
t483 = qJD(6) * t523 + t513;
t477 = rSges(6,1) * t524 - rSges(6,2) * t523 + rSges(6,3) * t544;
t476 = rSges(5,1) * t546 - rSges(5,2) * t521 + rSges(5,3) * t520;
t475 = rSges(5,1) * t545 - rSges(5,2) * t519 + rSges(5,3) * t518;
t474 = rSges(4,1) * t521 - rSges(4,2) * t520 + rSges(4,3) * t546;
t473 = rSges(4,1) * t519 - rSges(4,2) * t518 + rSges(4,3) * t545;
t472 = Icges(6,1) * t524 - Icges(6,4) * t523 + Icges(6,5) * t544;
t471 = Icges(6,4) * t524 - Icges(6,2) * t523 + Icges(6,6) * t544;
t470 = Icges(6,5) * t524 - Icges(6,6) * t523 + Icges(6,3) * t544;
t457 = t495 * t575 + t521 * t571;
t456 = -t495 * t571 + t521 * t575;
t455 = t493 * t575 + t519 * t571;
t454 = -t493 * t571 + t519 * t575;
t452 = pkin(5) * t495 + pkin(11) * t494;
t451 = pkin(5) * t493 + pkin(11) * t492;
t450 = qJD(6) * t492 + t488;
t449 = qJD(6) * t494 + t487;
t448 = rSges(6,1) * t495 - rSges(6,2) * t494 + rSges(6,3) * t521;
t447 = rSges(6,1) * t493 - rSges(6,2) * t492 + rSges(6,3) * t519;
t446 = rSges(7,1) * t491 + rSges(7,2) * t490 + rSges(7,3) * t523;
t445 = Icges(6,1) * t495 - Icges(6,4) * t494 + Icges(6,5) * t521;
t444 = Icges(6,1) * t493 - Icges(6,4) * t492 + Icges(6,5) * t519;
t443 = Icges(7,1) * t491 + Icges(7,4) * t490 + Icges(7,5) * t523;
t442 = Icges(6,4) * t495 - Icges(6,2) * t494 + Icges(6,6) * t521;
t441 = Icges(6,4) * t493 - Icges(6,2) * t492 + Icges(6,6) * t519;
t440 = Icges(7,4) * t491 + Icges(7,2) * t490 + Icges(7,6) * t523;
t439 = Icges(6,5) * t495 - Icges(6,6) * t494 + Icges(6,3) * t521;
t438 = Icges(6,5) * t493 - Icges(6,6) * t492 + Icges(6,3) * t519;
t437 = Icges(7,5) * t491 + Icges(7,6) * t490 + Icges(7,3) * t523;
t436 = -t473 * t548 + t504 * t537 + t581;
t435 = t474 * t548 - t504 * t536 + t583;
t434 = rSges(7,1) * t457 + rSges(7,2) * t456 + rSges(7,3) * t494;
t433 = rSges(7,1) * t455 + rSges(7,2) * t454 + rSges(7,3) * t492;
t432 = Icges(7,1) * t457 + Icges(7,4) * t456 + Icges(7,5) * t494;
t431 = Icges(7,1) * t455 + Icges(7,4) * t454 + Icges(7,5) * t492;
t430 = Icges(7,4) * t457 + Icges(7,2) * t456 + Icges(7,6) * t494;
t429 = Icges(7,4) * t455 + Icges(7,2) * t454 + Icges(7,6) * t492;
t428 = Icges(7,5) * t457 + Icges(7,6) * t456 + Icges(7,3) * t494;
t427 = Icges(7,5) * t455 + Icges(7,6) * t454 + Icges(7,3) * t492;
t426 = t473 * t536 - t474 * t537 + t592;
t425 = t505 * t537 + (-t475 - t480) * t548 + t579;
t424 = t476 * t548 + (-t505 - t508) * t536 + t582;
t423 = t475 * t536 + (-t476 - t481) * t537 + t585;
t422 = -t447 * t513 + t477 * t488 + t577;
t421 = t448 * t513 - t477 * t487 + t578;
t420 = t447 * t487 - t448 * t488 + t580;
t419 = -t433 * t483 + t446 * t450 - t451 * t513 + t484 * t488 + t577;
t418 = t434 * t483 - t446 * t449 + t452 * t513 - t484 * t487 + t578;
t417 = t433 * t449 - t434 * t450 + t451 * t487 - t452 * t488 + t580;
t1 = -t608 * ((-t529 * t599 + t531 * t557 + t533 * t558) * t600 - (-t528 * t599 + t530 * t557 + t532 * t558) * t599 + (-t551 * t599 + t552 * t557 + t553 * t558) * t570) * t599 / 0.2e1 + t450 * ((t428 * t492 + t430 * t454 + t432 * t455) * t449 + (t427 * t492 + t429 * t454 + t431 * t455) * t450 + (t437 * t492 + t440 * t454 + t443 * t455) * t483) / 0.2e1 + t483 * ((t428 * t523 + t430 * t490 + t432 * t491) * t449 + (t427 * t523 + t429 * t490 + t431 * t491) * t450 + (t437 * t523 + t440 * t490 + t443 * t491) * t483) / 0.2e1 + t449 * ((t428 * t494 + t430 * t456 + t432 * t457) * t449 + (t427 * t494 + t429 * t456 + t431 * t457) * t450 + (t437 * t494 + t440 * t456 + t443 * t457) * t483) / 0.2e1 + t487 * ((t439 * t521 - t442 * t494 + t445 * t495) * t487 + (t438 * t521 - t441 * t494 + t444 * t495) * t488 + (t470 * t521 - t471 * t494 + t472 * t495) * t513) / 0.2e1 + t488 * ((t439 * t519 - t442 * t492 + t445 * t493) * t487 + (t438 * t519 - t441 * t492 + t444 * t493) * t488 + (t470 * t519 - t471 * t492 + t472 * t493) * t513) / 0.2e1 + t513 * ((t439 * t544 - t442 * t523 + t445 * t524) * t487 + (t438 * t544 - t441 * t523 + t444 * t524) * t488 + (t470 * t544 - t471 * t523 + t472 * t524) * t513) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(6) * (t420 ^ 2 + t421 ^ 2 + t422 ^ 2) / 0.2e1 + m(7) * (t417 ^ 2 + t418 ^ 2 + t419 ^ 2) / 0.2e1 + m(5) * (t423 ^ 2 + t424 ^ 2 + t425 ^ 2) / 0.2e1 + m(4) * (t426 ^ 2 + t435 ^ 2 + t436 ^ 2) / 0.2e1 + m(3) * (t486 ^ 2 + t506 ^ 2 + t507 ^ 2) / 0.2e1 + ((t611 * t520 + t609 * t521 + t610 * t546) * t548 + (t617 * t520 + t613 * t521 + t615 * t546) * t537 + (t616 * t520 + t612 * t521 + t614 * t546) * t536) * t536 / 0.2e1 + ((t611 * t518 + t609 * t519 + t610 * t545) * t548 + (t617 * t518 + t613 * t519 + t615 * t545) * t537 + (t616 * t518 + t612 * t519 + t614 * t545) * t536) * t537 / 0.2e1 + ((t611 * t543 + t609 * t544 + t610 * t556) * t548 + (t617 * t543 + t613 * t544 + t615 * t556) * t537 + (t616 * t543 + t612 * t544 + t614 * t556) * t536) * t548 / 0.2e1 + (((t529 * t600 + t531 * t559 + t533 * t560) * t600 - (t528 * t600 + t530 * t559 + t532 * t560) * t599 + (t551 * t600 + t552 * t559 + t553 * t560) * t570) * t600 + t570 * (t570 ^ 2 * t551 + (((t531 * t576 + t533 * t574) * t567 - (t530 * t576 + t532 * t574) * t569) * t568 + (-t528 * t569 + t529 * t567 + t552 * t576 + t553 * t574) * t570) * t568)) * t608 / 0.2e1;
T  = t1;
