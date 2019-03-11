% Calculate kinetic energy for
% S6RPRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 04:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR9_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_energykin_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR9_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR9_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR9_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:01:14
% EndTime: 2019-03-09 04:01:17
% DurationCPUTime: 3.09s
% Computational Cost: add. (5488->352), mult. (14771->537), div. (0->0), fcn. (19266->16), ass. (0->165)
t623 = Icges(4,3) + Icges(5,3);
t571 = sin(pkin(12));
t574 = cos(pkin(12));
t583 = cos(qJ(1));
t576 = cos(pkin(6));
t580 = sin(qJ(1));
t603 = t576 * t580;
t560 = -t571 * t583 - t574 * t603;
t572 = sin(pkin(7));
t575 = cos(pkin(7));
t573 = sin(pkin(6));
t606 = t573 * t580;
t546 = -t560 * t572 + t575 * t606;
t622 = pkin(9) * t546;
t602 = t576 * t583;
t558 = -t571 * t580 + t574 * t602;
t605 = t573 * t583;
t594 = t558 * t572 + t575 * t605;
t621 = pkin(9) * t594;
t579 = sin(qJ(3));
t582 = cos(qJ(3));
t609 = sin(pkin(13));
t610 = cos(pkin(13));
t564 = -t579 * t609 + t582 * t610;
t551 = t564 * t575;
t559 = t571 * t602 + t574 * t580;
t563 = -t579 * t610 - t582 * t609;
t588 = t572 * t564;
t587 = t573 * t588;
t516 = t551 * t558 + t559 * t563 - t583 * t587;
t550 = t563 * t572;
t552 = t563 * t575;
t517 = t550 * t605 - t552 * t558 + t559 * t564;
t593 = t558 * t575 - t572 * t605;
t531 = -t559 * t579 + t582 * t593;
t532 = t559 * t582 + t579 * t593;
t620 = Icges(4,5) * t532 + Icges(5,5) * t517 + Icges(4,6) * t531 + Icges(5,6) * t516 - t623 * t594;
t561 = -t571 * t603 + t574 * t583;
t518 = t560 * t551 + t561 * t563 + t580 * t587;
t519 = -t550 * t606 - t552 * t560 + t561 * t564;
t592 = t560 * t575 + t572 * t606;
t533 = -t561 * t579 + t582 * t592;
t534 = t561 * t582 + t579 * t592;
t619 = Icges(4,5) * t534 + Icges(5,5) * t519 + Icges(4,6) * t533 + Icges(5,6) * t518 + t623 * t546;
t607 = t573 * t574;
t527 = t573 * t571 * t563 + t551 * t607 + t576 * t588;
t528 = -t550 * t576 + (-t552 * t574 + t564 * t571) * t573;
t543 = t572 * t576 * t582 + (t574 * t575 * t582 - t571 * t579) * t573;
t604 = t575 * t579;
t608 = t572 * t579;
t544 = t576 * t608 + (t571 * t582 + t574 * t604) * t573;
t557 = -t572 * t607 + t575 * t576;
t618 = Icges(4,5) * t544 + Icges(5,5) * t528 + Icges(4,6) * t543 + Icges(5,6) * t527 + t623 * t557;
t617 = t619 * t546 - t594 * t620;
t616 = t576 ^ 2;
t613 = cos(qJ(5));
t612 = pkin(3) * t582;
t611 = pkin(9) + qJ(4);
t541 = qJD(3) * t594;
t509 = -qJD(5) * t516 - t541;
t542 = qJD(3) * t546;
t510 = -qJD(5) * t518 + t542;
t601 = qJD(2) * t573;
t548 = qJD(3) * t557 + qJD(1);
t555 = pkin(3) * t608 + t575 * t611;
t556 = pkin(3) * t604 - t572 * t611;
t503 = -t555 * t605 + t556 * t558 + t559 * t612 + t621;
t570 = qJD(2) * t576;
t600 = qJD(4) * t557 + t503 * t542 + t570;
t520 = -qJD(5) * t527 + t548;
t597 = -t583 * t601 + qJD(1) * (pkin(1) * t583 + qJ(2) * t606);
t565 = pkin(1) * t580 - qJ(2) * t605;
t568 = t580 * t601;
t596 = t568 + (-pkin(2) * t559 - t565 + t621) * qJD(1);
t595 = qJD(1) * (pkin(2) * t561 + t622) + t597;
t525 = (-pkin(9) * t575 + t555) * t576 + ((pkin(9) * t572 + t556) * t574 + t612 * t571) * t573;
t591 = qJD(4) * t546 - t525 * t541 + t596;
t504 = t555 * t606 + t556 * t560 + t561 * t612 - t622;
t590 = -qJD(4) * t594 + t548 * t504 + t595;
t482 = pkin(4) * t517 - pkin(10) * t516;
t483 = pkin(4) * t519 - pkin(10) * t518;
t589 = t482 * t542 - (-t483 - t504) * t541 + t600;
t502 = pkin(4) * t528 - pkin(10) * t527;
t586 = -t502 * t541 + (-t482 - t503) * t548 + t591;
t585 = t548 * t483 + (-t502 - t525) * t542 + t590;
t581 = cos(qJ(6));
t578 = sin(qJ(5));
t577 = sin(qJ(6));
t567 = rSges(2,1) * t583 - rSges(2,2) * t580;
t566 = rSges(2,1) * t580 + rSges(2,2) * t583;
t530 = qJD(1) * (rSges(3,1) * t561 + rSges(3,2) * t560 + rSges(3,3) * t606) + t597;
t529 = t568 + (-rSges(3,1) * t559 - rSges(3,2) * t558 + rSges(3,3) * t605 - t565) * qJD(1);
t524 = rSges(4,1) * t544 + rSges(4,2) * t543 + rSges(4,3) * t557;
t523 = Icges(4,1) * t544 + Icges(4,4) * t543 + Icges(4,5) * t557;
t522 = Icges(4,4) * t544 + Icges(4,2) * t543 + Icges(4,6) * t557;
t515 = t528 * t613 + t557 * t578;
t514 = t528 * t578 - t557 * t613;
t508 = t519 * t613 + t546 * t578;
t507 = t519 * t578 - t546 * t613;
t506 = t517 * t613 - t578 * t594;
t505 = t517 * t578 + t594 * t613;
t501 = rSges(4,1) * t534 + rSges(4,2) * t533 + rSges(4,3) * t546;
t500 = rSges(4,1) * t532 + rSges(4,2) * t531 - rSges(4,3) * t594;
t499 = Icges(4,1) * t534 + Icges(4,4) * t533 + Icges(4,5) * t546;
t498 = Icges(4,1) * t532 + Icges(4,4) * t531 - Icges(4,5) * t594;
t497 = Icges(4,4) * t534 + Icges(4,2) * t533 + Icges(4,6) * t546;
t496 = Icges(4,4) * t532 + Icges(4,2) * t531 - Icges(4,6) * t594;
t492 = rSges(5,1) * t528 + rSges(5,2) * t527 + rSges(5,3) * t557;
t490 = Icges(5,1) * t528 + Icges(5,4) * t527 + Icges(5,5) * t557;
t489 = Icges(5,4) * t528 + Icges(5,2) * t527 + Icges(5,6) * t557;
t487 = t515 * t581 - t527 * t577;
t486 = -t515 * t577 - t527 * t581;
t484 = qJD(6) * t514 + t520;
t481 = pkin(5) * t515 + pkin(11) * t514;
t478 = rSges(5,1) * t519 + rSges(5,2) * t518 + rSges(5,3) * t546;
t477 = rSges(5,1) * t517 + rSges(5,2) * t516 - rSges(5,3) * t594;
t476 = Icges(5,1) * t519 + Icges(5,4) * t518 + Icges(5,5) * t546;
t475 = Icges(5,1) * t517 + Icges(5,4) * t516 - Icges(5,5) * t594;
t474 = Icges(5,4) * t519 + Icges(5,2) * t518 + Icges(5,6) * t546;
t473 = Icges(5,4) * t517 + Icges(5,2) * t516 - Icges(5,6) * t594;
t470 = t508 * t581 - t518 * t577;
t469 = -t508 * t577 - t518 * t581;
t468 = t506 * t581 - t516 * t577;
t467 = -t506 * t577 - t516 * t581;
t466 = qJD(6) * t507 + t510;
t465 = qJD(6) * t505 + t509;
t464 = rSges(6,1) * t515 - rSges(6,2) * t514 - rSges(6,3) * t527;
t463 = Icges(6,1) * t515 - Icges(6,4) * t514 - Icges(6,5) * t527;
t462 = Icges(6,4) * t515 - Icges(6,2) * t514 - Icges(6,6) * t527;
t461 = Icges(6,5) * t515 - Icges(6,6) * t514 - Icges(6,3) * t527;
t460 = pkin(5) * t508 + pkin(11) * t507;
t459 = pkin(5) * t506 + pkin(11) * t505;
t458 = rSges(6,1) * t508 - rSges(6,2) * t507 - rSges(6,3) * t518;
t457 = rSges(6,1) * t506 - rSges(6,2) * t505 - rSges(6,3) * t516;
t456 = Icges(6,1) * t508 - Icges(6,4) * t507 - Icges(6,5) * t518;
t455 = Icges(6,1) * t506 - Icges(6,4) * t505 - Icges(6,5) * t516;
t454 = Icges(6,4) * t508 - Icges(6,2) * t507 - Icges(6,6) * t518;
t453 = Icges(6,4) * t506 - Icges(6,2) * t505 - Icges(6,6) * t516;
t452 = Icges(6,5) * t508 - Icges(6,6) * t507 - Icges(6,3) * t518;
t451 = Icges(6,5) * t506 - Icges(6,6) * t505 - Icges(6,3) * t516;
t450 = t501 * t548 - t524 * t542 + t595;
t449 = -t500 * t548 - t524 * t541 + t596;
t448 = t570 + (t500 * t546 + t501 * t594) * qJD(3);
t447 = rSges(7,1) * t487 + rSges(7,2) * t486 + rSges(7,3) * t514;
t446 = Icges(7,1) * t487 + Icges(7,4) * t486 + Icges(7,5) * t514;
t445 = Icges(7,4) * t487 + Icges(7,2) * t486 + Icges(7,6) * t514;
t444 = Icges(7,5) * t487 + Icges(7,6) * t486 + Icges(7,3) * t514;
t443 = rSges(7,1) * t470 + rSges(7,2) * t469 + rSges(7,3) * t507;
t442 = rSges(7,1) * t468 + rSges(7,2) * t467 + rSges(7,3) * t505;
t441 = Icges(7,1) * t470 + Icges(7,4) * t469 + Icges(7,5) * t507;
t440 = Icges(7,1) * t468 + Icges(7,4) * t467 + Icges(7,5) * t505;
t439 = Icges(7,4) * t470 + Icges(7,2) * t469 + Icges(7,6) * t507;
t438 = Icges(7,4) * t468 + Icges(7,2) * t467 + Icges(7,6) * t505;
t437 = Icges(7,5) * t470 + Icges(7,6) * t469 + Icges(7,3) * t507;
t436 = Icges(7,5) * t468 + Icges(7,6) * t467 + Icges(7,3) * t505;
t435 = t478 * t548 + (-t492 - t525) * t542 + t590;
t434 = -t492 * t541 + (-t477 - t503) * t548 + t591;
t433 = (t477 * t546 - (-t478 - t504) * t594) * qJD(3) + t600;
t432 = t458 * t520 - t464 * t510 + t585;
t431 = -t457 * t520 + t464 * t509 + t586;
t430 = t457 * t510 - t458 * t509 + t589;
t429 = t443 * t484 - t447 * t466 + t460 * t520 - t481 * t510 + t585;
t428 = -t442 * t484 + t447 * t465 - t459 * t520 + t481 * t509 + t586;
t427 = t442 * t466 - t443 * t465 + t459 * t510 - t460 * t509 + t589;
t1 = m(4) * (t448 ^ 2 + t449 ^ 2 + t450 ^ 2) / 0.2e1 + m(5) * (t433 ^ 2 + t434 ^ 2 + t435 ^ 2) / 0.2e1 + m(6) * (t430 ^ 2 + t431 ^ 2 + t432 ^ 2) / 0.2e1 + m(7) * (t427 ^ 2 + t428 ^ 2 + t429 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 * t616 + t529 ^ 2 + t530 ^ 2) / 0.2e1 + t510 * ((-t518 * t452 - t507 * t454 + t508 * t456) * t510 + (-t451 * t518 - t453 * t507 + t455 * t508) * t509 + (-t461 * t518 - t462 * t507 + t463 * t508) * t520) / 0.2e1 + t509 * ((-t452 * t516 - t454 * t505 + t456 * t506) * t510 + (-t516 * t451 - t505 * t453 + t506 * t455) * t509 + (-t461 * t516 - t462 * t505 + t463 * t506) * t520) / 0.2e1 + t520 * ((-t452 * t527 - t454 * t514 + t456 * t515) * t510 + (-t451 * t527 - t453 * t514 + t455 * t515) * t509 + (-t527 * t461 - t514 * t462 + t515 * t463) * t520) / 0.2e1 + t466 * ((t507 * t437 + t469 * t439 + t470 * t441) * t466 + (t436 * t507 + t438 * t469 + t440 * t470) * t465 + (t444 * t507 + t445 * t469 + t446 * t470) * t484) / 0.2e1 + t465 * ((t437 * t505 + t439 * t467 + t441 * t468) * t466 + (t505 * t436 + t467 * t438 + t468 * t440) * t465 + (t444 * t505 + t445 * t467 + t446 * t468) * t484) / 0.2e1 + t484 * ((t437 * t514 + t439 * t486 + t441 * t487) * t466 + (t436 * t514 + t438 * t486 + t440 * t487) * t465 + (t514 * t444 + t486 * t445 + t487 * t446) * t484) / 0.2e1 + ((t527 * t489 + t528 * t490 + t543 * t522 + t544 * t523 + t618 * t557) * t548 + ((t474 * t527 + t476 * t528 + t497 * t543 + t499 * t544 + t619 * t557) * t546 - (t473 * t527 + t475 * t528 + t496 * t543 + t498 * t544 + t620 * t557) * t594) * qJD(3)) * t548 / 0.2e1 - ((t489 * t516 + t490 * t517 + t522 * t531 + t523 * t532 - t594 * t618) * t548 + ((t474 * t516 + t476 * t517 + t497 * t531 + t499 * t532) * t546 - (t473 * t516 + t475 * t517 + t496 * t531 + t498 * t532 + t617) * t594) * qJD(3)) * t541 / 0.2e1 + ((t489 * t518 + t490 * t519 + t522 * t533 + t523 * t534 + t618 * t546) * t548 + (-(t473 * t518 + t475 * t519 + t496 * t533 + t498 * t534) * t594 + (t474 * t518 + t476 * t519 + t497 * t533 + t499 * t534 + t617) * t546) * qJD(3)) * t542 / 0.2e1 + (m(2) * (t566 ^ 2 + t567 ^ 2) + Icges(2,3) + Icges(3,3) * t616 + ((Icges(3,2) * t574 ^ 2 + (Icges(3,1) * t571 + 0.2e1 * Icges(3,4) * t574) * t571) * t573 + 0.2e1 * t576 * (Icges(3,5) * t571 + Icges(3,6) * t574)) * t573) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
