% Calculate kinetic energy for
% S6RRRRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-10 00:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR15_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR15_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR15_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR15_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR15_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:34:11
% EndTime: 2019-03-10 00:34:14
% DurationCPUTime: 3.22s
% Computational Cost: add. (4647->346), mult. (12605->531), div. (0->0), fcn. (16124->14), ass. (0->162)
t628 = Icges(5,1) + Icges(6,2);
t627 = Icges(6,1) + Icges(5,3);
t626 = -Icges(5,4) - Icges(6,6);
t625 = -Icges(6,4) + Icges(5,5);
t624 = Icges(6,5) - Icges(5,6);
t623 = Icges(5,2) + Icges(6,3);
t569 = sin(pkin(6));
t570 = cos(pkin(6));
t577 = cos(qJ(2));
t578 = cos(qJ(1));
t598 = t577 * t578;
t574 = sin(qJ(2));
t575 = sin(qJ(1));
t601 = t574 * t575;
t555 = t570 * t598 - t601;
t599 = t575 * t577;
t600 = t574 * t578;
t556 = t570 * t600 + t599;
t557 = -t570 * t599 - t600;
t558 = -t570 * t601 + t598;
t602 = t569 * t578;
t603 = t569 * t575;
t589 = (Icges(3,5) * t556 + Icges(3,6) * t555 - Icges(3,3) * t602) * t578 - (Icges(3,5) * t558 + Icges(3,6) * t557 + Icges(3,3) * t603) * t575;
t622 = t569 * t589;
t573 = sin(qJ(3));
t605 = sin(pkin(7));
t592 = t569 * t605;
t606 = cos(pkin(7));
t608 = cos(qJ(3));
t520 = t556 * t608 + (t555 * t606 - t578 * t592) * t573;
t593 = t569 * t606;
t542 = -t555 * t605 - t578 * t593;
t572 = sin(qJ(4));
t607 = cos(qJ(4));
t500 = t520 * t572 - t542 * t607;
t501 = t520 * t607 + t542 * t572;
t590 = t608 * t605;
t588 = t569 * t590;
t591 = t606 * t608;
t519 = -t555 * t591 + t556 * t573 + t578 * t588;
t621 = t500 * t623 + t501 * t626 + t519 * t624;
t522 = t558 * t608 + (t557 * t606 + t575 * t592) * t573;
t543 = -t557 * t605 + t575 * t593;
t502 = t522 * t572 - t543 * t607;
t503 = t522 * t607 + t543 * t572;
t521 = -t557 * t591 + t558 * t573 - t575 * t588;
t620 = t502 * t623 + t503 * t626 + t521 * t624;
t619 = t500 * t624 + t501 * t625 + t519 * t627;
t618 = t502 * t624 + t503 * t625 + t521 * t627;
t617 = t626 * t500 + t501 * t628 + t625 * t519;
t616 = t626 * t502 + t503 * t628 + t625 * t521;
t541 = t570 * t605 * t573 + (t573 * t577 * t606 + t574 * t608) * t569;
t554 = t570 * t606 - t577 * t592;
t517 = t541 * t572 - t554 * t607;
t518 = t541 * t607 + t554 * t572;
t604 = t569 * t574;
t540 = -t569 * t577 * t591 - t570 * t590 + t573 * t604;
t615 = t517 * t623 + t518 * t626 + t540 * t624;
t614 = t517 * t624 + t518 * t625 + t540 * t627;
t613 = t626 * t517 + t518 * t628 + t625 * t540;
t523 = t556 * pkin(2) + pkin(10) * t542;
t524 = t558 * pkin(2) + pkin(10) * t543;
t595 = qJD(2) * t569;
t566 = t575 * t595;
t594 = t578 * t595;
t597 = t523 * t566 + t524 * t594;
t533 = qJD(3) * t543 + t566;
t596 = qJD(1) * (pkin(1) * t575 - pkin(9) * t602);
t567 = qJD(2) * t570 + qJD(1);
t492 = qJD(4) * t521 + t533;
t544 = qJD(3) * t554 + t567;
t510 = qJD(4) * t540 + t544;
t534 = qJD(3) * t542 - t594;
t488 = pkin(3) * t520 + pkin(11) * t519;
t489 = pkin(3) * t522 + pkin(11) * t521;
t587 = t533 * t488 - t489 * t534 + t597;
t493 = qJD(4) * t519 + t534;
t545 = pkin(2) * t604 + pkin(10) * t554;
t559 = qJD(1) * (pkin(1) * t578 + pkin(9) * t603);
t586 = t567 * t524 - t545 * t566 + t559;
t460 = pkin(4) * t501 + qJ(5) * t500;
t585 = qJD(5) * t517 + t492 * t460 + t587;
t584 = -t523 * t567 - t545 * t594 - t596;
t509 = pkin(3) * t541 + pkin(11) * t540;
t583 = t544 * t489 - t509 * t533 + t586;
t461 = pkin(4) * t503 + qJ(5) * t502;
t582 = qJD(5) * t500 + t510 * t461 + t583;
t581 = -t488 * t544 + t534 * t509 + t584;
t487 = pkin(4) * t518 + qJ(5) * t517;
t580 = qJD(5) * t502 + t493 * t487 + t581;
t576 = cos(qJ(6));
t571 = sin(qJ(6));
t564 = rSges(2,1) * t578 - rSges(2,2) * t575;
t563 = rSges(2,1) * t575 + rSges(2,2) * t578;
t551 = rSges(3,3) * t570 + (rSges(3,1) * t574 + rSges(3,2) * t577) * t569;
t550 = Icges(3,5) * t570 + (Icges(3,1) * t574 + Icges(3,4) * t577) * t569;
t549 = Icges(3,6) * t570 + (Icges(3,4) * t574 + Icges(3,2) * t577) * t569;
t548 = Icges(3,3) * t570 + (Icges(3,5) * t574 + Icges(3,6) * t577) * t569;
t532 = rSges(3,1) * t558 + rSges(3,2) * t557 + rSges(3,3) * t603;
t531 = rSges(3,1) * t556 + rSges(3,2) * t555 - rSges(3,3) * t602;
t530 = Icges(3,1) * t558 + Icges(3,4) * t557 + Icges(3,5) * t603;
t529 = Icges(3,1) * t556 + Icges(3,4) * t555 - Icges(3,5) * t602;
t528 = Icges(3,4) * t558 + Icges(3,2) * t557 + Icges(3,6) * t603;
t527 = Icges(3,4) * t556 + Icges(3,2) * t555 - Icges(3,6) * t602;
t508 = rSges(4,1) * t541 - rSges(4,2) * t540 + rSges(4,3) * t554;
t507 = Icges(4,1) * t541 - Icges(4,4) * t540 + Icges(4,5) * t554;
t506 = Icges(4,4) * t541 - Icges(4,2) * t540 + Icges(4,6) * t554;
t505 = Icges(4,5) * t541 - Icges(4,6) * t540 + Icges(4,3) * t554;
t504 = pkin(5) * t540 + pkin(12) * t518;
t497 = t517 * t571 + t540 * t576;
t496 = t517 * t576 - t540 * t571;
t495 = t532 * t567 - t551 * t566 + t559;
t494 = -t531 * t567 - t551 * t594 - t596;
t491 = (t531 * t575 + t532 * t578) * t595;
t486 = qJD(6) * t518 + t510;
t484 = rSges(4,1) * t522 - rSges(4,2) * t521 + rSges(4,3) * t543;
t483 = rSges(4,1) * t520 - rSges(4,2) * t519 + rSges(4,3) * t542;
t482 = Icges(4,1) * t522 - Icges(4,4) * t521 + Icges(4,5) * t543;
t481 = Icges(4,1) * t520 - Icges(4,4) * t519 + Icges(4,5) * t542;
t480 = Icges(4,4) * t522 - Icges(4,2) * t521 + Icges(4,6) * t543;
t479 = Icges(4,4) * t520 - Icges(4,2) * t519 + Icges(4,6) * t542;
t478 = Icges(4,5) * t522 - Icges(4,6) * t521 + Icges(4,3) * t543;
t477 = Icges(4,5) * t520 - Icges(4,6) * t519 + Icges(4,3) * t542;
t476 = rSges(5,1) * t518 - rSges(5,2) * t517 + rSges(5,3) * t540;
t475 = rSges(6,1) * t540 - rSges(6,2) * t518 + rSges(6,3) * t517;
t468 = pkin(5) * t521 + pkin(12) * t503;
t467 = pkin(5) * t519 + pkin(12) * t501;
t466 = t502 * t571 + t521 * t576;
t465 = t502 * t576 - t521 * t571;
t464 = t500 * t571 + t519 * t576;
t463 = t500 * t576 - t519 * t571;
t459 = qJD(6) * t501 + t493;
t458 = qJD(6) * t503 + t492;
t456 = rSges(5,1) * t503 - rSges(5,2) * t502 + rSges(5,3) * t521;
t455 = rSges(5,1) * t501 - rSges(5,2) * t500 + rSges(5,3) * t519;
t454 = rSges(6,1) * t521 - rSges(6,2) * t503 + rSges(6,3) * t502;
t453 = rSges(6,1) * t519 - rSges(6,2) * t501 + rSges(6,3) * t500;
t439 = rSges(7,1) * t497 + rSges(7,2) * t496 + rSges(7,3) * t518;
t438 = Icges(7,1) * t497 + Icges(7,4) * t496 + Icges(7,5) * t518;
t437 = Icges(7,4) * t497 + Icges(7,2) * t496 + Icges(7,6) * t518;
t436 = Icges(7,5) * t497 + Icges(7,6) * t496 + Icges(7,3) * t518;
t434 = rSges(7,1) * t466 + rSges(7,2) * t465 + rSges(7,3) * t503;
t433 = rSges(7,1) * t464 + rSges(7,2) * t463 + rSges(7,3) * t501;
t432 = t484 * t544 - t508 * t533 + t586;
t431 = -t483 * t544 + t508 * t534 + t584;
t430 = Icges(7,1) * t466 + Icges(7,4) * t465 + Icges(7,5) * t503;
t429 = Icges(7,1) * t464 + Icges(7,4) * t463 + Icges(7,5) * t501;
t428 = Icges(7,4) * t466 + Icges(7,2) * t465 + Icges(7,6) * t503;
t427 = Icges(7,4) * t464 + Icges(7,2) * t463 + Icges(7,6) * t501;
t426 = Icges(7,5) * t466 + Icges(7,6) * t465 + Icges(7,3) * t503;
t425 = Icges(7,5) * t464 + Icges(7,6) * t463 + Icges(7,3) * t501;
t424 = t483 * t533 - t484 * t534 + t597;
t423 = t456 * t510 - t476 * t492 + t583;
t422 = -t455 * t510 + t476 * t493 + t581;
t421 = t455 * t492 - t456 * t493 + t587;
t420 = t454 * t510 + (-t475 - t487) * t492 + t582;
t419 = t475 * t493 + (-t453 - t460) * t510 + t580;
t418 = t453 * t492 + (-t454 - t461) * t493 + t585;
t417 = t434 * t486 - t439 * t458 + t468 * t510 + (-t487 - t504) * t492 + t582;
t416 = -t433 * t486 + t439 * t459 + t493 * t504 + (-t460 - t467) * t510 + t580;
t415 = t433 * t458 - t434 * t459 + t467 * t492 + (-t461 - t468) * t493 + t585;
t1 = ((t548 * t603 + t549 * t557 + t550 * t558) * t567 + (-(t527 * t557 + t529 * t558) * t578 + (t557 * t528 + t558 * t530 - t622) * t575) * t595) * t566 / 0.2e1 - ((-t548 * t602 + t549 * t555 + t550 * t556) * t567 + ((t528 * t555 + t530 * t556) * t575 + (-t555 * t527 - t556 * t529 + t622) * t578) * t595) * t594 / 0.2e1 + t534 * ((t478 * t542 - t480 * t519 + t482 * t520) * t533 + (t542 * t477 - t519 * t479 + t520 * t481) * t534 + (t505 * t542 - t506 * t519 + t507 * t520) * t544) / 0.2e1 + t544 * ((t478 * t554 - t480 * t540 + t482 * t541) * t533 + (t477 * t554 - t479 * t540 + t481 * t541) * t534 + (t554 * t505 - t540 * t506 + t541 * t507) * t544) / 0.2e1 + t533 * ((t543 * t478 - t521 * t480 + t522 * t482) * t533 + (t477 * t543 - t479 * t521 + t481 * t522) * t534 + (t505 * t543 - t506 * t521 + t507 * t522) * t544) / 0.2e1 + t567 * ((t570 * t548 + (t549 * t577 + t550 * t574) * t569) * t567 + (((t528 * t577 + t530 * t574) * t575 - (t527 * t577 + t529 * t574) * t578) * t569 - t589 * t570) * t595) / 0.2e1 + m(7) * (t415 ^ 2 + t416 ^ 2 + t417 ^ 2) / 0.2e1 + m(6) * (t418 ^ 2 + t419 ^ 2 + t420 ^ 2) / 0.2e1 + m(5) * (t421 ^ 2 + t422 ^ 2 + t423 ^ 2) / 0.2e1 + m(4) * (t424 ^ 2 + t431 ^ 2 + t432 ^ 2) / 0.2e1 + m(3) * (t491 ^ 2 + t494 ^ 2 + t495 ^ 2) / 0.2e1 + t459 * ((t426 * t501 + t428 * t463 + t430 * t464) * t458 + (t425 * t501 + t427 * t463 + t429 * t464) * t459 + (t436 * t501 + t437 * t463 + t438 * t464) * t486) / 0.2e1 + t486 * ((t426 * t518 + t428 * t496 + t430 * t497) * t458 + (t425 * t518 + t427 * t496 + t429 * t497) * t459 + (t436 * t518 + t496 * t437 + t497 * t438) * t486) / 0.2e1 + t458 * ((t426 * t503 + t428 * t465 + t430 * t466) * t458 + (t425 * t503 + t427 * t465 + t429 * t466) * t459 + (t436 * t503 + t437 * t465 + t438 * t466) * t486) / 0.2e1 + ((t502 * t615 + t503 * t613 + t521 * t614) * t510 + (t502 * t621 + t503 * t617 + t521 * t619) * t493 + (t502 * t620 + t503 * t616 + t521 * t618) * t492) * t492 / 0.2e1 + ((t500 * t615 + t501 * t613 + t519 * t614) * t510 + (t500 * t621 + t501 * t617 + t519 * t619) * t493 + (t500 * t620 + t501 * t616 + t519 * t618) * t492) * t493 / 0.2e1 + ((t517 * t615 + t518 * t613 + t540 * t614) * t510 + (t517 * t621 + t518 * t617 + t540 * t619) * t493 + (t517 * t620 + t518 * t616 + t540 * t618) * t492) * t510 / 0.2e1 + (m(2) * (t563 ^ 2 + t564 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
