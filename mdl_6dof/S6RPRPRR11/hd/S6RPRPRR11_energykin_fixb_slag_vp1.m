% Calculate kinetic energy for
% S6RPRPRR11
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
% Datum: 2019-03-09 04:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR11_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR11_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_energykin_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR11_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR11_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR11_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:11:16
% EndTime: 2019-03-09 04:11:19
% DurationCPUTime: 2.79s
% Computational Cost: add. (4799->343), mult. (11781->516), div. (0->0), fcn. (15209->16), ass. (0->159)
t622 = Icges(4,2) + Icges(5,3);
t568 = cos(pkin(6));
t608 = sin(pkin(12));
t613 = sin(qJ(1));
t589 = t613 * t608;
t573 = cos(qJ(1));
t610 = cos(pkin(12));
t594 = t573 * t610;
t553 = -t568 * t589 + t594;
t571 = sin(qJ(3));
t614 = cos(qJ(3));
t590 = t613 * t610;
t593 = t573 * t608;
t580 = t568 * t590 + t593;
t566 = sin(pkin(6));
t609 = sin(pkin(7));
t596 = t566 * t609;
t611 = cos(pkin(7));
t617 = t580 * t611 - t613 * t596;
t531 = t553 * t571 + t617 * t614;
t532 = t553 * t614 - t617 * t571;
t597 = t566 * t611;
t545 = t580 * t609 + t597 * t613;
t488 = Icges(4,5) * t532 - Icges(4,6) * t531 + Icges(4,3) * t545;
t621 = t488 * t545;
t552 = t568 * t593 + t590;
t579 = -t568 * t594 + t589;
t576 = t579 * t611;
t530 = t552 * t614 + (-t573 * t596 - t576) * t571;
t544 = -t573 * t597 + t579 * t609;
t565 = sin(pkin(13));
t567 = cos(pkin(13));
t508 = -t530 * t565 + t544 * t567;
t607 = t544 * t565;
t509 = t530 * t567 + t607;
t591 = t614 * t609;
t604 = t566 * t573;
t529 = t552 * t571 + t576 * t614 + t591 * t604;
t620 = -Icges(4,4) * t530 + Icges(5,5) * t509 - Icges(4,6) * t544 + Icges(5,6) * t508 + t529 * t622;
t510 = -t532 * t565 + t545 * t567;
t606 = t545 * t565;
t511 = t532 * t567 + t606;
t619 = -Icges(4,4) * t532 + Icges(5,5) * t511 - Icges(4,6) * t545 + Icges(5,6) * t510 + t531 * t622;
t587 = t611 * t610;
t543 = t568 * t609 * t571 + (t571 * t587 + t608 * t614) * t566;
t551 = t568 * t611 - t596 * t610;
t527 = -t543 * t565 + t551 * t567;
t605 = t551 * t565;
t528 = t543 * t567 + t605;
t595 = t566 * t608;
t542 = -t566 * t587 * t614 - t568 * t591 + t571 * t595;
t618 = -Icges(4,4) * t543 + Icges(5,5) * t528 - Icges(4,6) * t551 + Icges(5,6) * t527 + t542 * t622;
t612 = pkin(4) * t567;
t540 = qJD(3) * t544;
t512 = qJD(5) * t529 + t540;
t541 = qJD(3) * t545;
t513 = qJD(5) * t531 + t541;
t549 = qJD(3) * t551 + qJD(1);
t602 = pkin(13) + qJ(5);
t499 = pkin(3) * t530 + qJ(4) * t529;
t563 = qJD(2) * t568;
t601 = qJD(4) * t542 + t499 * t541 + t563;
t600 = t566 * t613;
t533 = qJD(5) * t542 + t549;
t592 = cos(t602);
t588 = -qJD(2) * t604 + qJD(1) * (t573 * pkin(1) + qJ(2) * t600);
t555 = pkin(1) * t613 - qJ(2) * t604;
t561 = qJD(2) * t600;
t586 = t561 + (-t552 * pkin(2) - pkin(9) * t544 - t555) * qJD(1);
t584 = qJD(1) * (t553 * pkin(2) + pkin(9) * t545) + t588;
t520 = pkin(3) * t543 + qJ(4) * t542;
t583 = qJD(4) * t531 + t520 * t540 + t586;
t500 = pkin(3) * t532 + qJ(4) * t531;
t582 = qJD(4) * t529 + t549 * t500 + t584;
t447 = pkin(4) * t607 + pkin(10) * t529 + t530 * t612;
t448 = pkin(4) * t606 + pkin(10) * t531 + t532 * t612;
t581 = t447 * t541 + (-t448 - t500) * t540 + t601;
t474 = pkin(4) * t605 + pkin(10) * t542 + t543 * t612;
t578 = t474 * t540 + (-t447 - t499) * t549 + t583;
t575 = t549 * t448 + (-t474 - t520) * t541 + t582;
t572 = cos(qJ(6));
t570 = sin(qJ(6));
t564 = sin(t602);
t559 = t573 * rSges(2,1) - rSges(2,2) * t613;
t558 = rSges(2,1) * t613 + t573 * rSges(2,2);
t522 = t543 * t592 + t551 * t564;
t521 = t543 * t564 - t551 * t592;
t519 = qJD(1) * (t553 * rSges(3,1) - rSges(3,2) * t580 + rSges(3,3) * t600) + t588;
t518 = t561 + (-t552 * rSges(3,1) + rSges(3,2) * t579 + rSges(3,3) * t604 - t555) * qJD(1);
t517 = rSges(4,1) * t543 - rSges(4,2) * t542 + rSges(4,3) * t551;
t516 = Icges(4,1) * t543 - Icges(4,4) * t542 + Icges(4,5) * t551;
t514 = Icges(4,5) * t543 - Icges(4,6) * t542 + Icges(4,3) * t551;
t507 = t532 * t592 + t545 * t564;
t506 = t532 * t564 - t545 * t592;
t505 = t530 * t592 + t544 * t564;
t504 = t530 * t564 - t544 * t592;
t502 = t522 * t572 + t542 * t570;
t501 = -t522 * t570 + t542 * t572;
t498 = qJD(6) * t521 + t533;
t497 = pkin(5) * t522 + pkin(11) * t521;
t495 = rSges(4,1) * t532 - rSges(4,2) * t531 + rSges(4,3) * t545;
t494 = rSges(4,1) * t530 - rSges(4,2) * t529 + rSges(4,3) * t544;
t492 = Icges(4,1) * t532 - Icges(4,4) * t531 + Icges(4,5) * t545;
t491 = Icges(4,1) * t530 - Icges(4,4) * t529 + Icges(4,5) * t544;
t487 = Icges(4,5) * t530 - Icges(4,6) * t529 + Icges(4,3) * t544;
t486 = rSges(5,1) * t528 + rSges(5,2) * t527 + rSges(5,3) * t542;
t485 = Icges(5,1) * t528 + Icges(5,4) * t527 + Icges(5,5) * t542;
t484 = Icges(5,4) * t528 + Icges(5,2) * t527 + Icges(5,6) * t542;
t482 = rSges(6,1) * t522 - rSges(6,2) * t521 + rSges(6,3) * t542;
t481 = t507 * t572 + t531 * t570;
t480 = -t507 * t570 + t531 * t572;
t479 = t505 * t572 + t529 * t570;
t478 = -t505 * t570 + t529 * t572;
t477 = Icges(6,1) * t522 - Icges(6,4) * t521 + Icges(6,5) * t542;
t476 = Icges(6,4) * t522 - Icges(6,2) * t521 + Icges(6,6) * t542;
t475 = Icges(6,5) * t522 - Icges(6,6) * t521 + Icges(6,3) * t542;
t473 = qJD(6) * t506 + t513;
t472 = qJD(6) * t504 + t512;
t470 = pkin(5) * t507 + pkin(11) * t506;
t469 = pkin(5) * t505 + pkin(11) * t504;
t468 = rSges(5,1) * t511 + rSges(5,2) * t510 + rSges(5,3) * t531;
t467 = rSges(5,1) * t509 + rSges(5,2) * t508 + rSges(5,3) * t529;
t466 = Icges(5,1) * t511 + Icges(5,4) * t510 + Icges(5,5) * t531;
t465 = Icges(5,1) * t509 + Icges(5,4) * t508 + Icges(5,5) * t529;
t464 = Icges(5,4) * t511 + Icges(5,2) * t510 + Icges(5,6) * t531;
t463 = Icges(5,4) * t509 + Icges(5,2) * t508 + Icges(5,6) * t529;
t460 = rSges(6,1) * t507 - rSges(6,2) * t506 + rSges(6,3) * t531;
t459 = rSges(6,1) * t505 - rSges(6,2) * t504 + rSges(6,3) * t529;
t458 = Icges(6,1) * t507 - Icges(6,4) * t506 + Icges(6,5) * t531;
t457 = Icges(6,1) * t505 - Icges(6,4) * t504 + Icges(6,5) * t529;
t456 = Icges(6,4) * t507 - Icges(6,2) * t506 + Icges(6,6) * t531;
t455 = Icges(6,4) * t505 - Icges(6,2) * t504 + Icges(6,6) * t529;
t454 = Icges(6,5) * t507 - Icges(6,6) * t506 + Icges(6,3) * t531;
t453 = Icges(6,5) * t505 - Icges(6,6) * t504 + Icges(6,3) * t529;
t452 = rSges(7,1) * t502 + rSges(7,2) * t501 + rSges(7,3) * t521;
t451 = Icges(7,1) * t502 + Icges(7,4) * t501 + Icges(7,5) * t521;
t450 = Icges(7,4) * t502 + Icges(7,2) * t501 + Icges(7,6) * t521;
t449 = Icges(7,5) * t502 + Icges(7,6) * t501 + Icges(7,3) * t521;
t444 = t495 * t549 - t517 * t541 + t584;
t443 = -t494 * t549 + t517 * t540 + t586;
t442 = t563 + (t494 * t545 - t495 * t544) * qJD(3);
t441 = rSges(7,1) * t481 + rSges(7,2) * t480 + rSges(7,3) * t506;
t440 = rSges(7,1) * t479 + rSges(7,2) * t478 + rSges(7,3) * t504;
t439 = Icges(7,1) * t481 + Icges(7,4) * t480 + Icges(7,5) * t506;
t438 = Icges(7,1) * t479 + Icges(7,4) * t478 + Icges(7,5) * t504;
t437 = Icges(7,4) * t481 + Icges(7,2) * t480 + Icges(7,6) * t506;
t436 = Icges(7,4) * t479 + Icges(7,2) * t478 + Icges(7,6) * t504;
t435 = Icges(7,5) * t481 + Icges(7,6) * t480 + Icges(7,3) * t506;
t434 = Icges(7,5) * t479 + Icges(7,6) * t478 + Icges(7,3) * t504;
t433 = t468 * t549 + (-t486 - t520) * t541 + t582;
t432 = t486 * t540 + (-t467 - t499) * t549 + t583;
t431 = (t467 * t545 + (-t468 - t500) * t544) * qJD(3) + t601;
t430 = t460 * t533 - t482 * t513 + t575;
t429 = -t459 * t533 + t482 * t512 + t578;
t428 = t459 * t513 - t460 * t512 + t581;
t427 = t441 * t498 - t452 * t473 + t470 * t533 - t497 * t513 + t575;
t426 = -t440 * t498 + t452 * t472 - t469 * t533 + t497 * t512 + t578;
t425 = t440 * t473 - t441 * t472 + t469 * t513 - t470 * t512 + t581;
t1 = m(7) * (t425 ^ 2 + t426 ^ 2 + t427 ^ 2) / 0.2e1 + m(6) * (t428 ^ 2 + t429 ^ 2 + t430 ^ 2) / 0.2e1 + m(5) * (t431 ^ 2 + t432 ^ 2 + t433 ^ 2) / 0.2e1 + m(4) * (t442 ^ 2 + t443 ^ 2 + t444 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 * t568 ^ 2 + t518 ^ 2 + t519 ^ 2) / 0.2e1 + t498 * ((t435 * t521 + t437 * t501 + t439 * t502) * t473 + (t434 * t521 + t436 * t501 + t438 * t502) * t472 + (t449 * t521 + t450 * t501 + t451 * t502) * t498) / 0.2e1 + t533 * ((t454 * t542 - t456 * t521 + t458 * t522) * t513 + (t453 * t542 - t455 * t521 + t457 * t522) * t512 + (t542 * t475 - t521 * t476 + t522 * t477) * t533) / 0.2e1 + t473 * ((t435 * t506 + t437 * t480 + t439 * t481) * t473 + (t434 * t506 + t436 * t480 + t438 * t481) * t472 + (t449 * t506 + t450 * t480 + t451 * t481) * t498) / 0.2e1 + t472 * ((t435 * t504 + t437 * t478 + t439 * t479) * t473 + (t434 * t504 + t436 * t478 + t438 * t479) * t472 + (t449 * t504 + t450 * t478 + t451 * t479) * t498) / 0.2e1 + t513 * ((t454 * t531 - t456 * t506 + t458 * t507) * t513 + (t453 * t531 - t455 * t506 + t457 * t507) * t512 + (t475 * t531 - t476 * t506 + t477 * t507) * t533) / 0.2e1 + t512 * ((t454 * t529 - t456 * t504 + t458 * t505) * t513 + (t453 * t529 - t455 * t504 + t457 * t505) * t512 + (t475 * t529 - t476 * t504 + t477 * t505) * t533) / 0.2e1 + ((t484 * t527 + t485 * t528 + t514 * t551 + t516 * t543 + t618 * t542) * t549 + ((t464 * t527 + t466 * t528 + t488 * t551 + t492 * t543 + t619 * t542) * t545 + (t463 * t527 + t465 * t528 + t487 * t551 + t491 * t543 + t620 * t542) * t544) * qJD(3)) * t549 / 0.2e1 + ((t484 * t508 + t485 * t509 + t514 * t544 + t516 * t530 + t618 * t529) * t549 + ((t464 * t508 + t466 * t509 + t492 * t530 + t619 * t529) * t545 + (t463 * t508 + t465 * t509 + t487 * t544 + t491 * t530 + t620 * t529 + t621) * t544) * qJD(3)) * t540 / 0.2e1 + ((t484 * t510 + t485 * t511 + t514 * t545 + t516 * t532 + t618 * t531) * t549 + ((t464 * t510 + t466 * t511 + t492 * t532 + t619 * t531 + t621) * t545 + (t463 * t510 + t465 * t511 + t487 * t545 + t491 * t532 + t620 * t531) * t544) * qJD(3)) * t541 / 0.2e1 + ((Icges(3,5) * t568 + (Icges(3,1) * t608 + Icges(3,4) * t610) * t566) * t595 + t566 * t610 * (Icges(3,6) * t568 + (Icges(3,4) * t608 + Icges(3,2) * t610) * t566) + t568 * (Icges(3,3) * t568 + (Icges(3,5) * t608 + Icges(3,6) * t610) * t566) + Icges(2,3) + m(2) * (t558 ^ 2 + t559 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
