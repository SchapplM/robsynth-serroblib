% Calculate kinetic energy for
% S6RPRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_energykin_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR9_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR9_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR9_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:26:49
% EndTime: 2019-03-09 05:26:52
% DurationCPUTime: 2.74s
% Computational Cost: add. (4916->343), mult. (12087->516), div. (0->0), fcn. (15605->16), ass. (0->158)
t621 = Icges(5,3) + Icges(6,3);
t566 = cos(pkin(6));
t607 = cos(pkin(12));
t611 = sin(qJ(1));
t591 = t611 * t607;
t573 = cos(qJ(1));
t605 = sin(pkin(12));
t594 = t573 * t605;
t552 = t566 * t594 + t591;
t570 = sin(qJ(3));
t590 = t611 * t605;
t595 = t573 * t607;
t580 = -t566 * t595 + t590;
t608 = cos(pkin(7));
t575 = t580 * t608;
t565 = sin(pkin(6));
t606 = sin(pkin(7));
t597 = t565 * t606;
t612 = cos(qJ(3));
t530 = t552 * t612 + (-t573 * t597 - t575) * t570;
t598 = t565 * t608;
t544 = -t573 * t598 + t580 * t606;
t600 = qJ(4) + pkin(13);
t564 = sin(t600);
t593 = cos(t600);
t504 = t530 * t564 - t544 * t593;
t505 = t530 * t593 + t544 * t564;
t569 = sin(qJ(4));
t572 = cos(qJ(4));
t508 = -t530 * t569 + t544 * t572;
t604 = t544 * t569;
t509 = t530 * t572 + t604;
t592 = t612 * t606;
t601 = t565 * t573;
t529 = t552 * t570 + t575 * t612 + t592 * t601;
t620 = Icges(5,5) * t509 + Icges(6,5) * t505 + Icges(5,6) * t508 - Icges(6,6) * t504 + t621 * t529;
t553 = -t566 * t590 + t595;
t581 = t566 * t591 + t594;
t617 = t581 * t608 - t611 * t597;
t532 = t553 * t612 - t570 * t617;
t545 = t581 * t606 + t598 * t611;
t506 = t532 * t564 - t545 * t593;
t507 = t532 * t593 + t545 * t564;
t510 = -t532 * t569 + t545 * t572;
t603 = t545 * t569;
t511 = t532 * t572 + t603;
t531 = t553 * t570 + t612 * t617;
t619 = Icges(5,5) * t511 + Icges(6,5) * t507 + Icges(5,6) * t510 - Icges(6,6) * t506 + t621 * t531;
t588 = t608 * t607;
t543 = t566 * t606 * t570 + (t570 * t588 + t605 * t612) * t565;
t551 = t566 * t608 - t597 * t607;
t521 = t543 * t564 - t551 * t593;
t522 = t543 * t593 + t551 * t564;
t527 = -t543 * t569 + t551 * t572;
t602 = t551 * t569;
t528 = t543 * t572 + t602;
t596 = t565 * t605;
t542 = -t565 * t588 * t612 - t566 * t592 + t570 * t596;
t618 = Icges(5,5) * t528 + Icges(6,5) * t522 + Icges(5,6) * t527 - Icges(6,6) * t521 + t621 * t542;
t610 = pkin(4) * t572;
t540 = qJD(3) * t544;
t512 = qJD(4) * t529 + t540;
t541 = qJD(3) * t545;
t513 = qJD(4) * t531 + t541;
t549 = qJD(3) * t551 + qJD(1);
t599 = t565 * t611;
t533 = qJD(4) * t542 + t549;
t589 = -qJD(2) * t601 + qJD(1) * (t573 * pkin(1) + qJ(2) * t599);
t555 = pkin(1) * t611 - qJ(2) * t601;
t561 = qJD(2) * t599;
t587 = t561 + (-t552 * pkin(2) - pkin(9) * t544 - t555) * qJD(1);
t585 = qJD(1) * (t553 * pkin(2) + pkin(9) * t545) + t589;
t499 = pkin(3) * t530 + pkin(10) * t529;
t500 = pkin(3) * t532 + pkin(10) * t531;
t563 = qJD(2) * t566;
t584 = t499 * t541 - t500 * t540 + t563;
t447 = pkin(4) * t604 + qJ(5) * t529 + t530 * t610;
t583 = qJD(5) * t542 + t513 * t447 + t584;
t520 = pkin(3) * t543 + pkin(10) * t542;
t582 = -t499 * t549 + t520 * t540 + t587;
t579 = t549 * t500 - t520 * t541 + t585;
t474 = pkin(4) * t602 + qJ(5) * t542 + t543 * t610;
t578 = qJD(5) * t531 + t512 * t474 + t582;
t448 = pkin(4) * t603 + qJ(5) * t531 + t532 * t610;
t577 = qJD(5) * t529 + t533 * t448 + t579;
t571 = cos(qJ(6));
t568 = sin(qJ(6));
t559 = t573 * rSges(2,1) - rSges(2,2) * t611;
t558 = rSges(2,1) * t611 + t573 * rSges(2,2);
t519 = qJD(1) * (t553 * rSges(3,1) - rSges(3,2) * t581 + rSges(3,3) * t599) + t589;
t518 = t561 + (-t552 * rSges(3,1) + rSges(3,2) * t580 + rSges(3,3) * t601 - t555) * qJD(1);
t517 = rSges(4,1) * t543 - rSges(4,2) * t542 + rSges(4,3) * t551;
t516 = Icges(4,1) * t543 - Icges(4,4) * t542 + Icges(4,5) * t551;
t515 = Icges(4,4) * t543 - Icges(4,2) * t542 + Icges(4,6) * t551;
t514 = Icges(4,5) * t543 - Icges(4,6) * t542 + Icges(4,3) * t551;
t502 = t522 * t571 + t542 * t568;
t501 = -t522 * t568 + t542 * t571;
t498 = qJD(6) * t521 + t533;
t497 = pkin(5) * t522 + pkin(11) * t521;
t494 = rSges(4,1) * t532 - rSges(4,2) * t531 + rSges(4,3) * t545;
t493 = rSges(4,1) * t530 - rSges(4,2) * t529 + rSges(4,3) * t544;
t492 = Icges(4,1) * t532 - Icges(4,4) * t531 + Icges(4,5) * t545;
t491 = Icges(4,1) * t530 - Icges(4,4) * t529 + Icges(4,5) * t544;
t490 = Icges(4,4) * t532 - Icges(4,2) * t531 + Icges(4,6) * t545;
t489 = Icges(4,4) * t530 - Icges(4,2) * t529 + Icges(4,6) * t544;
t488 = Icges(4,5) * t532 - Icges(4,6) * t531 + Icges(4,3) * t545;
t487 = Icges(4,5) * t530 - Icges(4,6) * t529 + Icges(4,3) * t544;
t486 = rSges(5,1) * t528 + rSges(5,2) * t527 + rSges(5,3) * t542;
t485 = Icges(5,1) * t528 + Icges(5,4) * t527 + Icges(5,5) * t542;
t484 = Icges(5,4) * t528 + Icges(5,2) * t527 + Icges(5,6) * t542;
t482 = rSges(6,1) * t522 - rSges(6,2) * t521 + rSges(6,3) * t542;
t481 = t507 * t571 + t531 * t568;
t480 = -t507 * t568 + t531 * t571;
t479 = t505 * t571 + t529 * t568;
t478 = -t505 * t568 + t529 * t571;
t477 = Icges(6,1) * t522 - Icges(6,4) * t521 + Icges(6,5) * t542;
t476 = Icges(6,4) * t522 - Icges(6,2) * t521 + Icges(6,6) * t542;
t473 = qJD(6) * t506 + t513;
t472 = qJD(6) * t504 + t512;
t471 = pkin(5) * t507 + pkin(11) * t506;
t470 = pkin(5) * t505 + pkin(11) * t504;
t469 = rSges(5,1) * t511 + rSges(5,2) * t510 + rSges(5,3) * t531;
t468 = rSges(5,1) * t509 + rSges(5,2) * t508 + rSges(5,3) * t529;
t467 = Icges(5,1) * t511 + Icges(5,4) * t510 + Icges(5,5) * t531;
t466 = Icges(5,1) * t509 + Icges(5,4) * t508 + Icges(5,5) * t529;
t465 = Icges(5,4) * t511 + Icges(5,2) * t510 + Icges(5,6) * t531;
t464 = Icges(5,4) * t509 + Icges(5,2) * t508 + Icges(5,6) * t529;
t460 = rSges(6,1) * t507 - rSges(6,2) * t506 + rSges(6,3) * t531;
t459 = rSges(6,1) * t505 - rSges(6,2) * t504 + rSges(6,3) * t529;
t458 = Icges(6,1) * t507 - Icges(6,4) * t506 + Icges(6,5) * t531;
t457 = Icges(6,1) * t505 - Icges(6,4) * t504 + Icges(6,5) * t529;
t456 = Icges(6,4) * t507 - Icges(6,2) * t506 + Icges(6,6) * t531;
t455 = Icges(6,4) * t505 - Icges(6,2) * t504 + Icges(6,6) * t529;
t452 = rSges(7,1) * t502 + rSges(7,2) * t501 + rSges(7,3) * t521;
t451 = Icges(7,1) * t502 + Icges(7,4) * t501 + Icges(7,5) * t521;
t450 = Icges(7,4) * t502 + Icges(7,2) * t501 + Icges(7,6) * t521;
t449 = Icges(7,5) * t502 + Icges(7,6) * t501 + Icges(7,3) * t521;
t444 = t494 * t549 - t517 * t541 + t585;
t443 = -t493 * t549 + t517 * t540 + t587;
t442 = t563 + (t493 * t545 - t494 * t544) * qJD(3);
t441 = rSges(7,1) * t481 + rSges(7,2) * t480 + rSges(7,3) * t506;
t440 = rSges(7,1) * t479 + rSges(7,2) * t478 + rSges(7,3) * t504;
t439 = Icges(7,1) * t481 + Icges(7,4) * t480 + Icges(7,5) * t506;
t438 = Icges(7,1) * t479 + Icges(7,4) * t478 + Icges(7,5) * t504;
t437 = Icges(7,4) * t481 + Icges(7,2) * t480 + Icges(7,6) * t506;
t436 = Icges(7,4) * t479 + Icges(7,2) * t478 + Icges(7,6) * t504;
t435 = Icges(7,5) * t481 + Icges(7,6) * t480 + Icges(7,3) * t506;
t434 = Icges(7,5) * t479 + Icges(7,6) * t478 + Icges(7,3) * t504;
t433 = t469 * t533 - t486 * t513 + t579;
t432 = -t468 * t533 + t486 * t512 + t582;
t431 = t468 * t513 - t469 * t512 + t584;
t430 = t460 * t533 + (-t474 - t482) * t513 + t577;
t429 = t482 * t512 + (-t447 - t459) * t533 + t578;
t428 = t459 * t513 + (-t448 - t460) * t512 + t583;
t427 = t441 * t498 - t452 * t473 + t471 * t533 + (-t474 - t497) * t513 + t577;
t426 = -t440 * t498 + t452 * t472 + t497 * t512 + (-t447 - t470) * t533 + t578;
t425 = t440 * t473 - t441 * t472 + t470 * t513 + (-t448 - t471) * t512 + t583;
t1 = ((t514 * t544 - t515 * t529 + t516 * t530) * t549 + ((t488 * t544 - t490 * t529 + t492 * t530) * t545 + (t487 * t544 - t489 * t529 + t491 * t530) * t544) * qJD(3)) * t540 / 0.2e1 + t549 * ((t551 * t514 - t542 * t515 + t543 * t516) * t549 + ((t488 * t551 - t490 * t542 + t492 * t543) * t545 + (t487 * t551 - t489 * t542 + t491 * t543) * t544) * qJD(3)) / 0.2e1 + m(6) * (t428 ^ 2 + t429 ^ 2 + t430 ^ 2) / 0.2e1 + m(7) * (t425 ^ 2 + t426 ^ 2 + t427 ^ 2) / 0.2e1 + m(5) * (t431 ^ 2 + t432 ^ 2 + t433 ^ 2) / 0.2e1 + m(4) * (t442 ^ 2 + t443 ^ 2 + t444 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 * t566 ^ 2 + t518 ^ 2 + t519 ^ 2) / 0.2e1 + t473 * ((t506 * t435 + t480 * t437 + t481 * t439) * t473 + (t434 * t506 + t436 * t480 + t438 * t481) * t472 + (t449 * t506 + t450 * t480 + t451 * t481) * t498) / 0.2e1 + t498 * ((t435 * t521 + t437 * t501 + t439 * t502) * t473 + (t434 * t521 + t436 * t501 + t438 * t502) * t472 + (t521 * t449 + t501 * t450 + t502 * t451) * t498) / 0.2e1 + t472 * ((t435 * t504 + t437 * t478 + t439 * t479) * t473 + (t504 * t434 + t478 * t436 + t479 * t438) * t472 + (t449 * t504 + t450 * t478 + t451 * t479) * t498) / 0.2e1 + ((t514 * t545 - t515 * t531 + t516 * t532) * t549 + ((t488 * t545 - t490 * t531 + t492 * t532) * t545 + (t487 * t545 - t489 * t531 + t491 * t532) * t544) * qJD(3)) * t541 / 0.2e1 + ((-t476 * t504 + t477 * t505 + t484 * t508 + t485 * t509 + t529 * t618) * t533 + (-t456 * t504 + t458 * t505 + t465 * t508 + t467 * t509 + t529 * t619) * t513 + (-t504 * t455 + t505 * t457 + t508 * t464 + t509 * t466 + t529 * t620) * t512) * t512 / 0.2e1 + ((-t476 * t506 + t477 * t507 + t484 * t510 + t485 * t511 + t531 * t618) * t533 + (-t506 * t456 + t507 * t458 + t510 * t465 + t511 * t467 + t531 * t619) * t513 + (-t455 * t506 + t457 * t507 + t464 * t510 + t466 * t511 + t531 * t620) * t512) * t513 / 0.2e1 + ((-t521 * t476 + t522 * t477 + t527 * t484 + t528 * t485 + t542 * t618) * t533 + (-t456 * t521 + t458 * t522 + t465 * t527 + t467 * t528 + t542 * t619) * t513 + (-t455 * t521 + t457 * t522 + t464 * t527 + t466 * t528 + t542 * t620) * t512) * t533 / 0.2e1 + ((Icges(3,5) * t566 + (Icges(3,1) * t605 + Icges(3,4) * t607) * t565) * t596 + t565 * t607 * (Icges(3,6) * t566 + (Icges(3,4) * t605 + Icges(3,2) * t607) * t565) + t566 * (Icges(3,3) * t566 + (Icges(3,5) * t605 + Icges(3,6) * t607) * t565) + Icges(2,3) + m(2) * (t558 ^ 2 + t559 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
