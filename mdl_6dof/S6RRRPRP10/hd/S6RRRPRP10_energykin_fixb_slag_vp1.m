% Calculate kinetic energy for
% S6RRRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 17:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP10_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP10_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP10_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP10_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:29:27
% EndTime: 2019-03-09 17:29:29
% DurationCPUTime: 2.86s
% Computational Cost: add. (2938->327), mult. (6398->479), div. (0->0), fcn. (7848->12), ass. (0->151)
t604 = Icges(6,1) + Icges(7,1);
t603 = -Icges(6,4) + Icges(7,5);
t602 = Icges(7,4) + Icges(6,5);
t601 = Icges(4,2) + Icges(5,3);
t600 = Icges(6,2) + Icges(7,3);
t599 = Icges(7,2) + Icges(6,3);
t598 = -Icges(6,6) + Icges(7,6);
t597 = rSges(7,1) + pkin(5);
t596 = rSges(7,3) + qJ(6);
t540 = sin(qJ(2));
t541 = sin(qJ(1));
t542 = cos(qJ(2));
t543 = cos(qJ(1));
t573 = cos(pkin(6));
t555 = t543 * t573;
t517 = t540 * t541 - t542 * t555;
t518 = t540 * t555 + t541 * t542;
t536 = sin(pkin(6));
t568 = t536 * t543;
t479 = Icges(3,5) * t518 - Icges(3,6) * t517 - Icges(3,3) * t568;
t556 = t541 * t573;
t519 = t543 * t540 + t542 * t556;
t520 = -t540 * t556 + t543 * t542;
t570 = t536 * t541;
t480 = Icges(3,5) * t520 - Icges(3,6) * t519 + Icges(3,3) * t570;
t595 = t536 * (t479 * t543 - t480 * t541);
t539 = sin(qJ(3));
t575 = cos(qJ(3));
t501 = t518 * t575 - t539 * t568;
t560 = pkin(11) + qJ(5);
t534 = sin(t560);
t554 = cos(t560);
t467 = t501 * t534 - t517 * t554;
t468 = t501 * t554 + t517 * t534;
t558 = t536 * t575;
t500 = t518 * t539 + t543 * t558;
t594 = t600 * t467 + t603 * t468 + t598 * t500;
t503 = t520 * t575 + t539 * t570;
t469 = t503 * t534 - t519 * t554;
t470 = t503 * t554 + t519 * t534;
t502 = t520 * t539 - t541 * t558;
t593 = t600 * t469 + t603 * t470 + t598 * t502;
t592 = t598 * t467 + t602 * t468 + t599 * t500;
t591 = t598 * t469 + t602 * t470 + t599 * t502;
t590 = t603 * t467 + t604 * t468 + t602 * t500;
t589 = t603 * t469 + t604 * t470 + t602 * t502;
t535 = sin(pkin(11));
t537 = cos(pkin(11));
t471 = -t501 * t535 + t517 * t537;
t572 = t517 * t535;
t472 = t501 * t537 + t572;
t588 = -Icges(4,4) * t501 + Icges(5,5) * t472 - Icges(4,6) * t517 + Icges(5,6) * t471 + t601 * t500;
t473 = -t503 * t535 + t519 * t537;
t571 = t519 * t535;
t474 = t503 * t537 + t571;
t587 = -Icges(4,4) * t503 + Icges(5,5) * t474 - Icges(4,6) * t519 + Icges(5,6) * t473 + t601 * t502;
t516 = t539 * t573 + t540 * t558;
t569 = t536 * t542;
t493 = t516 * t534 + t554 * t569;
t494 = t516 * t554 - t534 * t569;
t515 = t536 * t539 * t540 - t573 * t575;
t586 = t600 * t493 + t603 * t494 + t598 * t515;
t585 = t598 * t493 + t602 * t494 + t599 * t515;
t584 = t603 * t493 + t604 * t494 + t602 * t515;
t498 = -t516 * t535 - t537 * t569;
t559 = t535 * t569;
t499 = t516 * t537 - t559;
t583 = -Icges(4,4) * t516 + Icges(5,5) * t499 + Icges(4,6) * t569 + Icges(5,6) * t498 + t601 * t515;
t574 = pkin(4) * t537;
t566 = rSges(7,2) * t500 + t596 * t467 + t597 * t468;
t565 = rSges(7,2) * t502 + t596 * t469 + t597 * t470;
t564 = rSges(7,2) * t515 + t596 * t493 + t597 * t494;
t491 = pkin(2) * t518 + pkin(9) * t517;
t492 = pkin(2) * t520 + pkin(9) * t519;
t561 = qJD(2) * t536;
t529 = t541 * t561;
t557 = t543 * t561;
t563 = t491 * t529 + t492 * t557;
t504 = qJD(3) * t519 + t529;
t562 = qJD(1) * (pkin(1) * t541 - pkin(8) * t568);
t530 = qJD(2) * t573 + qJD(1);
t463 = pkin(3) * t501 + qJ(4) * t500;
t553 = qJD(4) * t515 + t504 * t463 + t563;
t505 = qJD(3) * t517 - t557;
t522 = -qJD(3) * t569 + t530;
t521 = (pkin(2) * t540 - pkin(9) * t542) * t536;
t523 = qJD(1) * (pkin(1) * t543 + pkin(8) * t570);
t551 = t530 * t492 - t521 * t529 + t523;
t464 = pkin(3) * t503 + qJ(4) * t502;
t550 = qJD(4) * t500 + t522 * t464 + t551;
t406 = pkin(4) * t572 + pkin(10) * t500 + t501 * t574;
t407 = pkin(4) * t571 + pkin(10) * t502 + t503 * t574;
t549 = t504 * t406 + (-t407 - t464) * t505 + t553;
t548 = -t491 * t530 - t521 * t557 - t562;
t490 = pkin(3) * t516 + qJ(4) * t515;
t547 = qJD(4) * t502 + t505 * t490 + t548;
t446 = -pkin(4) * t559 + pkin(10) * t515 + t516 * t574;
t546 = t522 * t407 + (-t446 - t490) * t504 + t550;
t545 = t505 * t446 + (-t406 - t463) * t522 + t547;
t526 = rSges(2,1) * t543 - rSges(2,2) * t541;
t525 = rSges(2,1) * t541 + rSges(2,2) * t543;
t511 = t573 * rSges(3,3) + (rSges(3,1) * t540 + rSges(3,2) * t542) * t536;
t510 = Icges(3,5) * t573 + (Icges(3,1) * t540 + Icges(3,4) * t542) * t536;
t509 = Icges(3,6) * t573 + (Icges(3,4) * t540 + Icges(3,2) * t542) * t536;
t508 = Icges(3,3) * t573 + (Icges(3,5) * t540 + Icges(3,6) * t542) * t536;
t495 = qJD(5) * t515 + t522;
t487 = rSges(3,1) * t520 - rSges(3,2) * t519 + rSges(3,3) * t570;
t486 = rSges(3,1) * t518 - rSges(3,2) * t517 - rSges(3,3) * t568;
t484 = Icges(3,1) * t520 - Icges(3,4) * t519 + Icges(3,5) * t570;
t483 = Icges(3,1) * t518 - Icges(3,4) * t517 - Icges(3,5) * t568;
t482 = Icges(3,4) * t520 - Icges(3,2) * t519 + Icges(3,6) * t570;
t481 = Icges(3,4) * t518 - Icges(3,2) * t517 - Icges(3,6) * t568;
t478 = rSges(4,1) * t516 - rSges(4,2) * t515 - rSges(4,3) * t569;
t477 = Icges(4,1) * t516 - Icges(4,4) * t515 - Icges(4,5) * t569;
t475 = Icges(4,5) * t516 - Icges(4,6) * t515 - Icges(4,3) * t569;
t466 = qJD(5) * t500 + t505;
t465 = qJD(5) * t502 + t504;
t459 = rSges(4,1) * t503 - rSges(4,2) * t502 + rSges(4,3) * t519;
t458 = rSges(4,1) * t501 - rSges(4,2) * t500 + rSges(4,3) * t517;
t457 = Icges(4,1) * t503 - Icges(4,4) * t502 + Icges(4,5) * t519;
t456 = Icges(4,1) * t501 - Icges(4,4) * t500 + Icges(4,5) * t517;
t453 = Icges(4,5) * t503 - Icges(4,6) * t502 + Icges(4,3) * t519;
t452 = Icges(4,5) * t501 - Icges(4,6) * t500 + Icges(4,3) * t517;
t451 = rSges(5,1) * t499 + rSges(5,2) * t498 + rSges(5,3) * t515;
t450 = Icges(5,1) * t499 + Icges(5,4) * t498 + Icges(5,5) * t515;
t449 = Icges(5,4) * t499 + Icges(5,2) * t498 + Icges(5,6) * t515;
t445 = rSges(6,1) * t494 - rSges(6,2) * t493 + rSges(6,3) * t515;
t437 = t487 * t530 - t511 * t529 + t523;
t436 = -t486 * t530 - t511 * t557 - t562;
t435 = (t486 * t541 + t487 * t543) * t561;
t431 = rSges(5,1) * t474 + rSges(5,2) * t473 + rSges(5,3) * t502;
t430 = rSges(5,1) * t472 + rSges(5,2) * t471 + rSges(5,3) * t500;
t429 = Icges(5,1) * t474 + Icges(5,4) * t473 + Icges(5,5) * t502;
t428 = Icges(5,1) * t472 + Icges(5,4) * t471 + Icges(5,5) * t500;
t427 = Icges(5,4) * t474 + Icges(5,2) * t473 + Icges(5,6) * t502;
t426 = Icges(5,4) * t472 + Icges(5,2) * t471 + Icges(5,6) * t500;
t423 = rSges(6,1) * t470 - rSges(6,2) * t469 + rSges(6,3) * t502;
t421 = rSges(6,1) * t468 - rSges(6,2) * t467 + rSges(6,3) * t500;
t403 = t459 * t522 - t478 * t504 + t551;
t402 = -t458 * t522 + t478 * t505 + t548;
t401 = t458 * t504 - t459 * t505 + t563;
t400 = t431 * t522 + (-t451 - t490) * t504 + t550;
t399 = t451 * t505 + (-t430 - t463) * t522 + t547;
t398 = t430 * t504 + (-t431 - t464) * t505 + t553;
t397 = t423 * t495 - t445 * t465 + t546;
t396 = -t421 * t495 + t445 * t466 + t545;
t395 = t421 * t465 - t423 * t466 + t549;
t394 = qJD(6) * t467 - t465 * t564 + t495 * t565 + t546;
t393 = qJD(6) * t469 + t466 * t564 - t495 * t566 + t545;
t392 = qJD(6) * t493 + t465 * t566 - t466 * t565 + t549;
t1 = ((t508 * t570 - t509 * t519 + t510 * t520) * t530 + (-(-t481 * t519 + t483 * t520) * t543 + (-t519 * t482 + t520 * t484 - t595) * t541) * t561) * t529 / 0.2e1 - ((-t508 * t568 - t509 * t517 + t510 * t518) * t530 + ((-t482 * t517 + t484 * t518) * t541 + (t517 * t481 - t518 * t483 + t595) * t543) * t561) * t557 / 0.2e1 + t530 * ((t573 * t480 + (t482 * t542 + t484 * t540) * t536) * t529 - (t573 * t479 + (t481 * t542 + t483 * t540) * t536) * t557 + (t573 * t508 + (t509 * t542 + t510 * t540) * t536) * t530) / 0.2e1 + m(3) * (t435 ^ 2 + t436 ^ 2 + t437 ^ 2) / 0.2e1 + m(7) * (t392 ^ 2 + t393 ^ 2 + t394 ^ 2) / 0.2e1 + m(6) * (t395 ^ 2 + t396 ^ 2 + t397 ^ 2) / 0.2e1 + m(5) * (t398 ^ 2 + t399 ^ 2 + t400 ^ 2) / 0.2e1 + m(4) * (t401 ^ 2 + t402 ^ 2 + t403 ^ 2) / 0.2e1 + ((t469 * t586 + t470 * t584 + t502 * t585) * t495 + (t469 * t594 + t470 * t590 + t502 * t592) * t466 + (t593 * t469 + t589 * t470 + t591 * t502) * t465) * t465 / 0.2e1 + ((t467 * t586 + t468 * t584 + t500 * t585) * t495 + (t594 * t467 + t590 * t468 + t592 * t500) * t466 + (t467 * t593 + t468 * t589 + t500 * t591) * t465) * t466 / 0.2e1 + ((t586 * t493 + t584 * t494 + t585 * t515) * t495 + (t493 * t594 + t494 * t590 + t515 * t592) * t466 + (t493 * t593 + t494 * t589 + t515 * t591) * t465) * t495 / 0.2e1 + ((t449 * t473 + t450 * t474 + t475 * t519 + t477 * t503 + t502 * t583) * t522 + (t426 * t473 + t428 * t474 + t452 * t519 + t456 * t503 + t502 * t588) * t505 + (t427 * t473 + t429 * t474 + t453 * t519 + t457 * t503 + t587 * t502) * t504) * t504 / 0.2e1 + ((t449 * t471 + t450 * t472 + t475 * t517 + t477 * t501 + t500 * t583) * t522 + (t426 * t471 + t428 * t472 + t452 * t517 + t456 * t501 + t588 * t500) * t505 + (t427 * t471 + t429 * t472 + t453 * t517 + t457 * t501 + t500 * t587) * t504) * t505 / 0.2e1 + ((t449 * t498 + t450 * t499 - t475 * t569 + t477 * t516 + t583 * t515) * t522 + (t426 * t498 + t428 * t499 - t452 * t569 + t456 * t516 + t515 * t588) * t505 + (t427 * t498 + t429 * t499 - t453 * t569 + t457 * t516 + t515 * t587) * t504) * t522 / 0.2e1 + (Icges(2,3) + m(2) * (t525 ^ 2 + t526 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
