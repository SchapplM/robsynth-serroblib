% Calculate kinetic energy for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
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
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPPRRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_energykin_fixb_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPPRRR1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPPRRR1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPPRRR1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:39:08
% EndTime: 2019-03-08 18:39:10
% DurationCPUTime: 1.77s
% Computational Cost: add. (8225->257), mult. (23364->420), div. (0->0), fcn. (31051->18), ass. (0->134)
t565 = sin(pkin(12));
t567 = cos(pkin(12));
t609 = cos(pkin(13));
t612 = cos(pkin(6));
t600 = t612 * t609;
t605 = sin(pkin(13));
t587 = -t565 * t600 - t567 * t605;
t611 = cos(pkin(7));
t584 = t587 * t611;
t599 = t612 * t605;
t588 = -t565 * t599 + t567 * t609;
t566 = sin(pkin(6));
t607 = sin(pkin(7));
t608 = cos(pkin(14));
t597 = t607 * t608;
t595 = t566 * t597;
t604 = sin(pkin(14));
t575 = -t565 * t595 - t584 * t608 + t588 * t604;
t601 = t566 * t611;
t581 = t565 * t601 - t587 * t607;
t606 = sin(pkin(8));
t610 = cos(pkin(8));
t617 = t575 * t610 - t581 * t606;
t589 = -t565 * t605 + t567 * t600;
t585 = t589 * t611;
t590 = t565 * t609 + t567 * t599;
t576 = t567 * t595 - t585 * t608 + t590 * t604;
t582 = -t567 * t601 - t589 * t607;
t616 = t576 * t610 - t582 * t606;
t598 = t611 * t609;
t578 = t612 * t597 + (t598 * t608 - t604 * t605) * t566;
t586 = -t566 * t607 * t609 + t611 * t612;
t615 = t578 * t610 + t586 * t606;
t614 = cos(qJ(4));
t613 = cos(qJ(5));
t596 = t607 * t604;
t594 = t566 * t596;
t547 = -t567 * t594 + t585 * t604 + t590 * t608;
t570 = sin(qJ(4));
t529 = t547 * t570 + t614 * t616;
t541 = t576 * t606 + t582 * t610;
t539 = qJD(4) * t541;
t518 = qJD(5) * t529 + t539;
t548 = t565 * t594 + t584 * t604 + t588 * t608;
t531 = t548 * t570 + t614 * t617;
t542 = t575 * t606 + t581 * t610;
t540 = qJD(4) * t542;
t519 = qJD(5) * t531 + t540;
t557 = t612 * t596 + (t598 * t604 + t605 * t608) * t566;
t537 = t557 * t570 - t614 * t615;
t549 = -t578 * t606 + t586 * t610;
t546 = qJD(4) * t549;
t533 = qJD(5) * t537 + t546;
t603 = qJD(2) * t566;
t550 = qJD(3) * t581 + t565 * t603;
t562 = qJD(2) * t612 + qJD(1);
t558 = qJD(3) * t586 + t562;
t551 = qJD(3) * t582 - t567 * t603;
t530 = t547 * t614 - t570 * t616;
t508 = pkin(4) * t530 + pkin(10) * t529;
t538 = t557 * t614 + t570 * t615;
t524 = pkin(4) * t538 + pkin(10) * t537;
t593 = -t508 * t546 + t524 * t539 + t550;
t532 = t548 * t614 - t570 * t617;
t509 = pkin(4) * t532 + pkin(10) * t531;
t592 = t508 * t540 - t509 * t539 + t558;
t591 = t509 * t546 - t524 * t540 + t551;
t571 = cos(qJ(6));
t569 = sin(qJ(5));
t568 = sin(qJ(6));
t528 = t538 * t613 + t549 * t569;
t527 = t538 * t569 - t549 * t613;
t523 = rSges(5,1) * t538 - rSges(5,2) * t537 + rSges(5,3) * t549;
t522 = Icges(5,1) * t538 - Icges(5,4) * t537 + Icges(5,5) * t549;
t521 = Icges(5,4) * t538 - Icges(5,2) * t537 + Icges(5,6) * t549;
t520 = Icges(5,5) * t538 - Icges(5,6) * t537 + Icges(5,3) * t549;
t517 = t532 * t613 + t542 * t569;
t516 = t532 * t569 - t542 * t613;
t515 = t530 * t613 + t541 * t569;
t514 = t530 * t569 - t541 * t613;
t513 = t528 * t571 + t537 * t568;
t512 = -t528 * t568 + t537 * t571;
t510 = qJD(6) * t527 + t533;
t507 = pkin(5) * t528 + pkin(11) * t527;
t504 = rSges(5,1) * t532 - rSges(5,2) * t531 + rSges(5,3) * t542;
t503 = rSges(5,1) * t530 - rSges(5,2) * t529 + rSges(5,3) * t541;
t502 = Icges(5,1) * t532 - Icges(5,4) * t531 + Icges(5,5) * t542;
t501 = Icges(5,1) * t530 - Icges(5,4) * t529 + Icges(5,5) * t541;
t500 = Icges(5,4) * t532 - Icges(5,2) * t531 + Icges(5,6) * t542;
t499 = Icges(5,4) * t530 - Icges(5,2) * t529 + Icges(5,6) * t541;
t498 = Icges(5,5) * t532 - Icges(5,6) * t531 + Icges(5,3) * t542;
t497 = Icges(5,5) * t530 - Icges(5,6) * t529 + Icges(5,3) * t541;
t496 = rSges(6,1) * t528 - rSges(6,2) * t527 + rSges(6,3) * t537;
t495 = Icges(6,1) * t528 - Icges(6,4) * t527 + Icges(6,5) * t537;
t494 = Icges(6,4) * t528 - Icges(6,2) * t527 + Icges(6,6) * t537;
t493 = Icges(6,5) * t528 - Icges(6,6) * t527 + Icges(6,3) * t537;
t492 = t517 * t571 + t531 * t568;
t491 = -t517 * t568 + t531 * t571;
t490 = t515 * t571 + t529 * t568;
t489 = -t515 * t568 + t529 * t571;
t488 = qJD(6) * t516 + t519;
t487 = qJD(6) * t514 + t518;
t486 = pkin(5) * t517 + pkin(11) * t516;
t485 = pkin(5) * t515 + pkin(11) * t514;
t484 = rSges(6,1) * t517 - rSges(6,2) * t516 + rSges(6,3) * t531;
t483 = rSges(6,1) * t515 - rSges(6,2) * t514 + rSges(6,3) * t529;
t482 = Icges(6,1) * t517 - Icges(6,4) * t516 + Icges(6,5) * t531;
t481 = Icges(6,1) * t515 - Icges(6,4) * t514 + Icges(6,5) * t529;
t480 = Icges(6,4) * t517 - Icges(6,2) * t516 + Icges(6,6) * t531;
t479 = Icges(6,4) * t515 - Icges(6,2) * t514 + Icges(6,6) * t529;
t478 = Icges(6,5) * t517 - Icges(6,6) * t516 + Icges(6,3) * t531;
t477 = Icges(6,5) * t515 - Icges(6,6) * t514 + Icges(6,3) * t529;
t476 = rSges(7,1) * t513 + rSges(7,2) * t512 + rSges(7,3) * t527;
t475 = Icges(7,1) * t513 + Icges(7,4) * t512 + Icges(7,5) * t527;
t474 = Icges(7,4) * t513 + Icges(7,2) * t512 + Icges(7,6) * t527;
t473 = Icges(7,5) * t513 + Icges(7,6) * t512 + Icges(7,3) * t527;
t472 = (t504 * t549 - t523 * t542) * qJD(4) + t551;
t471 = (-t503 * t549 + t523 * t541) * qJD(4) + t550;
t470 = (t503 * t542 - t504 * t541) * qJD(4) + t558;
t469 = t492 * rSges(7,1) + t491 * rSges(7,2) + rSges(7,3) * t516;
t468 = t490 * rSges(7,1) + t489 * rSges(7,2) + rSges(7,3) * t514;
t467 = Icges(7,1) * t492 + Icges(7,4) * t491 + Icges(7,5) * t516;
t466 = Icges(7,1) * t490 + Icges(7,4) * t489 + Icges(7,5) * t514;
t465 = Icges(7,4) * t492 + Icges(7,2) * t491 + Icges(7,6) * t516;
t464 = Icges(7,4) * t490 + Icges(7,2) * t489 + Icges(7,6) * t514;
t463 = Icges(7,5) * t492 + Icges(7,6) * t491 + Icges(7,3) * t516;
t462 = Icges(7,5) * t490 + Icges(7,6) * t489 + Icges(7,3) * t514;
t461 = t484 * t533 - t496 * t519 + t591;
t460 = -t483 * t533 + t496 * t518 + t593;
t459 = t483 * t519 - t484 * t518 + t592;
t458 = t469 * t510 - t488 * t476 + t486 * t533 - t507 * t519 + t591;
t457 = -t468 * t510 + t487 * t476 - t485 * t533 + t507 * t518 + t593;
t456 = t488 * t468 - t487 * t469 + t485 * t519 - t486 * t518 + t592;
t1 = m(5) * (t470 ^ 2 + t471 ^ 2 + t472 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t562 ^ 2 + (t565 ^ 2 + t567 ^ 2) * qJD(2) ^ 2 * t566 ^ 2) / 0.2e1 + m(4) * (t550 ^ 2 + t551 ^ 2 + t558 ^ 2) / 0.2e1 + t488 * ((t516 * t463 + t491 * t465 + t492 * t467) * t488 + (t462 * t516 + t491 * t464 + t492 * t466) * t487 + (t473 * t516 + t491 * t474 + t492 * t475) * t510) / 0.2e1 + t510 * ((t463 * t527 + t465 * t512 + t467 * t513) * t488 + (t462 * t527 + t464 * t512 + t466 * t513) * t487 + (t527 * t473 + t512 * t474 + t513 * t475) * t510) / 0.2e1 + t487 * ((t463 * t514 + t489 * t465 + t490 * t467) * t488 + (t514 * t462 + t489 * t464 + t490 * t466) * t487 + (t473 * t514 + t489 * t474 + t490 * t475) * t510) / 0.2e1 + t518 * ((t478 * t529 - t480 * t514 + t482 * t515) * t519 + (t529 * t477 - t514 * t479 + t515 * t481) * t518 + (t493 * t529 - t494 * t514 + t495 * t515) * t533) / 0.2e1 + t533 * ((t478 * t537 - t480 * t527 + t482 * t528) * t519 + (t477 * t537 - t479 * t527 + t481 * t528) * t518 + (t493 * t537 - t494 * t527 + t495 * t528) * t533) / 0.2e1 + t519 * ((t478 * t531 - t480 * t516 + t482 * t517) * t519 + (t477 * t531 - t479 * t516 + t481 * t517) * t518 + (t493 * t531 - t494 * t516 + t495 * t517) * t533) / 0.2e1 + m(7) * (t456 ^ 2 + t457 ^ 2 + t458 ^ 2) / 0.2e1 + m(6) * (t459 ^ 2 + t460 ^ 2 + t461 ^ 2) / 0.2e1 + (t542 * ((t498 * t542 - t531 * t500 + t532 * t502) * t542 + (t497 * t542 - t499 * t531 + t501 * t532) * t541 + (t520 * t542 - t521 * t531 + t522 * t532) * t549) + t541 * ((t498 * t541 - t500 * t529 + t502 * t530) * t542 + (t497 * t541 - t499 * t529 + t501 * t530) * t541 + (t520 * t541 - t521 * t529 + t522 * t530) * t549) + t549 * ((t498 * t549 - t500 * t537 + t502 * t538) * t542 + (t497 * t549 - t499 * t537 + t501 * t538) * t541 + (t520 * t549 - t521 * t537 + t522 * t538) * t549)) * qJD(4) ^ 2 / 0.2e1;
T  = t1;
