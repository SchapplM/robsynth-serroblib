% Calculate kinetic energy for
% S6PRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 01:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_energykin_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRR4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRR4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:54:49
% EndTime: 2019-03-09 00:54:51
% DurationCPUTime: 2.26s
% Computational Cost: add. (5310->385), mult. (13161->613), div. (0->0), fcn. (16854->16), ass. (0->176)
t627 = qJD(2) ^ 2;
t626 = cos(qJ(3));
t592 = cos(qJ(4));
t625 = pkin(4) * t592;
t623 = cos(pkin(7));
t622 = sin(pkin(7));
t583 = sin(pkin(13));
t585 = cos(pkin(13));
t590 = sin(qJ(2));
t586 = cos(pkin(6));
t593 = cos(qJ(2));
t614 = t586 * t593;
t572 = -t583 * t590 + t585 * t614;
t584 = sin(pkin(6));
t607 = t584 * t623;
t559 = -t572 * t622 - t585 * t607;
t588 = sin(qJ(4));
t621 = t559 * t588;
t574 = -t583 * t614 - t585 * t590;
t560 = -t574 * t622 + t583 * t607;
t620 = t560 * t588;
t606 = t584 * t622;
t571 = t586 * t623 - t593 * t606;
t619 = t571 * t588;
t618 = t583 * t584;
t617 = t584 * t590;
t616 = t585 * t584;
t615 = t586 * t590;
t613 = qJ(4) + qJ(5);
t612 = qJD(2) * t584;
t579 = t583 * t612;
t549 = qJD(3) * t560 + t579;
t581 = qJD(2) * t586;
t562 = qJD(3) * t571 + t581;
t575 = -t583 * t615 + t585 * t593;
t589 = sin(qJ(3));
t604 = t626 * t622;
t603 = t584 * t604;
t605 = t623 * t626;
t534 = -t574 * t605 + t575 * t589 - t583 * t603;
t503 = qJD(4) * t534 + t549;
t557 = -t584 * t593 * t605 - t586 * t604 + t589 * t617;
t527 = qJD(4) * t557 + t562;
t610 = t585 * t612;
t573 = t583 * t593 + t585 * t615;
t539 = t573 * pkin(2) + pkin(9) * t559;
t540 = t575 * pkin(2) + pkin(9) * t560;
t609 = t539 * t579 + t540 * t610 + qJD(1);
t608 = cos(t613);
t478 = qJD(5) * t534 + t503;
t515 = qJD(5) * t557 + t527;
t550 = qJD(3) * t559 - t610;
t561 = pkin(2) * t617 + pkin(9) * t571;
t602 = t540 * t581 - t561 * t579;
t532 = -t572 * t605 + t573 * t589 + t585 * t603;
t504 = qJD(4) * t532 + t550;
t479 = qJD(5) * t532 + t504;
t533 = t573 * t626 + (t572 * t623 - t585 * t606) * t589;
t499 = t533 * pkin(3) + t532 * pkin(10);
t535 = t575 * t626 + (t574 * t623 + t583 * t606) * t589;
t500 = t535 * pkin(3) + t534 * pkin(10);
t601 = t549 * t499 - t500 * t550 + t609;
t600 = (-t539 * t586 - t561 * t616) * qJD(2);
t558 = t586 * t622 * t589 + (t589 * t593 * t623 + t590 * t626) * t584;
t522 = t558 * pkin(3) + t557 * pkin(10);
t599 = t562 * t500 - t522 * t549 + t602;
t445 = pkin(4) * t621 + pkin(11) * t532 + t533 * t625;
t446 = pkin(4) * t620 + pkin(11) * t534 + t535 * t625;
t598 = t503 * t445 - t446 * t504 + t601;
t597 = -t499 * t562 + t550 * t522 + t600;
t476 = pkin(4) * t619 + pkin(11) * t557 + t558 * t625;
t596 = t527 * t446 - t476 * t503 + t599;
t595 = -t445 * t527 + t504 * t476 + t597;
t591 = cos(qJ(6));
t587 = sin(qJ(6));
t582 = sin(t613);
t569 = rSges(3,3) * t586 + (rSges(3,1) * t590 + rSges(3,2) * t593) * t584;
t568 = Icges(3,5) * t586 + (Icges(3,1) * t590 + Icges(3,4) * t593) * t584;
t567 = Icges(3,6) * t586 + (Icges(3,4) * t590 + Icges(3,2) * t593) * t584;
t566 = Icges(3,3) * t586 + (Icges(3,5) * t590 + Icges(3,6) * t593) * t584;
t548 = rSges(3,1) * t575 + rSges(3,2) * t574 + rSges(3,3) * t618;
t547 = rSges(3,1) * t573 + rSges(3,2) * t572 - rSges(3,3) * t616;
t546 = Icges(3,1) * t575 + Icges(3,4) * t574 + Icges(3,5) * t618;
t545 = Icges(3,1) * t573 + Icges(3,4) * t572 - Icges(3,5) * t616;
t544 = Icges(3,4) * t575 + Icges(3,2) * t574 + Icges(3,6) * t618;
t543 = Icges(3,4) * t573 + Icges(3,2) * t572 - Icges(3,6) * t616;
t542 = Icges(3,5) * t575 + Icges(3,6) * t574 + Icges(3,3) * t618;
t541 = Icges(3,5) * t573 + Icges(3,6) * t572 - Icges(3,3) * t616;
t538 = t558 * t592 + t619;
t537 = -t558 * t588 + t571 * t592;
t524 = t558 * t608 + t571 * t582;
t523 = t558 * t582 - t571 * t608;
t521 = (-t547 * t586 - t569 * t616) * qJD(2);
t520 = (t548 * t586 - t569 * t618) * qJD(2);
t519 = rSges(4,1) * t558 - rSges(4,2) * t557 + rSges(4,3) * t571;
t518 = Icges(4,1) * t558 - Icges(4,4) * t557 + Icges(4,5) * t571;
t517 = Icges(4,4) * t558 - Icges(4,2) * t557 + Icges(4,6) * t571;
t516 = Icges(4,5) * t558 - Icges(4,6) * t557 + Icges(4,3) * t571;
t514 = t535 * t592 + t620;
t513 = -t535 * t588 + t560 * t592;
t512 = t533 * t592 + t621;
t511 = -t533 * t588 + t559 * t592;
t510 = t535 * t608 + t560 * t582;
t509 = t535 * t582 - t560 * t608;
t508 = t533 * t608 + t559 * t582;
t507 = t533 * t582 - t559 * t608;
t506 = t524 * t591 + t557 * t587;
t505 = -t524 * t587 + t557 * t591;
t502 = qJD(1) + (t547 * t583 + t548 * t585) * t612;
t498 = pkin(5) * t524 + pkin(12) * t523;
t496 = rSges(5,1) * t538 + rSges(5,2) * t537 + rSges(5,3) * t557;
t495 = rSges(4,1) * t535 - rSges(4,2) * t534 + rSges(4,3) * t560;
t494 = rSges(4,1) * t533 - rSges(4,2) * t532 + rSges(4,3) * t559;
t493 = Icges(5,1) * t538 + Icges(5,4) * t537 + Icges(5,5) * t557;
t492 = Icges(5,4) * t538 + Icges(5,2) * t537 + Icges(5,6) * t557;
t491 = Icges(5,5) * t538 + Icges(5,6) * t537 + Icges(5,3) * t557;
t490 = Icges(4,1) * t535 - Icges(4,4) * t534 + Icges(4,5) * t560;
t489 = Icges(4,1) * t533 - Icges(4,4) * t532 + Icges(4,5) * t559;
t488 = Icges(4,4) * t535 - Icges(4,2) * t534 + Icges(4,6) * t560;
t487 = Icges(4,4) * t533 - Icges(4,2) * t532 + Icges(4,6) * t559;
t486 = Icges(4,5) * t535 - Icges(4,6) * t534 + Icges(4,3) * t560;
t485 = Icges(4,5) * t533 - Icges(4,6) * t532 + Icges(4,3) * t559;
t484 = qJD(6) * t523 + t515;
t483 = rSges(6,1) * t524 - rSges(6,2) * t523 + rSges(6,3) * t557;
t482 = Icges(6,1) * t524 - Icges(6,4) * t523 + Icges(6,5) * t557;
t481 = Icges(6,4) * t524 - Icges(6,2) * t523 + Icges(6,6) * t557;
t480 = Icges(6,5) * t524 - Icges(6,6) * t523 + Icges(6,3) * t557;
t475 = t510 * t591 + t534 * t587;
t474 = -t510 * t587 + t534 * t591;
t473 = t508 * t591 + t532 * t587;
t472 = -t508 * t587 + t532 * t591;
t471 = pkin(5) * t510 + pkin(12) * t509;
t470 = pkin(5) * t508 + pkin(12) * t507;
t469 = rSges(5,1) * t514 + rSges(5,2) * t513 + rSges(5,3) * t534;
t468 = rSges(5,1) * t512 + rSges(5,2) * t511 + rSges(5,3) * t532;
t467 = Icges(5,1) * t514 + Icges(5,4) * t513 + Icges(5,5) * t534;
t466 = Icges(5,1) * t512 + Icges(5,4) * t511 + Icges(5,5) * t532;
t465 = Icges(5,4) * t514 + Icges(5,2) * t513 + Icges(5,6) * t534;
t464 = Icges(5,4) * t512 + Icges(5,2) * t511 + Icges(5,6) * t532;
t463 = Icges(5,5) * t514 + Icges(5,6) * t513 + Icges(5,3) * t534;
t462 = Icges(5,5) * t512 + Icges(5,6) * t511 + Icges(5,3) * t532;
t461 = rSges(6,1) * t510 - rSges(6,2) * t509 + rSges(6,3) * t534;
t460 = rSges(6,1) * t508 - rSges(6,2) * t507 + rSges(6,3) * t532;
t459 = qJD(6) * t507 + t479;
t458 = qJD(6) * t509 + t478;
t457 = Icges(6,1) * t510 - Icges(6,4) * t509 + Icges(6,5) * t534;
t456 = Icges(6,1) * t508 - Icges(6,4) * t507 + Icges(6,5) * t532;
t455 = Icges(6,4) * t510 - Icges(6,2) * t509 + Icges(6,6) * t534;
t454 = Icges(6,4) * t508 - Icges(6,2) * t507 + Icges(6,6) * t532;
t453 = Icges(6,5) * t510 - Icges(6,6) * t509 + Icges(6,3) * t534;
t452 = Icges(6,5) * t508 - Icges(6,6) * t507 + Icges(6,3) * t532;
t450 = rSges(7,1) * t506 + rSges(7,2) * t505 + rSges(7,3) * t523;
t449 = Icges(7,1) * t506 + Icges(7,4) * t505 + Icges(7,5) * t523;
t448 = Icges(7,4) * t506 + Icges(7,2) * t505 + Icges(7,6) * t523;
t447 = Icges(7,5) * t506 + Icges(7,6) * t505 + Icges(7,3) * t523;
t442 = -t494 * t562 + t519 * t550 + t600;
t441 = t495 * t562 - t519 * t549 + t602;
t440 = rSges(7,1) * t475 + rSges(7,2) * t474 + rSges(7,3) * t509;
t439 = rSges(7,1) * t473 + rSges(7,2) * t472 + rSges(7,3) * t507;
t438 = Icges(7,1) * t475 + Icges(7,4) * t474 + Icges(7,5) * t509;
t437 = Icges(7,1) * t473 + Icges(7,4) * t472 + Icges(7,5) * t507;
t436 = Icges(7,4) * t475 + Icges(7,2) * t474 + Icges(7,6) * t509;
t435 = Icges(7,4) * t473 + Icges(7,2) * t472 + Icges(7,6) * t507;
t434 = Icges(7,5) * t475 + Icges(7,6) * t474 + Icges(7,3) * t509;
t433 = Icges(7,5) * t473 + Icges(7,6) * t472 + Icges(7,3) * t507;
t432 = t494 * t549 - t495 * t550 + t609;
t431 = -t468 * t527 + t496 * t504 + t597;
t430 = t469 * t527 - t496 * t503 + t599;
t429 = t468 * t503 - t469 * t504 + t601;
t428 = -t460 * t515 + t479 * t483 + t595;
t427 = t461 * t515 - t478 * t483 + t596;
t426 = t460 * t478 - t461 * t479 + t598;
t425 = -t439 * t484 + t450 * t459 - t470 * t515 + t479 * t498 + t595;
t424 = t440 * t484 - t450 * t458 + t471 * t515 - t478 * t498 + t596;
t423 = t439 * t458 - t440 * t459 + t470 * t478 - t471 * t479 + t598;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + t459 * ((t434 * t507 + t436 * t472 + t438 * t473) * t458 + (t507 * t433 + t472 * t435 + t473 * t437) * t459 + (t447 * t507 + t448 * t472 + t449 * t473) * t484) / 0.2e1 + t458 * ((t509 * t434 + t474 * t436 + t475 * t438) * t458 + (t433 * t509 + t435 * t474 + t437 * t475) * t459 + (t447 * t509 + t448 * t474 + t449 * t475) * t484) / 0.2e1 + t484 * ((t434 * t523 + t436 * t505 + t438 * t506) * t458 + (t433 * t523 + t435 * t505 + t437 * t506) * t459 + (t523 * t447 + t505 * t448 + t506 * t449) * t484) / 0.2e1 + t478 * ((t534 * t453 - t509 * t455 + t510 * t457) * t478 + (t452 * t534 - t454 * t509 + t456 * t510) * t479 + (t480 * t534 - t481 * t509 + t482 * t510) * t515) / 0.2e1 + t479 * ((t453 * t532 - t455 * t507 + t457 * t508) * t478 + (t532 * t452 - t507 * t454 + t508 * t456) * t479 + (t480 * t532 - t481 * t507 + t482 * t508) * t515) / 0.2e1 + t515 * ((t453 * t557 - t455 * t523 + t457 * t524) * t478 + (t452 * t557 - t454 * t523 + t456 * t524) * t479 + (t557 * t480 - t523 * t481 + t524 * t482) * t515) / 0.2e1 + t504 * ((t463 * t532 + t465 * t511 + t467 * t512) * t503 + (t532 * t462 + t511 * t464 + t512 * t466) * t504 + (t491 * t532 + t492 * t511 + t493 * t512) * t527) / 0.2e1 + t503 * ((t534 * t463 + t513 * t465 + t514 * t467) * t503 + (t462 * t534 + t464 * t513 + t466 * t514) * t504 + (t491 * t534 + t492 * t513 + t493 * t514) * t527) / 0.2e1 + t527 * ((t463 * t557 + t465 * t537 + t467 * t538) * t503 + (t462 * t557 + t464 * t537 + t466 * t538) * t504 + (t557 * t491 + t537 * t492 + t538 * t493) * t527) / 0.2e1 + t562 * ((t486 * t571 - t488 * t557 + t490 * t558) * t549 + (t485 * t571 - t487 * t557 + t489 * t558) * t550 + (t571 * t516 - t557 * t517 + t558 * t518) * t562) / 0.2e1 + t549 * ((t560 * t486 - t534 * t488 + t535 * t490) * t549 + (t485 * t560 - t487 * t534 + t489 * t535) * t550 + (t516 * t560 - t517 * t534 + t518 * t535) * t562) / 0.2e1 + t550 * ((t486 * t559 - t488 * t532 + t490 * t533) * t549 + (t559 * t485 - t532 * t487 + t533 * t489) * t550 + (t516 * t559 - t517 * t532 + t518 * t533) * t562) / 0.2e1 + m(7) * (t423 ^ 2 + t424 ^ 2 + t425 ^ 2) / 0.2e1 + m(6) * (t426 ^ 2 + t427 ^ 2 + t428 ^ 2) / 0.2e1 + m(5) * (t429 ^ 2 + t430 ^ 2 + t431 ^ 2) / 0.2e1 + m(4) * (t432 ^ 2 + t441 ^ 2 + t442 ^ 2) / 0.2e1 + m(3) * (t502 ^ 2 + t520 ^ 2 + t521 ^ 2) / 0.2e1 - t627 * ((-t542 * t616 + t544 * t572 + t546 * t573) * t618 - (-t541 * t616 + t543 * t572 + t545 * t573) * t616 + (-t566 * t616 + t567 * t572 + t568 * t573) * t586) * t616 / 0.2e1 + (t586 * (t586 ^ 2 * t566 + (((t544 * t593 + t546 * t590) * t583 - (t543 * t593 + t545 * t590) * t585) * t584 + (-t541 * t585 + t542 * t583 + t567 * t593 + t568 * t590) * t586) * t584) + ((t542 * t618 + t544 * t574 + t546 * t575) * t618 - (t541 * t618 + t543 * t574 + t545 * t575) * t616 + (t566 * t618 + t567 * t574 + t568 * t575) * t586) * t618) * t627 / 0.2e1;
T  = t1;
