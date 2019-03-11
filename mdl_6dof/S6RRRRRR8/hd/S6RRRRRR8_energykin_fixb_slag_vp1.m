% Calculate kinetic energy for
% S6RRRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 05:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR8_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_energykin_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR8_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR8_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR8_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:48:16
% EndTime: 2019-03-10 04:48:19
% DurationCPUTime: 2.68s
% Computational Cost: add. (5379->392), mult. (13213->619), div. (0->0), fcn. (16888->16), ass. (0->181)
t589 = sin(pkin(6));
t590 = cos(pkin(6));
t598 = cos(qJ(2));
t599 = cos(qJ(1));
t622 = t598 * t599;
t594 = sin(qJ(2));
t595 = sin(qJ(1));
t625 = t594 * t595;
t573 = t590 * t622 - t625;
t623 = t595 * t598;
t624 = t594 * t599;
t574 = t590 * t624 + t623;
t575 = -t590 * t623 - t624;
t576 = -t590 * t625 + t622;
t626 = t589 * t599;
t627 = t589 * t595;
t611 = (Icges(3,5) * t574 + Icges(3,6) * t573 - Icges(3,3) * t626) * t599 - (Icges(3,5) * t576 + Icges(3,6) * t575 + Icges(3,3) * t627) * t595;
t638 = t589 * t611;
t636 = cos(qJ(3));
t597 = cos(qJ(4));
t635 = pkin(4) * t597;
t633 = cos(pkin(7));
t632 = sin(pkin(7));
t615 = t589 * t633;
t560 = -t573 * t632 - t599 * t615;
t592 = sin(qJ(4));
t631 = t560 * t592;
t561 = -t575 * t632 + t595 * t615;
t630 = t561 * t592;
t614 = t589 * t632;
t572 = t590 * t633 - t598 * t614;
t629 = t572 * t592;
t628 = t589 * t594;
t621 = qJ(4) + qJ(5);
t540 = t574 * pkin(2) + pkin(10) * t560;
t541 = t576 * pkin(2) + pkin(10) * t561;
t618 = qJD(2) * t589;
t584 = t595 * t618;
t617 = t599 * t618;
t620 = t540 * t584 + t541 * t617;
t550 = qJD(3) * t561 + t584;
t619 = qJD(1) * (pkin(1) * t595 - pkin(9) * t626);
t585 = qJD(2) * t590 + qJD(1);
t593 = sin(qJ(3));
t612 = t636 * t632;
t610 = t589 * t612;
t613 = t633 * t636;
t538 = -t575 * t613 + t576 * t593 - t595 * t610;
t506 = qJD(4) * t538 + t550;
t562 = qJD(3) * t572 + t585;
t616 = cos(t621);
t484 = qJD(5) * t538 + t506;
t558 = -t589 * t598 * t613 - t590 * t612 + t593 * t628;
t524 = qJD(4) * t558 + t562;
t551 = qJD(3) * t560 - t617;
t514 = qJD(5) * t558 + t524;
t536 = -t573 * t613 + t574 * t593 + t599 * t610;
t537 = t574 * t636 + (t573 * t633 - t599 * t614) * t593;
t500 = pkin(3) * t537 + pkin(11) * t536;
t539 = t576 * t636 + (t575 * t633 + t595 * t614) * t593;
t501 = pkin(3) * t539 + pkin(11) * t538;
t609 = t550 * t500 - t501 * t551 + t620;
t507 = qJD(4) * t536 + t551;
t563 = pkin(2) * t628 + pkin(10) * t572;
t577 = qJD(1) * (pkin(1) * t599 + pkin(9) * t627);
t608 = t585 * t541 - t563 * t584 + t577;
t485 = qJD(5) * t536 + t507;
t446 = pkin(4) * t631 + pkin(12) * t536 + t537 * t635;
t447 = pkin(4) * t630 + pkin(12) * t538 + t539 * t635;
t607 = t506 * t446 - t447 * t507 + t609;
t606 = -t540 * t585 - t563 * t617 - t619;
t559 = t590 * t632 * t593 + (t593 * t598 * t633 + t594 * t636) * t589;
t523 = pkin(3) * t559 + pkin(11) * t558;
t605 = t562 * t501 - t523 * t550 + t608;
t604 = -t500 * t562 + t551 * t523 + t606;
t473 = pkin(4) * t629 + pkin(12) * t558 + t559 * t635;
t603 = t524 * t447 - t473 * t506 + t605;
t602 = -t446 * t524 + t507 * t473 + t604;
t596 = cos(qJ(6));
t591 = sin(qJ(6));
t588 = sin(t621);
t582 = rSges(2,1) * t599 - rSges(2,2) * t595;
t581 = rSges(2,1) * t595 + rSges(2,2) * t599;
t570 = rSges(3,3) * t590 + (rSges(3,1) * t594 + rSges(3,2) * t598) * t589;
t569 = Icges(3,5) * t590 + (Icges(3,1) * t594 + Icges(3,4) * t598) * t589;
t568 = Icges(3,6) * t590 + (Icges(3,4) * t594 + Icges(3,2) * t598) * t589;
t567 = Icges(3,3) * t590 + (Icges(3,5) * t594 + Icges(3,6) * t598) * t589;
t549 = rSges(3,1) * t576 + rSges(3,2) * t575 + rSges(3,3) * t627;
t548 = rSges(3,1) * t574 + rSges(3,2) * t573 - rSges(3,3) * t626;
t547 = Icges(3,1) * t576 + Icges(3,4) * t575 + Icges(3,5) * t627;
t546 = Icges(3,1) * t574 + Icges(3,4) * t573 - Icges(3,5) * t626;
t545 = Icges(3,4) * t576 + Icges(3,2) * t575 + Icges(3,6) * t627;
t544 = Icges(3,4) * t574 + Icges(3,2) * t573 - Icges(3,6) * t626;
t535 = t559 * t597 + t629;
t534 = -t559 * t592 + t572 * t597;
t526 = t559 * t616 + t572 * t588;
t525 = t559 * t588 - t572 * t616;
t522 = rSges(4,1) * t559 - rSges(4,2) * t558 + rSges(4,3) * t572;
t521 = Icges(4,1) * t559 - Icges(4,4) * t558 + Icges(4,5) * t572;
t520 = Icges(4,4) * t559 - Icges(4,2) * t558 + Icges(4,6) * t572;
t519 = Icges(4,5) * t559 - Icges(4,6) * t558 + Icges(4,3) * t572;
t518 = t539 * t597 + t630;
t517 = -t539 * t592 + t561 * t597;
t516 = t537 * t597 + t631;
t515 = -t537 * t592 + t560 * t597;
t513 = t549 * t585 - t570 * t584 + t577;
t512 = -t548 * t585 - t570 * t617 - t619;
t511 = t539 * t616 + t561 * t588;
t510 = t539 * t588 - t561 * t616;
t509 = t537 * t616 + t560 * t588;
t508 = t537 * t588 - t560 * t616;
t505 = t526 * t596 + t558 * t591;
t504 = -t526 * t591 + t558 * t596;
t503 = (t548 * t595 + t549 * t599) * t618;
t499 = pkin(5) * t526 + pkin(13) * t525;
t497 = rSges(4,1) * t539 - rSges(4,2) * t538 + rSges(4,3) * t561;
t496 = rSges(4,1) * t537 - rSges(4,2) * t536 + rSges(4,3) * t560;
t495 = Icges(4,1) * t539 - Icges(4,4) * t538 + Icges(4,5) * t561;
t494 = Icges(4,1) * t537 - Icges(4,4) * t536 + Icges(4,5) * t560;
t493 = Icges(4,4) * t539 - Icges(4,2) * t538 + Icges(4,6) * t561;
t492 = Icges(4,4) * t537 - Icges(4,2) * t536 + Icges(4,6) * t560;
t491 = Icges(4,5) * t539 - Icges(4,6) * t538 + Icges(4,3) * t561;
t490 = Icges(4,5) * t537 - Icges(4,6) * t536 + Icges(4,3) * t560;
t489 = rSges(5,1) * t535 + rSges(5,2) * t534 + rSges(5,3) * t558;
t488 = Icges(5,1) * t535 + Icges(5,4) * t534 + Icges(5,5) * t558;
t487 = Icges(5,4) * t535 + Icges(5,2) * t534 + Icges(5,6) * t558;
t486 = Icges(5,5) * t535 + Icges(5,6) * t534 + Icges(5,3) * t558;
t482 = qJD(6) * t525 + t514;
t481 = rSges(6,1) * t526 - rSges(6,2) * t525 + rSges(6,3) * t558;
t480 = t511 * t596 + t538 * t591;
t479 = -t511 * t591 + t538 * t596;
t478 = t509 * t596 + t536 * t591;
t477 = -t509 * t591 + t536 * t596;
t476 = Icges(6,1) * t526 - Icges(6,4) * t525 + Icges(6,5) * t558;
t475 = Icges(6,4) * t526 - Icges(6,2) * t525 + Icges(6,6) * t558;
t474 = Icges(6,5) * t526 - Icges(6,6) * t525 + Icges(6,3) * t558;
t472 = pkin(5) * t511 + pkin(13) * t510;
t471 = pkin(5) * t509 + pkin(13) * t508;
t470 = rSges(5,1) * t518 + rSges(5,2) * t517 + rSges(5,3) * t538;
t469 = rSges(5,1) * t516 + rSges(5,2) * t515 + rSges(5,3) * t536;
t468 = Icges(5,1) * t518 + Icges(5,4) * t517 + Icges(5,5) * t538;
t467 = Icges(5,1) * t516 + Icges(5,4) * t515 + Icges(5,5) * t536;
t466 = Icges(5,4) * t518 + Icges(5,2) * t517 + Icges(5,6) * t538;
t465 = Icges(5,4) * t516 + Icges(5,2) * t515 + Icges(5,6) * t536;
t464 = Icges(5,5) * t518 + Icges(5,6) * t517 + Icges(5,3) * t538;
t463 = Icges(5,5) * t516 + Icges(5,6) * t515 + Icges(5,3) * t536;
t462 = qJD(6) * t508 + t485;
t461 = qJD(6) * t510 + t484;
t460 = rSges(6,1) * t511 - rSges(6,2) * t510 + rSges(6,3) * t538;
t459 = rSges(6,1) * t509 - rSges(6,2) * t508 + rSges(6,3) * t536;
t458 = Icges(6,1) * t511 - Icges(6,4) * t510 + Icges(6,5) * t538;
t457 = Icges(6,1) * t509 - Icges(6,4) * t508 + Icges(6,5) * t536;
t456 = Icges(6,4) * t511 - Icges(6,2) * t510 + Icges(6,6) * t538;
t455 = Icges(6,4) * t509 - Icges(6,2) * t508 + Icges(6,6) * t536;
t454 = Icges(6,5) * t511 - Icges(6,6) * t510 + Icges(6,3) * t538;
t453 = Icges(6,5) * t509 - Icges(6,6) * t508 + Icges(6,3) * t536;
t451 = rSges(7,1) * t505 + rSges(7,2) * t504 + rSges(7,3) * t525;
t450 = Icges(7,1) * t505 + Icges(7,4) * t504 + Icges(7,5) * t525;
t449 = Icges(7,4) * t505 + Icges(7,2) * t504 + Icges(7,6) * t525;
t448 = Icges(7,5) * t505 + Icges(7,6) * t504 + Icges(7,3) * t525;
t443 = t497 * t562 - t522 * t550 + t608;
t442 = -t496 * t562 + t522 * t551 + t606;
t441 = rSges(7,1) * t480 + rSges(7,2) * t479 + rSges(7,3) * t510;
t440 = rSges(7,1) * t478 + rSges(7,2) * t477 + rSges(7,3) * t508;
t439 = Icges(7,1) * t480 + Icges(7,4) * t479 + Icges(7,5) * t510;
t438 = Icges(7,1) * t478 + Icges(7,4) * t477 + Icges(7,5) * t508;
t437 = Icges(7,4) * t480 + Icges(7,2) * t479 + Icges(7,6) * t510;
t436 = Icges(7,4) * t478 + Icges(7,2) * t477 + Icges(7,6) * t508;
t435 = Icges(7,5) * t480 + Icges(7,6) * t479 + Icges(7,3) * t510;
t434 = Icges(7,5) * t478 + Icges(7,6) * t477 + Icges(7,3) * t508;
t433 = t496 * t550 - t497 * t551 + t620;
t432 = t470 * t524 - t489 * t506 + t605;
t431 = -t469 * t524 + t489 * t507 + t604;
t430 = t469 * t506 - t470 * t507 + t609;
t429 = t460 * t514 - t481 * t484 + t603;
t428 = -t459 * t514 + t481 * t485 + t602;
t427 = t459 * t484 - t460 * t485 + t607;
t426 = t441 * t482 - t451 * t461 + t472 * t514 - t484 * t499 + t603;
t425 = -t440 * t482 + t451 * t462 - t471 * t514 + t485 * t499 + t602;
t424 = t440 * t461 - t441 * t462 + t471 * t484 - t472 * t485 + t607;
t1 = t524 * ((t464 * t558 + t466 * t534 + t468 * t535) * t506 + (t463 * t558 + t465 * t534 + t467 * t535) * t507 + (t558 * t486 + t534 * t487 + t535 * t488) * t524) / 0.2e1 + t506 * ((t538 * t464 + t517 * t466 + t518 * t468) * t506 + (t463 * t538 + t465 * t517 + t467 * t518) * t507 + (t486 * t538 + t487 * t517 + t488 * t518) * t524) / 0.2e1 + t551 * ((t491 * t560 - t493 * t536 + t495 * t537) * t550 + (t560 * t490 - t536 * t492 + t537 * t494) * t551 + (t519 * t560 - t520 * t536 + t521 * t537) * t562) / 0.2e1 + t562 * ((t491 * t572 - t493 * t558 + t495 * t559) * t550 + (t490 * t572 - t492 * t558 + t494 * t559) * t551 + (t572 * t519 - t558 * t520 + t559 * t521) * t562) / 0.2e1 + t550 * ((t561 * t491 - t538 * t493 + t539 * t495) * t550 + (t490 * t561 - t492 * t538 + t494 * t539) * t551 + (t519 * t561 - t520 * t538 + t521 * t539) * t562) / 0.2e1 + t585 * ((t590 * t567 + (t568 * t598 + t569 * t594) * t589) * t585 + (((t545 * t598 + t547 * t594) * t595 - (t544 * t598 + t546 * t594) * t599) * t589 - t611 * t590) * t618) / 0.2e1 + m(7) * (t424 ^ 2 + t425 ^ 2 + t426 ^ 2) / 0.2e1 + m(6) * (t427 ^ 2 + t428 ^ 2 + t429 ^ 2) / 0.2e1 + m(5) * (t430 ^ 2 + t431 ^ 2 + t432 ^ 2) / 0.2e1 + m(4) * (t433 ^ 2 + t442 ^ 2 + t443 ^ 2) / 0.2e1 + m(3) * (t503 ^ 2 + t512 ^ 2 + t513 ^ 2) / 0.2e1 + t461 * ((t510 * t435 + t479 * t437 + t480 * t439) * t461 + (t434 * t510 + t436 * t479 + t438 * t480) * t462 + (t448 * t510 + t449 * t479 + t450 * t480) * t482) / 0.2e1 + t482 * ((t435 * t525 + t437 * t504 + t439 * t505) * t461 + (t434 * t525 + t436 * t504 + t438 * t505) * t462 + (t525 * t448 + t504 * t449 + t505 * t450) * t482) / 0.2e1 + t462 * ((t435 * t508 + t437 * t477 + t439 * t478) * t461 + (t508 * t434 + t477 * t436 + t478 * t438) * t462 + (t448 * t508 + t449 * t477 + t450 * t478) * t482) / 0.2e1 + t514 * ((t454 * t558 - t525 * t456 + t526 * t458) * t484 + (t453 * t558 - t455 * t525 + t457 * t526) * t485 + (t558 * t474 - t525 * t475 + t526 * t476) * t514) / 0.2e1 + t485 * ((t454 * t536 - t456 * t508 + t458 * t509) * t484 + (t536 * t453 - t508 * t455 + t509 * t457) * t485 + (t474 * t536 - t475 * t508 + t476 * t509) * t514) / 0.2e1 + t484 * ((t538 * t454 - t510 * t456 + t511 * t458) * t484 + (t453 * t538 - t455 * t510 + t457 * t511) * t485 + (t474 * t538 - t475 * t510 + t476 * t511) * t514) / 0.2e1 + t507 * ((t464 * t536 + t466 * t515 + t468 * t516) * t506 + (t536 * t463 + t515 * t465 + t516 * t467) * t507 + (t486 * t536 + t487 * t515 + t488 * t516) * t524) / 0.2e1 - ((-t567 * t626 + t568 * t573 + t569 * t574) * t585 + ((t545 * t573 + t547 * t574) * t595 + (-t573 * t544 - t574 * t546 + t638) * t599) * t618) * t617 / 0.2e1 + ((t567 * t627 + t568 * t575 + t569 * t576) * t585 + (-(t544 * t575 + t546 * t576) * t599 + (t575 * t545 + t576 * t547 - t638) * t595) * t618) * t584 / 0.2e1 + (m(2) * (t581 ^ 2 + t582 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
