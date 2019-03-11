% Calculate kinetic energy for
% S6RRRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 02:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP11_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP11_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP11_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP11_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP11_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:32:34
% EndTime: 2019-03-10 02:32:38
% DurationCPUTime: 3.41s
% Computational Cost: add. (5408->341), mult. (14662->529), div. (0->0), fcn. (18891->14), ass. (0->164)
t644 = Icges(6,1) + Icges(7,1);
t643 = Icges(6,4) + Icges(7,4);
t642 = Icges(6,5) + Icges(7,5);
t641 = Icges(6,2) + Icges(7,2);
t640 = Icges(6,6) + Icges(7,6);
t639 = Icges(6,3) + Icges(7,3);
t638 = rSges(7,3) + qJ(6);
t572 = sin(pkin(6));
t573 = cos(pkin(6));
t581 = cos(qJ(2));
t582 = cos(qJ(1));
t608 = t581 * t582;
t578 = sin(qJ(2));
t579 = sin(qJ(1));
t611 = t578 * t579;
t557 = t573 * t608 - t611;
t609 = t579 * t581;
t610 = t578 * t582;
t558 = t573 * t610 + t609;
t559 = -t573 * t609 - t610;
t560 = -t573 * t611 + t608;
t612 = t572 * t582;
t613 = t572 * t579;
t596 = (Icges(3,5) * t558 + Icges(3,6) * t557 - Icges(3,3) * t612) * t582 - (Icges(3,5) * t560 + Icges(3,6) * t559 + Icges(3,3) * t613) * t579;
t637 = t572 * t596;
t577 = sin(qJ(3));
t618 = sin(pkin(7));
t599 = t572 * t618;
t619 = cos(pkin(7));
t623 = cos(qJ(3));
t525 = t558 * t623 + (t557 * t619 - t582 * t599) * t577;
t576 = sin(qJ(4));
t600 = t572 * t619;
t592 = -t557 * t618 - t582 * t600;
t622 = cos(qJ(4));
t508 = t525 * t622 + t576 * t592;
t597 = t623 * t618;
t595 = t572 * t597;
t598 = t619 * t623;
t524 = -t557 * t598 + t558 * t577 + t582 * t595;
t575 = sin(qJ(5));
t580 = cos(qJ(5));
t478 = -t508 * t575 + t524 * t580;
t617 = t524 * t575;
t479 = t508 * t580 + t617;
t507 = t525 * t576 - t592 * t622;
t636 = t640 * t478 + t642 * t479 + t639 * t507;
t527 = t560 * t623 + (t559 * t619 + t579 * t599) * t577;
t591 = -t559 * t618 + t579 * t600;
t510 = t527 * t622 + t576 * t591;
t526 = -t559 * t598 + t560 * t577 - t579 * t595;
t480 = -t510 * t575 + t526 * t580;
t616 = t526 * t575;
t481 = t510 * t580 + t616;
t509 = t527 * t576 - t591 * t622;
t635 = t640 * t480 + t642 * t481 + t639 * t509;
t634 = t641 * t478 + t643 * t479 + t640 * t507;
t633 = t641 * t480 + t643 * t481 + t640 * t509;
t632 = t643 * t478 + t644 * t479 + t642 * t507;
t631 = t643 * t480 + t644 * t481 + t642 * t509;
t546 = t573 * t618 * t577 + (t577 * t581 * t619 + t578 * t623) * t572;
t590 = t573 * t619 - t581 * t599;
t523 = t546 * t622 + t576 * t590;
t614 = t572 * t578;
t545 = -t572 * t581 * t598 - t573 * t597 + t577 * t614;
t505 = -t523 * t575 + t545 * t580;
t615 = t545 * t575;
t506 = t523 * t580 + t615;
t522 = t546 * t576 - t590 * t622;
t630 = t640 * t505 + t642 * t506 + t639 * t522;
t629 = t641 * t505 + t643 * t506 + t640 * t522;
t628 = t643 * t505 + t644 * t506 + t642 * t522;
t621 = pkin(5) * t580;
t607 = rSges(7,1) * t479 + rSges(7,2) * t478 + pkin(5) * t617 + t638 * t507 + t508 * t621;
t606 = rSges(7,1) * t481 + rSges(7,2) * t480 + pkin(5) * t616 + t638 * t509 + t510 * t621;
t605 = rSges(7,1) * t506 + rSges(7,2) * t505 + pkin(5) * t615 + t638 * t522 + t523 * t621;
t528 = t558 * pkin(2) + pkin(10) * t592;
t529 = t560 * pkin(2) + pkin(10) * t591;
t602 = qJD(2) * t572;
t568 = t579 * t602;
t601 = t582 * t602;
t604 = t528 * t568 + t529 * t601;
t538 = qJD(3) * t591 + t568;
t603 = qJD(1) * (pkin(1) * t579 - pkin(9) * t612);
t569 = qJD(2) * t573 + qJD(1);
t501 = qJD(4) * t526 + t538;
t547 = qJD(3) * t590 + t569;
t516 = qJD(4) * t545 + t547;
t539 = qJD(3) * t592 - t601;
t497 = pkin(3) * t525 + pkin(11) * t524;
t498 = pkin(3) * t527 + pkin(11) * t526;
t594 = t538 * t497 - t498 * t539 + t604;
t502 = qJD(4) * t524 + t539;
t548 = pkin(2) * t614 + pkin(10) * t590;
t561 = qJD(1) * (pkin(1) * t582 + pkin(9) * t613);
t593 = t569 * t529 - t548 * t568 + t561;
t475 = pkin(4) * t508 + pkin(12) * t507;
t476 = pkin(4) * t510 + pkin(12) * t509;
t589 = t501 * t475 - t476 * t502 + t594;
t588 = -t528 * t569 - t548 * t601 - t603;
t515 = pkin(3) * t546 + pkin(11) * t545;
t587 = t547 * t498 - t515 * t538 + t593;
t586 = -t497 * t547 + t539 * t515 + t588;
t496 = pkin(4) * t523 + pkin(12) * t522;
t585 = t516 * t476 - t496 * t501 + t587;
t584 = -t475 * t516 + t502 * t496 + t586;
t566 = rSges(2,1) * t582 - rSges(2,2) * t579;
t565 = rSges(2,1) * t579 + rSges(2,2) * t582;
t554 = rSges(3,3) * t573 + (rSges(3,1) * t578 + rSges(3,2) * t581) * t572;
t553 = Icges(3,5) * t573 + (Icges(3,1) * t578 + Icges(3,4) * t581) * t572;
t552 = Icges(3,6) * t573 + (Icges(3,4) * t578 + Icges(3,2) * t581) * t572;
t551 = Icges(3,3) * t573 + (Icges(3,5) * t578 + Icges(3,6) * t581) * t572;
t537 = rSges(3,1) * t560 + rSges(3,2) * t559 + rSges(3,3) * t613;
t536 = rSges(3,1) * t558 + rSges(3,2) * t557 - rSges(3,3) * t612;
t535 = Icges(3,1) * t560 + Icges(3,4) * t559 + Icges(3,5) * t613;
t534 = Icges(3,1) * t558 + Icges(3,4) * t557 - Icges(3,5) * t612;
t533 = Icges(3,4) * t560 + Icges(3,2) * t559 + Icges(3,6) * t613;
t532 = Icges(3,4) * t558 + Icges(3,2) * t557 - Icges(3,6) * t612;
t514 = t546 * rSges(4,1) - t545 * rSges(4,2) + rSges(4,3) * t590;
t513 = Icges(4,1) * t546 - Icges(4,4) * t545 + Icges(4,5) * t590;
t512 = Icges(4,4) * t546 - Icges(4,2) * t545 + Icges(4,6) * t590;
t511 = Icges(4,5) * t546 - Icges(4,6) * t545 + Icges(4,3) * t590;
t504 = t537 * t569 - t554 * t568 + t561;
t503 = -t536 * t569 - t554 * t601 - t603;
t500 = (t536 * t579 + t537 * t582) * t602;
t495 = qJD(5) * t522 + t516;
t493 = t527 * rSges(4,1) - t526 * rSges(4,2) + rSges(4,3) * t591;
t492 = t525 * rSges(4,1) - t524 * rSges(4,2) + rSges(4,3) * t592;
t491 = Icges(4,1) * t527 - Icges(4,4) * t526 + Icges(4,5) * t591;
t490 = Icges(4,1) * t525 - Icges(4,4) * t524 + Icges(4,5) * t592;
t489 = Icges(4,4) * t527 - Icges(4,2) * t526 + Icges(4,6) * t591;
t488 = Icges(4,4) * t525 - Icges(4,2) * t524 + Icges(4,6) * t592;
t487 = Icges(4,5) * t527 - Icges(4,6) * t526 + Icges(4,3) * t591;
t486 = Icges(4,5) * t525 - Icges(4,6) * t524 + Icges(4,3) * t592;
t485 = rSges(5,1) * t523 - rSges(5,2) * t522 + rSges(5,3) * t545;
t484 = Icges(5,1) * t523 - Icges(5,4) * t522 + Icges(5,5) * t545;
t483 = Icges(5,4) * t523 - Icges(5,2) * t522 + Icges(5,6) * t545;
t482 = Icges(5,5) * t523 - Icges(5,6) * t522 + Icges(5,3) * t545;
t474 = qJD(5) * t507 + t502;
t473 = qJD(5) * t509 + t501;
t471 = rSges(5,1) * t510 - rSges(5,2) * t509 + rSges(5,3) * t526;
t470 = rSges(5,1) * t508 - rSges(5,2) * t507 + rSges(5,3) * t524;
t469 = Icges(5,1) * t510 - Icges(5,4) * t509 + Icges(5,5) * t526;
t468 = Icges(5,1) * t508 - Icges(5,4) * t507 + Icges(5,5) * t524;
t467 = Icges(5,4) * t510 - Icges(5,2) * t509 + Icges(5,6) * t526;
t466 = Icges(5,4) * t508 - Icges(5,2) * t507 + Icges(5,6) * t524;
t465 = Icges(5,5) * t510 - Icges(5,6) * t509 + Icges(5,3) * t526;
t464 = Icges(5,5) * t508 - Icges(5,6) * t507 + Icges(5,3) * t524;
t462 = rSges(6,1) * t506 + rSges(6,2) * t505 + rSges(6,3) * t522;
t452 = rSges(6,1) * t481 + rSges(6,2) * t480 + rSges(6,3) * t509;
t450 = rSges(6,1) * t479 + rSges(6,2) * t478 + rSges(6,3) * t507;
t448 = t493 * t547 - t514 * t538 + t593;
t447 = -t492 * t547 + t514 * t539 + t588;
t432 = t492 * t538 - t493 * t539 + t604;
t431 = t471 * t516 - t485 * t501 + t587;
t430 = -t470 * t516 + t485 * t502 + t586;
t429 = t470 * t501 - t471 * t502 + t594;
t428 = t452 * t495 - t462 * t473 + t585;
t427 = -t450 * t495 + t462 * t474 + t584;
t426 = t450 * t473 - t452 * t474 + t589;
t425 = qJD(6) * t507 - t473 * t605 + t495 * t606 + t585;
t424 = qJD(6) * t509 + t474 * t605 - t495 * t607 + t584;
t423 = qJD(6) * t522 + t473 * t607 - t474 * t606 + t589;
t1 = ((t551 * t613 + t552 * t559 + t553 * t560) * t569 + (-(t532 * t559 + t534 * t560) * t582 + (t559 * t533 + t560 * t535 - t637) * t579) * t602) * t568 / 0.2e1 - ((-t551 * t612 + t552 * t557 + t553 * t558) * t569 + ((t533 * t557 + t535 * t558) * t579 + (-t557 * t532 - t558 * t534 + t637) * t582) * t602) * t601 / 0.2e1 + t539 * ((t487 * t592 - t524 * t489 + t525 * t491) * t538 + (t592 * t486 - t524 * t488 + t525 * t490) * t539 + (t511 * t592 - t524 * t512 + t525 * t513) * t547) / 0.2e1 + t547 * ((t487 * t590 - t545 * t489 + t546 * t491) * t538 + (t486 * t590 - t545 * t488 + t546 * t490) * t539 + (t590 * t511 - t545 * t512 + t546 * t513) * t547) / 0.2e1 + t569 * ((t573 * t551 + (t552 * t581 + t553 * t578) * t572) * t569 + (((t533 * t581 + t535 * t578) * t579 - (t532 * t581 + t534 * t578) * t582) * t572 - t596 * t573) * t602) / 0.2e1 + m(6) * (t426 ^ 2 + t427 ^ 2 + t428 ^ 2) / 0.2e1 + m(7) * (t423 ^ 2 + t424 ^ 2 + t425 ^ 2) / 0.2e1 + m(5) * (t429 ^ 2 + t430 ^ 2 + t431 ^ 2) / 0.2e1 + m(4) * (t432 ^ 2 + t447 ^ 2 + t448 ^ 2) / 0.2e1 + m(3) * (t500 ^ 2 + t503 ^ 2 + t504 ^ 2) / 0.2e1 + t501 * ((t526 * t465 - t509 * t467 + t510 * t469) * t501 + (t464 * t526 - t466 * t509 + t468 * t510) * t502 + (t482 * t526 - t483 * t509 + t484 * t510) * t516) / 0.2e1 + t502 * ((t465 * t524 - t467 * t507 + t469 * t508) * t501 + (t524 * t464 - t507 * t466 + t508 * t468) * t502 + (t482 * t524 - t483 * t507 + t484 * t508) * t516) / 0.2e1 + t516 * ((t465 * t545 - t467 * t522 + t469 * t523) * t501 + (t464 * t545 - t466 * t522 + t468 * t523) * t502 + (t545 * t482 - t522 * t483 + t523 * t484) * t516) / 0.2e1 + t538 * ((t591 * t487 - t526 * t489 + t527 * t491) * t538 + (t486 * t591 - t526 * t488 + t527 * t490) * t539 + (t511 * t591 - t526 * t512 + t527 * t513) * t547) / 0.2e1 + ((t480 * t629 + t481 * t628 + t509 * t630) * t495 + (t480 * t634 + t481 * t632 + t509 * t636) * t474 + (t633 * t480 + t631 * t481 + t635 * t509) * t473) * t473 / 0.2e1 + ((t478 * t629 + t479 * t628 + t507 * t630) * t495 + (t634 * t478 + t632 * t479 + t636 * t507) * t474 + (t478 * t633 + t479 * t631 + t507 * t635) * t473) * t474 / 0.2e1 + ((t629 * t505 + t628 * t506 + t630 * t522) * t495 + (t505 * t634 + t506 * t632 + t522 * t636) * t474 + (t505 * t633 + t506 * t631 + t522 * t635) * t473) * t495 / 0.2e1 + (Icges(2,3) + m(2) * (t565 ^ 2 + t566 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
