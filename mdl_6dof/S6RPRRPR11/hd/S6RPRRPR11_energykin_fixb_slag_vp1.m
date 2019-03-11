% Calculate kinetic energy for
% S6RPRRPR11
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
% Datum: 2019-03-09 05:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR11_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR11_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_energykin_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR11_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR11_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR11_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:39:05
% EndTime: 2019-03-09 05:39:07
% DurationCPUTime: 2.83s
% Computational Cost: add. (5122->343), mult. (13559->516), div. (0->0), fcn. (17625->16), ass. (0->158)
t609 = Icges(5,2) + Icges(6,3);
t558 = cos(pkin(6));
t595 = cos(pkin(12));
t598 = sin(qJ(1));
t580 = t598 * t595;
t555 = sin(pkin(12));
t562 = cos(qJ(1));
t588 = t562 * t555;
t539 = t558 * t588 + t580;
t561 = sin(qJ(3));
t582 = t562 * t595;
t585 = t598 * t555;
t571 = -t558 * t582 + t585;
t596 = cos(pkin(7));
t565 = t571 * t596;
t556 = sin(pkin(6));
t594 = sin(pkin(7));
t583 = t556 * t594;
t600 = cos(qJ(3));
t519 = t539 * t600 + (-t562 * t583 - t565) * t561;
t584 = t556 * t596;
t532 = -t562 * t584 + t571 * t594;
t560 = sin(qJ(4));
t599 = cos(qJ(4));
t501 = t519 * t599 + t532 * t560;
t581 = t600 * t594;
t589 = t556 * t562;
t518 = t539 * t561 + t565 * t600 + t581 * t589;
t554 = sin(pkin(13));
t557 = cos(pkin(13));
t471 = -t501 * t554 + t518 * t557;
t593 = t518 * t554;
t472 = t501 * t557 + t593;
t500 = t519 * t560 - t532 * t599;
t608 = -Icges(5,4) * t501 + Icges(6,5) * t472 - Icges(5,6) * t518 + Icges(6,6) * t471 + t609 * t500;
t540 = -t558 * t585 + t582;
t572 = t558 * t580 + t588;
t605 = t572 * t596 - t598 * t583;
t521 = t540 * t600 - t561 * t605;
t533 = t572 * t594 + t584 * t598;
t503 = t521 * t599 + t533 * t560;
t520 = t540 * t561 + t600 * t605;
t473 = -t503 * t554 + t520 * t557;
t592 = t520 * t554;
t474 = t503 * t557 + t592;
t502 = t521 * t560 - t533 * t599;
t607 = -Icges(5,4) * t503 + Icges(6,5) * t474 - Icges(5,6) * t520 + Icges(6,6) * t473 + t609 * t502;
t578 = t596 * t595;
t531 = t558 * t594 * t561 + (t555 * t600 + t561 * t578) * t556;
t569 = t558 * t596 - t583 * t595;
t517 = t531 * t599 + t560 * t569;
t590 = t556 * t555;
t530 = -t556 * t578 * t600 - t558 * t581 + t561 * t590;
t496 = -t517 * t554 + t530 * t557;
t591 = t530 * t554;
t497 = t517 * t557 + t591;
t516 = t531 * t560 - t569 * t599;
t606 = -Icges(5,4) * t517 + Icges(6,5) * t497 - Icges(5,6) * t530 + Icges(6,6) * t496 + t609 * t516;
t597 = pkin(5) * t557;
t528 = qJD(3) * t532;
t504 = qJD(4) * t518 + t528;
t529 = qJD(3) * t533;
t505 = qJD(4) * t520 + t529;
t536 = qJD(3) * t569 + qJD(1);
t586 = t556 * t598;
t522 = qJD(4) * t530 + t536;
t579 = -qJD(2) * t589 + qJD(1) * (t562 * pkin(1) + qJ(2) * t586);
t542 = pkin(1) * t598 - qJ(2) * t589;
t548 = qJD(2) * t586;
t577 = t548 + (-t539 * pkin(2) - pkin(9) * t532 - t542) * qJD(1);
t575 = qJD(1) * (t540 * pkin(2) + pkin(9) * t533) + t579;
t490 = pkin(3) * t519 + pkin(10) * t518;
t491 = pkin(3) * t521 + pkin(10) * t520;
t550 = qJD(2) * t558;
t574 = t490 * t529 - t491 * t528 + t550;
t463 = pkin(4) * t501 + qJ(5) * t500;
t573 = qJD(5) * t516 + t505 * t463 + t574;
t512 = pkin(3) * t531 + pkin(10) * t530;
t570 = -t490 * t536 + t512 * t528 + t577;
t568 = t536 * t491 - t512 * t529 + t575;
t489 = pkin(4) * t517 + qJ(5) * t516;
t567 = qJD(5) * t502 + t504 * t489 + t570;
t464 = pkin(4) * t503 + qJ(5) * t502;
t564 = qJD(5) * t500 + t522 * t464 + t568;
t553 = pkin(13) + qJ(6);
t552 = cos(t553);
t551 = sin(t553);
t546 = t562 * rSges(2,1) - rSges(2,2) * t598;
t545 = rSges(2,1) * t598 + t562 * rSges(2,2);
t511 = qJD(1) * (t540 * rSges(3,1) - rSges(3,2) * t572 + rSges(3,3) * t586) + t579;
t510 = t548 + (-t539 * rSges(3,1) + rSges(3,2) * t571 + rSges(3,3) * t589 - t542) * qJD(1);
t509 = t531 * rSges(4,1) - t530 * rSges(4,2) + rSges(4,3) * t569;
t508 = Icges(4,1) * t531 - Icges(4,4) * t530 + Icges(4,5) * t569;
t507 = Icges(4,4) * t531 - Icges(4,2) * t530 + Icges(4,6) * t569;
t506 = Icges(4,5) * t531 - Icges(4,6) * t530 + Icges(4,3) * t569;
t494 = t517 * t552 + t530 * t551;
t493 = -t517 * t551 + t530 * t552;
t492 = qJD(6) * t516 + t522;
t486 = rSges(4,1) * t521 - rSges(4,2) * t520 + rSges(4,3) * t533;
t485 = rSges(4,1) * t519 - rSges(4,2) * t518 + rSges(4,3) * t532;
t484 = Icges(4,1) * t521 - Icges(4,4) * t520 + Icges(4,5) * t533;
t483 = Icges(4,1) * t519 - Icges(4,4) * t518 + Icges(4,5) * t532;
t482 = Icges(4,4) * t521 - Icges(4,2) * t520 + Icges(4,6) * t533;
t481 = Icges(4,4) * t519 - Icges(4,2) * t518 + Icges(4,6) * t532;
t480 = Icges(4,5) * t521 - Icges(4,6) * t520 + Icges(4,3) * t533;
t479 = Icges(4,5) * t519 - Icges(4,6) * t518 + Icges(4,3) * t532;
t478 = rSges(5,1) * t517 - rSges(5,2) * t516 + rSges(5,3) * t530;
t477 = Icges(5,1) * t517 - Icges(5,4) * t516 + Icges(5,5) * t530;
t475 = Icges(5,5) * t517 - Icges(5,6) * t516 + Icges(5,3) * t530;
t470 = t503 * t552 + t520 * t551;
t469 = -t503 * t551 + t520 * t552;
t468 = t501 * t552 + t518 * t551;
t467 = -t501 * t551 + t518 * t552;
t466 = qJD(6) * t502 + t505;
t465 = qJD(6) * t500 + t504;
t460 = rSges(5,1) * t503 - rSges(5,2) * t502 + rSges(5,3) * t520;
t459 = rSges(5,1) * t501 - rSges(5,2) * t500 + rSges(5,3) * t518;
t458 = Icges(5,1) * t503 - Icges(5,4) * t502 + Icges(5,5) * t520;
t457 = Icges(5,1) * t501 - Icges(5,4) * t500 + Icges(5,5) * t518;
t454 = Icges(5,5) * t503 - Icges(5,6) * t502 + Icges(5,3) * t520;
t453 = Icges(5,5) * t501 - Icges(5,6) * t500 + Icges(5,3) * t518;
t452 = rSges(6,1) * t497 + rSges(6,2) * t496 + rSges(6,3) * t516;
t451 = Icges(6,1) * t497 + Icges(6,4) * t496 + Icges(6,5) * t516;
t450 = Icges(6,4) * t497 + Icges(6,2) * t496 + Icges(6,6) * t516;
t448 = rSges(7,1) * t494 + rSges(7,2) * t493 + rSges(7,3) * t516;
t447 = Icges(7,1) * t494 + Icges(7,4) * t493 + Icges(7,5) * t516;
t446 = Icges(7,4) * t494 + Icges(7,2) * t493 + Icges(7,6) * t516;
t445 = Icges(7,5) * t494 + Icges(7,6) * t493 + Icges(7,3) * t516;
t444 = pkin(5) * t591 + pkin(11) * t516 + t517 * t597;
t442 = t486 * t536 - t509 * t529 + t575;
t441 = -t485 * t536 + t509 * t528 + t577;
t440 = t550 + (t485 * t533 - t486 * t532) * qJD(3);
t439 = rSges(6,1) * t474 + rSges(6,2) * t473 + rSges(6,3) * t502;
t438 = rSges(6,1) * t472 + rSges(6,2) * t471 + rSges(6,3) * t500;
t437 = Icges(6,1) * t474 + Icges(6,4) * t473 + Icges(6,5) * t502;
t436 = Icges(6,1) * t472 + Icges(6,4) * t471 + Icges(6,5) * t500;
t435 = Icges(6,4) * t474 + Icges(6,2) * t473 + Icges(6,6) * t502;
t434 = Icges(6,4) * t472 + Icges(6,2) * t471 + Icges(6,6) * t500;
t431 = rSges(7,1) * t470 + rSges(7,2) * t469 + rSges(7,3) * t502;
t430 = rSges(7,1) * t468 + rSges(7,2) * t467 + rSges(7,3) * t500;
t429 = Icges(7,1) * t470 + Icges(7,4) * t469 + Icges(7,5) * t502;
t428 = Icges(7,1) * t468 + Icges(7,4) * t467 + Icges(7,5) * t500;
t427 = Icges(7,4) * t470 + Icges(7,2) * t469 + Icges(7,6) * t502;
t426 = Icges(7,4) * t468 + Icges(7,2) * t467 + Icges(7,6) * t500;
t425 = Icges(7,5) * t470 + Icges(7,6) * t469 + Icges(7,3) * t502;
t424 = Icges(7,5) * t468 + Icges(7,6) * t467 + Icges(7,3) * t500;
t423 = pkin(5) * t592 + pkin(11) * t502 + t503 * t597;
t422 = pkin(5) * t593 + pkin(11) * t500 + t501 * t597;
t421 = t460 * t522 - t478 * t505 + t568;
t420 = -t459 * t522 + t478 * t504 + t570;
t419 = t459 * t505 - t460 * t504 + t574;
t418 = t439 * t522 + (-t452 - t489) * t505 + t564;
t417 = t452 * t504 + (-t438 - t463) * t522 + t567;
t416 = t438 * t505 + (-t439 - t464) * t504 + t573;
t415 = t423 * t522 + t431 * t492 - t448 * t466 + (-t444 - t489) * t505 + t564;
t414 = -t430 * t492 + t444 * t504 + t448 * t465 + (-t422 - t463) * t522 + t567;
t413 = t422 * t505 + t430 * t466 - t431 * t465 + (-t423 - t464) * t504 + t573;
t1 = ((t506 * t532 - t507 * t518 + t508 * t519) * t536 + ((t480 * t532 - t518 * t482 + t519 * t484) * t533 + (t479 * t532 - t481 * t518 + t483 * t519) * t532) * qJD(3)) * t528 / 0.2e1 + t536 * ((t506 * t569 - t530 * t507 + t531 * t508) * t536 + ((t480 * t569 - t530 * t482 + t531 * t484) * t533 + (t479 * t569 - t530 * t481 + t531 * t483) * t532) * qJD(3)) / 0.2e1 + m(7) * (t413 ^ 2 + t414 ^ 2 + t415 ^ 2) / 0.2e1 + m(6) * (t416 ^ 2 + t417 ^ 2 + t418 ^ 2) / 0.2e1 + m(5) * (t419 ^ 2 + t420 ^ 2 + t421 ^ 2) / 0.2e1 + m(4) * (t440 ^ 2 + t441 ^ 2 + t442 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 * t558 ^ 2 + t510 ^ 2 + t511 ^ 2) / 0.2e1 + t466 * ((t502 * t425 + t469 * t427 + t470 * t429) * t466 + (t424 * t502 + t426 * t469 + t428 * t470) * t465 + (t445 * t502 + t446 * t469 + t447 * t470) * t492) / 0.2e1 + t465 * ((t425 * t500 + t427 * t467 + t429 * t468) * t466 + (t500 * t424 + t467 * t426 + t468 * t428) * t465 + (t445 * t500 + t446 * t467 + t447 * t468) * t492) / 0.2e1 + t492 * ((t425 * t516 + t427 * t493 + t429 * t494) * t466 + (t424 * t516 + t426 * t493 + t428 * t494) * t465 + (t516 * t445 + t493 * t446 + t494 * t447) * t492) / 0.2e1 + ((t506 * t533 - t507 * t520 + t508 * t521) * t536 + ((t480 * t533 - t482 * t520 + t484 * t521) * t533 + (t479 * t533 - t481 * t520 + t483 * t521) * t532) * qJD(3)) * t529 / 0.2e1 + ((t450 * t471 + t451 * t472 + t475 * t518 + t477 * t501 + t500 * t606) * t522 + (t435 * t471 + t437 * t472 + t454 * t518 + t458 * t501 + t500 * t607) * t505 + (t471 * t434 + t472 * t436 + t518 * t453 + t501 * t457 + t500 * t608) * t504) * t504 / 0.2e1 + ((t450 * t473 + t451 * t474 + t475 * t520 + t477 * t503 + t502 * t606) * t522 + (t473 * t435 + t474 * t437 + t520 * t454 + t503 * t458 + t502 * t607) * t505 + (t434 * t473 + t436 * t474 + t453 * t520 + t457 * t503 + t502 * t608) * t504) * t505 / 0.2e1 + ((t496 * t450 + t497 * t451 + t530 * t475 + t517 * t477 + t516 * t606) * t522 + (t435 * t496 + t437 * t497 + t454 * t530 + t458 * t517 + t516 * t607) * t505 + (t434 * t496 + t436 * t497 + t453 * t530 + t457 * t517 + t516 * t608) * t504) * t522 / 0.2e1 + ((Icges(3,5) * t558 + (Icges(3,1) * t555 + Icges(3,4) * t595) * t556) * t590 + t556 * t595 * (Icges(3,6) * t558 + (Icges(3,4) * t555 + Icges(3,2) * t595) * t556) + t558 * (Icges(3,3) * t558 + (Icges(3,5) * t555 + Icges(3,6) * t595) * t556) + Icges(2,3) + m(2) * (t545 ^ 2 + t546 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
