% Calculate kinetic energy for
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPRRRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_energykin_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRR1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRRRR1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:59:42
% EndTime: 2019-03-08 18:59:43
% DurationCPUTime: 1.65s
% Computational Cost: add. (4958->321), mult. (12297->506), div. (0->0), fcn. (15943->16), ass. (0->153)
t599 = cos(qJ(3));
t564 = cos(qJ(4));
t598 = t564 * pkin(4);
t596 = cos(pkin(6));
t595 = cos(pkin(7));
t594 = cos(pkin(13));
t593 = sin(pkin(7));
t592 = sin(pkin(13));
t559 = cos(pkin(12));
t557 = sin(pkin(12));
t581 = t596 * t594;
t572 = t557 * t592 - t559 * t581;
t558 = sin(pkin(6));
t583 = t558 * t595;
t539 = -t559 * t583 + t572 * t593;
t561 = sin(qJ(4));
t591 = t539 * t561;
t571 = t557 * t581 + t559 * t592;
t540 = t557 * t583 + t571 * t593;
t590 = t540 * t561;
t582 = t558 * t593;
t545 = -t582 * t594 + t595 * t596;
t589 = t545 * t561;
t588 = qJ(4) + qJ(5);
t579 = t596 * t592;
t546 = t557 * t594 + t559 * t579;
t562 = sin(qJ(3));
t568 = t572 * t595;
t577 = t599 * t582;
t524 = t546 * t562 + t559 * t577 + t568 * t599;
t535 = qJD(3) * t539;
t510 = qJD(4) * t524 + t535;
t547 = -t557 * t579 + t559 * t594;
t567 = t571 * t595;
t526 = t547 * t562 - t557 * t577 + t567 * t599;
t536 = qJD(3) * t540;
t511 = qJD(4) * t526 + t536;
t578 = t595 * t594;
t580 = t596 * t593;
t537 = -t580 * t599 + (t562 * t592 - t578 * t599) * t558;
t544 = qJD(3) * t545;
t530 = qJD(4) * t537 + t544;
t587 = qJD(2) * t558;
t552 = qJD(2) * t596 + qJD(1);
t479 = qJD(5) * t524 + t510;
t480 = qJD(5) * t526 + t511;
t512 = qJD(5) * t537 + t530;
t585 = t559 * t587;
t584 = cos(t588);
t525 = t546 * t599 + (-t559 * t582 - t568) * t562;
t497 = t525 * pkin(3) + t524 * pkin(9);
t538 = t562 * t580 + (t562 * t578 + t592 * t599) * t558;
t517 = t538 * pkin(3) + t537 * pkin(9);
t551 = t557 * t587;
t576 = -t497 * t544 + t517 * t535 + t551;
t527 = t547 * t599 + (t557 * t582 - t567) * t562;
t498 = t527 * pkin(3) + t526 * pkin(9);
t575 = t497 * t536 - t498 * t535 + t552;
t574 = t498 * t544 - t517 * t536 - t585;
t443 = pkin(4) * t591 + pkin(10) * t524 + t525 * t598;
t474 = pkin(4) * t589 + pkin(10) * t537 + t538 * t598;
t573 = -t530 * t443 + t510 * t474 + t576;
t444 = pkin(4) * t590 + pkin(10) * t526 + t527 * t598;
t570 = t511 * t443 - t510 * t444 + t575;
t569 = t530 * t444 - t511 * t474 + t574;
t563 = cos(qJ(6));
t560 = sin(qJ(6));
t555 = sin(t588);
t529 = t538 * t564 + t589;
t528 = -t538 * t561 + t545 * t564;
t519 = t538 * t584 + t545 * t555;
t518 = t538 * t555 - t545 * t584;
t516 = t538 * rSges(4,1) - t537 * rSges(4,2) + t545 * rSges(4,3);
t515 = Icges(4,1) * t538 - Icges(4,4) * t537 + Icges(4,5) * t545;
t514 = Icges(4,4) * t538 - Icges(4,2) * t537 + Icges(4,6) * t545;
t513 = Icges(4,5) * t538 - Icges(4,6) * t537 + Icges(4,3) * t545;
t509 = t527 * t564 + t590;
t508 = -t527 * t561 + t540 * t564;
t507 = t525 * t564 + t591;
t506 = -t525 * t561 + t539 * t564;
t505 = t519 * t563 + t537 * t560;
t504 = -t519 * t560 + t537 * t563;
t503 = t527 * t584 + t540 * t555;
t502 = t527 * t555 - t540 * t584;
t501 = t525 * t584 + t539 * t555;
t500 = t525 * t555 - t539 * t584;
t496 = t519 * pkin(5) + t518 * pkin(11);
t494 = qJD(6) * t518 + t512;
t493 = t529 * rSges(5,1) + t528 * rSges(5,2) + t537 * rSges(5,3);
t492 = Icges(5,1) * t529 + Icges(5,4) * t528 + Icges(5,5) * t537;
t491 = Icges(5,4) * t529 + Icges(5,2) * t528 + Icges(5,6) * t537;
t490 = Icges(5,5) * t529 + Icges(5,6) * t528 + Icges(5,3) * t537;
t488 = t527 * rSges(4,1) - t526 * rSges(4,2) + t540 * rSges(4,3);
t487 = t525 * rSges(4,1) - t524 * rSges(4,2) + t539 * rSges(4,3);
t486 = Icges(4,1) * t527 - Icges(4,4) * t526 + Icges(4,5) * t540;
t485 = Icges(4,1) * t525 - Icges(4,4) * t524 + Icges(4,5) * t539;
t484 = Icges(4,4) * t527 - Icges(4,2) * t526 + Icges(4,6) * t540;
t483 = Icges(4,4) * t525 - Icges(4,2) * t524 + Icges(4,6) * t539;
t482 = Icges(4,5) * t527 - Icges(4,6) * t526 + Icges(4,3) * t540;
t481 = Icges(4,5) * t525 - Icges(4,6) * t524 + Icges(4,3) * t539;
t478 = t519 * rSges(6,1) - t518 * rSges(6,2) + t537 * rSges(6,3);
t477 = Icges(6,1) * t519 - Icges(6,4) * t518 + Icges(6,5) * t537;
t476 = Icges(6,4) * t519 - Icges(6,2) * t518 + Icges(6,6) * t537;
t475 = Icges(6,5) * t519 - Icges(6,6) * t518 + Icges(6,3) * t537;
t473 = t503 * t563 + t526 * t560;
t472 = -t503 * t560 + t526 * t563;
t471 = t501 * t563 + t524 * t560;
t470 = -t501 * t560 + t524 * t563;
t469 = pkin(5) * t503 + pkin(11) * t502;
t468 = pkin(5) * t501 + pkin(11) * t500;
t467 = qJD(6) * t502 + t480;
t466 = qJD(6) * t500 + t479;
t465 = t509 * rSges(5,1) + t508 * rSges(5,2) + t526 * rSges(5,3);
t464 = t507 * rSges(5,1) + t506 * rSges(5,2) + t524 * rSges(5,3);
t463 = Icges(5,1) * t509 + Icges(5,4) * t508 + Icges(5,5) * t526;
t462 = Icges(5,1) * t507 + Icges(5,4) * t506 + Icges(5,5) * t524;
t461 = Icges(5,4) * t509 + Icges(5,2) * t508 + Icges(5,6) * t526;
t460 = Icges(5,4) * t507 + Icges(5,2) * t506 + Icges(5,6) * t524;
t459 = Icges(5,5) * t509 + Icges(5,6) * t508 + Icges(5,3) * t526;
t458 = Icges(5,5) * t507 + Icges(5,6) * t506 + Icges(5,3) * t524;
t456 = t503 * rSges(6,1) - t502 * rSges(6,2) + t526 * rSges(6,3);
t455 = t501 * rSges(6,1) - t500 * rSges(6,2) + t524 * rSges(6,3);
t454 = Icges(6,1) * t503 - Icges(6,4) * t502 + Icges(6,5) * t526;
t453 = Icges(6,1) * t501 - Icges(6,4) * t500 + Icges(6,5) * t524;
t452 = Icges(6,4) * t503 - Icges(6,2) * t502 + Icges(6,6) * t526;
t451 = Icges(6,4) * t501 - Icges(6,2) * t500 + Icges(6,6) * t524;
t450 = Icges(6,5) * t503 - Icges(6,6) * t502 + Icges(6,3) * t526;
t449 = Icges(6,5) * t501 - Icges(6,6) * t500 + Icges(6,3) * t524;
t448 = t505 * rSges(7,1) + t504 * rSges(7,2) + t518 * rSges(7,3);
t447 = Icges(7,1) * t505 + Icges(7,4) * t504 + Icges(7,5) * t518;
t446 = Icges(7,4) * t505 + Icges(7,2) * t504 + Icges(7,6) * t518;
t445 = Icges(7,5) * t505 + Icges(7,6) * t504 + Icges(7,3) * t518;
t442 = -t585 + (t488 * t545 - t516 * t540) * qJD(3);
t441 = t551 + (-t487 * t545 + t516 * t539) * qJD(3);
t438 = (t487 * t540 - t488 * t539) * qJD(3) + t552;
t437 = rSges(7,1) * t473 + rSges(7,2) * t472 + rSges(7,3) * t502;
t436 = rSges(7,1) * t471 + rSges(7,2) * t470 + rSges(7,3) * t500;
t435 = Icges(7,1) * t473 + Icges(7,4) * t472 + Icges(7,5) * t502;
t434 = Icges(7,1) * t471 + Icges(7,4) * t470 + Icges(7,5) * t500;
t433 = Icges(7,4) * t473 + Icges(7,2) * t472 + Icges(7,6) * t502;
t432 = Icges(7,4) * t471 + Icges(7,2) * t470 + Icges(7,6) * t500;
t431 = Icges(7,5) * t473 + Icges(7,6) * t472 + Icges(7,3) * t502;
t430 = Icges(7,5) * t471 + Icges(7,6) * t470 + Icges(7,3) * t500;
t429 = t530 * t465 - t511 * t493 + t574;
t428 = -t530 * t464 + t510 * t493 + t576;
t427 = t511 * t464 - t510 * t465 + t575;
t426 = t512 * t456 - t480 * t478 + t569;
t425 = -t512 * t455 + t479 * t478 + t573;
t424 = t480 * t455 - t479 * t456 + t570;
t423 = t494 * t437 - t467 * t448 + t512 * t469 - t480 * t496 + t569;
t422 = -t494 * t436 + t466 * t448 - t512 * t468 + t479 * t496 + t573;
t421 = t467 * t436 - t466 * t437 + t480 * t468 - t479 * t469 + t570;
t1 = t466 * ((t431 * t500 + t433 * t470 + t435 * t471) * t467 + (t430 * t500 + t432 * t470 + t434 * t471) * t466 + (t445 * t500 + t446 * t470 + t447 * t471) * t494) / 0.2e1 + t494 * ((t518 * t431 + t504 * t433 + t505 * t435) * t467 + (t518 * t430 + t504 * t432 + t505 * t434) * t466 + (t518 * t445 + t504 * t446 + t505 * t447) * t494) / 0.2e1 + t467 * ((t431 * t502 + t433 * t472 + t435 * t473) * t467 + (t430 * t502 + t432 * t472 + t434 * t473) * t466 + (t445 * t502 + t446 * t472 + t447 * t473) * t494) / 0.2e1 + t480 * ((t526 * t450 - t502 * t452 + t503 * t454) * t480 + (t526 * t449 - t502 * t451 + t503 * t453) * t479 + (t526 * t475 - t502 * t476 + t503 * t477) * t512) / 0.2e1 + t479 * ((t524 * t450 - t500 * t452 + t501 * t454) * t480 + (t524 * t449 - t500 * t451 + t501 * t453) * t479 + (t524 * t475 - t500 * t476 + t501 * t477) * t512) / 0.2e1 + t512 * ((t537 * t450 - t518 * t452 + t519 * t454) * t480 + (t537 * t449 - t518 * t451 + t519 * t453) * t479 + (t537 * t475 - t518 * t476 + t519 * t477) * t512) / 0.2e1 + t511 * ((t526 * t459 + t508 * t461 + t509 * t463) * t511 + (t526 * t458 + t508 * t460 + t509 * t462) * t510 + (t526 * t490 + t508 * t491 + t509 * t492) * t530) / 0.2e1 + t530 * ((t537 * t459 + t528 * t461 + t529 * t463) * t511 + (t537 * t458 + t528 * t460 + t529 * t462) * t510 + (t537 * t490 + t528 * t491 + t529 * t492) * t530) / 0.2e1 + t510 * ((t524 * t459 + t506 * t461 + t507 * t463) * t511 + (t524 * t458 + t506 * t460 + t507 * t462) * t510 + (t524 * t490 + t506 * t491 + t507 * t492) * t530) / 0.2e1 + m(7) * (t421 ^ 2 + t422 ^ 2 + t423 ^ 2) / 0.2e1 + m(6) * (t424 ^ 2 + t425 ^ 2 + t426 ^ 2) / 0.2e1 + m(5) * (t427 ^ 2 + t428 ^ 2 + t429 ^ 2) / 0.2e1 + m(4) * (t438 ^ 2 + t441 ^ 2 + t442 ^ 2) / 0.2e1 + m(3) * (t552 ^ 2 + (t557 ^ 2 + t559 ^ 2) * qJD(2) ^ 2 * t558 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + (t540 * ((t540 * t482 - t526 * t484 + t527 * t486) * t540 + (t540 * t481 - t526 * t483 + t527 * t485) * t539 + (t540 * t513 - t526 * t514 + t527 * t515) * t545) + t539 * ((t539 * t482 - t524 * t484 + t525 * t486) * t540 + (t539 * t481 - t524 * t483 + t525 * t485) * t539 + (t539 * t513 - t524 * t514 + t525 * t515) * t545) + t545 * ((t545 * t482 - t537 * t484 + t538 * t486) * t540 + (t545 * t481 - t537 * t483 + t538 * t485) * t539 + (t545 * t513 - t537 * t514 + t538 * t515) * t545)) * qJD(3) ^ 2 / 0.2e1;
T  = t1;
