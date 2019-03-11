% Calculate kinetic energy for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR13_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR13_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR13_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR13_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:21:41
% EndTime: 2019-03-09 04:21:43
% DurationCPUTime: 2.48s
% Computational Cost: add. (3845->301), mult. (10537->457), div. (0->0), fcn. (13551->14), ass. (0->146)
t610 = Icges(4,1) + Icges(5,2);
t609 = Icges(5,1) + Icges(4,3);
t608 = -Icges(4,4) - Icges(5,6);
t607 = -Icges(5,4) + Icges(4,5);
t606 = Icges(5,5) - Icges(4,6);
t605 = Icges(4,2) + Icges(5,3);
t551 = cos(pkin(6));
t587 = cos(pkin(12));
t589 = sin(qJ(1));
t572 = t589 * t587;
t549 = sin(pkin(12));
t556 = cos(qJ(1));
t583 = t556 * t549;
t537 = t551 * t583 + t572;
t554 = sin(qJ(3));
t574 = t556 * t587;
t580 = t589 * t549;
t563 = -t551 * t574 + t580;
t588 = cos(pkin(7));
t560 = t563 * t588;
t586 = sin(pkin(7));
t591 = cos(qJ(3));
t573 = t591 * t586;
t550 = sin(pkin(6));
t584 = t550 * t556;
t515 = t537 * t554 + t560 * t591 + t573 * t584;
t575 = t554 * t586;
t516 = t537 * t591 - t554 * t560 - t575 * t584;
t577 = t550 * t588;
t530 = -t556 * t577 + t563 * t586;
t604 = t515 * t605 + t516 * t608 + t530 * t606;
t538 = -t551 * t580 + t574;
t564 = t551 * t572 + t583;
t576 = t550 * t586;
t595 = t564 * t588 - t589 * t576;
t517 = t538 * t554 + t591 * t595;
t518 = t538 * t591 - t554 * t595;
t531 = t564 * t586 + t577 * t589;
t603 = t517 * t605 + t518 * t608 + t531 * t606;
t602 = t515 * t606 + t516 * t607 + t530 * t609;
t601 = t517 * t606 + t518 * t607 + t531 * t609;
t600 = t608 * t515 + t516 * t610 + t607 * t530;
t599 = t608 * t517 + t518 * t610 + t607 * t531;
t570 = t588 * t587;
t585 = t549 * t550;
t528 = -t550 * t570 * t591 - t551 * t573 + t554 * t585;
t529 = t551 * t575 + (t549 * t591 + t554 * t570) * t550;
t536 = t551 * t588 - t576 * t587;
t598 = t528 * t605 + t529 * t608 + t536 * t606;
t597 = t528 * t606 + t529 * t607 + t536 * t609;
t596 = t608 * t528 + t529 * t610 + t607 * t536;
t594 = t530 * t602 + t531 * t601;
t590 = cos(qJ(5));
t526 = qJD(3) * t530;
t491 = qJD(5) * t516 + t526;
t527 = qJD(3) * t531;
t492 = qJD(5) * t518 + t527;
t534 = qJD(3) * t536 + qJD(1);
t480 = pkin(3) * t516 + qJ(4) * t515;
t548 = qJD(2) * t551;
t582 = qJD(4) * t528 + t480 * t527 + t548;
t581 = t550 * t589;
t519 = qJD(5) * t529 + t534;
t571 = -qJD(2) * t584 + qJD(1) * (t556 * pkin(1) + qJ(2) * t581);
t540 = pkin(1) * t589 - qJ(2) * t584;
t547 = qJD(2) * t581;
t569 = t547 + (-t537 * pkin(2) - pkin(9) * t530 - t540) * qJD(1);
t567 = qJD(1) * (t538 * pkin(2) + pkin(9) * t531) + t571;
t506 = pkin(3) * t529 + qJ(4) * t528;
t566 = qJD(4) * t517 + t506 * t526 + t569;
t481 = pkin(3) * t518 + qJ(4) * t517;
t565 = qJD(4) * t515 + t534 * t481 + t567;
t493 = pkin(4) * t530 + pkin(10) * t516;
t494 = pkin(4) * t531 + pkin(10) * t518;
t562 = t493 * t527 + (-t481 - t494) * t526 + t582;
t520 = pkin(4) * t536 + pkin(10) * t529;
t559 = t520 * t526 + (-t480 - t493) * t534 + t566;
t558 = t534 * t494 + (-t506 - t520) * t527 + t565;
t555 = cos(qJ(6));
t553 = sin(qJ(5));
t552 = sin(qJ(6));
t545 = t556 * rSges(2,1) - rSges(2,2) * t589;
t544 = rSges(2,1) * t589 + t556 * rSges(2,2);
t514 = t528 * t553 + t536 * t590;
t513 = -t528 * t590 + t536 * t553;
t505 = qJD(1) * (t538 * rSges(3,1) - rSges(3,2) * t564 + rSges(3,3) * t581) + t571;
t504 = t547 + (-t537 * rSges(3,1) + rSges(3,2) * t563 + rSges(3,3) * t584 - t540) * qJD(1);
t502 = rSges(5,1) * t536 - rSges(5,2) * t529 + rSges(5,3) * t528;
t501 = rSges(4,1) * t529 - rSges(4,2) * t528 + rSges(4,3) * t536;
t490 = t517 * t553 + t531 * t590;
t489 = -t517 * t590 + t531 * t553;
t488 = t515 * t553 + t530 * t590;
t487 = -t515 * t590 + t530 * t553;
t486 = t514 * t555 + t529 * t552;
t485 = -t514 * t552 + t529 * t555;
t482 = qJD(6) * t513 + t519;
t479 = pkin(5) * t514 + pkin(11) * t513;
t476 = rSges(5,1) * t531 - rSges(5,2) * t518 + rSges(5,3) * t517;
t475 = rSges(5,1) * t530 - rSges(5,2) * t516 + rSges(5,3) * t515;
t474 = rSges(4,1) * t518 - rSges(4,2) * t517 + rSges(4,3) * t531;
t473 = rSges(4,1) * t516 - rSges(4,2) * t515 + rSges(4,3) * t530;
t459 = rSges(6,1) * t514 - rSges(6,2) * t513 + rSges(6,3) * t529;
t458 = Icges(6,1) * t514 - Icges(6,4) * t513 + Icges(6,5) * t529;
t457 = Icges(6,4) * t514 - Icges(6,2) * t513 + Icges(6,6) * t529;
t456 = Icges(6,5) * t514 - Icges(6,6) * t513 + Icges(6,3) * t529;
t455 = t490 * t555 + t518 * t552;
t454 = -t490 * t552 + t518 * t555;
t453 = t488 * t555 + t516 * t552;
t452 = -t488 * t552 + t516 * t555;
t451 = qJD(6) * t489 + t492;
t450 = qJD(6) * t487 + t491;
t449 = pkin(5) * t490 + pkin(11) * t489;
t448 = pkin(5) * t488 + pkin(11) * t487;
t447 = rSges(6,1) * t490 - rSges(6,2) * t489 + rSges(6,3) * t518;
t446 = rSges(6,1) * t488 - rSges(6,2) * t487 + rSges(6,3) * t516;
t445 = Icges(6,1) * t490 - Icges(6,4) * t489 + Icges(6,5) * t518;
t444 = Icges(6,1) * t488 - Icges(6,4) * t487 + Icges(6,5) * t516;
t443 = Icges(6,4) * t490 - Icges(6,2) * t489 + Icges(6,6) * t518;
t442 = Icges(6,4) * t488 - Icges(6,2) * t487 + Icges(6,6) * t516;
t441 = Icges(6,5) * t490 - Icges(6,6) * t489 + Icges(6,3) * t518;
t440 = Icges(6,5) * t488 - Icges(6,6) * t487 + Icges(6,3) * t516;
t439 = rSges(7,1) * t486 + rSges(7,2) * t485 + rSges(7,3) * t513;
t438 = Icges(7,1) * t486 + Icges(7,4) * t485 + Icges(7,5) * t513;
t437 = Icges(7,4) * t486 + Icges(7,2) * t485 + Icges(7,6) * t513;
t436 = Icges(7,5) * t486 + Icges(7,6) * t485 + Icges(7,3) * t513;
t435 = t474 * t534 - t501 * t527 + t567;
t434 = -t473 * t534 + t501 * t526 + t569;
t433 = t548 + (t473 * t531 - t474 * t530) * qJD(3);
t432 = rSges(7,1) * t455 + rSges(7,2) * t454 + rSges(7,3) * t489;
t431 = rSges(7,1) * t453 + rSges(7,2) * t452 + rSges(7,3) * t487;
t430 = Icges(7,1) * t455 + Icges(7,4) * t454 + Icges(7,5) * t489;
t429 = Icges(7,1) * t453 + Icges(7,4) * t452 + Icges(7,5) * t487;
t428 = Icges(7,4) * t455 + Icges(7,2) * t454 + Icges(7,6) * t489;
t427 = Icges(7,4) * t453 + Icges(7,2) * t452 + Icges(7,6) * t487;
t426 = Icges(7,5) * t455 + Icges(7,6) * t454 + Icges(7,3) * t489;
t425 = Icges(7,5) * t453 + Icges(7,6) * t452 + Icges(7,3) * t487;
t424 = t476 * t534 + (-t502 - t506) * t527 + t565;
t423 = t502 * t526 + (-t475 - t480) * t534 + t566;
t422 = (t475 * t531 + (-t476 - t481) * t530) * qJD(3) + t582;
t421 = t447 * t519 - t459 * t492 + t558;
t420 = -t446 * t519 + t459 * t491 + t559;
t419 = t446 * t492 - t447 * t491 + t562;
t418 = t432 * t482 - t439 * t451 + t449 * t519 - t479 * t492 + t558;
t417 = -t431 * t482 + t439 * t450 - t448 * t519 + t479 * t491 + t559;
t416 = t431 * t451 - t432 * t450 + t448 * t492 - t449 * t491 + t562;
t1 = m(7) * (t416 ^ 2 + t417 ^ 2 + t418 ^ 2) / 0.2e1 + m(6) * (t419 ^ 2 + t420 ^ 2 + t421 ^ 2) / 0.2e1 + m(5) * (t422 ^ 2 + t423 ^ 2 + t424 ^ 2) / 0.2e1 + m(4) * (t433 ^ 2 + t434 ^ 2 + t435 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 * t551 ^ 2 + t504 ^ 2 + t505 ^ 2) / 0.2e1 + t450 * ((t426 * t487 + t428 * t452 + t430 * t453) * t451 + (t425 * t487 + t427 * t452 + t429 * t453) * t450 + (t436 * t487 + t437 * t452 + t438 * t453) * t482) / 0.2e1 + t451 * ((t426 * t489 + t428 * t454 + t430 * t455) * t451 + (t425 * t489 + t427 * t454 + t429 * t455) * t450 + (t436 * t489 + t437 * t454 + t438 * t455) * t482) / 0.2e1 + t482 * ((t426 * t513 + t428 * t485 + t430 * t486) * t451 + (t425 * t513 + t427 * t485 + t429 * t486) * t450 + (t436 * t513 + t437 * t485 + t438 * t486) * t482) / 0.2e1 + t492 * ((t441 * t518 - t443 * t489 + t445 * t490) * t492 + (t440 * t518 - t442 * t489 + t444 * t490) * t491 + (t456 * t518 - t457 * t489 + t458 * t490) * t519) / 0.2e1 + t491 * ((t441 * t516 - t443 * t487 + t445 * t488) * t492 + (t440 * t516 - t442 * t487 + t444 * t488) * t491 + (t456 * t516 - t457 * t487 + t458 * t488) * t519) / 0.2e1 + t519 * ((t441 * t529 - t443 * t513 + t445 * t514) * t492 + (t440 * t529 - t442 * t513 + t444 * t514) * t491 + (t456 * t529 - t457 * t513 + t458 * t514) * t519) / 0.2e1 + ((t598 * t528 + t596 * t529 + t597 * t536) * t534 + ((t528 * t603 + t529 * t599 + t536 * t601) * t531 + (t528 * t604 + t600 * t529 + t602 * t536) * t530) * qJD(3)) * t534 / 0.2e1 + ((t515 * t598 + t516 * t596 + t530 * t597) * t534 + ((t515 * t603 + t516 * t599) * t531 + (t515 * t604 + t600 * t516 + t594) * t530) * qJD(3)) * t526 / 0.2e1 + ((t517 * t598 + t518 * t596 + t531 * t597) * t534 + ((t517 * t604 + t600 * t518) * t530 + (t517 * t603 + t518 * t599 + t594) * t531) * qJD(3)) * t527 / 0.2e1 + ((Icges(3,5) * t551 + (Icges(3,1) * t549 + Icges(3,4) * t587) * t550) * t585 + t550 * t587 * (Icges(3,6) * t551 + (Icges(3,4) * t549 + Icges(3,2) * t587) * t550) + t551 * (Icges(3,3) * t551 + (Icges(3,5) * t549 + Icges(3,6) * t587) * t550) + Icges(2,3) + m(2) * (t544 ^ 2 + t545 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
