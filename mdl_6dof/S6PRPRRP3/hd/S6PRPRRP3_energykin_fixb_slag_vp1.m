% Calculate kinetic energy for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRP3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:04:41
% EndTime: 2019-03-08 20:04:43
% DurationCPUTime: 2.95s
% Computational Cost: add. (3114->330), mult. (5657->504), div. (0->0), fcn. (6851->12), ass. (0->155)
t585 = Icges(6,1) + Icges(7,1);
t584 = Icges(6,4) + Icges(7,4);
t583 = Icges(6,5) + Icges(7,5);
t582 = Icges(6,2) + Icges(7,2);
t581 = Icges(6,6) + Icges(7,6);
t580 = Icges(6,3) + Icges(7,3);
t579 = rSges(7,3) + qJ(6);
t515 = sin(pkin(10));
t518 = cos(pkin(10));
t525 = cos(qJ(2));
t519 = cos(pkin(6));
t523 = sin(qJ(2));
t554 = t519 * t523;
t501 = t515 * t525 + t518 * t554;
t544 = pkin(11) + qJ(4);
t513 = sin(t544);
t537 = cos(t544);
t516 = sin(pkin(6));
t557 = t516 * t518;
t477 = t501 * t537 - t513 * t557;
t553 = t519 * t525;
t500 = t515 * t523 - t518 * t553;
t522 = sin(qJ(5));
t524 = cos(qJ(5));
t448 = -t477 * t522 + t500 * t524;
t561 = t500 * t522;
t449 = t477 * t524 + t561;
t532 = t516 * t537;
t476 = t501 * t513 + t518 * t532;
t578 = t448 * t581 + t449 * t583 + t476 * t580;
t503 = -t515 * t554 + t518 * t525;
t558 = t515 * t516;
t479 = t503 * t537 + t513 * t558;
t502 = t515 * t553 + t518 * t523;
t450 = -t479 * t522 + t502 * t524;
t560 = t502 * t522;
t451 = t479 * t524 + t560;
t478 = t503 * t513 - t515 * t532;
t577 = t450 * t581 + t451 * t583 + t478 * t580;
t576 = t448 * t582 + t449 * t584 + t476 * t581;
t575 = t450 * t582 + t451 * t584 + t478 * t581;
t574 = t448 * t584 + t449 * t585 + t476 * t583;
t573 = t450 * t584 + t451 * t585 + t478 * t583;
t492 = t513 * t519 + t523 * t532;
t555 = t516 * t525;
t480 = -t492 * t522 - t524 * t555;
t541 = t522 * t555;
t481 = t492 * t524 - t541;
t556 = t516 * t523;
t491 = t513 * t556 - t519 * t537;
t572 = t480 * t581 + t481 * t583 + t491 * t580;
t571 = t480 * t582 + t481 * t584 + t491 * t581;
t570 = t480 * t584 + t481 * t585 + t491 * t583;
t514 = sin(pkin(11));
t517 = cos(pkin(11));
t498 = -t514 * t556 + t517 * t519;
t559 = t514 * t519;
t499 = t517 * t556 + t559;
t456 = Icges(4,5) * t499 + Icges(4,6) * t498 - Icges(4,3) * t555;
t489 = Icges(3,6) * t519 + (Icges(3,4) * t523 + Icges(3,2) * t525) * t516;
t569 = t456 - t489;
t568 = qJD(2) ^ 2;
t564 = pkin(3) * t517;
t563 = pkin(5) * t524;
t551 = rSges(7,1) * t449 + rSges(7,2) * t448 + pkin(5) * t561 + t476 * t579 + t477 * t563;
t550 = rSges(7,1) * t451 + rSges(7,2) * t450 + pkin(5) * t560 + t478 * t579 + t479 * t563;
t549 = rSges(7,1) * t481 + rSges(7,2) * t480 - pkin(5) * t541 + t491 * t579 + t492 * t563;
t474 = pkin(2) * t503 + qJ(3) * t502;
t512 = qJD(2) * t519;
t548 = qJD(3) * t500 + t474 * t512;
t547 = qJD(2) * t516;
t509 = t515 * t547;
t486 = qJD(4) * t502 + t509;
t546 = qJD(3) * t525;
t543 = t514 * t558;
t542 = t514 * t557;
t540 = t518 * t547;
t473 = pkin(2) * t501 + qJ(3) * t500;
t539 = t473 * t509 + t474 * t540 + qJD(1);
t504 = (pkin(2) * t523 - qJ(3) * t525) * t516;
t536 = (-rSges(4,1) * t499 - rSges(4,2) * t498 + rSges(4,3) * t555 - t504) * t516;
t535 = (-pkin(3) * t559 - (-pkin(8) * t525 + t523 * t564) * t516 - t504) * t516;
t487 = qJD(4) * t500 - t540;
t505 = -qJD(4) * t555 + t512;
t430 = -pkin(3) * t542 + pkin(8) * t500 + t501 * t564;
t431 = pkin(3) * t543 + pkin(8) * t502 + t503 * t564;
t531 = t430 * t509 + t431 * t540 - t516 * t546 + t539;
t530 = qJD(2) * t515 * t535 + t431 * t512 + t548;
t444 = pkin(4) * t477 + pkin(9) * t476;
t445 = pkin(4) * t479 + pkin(9) * t478;
t529 = t444 * t486 - t487 * t445 + t531;
t497 = qJD(3) * t502;
t528 = t497 + ((-t430 - t473) * t519 + t518 * t535) * qJD(2);
t466 = pkin(4) * t492 + pkin(9) * t491;
t527 = t445 * t505 - t466 * t486 + t530;
t526 = -t444 * t505 + t487 * t466 + t528;
t493 = t519 * rSges(3,3) + (rSges(3,1) * t523 + rSges(3,2) * t525) * t516;
t490 = Icges(3,5) * t519 + (Icges(3,1) * t523 + Icges(3,4) * t525) * t516;
t488 = Icges(3,3) * t519 + (Icges(3,5) * t523 + Icges(3,6) * t525) * t516;
t485 = t503 * t517 + t543;
t484 = -t503 * t514 + t517 * t558;
t483 = t501 * t517 - t542;
t482 = -t501 * t514 - t517 * t557;
t475 = qJD(5) * t491 + t505;
t471 = rSges(3,1) * t503 - rSges(3,2) * t502 + rSges(3,3) * t558;
t470 = rSges(3,1) * t501 - rSges(3,2) * t500 - rSges(3,3) * t557;
t464 = Icges(3,1) * t503 - Icges(3,4) * t502 + Icges(3,5) * t558;
t463 = Icges(3,1) * t501 - Icges(3,4) * t500 - Icges(3,5) * t557;
t462 = Icges(3,4) * t503 - Icges(3,2) * t502 + Icges(3,6) * t558;
t461 = Icges(3,4) * t501 - Icges(3,2) * t500 - Icges(3,6) * t557;
t460 = Icges(3,5) * t503 - Icges(3,6) * t502 + Icges(3,3) * t558;
t459 = Icges(3,5) * t501 - Icges(3,6) * t500 - Icges(3,3) * t557;
t458 = Icges(4,1) * t499 + Icges(4,4) * t498 - Icges(4,5) * t555;
t457 = Icges(4,4) * t499 + Icges(4,2) * t498 - Icges(4,6) * t555;
t455 = rSges(5,1) * t492 - rSges(5,2) * t491 - rSges(5,3) * t555;
t454 = Icges(5,1) * t492 - Icges(5,4) * t491 - Icges(5,5) * t555;
t453 = Icges(5,4) * t492 - Icges(5,2) * t491 - Icges(5,6) * t555;
t452 = Icges(5,5) * t492 - Icges(5,6) * t491 - Icges(5,3) * t555;
t447 = qJD(5) * t476 + t487;
t446 = qJD(5) * t478 + t486;
t442 = (-t470 * t519 - t493 * t557) * qJD(2);
t441 = (t471 * t519 - t493 * t558) * qJD(2);
t440 = rSges(4,1) * t485 + rSges(4,2) * t484 + rSges(4,3) * t502;
t439 = rSges(4,1) * t483 + rSges(4,2) * t482 + rSges(4,3) * t500;
t438 = Icges(4,1) * t485 + Icges(4,4) * t484 + Icges(4,5) * t502;
t437 = Icges(4,1) * t483 + Icges(4,4) * t482 + Icges(4,5) * t500;
t436 = Icges(4,4) * t485 + Icges(4,2) * t484 + Icges(4,6) * t502;
t435 = Icges(4,4) * t483 + Icges(4,2) * t482 + Icges(4,6) * t500;
t434 = Icges(4,5) * t485 + Icges(4,6) * t484 + Icges(4,3) * t502;
t433 = Icges(4,5) * t483 + Icges(4,6) * t482 + Icges(4,3) * t500;
t429 = rSges(5,1) * t479 - rSges(5,2) * t478 + rSges(5,3) * t502;
t428 = rSges(5,1) * t477 - rSges(5,2) * t476 + rSges(5,3) * t500;
t427 = Icges(5,1) * t479 - Icges(5,4) * t478 + Icges(5,5) * t502;
t426 = Icges(5,1) * t477 - Icges(5,4) * t476 + Icges(5,5) * t500;
t425 = Icges(5,4) * t479 - Icges(5,2) * t478 + Icges(5,6) * t502;
t424 = Icges(5,4) * t477 - Icges(5,2) * t476 + Icges(5,6) * t500;
t423 = Icges(5,5) * t479 - Icges(5,6) * t478 + Icges(5,3) * t502;
t422 = Icges(5,5) * t477 - Icges(5,6) * t476 + Icges(5,3) * t500;
t421 = rSges(6,1) * t481 + rSges(6,2) * t480 + rSges(6,3) * t491;
t408 = qJD(1) + (t470 * t515 + t471 * t518) * t547;
t407 = rSges(6,1) * t451 + rSges(6,2) * t450 + rSges(6,3) * t478;
t405 = rSges(6,1) * t449 + rSges(6,2) * t448 + rSges(6,3) * t476;
t389 = t497 + ((-t439 - t473) * t519 + t518 * t536) * qJD(2);
t388 = (t440 * t519 + t515 * t536) * qJD(2) + t548;
t387 = (-t546 + (t439 * t515 + t440 * t518) * qJD(2)) * t516 + t539;
t386 = -t428 * t505 + t455 * t487 + t528;
t385 = t429 * t505 - t455 * t486 + t530;
t384 = t428 * t486 - t429 * t487 + t531;
t383 = -t405 * t475 + t421 * t447 + t526;
t382 = t407 * t475 - t421 * t446 + t527;
t381 = t405 * t446 - t407 * t447 + t529;
t380 = qJD(6) * t478 + t447 * t549 - t475 * t551 + t526;
t379 = qJD(6) * t476 - t446 * t549 + t475 * t550 + t527;
t378 = qJD(6) * t491 + t446 * t551 - t447 * t550 + t529;
t1 = t486 * ((t502 * t423 - t478 * t425 + t479 * t427) * t486 + (t422 * t502 - t424 * t478 + t426 * t479) * t487 + (t452 * t502 - t453 * t478 + t454 * t479) * t505) / 0.2e1 + t487 * ((t423 * t500 - t425 * t476 + t427 * t477) * t486 + (t422 * t500 - t424 * t476 + t426 * t477) * t487 + (t452 * t500 - t453 * t476 + t454 * t477) * t505) / 0.2e1 + t505 * ((-t423 * t555 - t425 * t491 + t427 * t492) * t486 + (-t422 * t555 - t424 * t491 + t426 * t492) * t487 + (-t452 * t555 - t491 * t453 + t492 * t454) * t505) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t408 ^ 2 + t441 ^ 2 + t442 ^ 2) / 0.2e1 + m(4) * (t387 ^ 2 + t388 ^ 2 + t389 ^ 2) / 0.2e1 + m(5) * (t384 ^ 2 + t385 ^ 2 + t386 ^ 2) / 0.2e1 + m(6) * (t381 ^ 2 + t382 ^ 2 + t383 ^ 2) / 0.2e1 + m(7) * (t378 ^ 2 + t379 ^ 2 + t380 ^ 2) / 0.2e1 + ((t450 * t571 + t451 * t570 + t478 * t572) * t475 + (t450 * t576 + t451 * t574 + t478 * t578) * t447 + (t575 * t450 + t573 * t451 + t577 * t478) * t446) * t446 / 0.2e1 + ((t448 * t571 + t449 * t570 + t476 * t572) * t475 + (t576 * t448 + t574 * t449 + t578 * t476) * t447 + (t448 * t575 + t449 * t573 + t476 * t577) * t446) * t447 / 0.2e1 + ((t571 * t480 + t570 * t481 + t572 * t491) * t475 + (t480 * t576 + t481 * t574 + t491 * t578) * t447 + (t480 * t575 + t481 * t573 + t491 * t577) * t446) * t475 / 0.2e1 - (((t434 * t500 + t436 * t482 + t438 * t483) * t515 - (t433 * t500 + t435 * t482 + t437 * t483) * t518) * t516 + (-t460 * t557 - t462 * t500 + t464 * t501) * t558 - (-t459 * t557 - t461 * t500 + t463 * t501) * t557 + (t457 * t482 + t458 * t483 - t488 * t557 + t490 * t501 + t500 * t569) * t519) * t568 * t557 / 0.2e1 + ((((t462 * t525 + t464 * t523) * t515 - (t461 * t525 + t463 * t523) * t518) * t516 ^ 2 + (-t434 * t555 + t436 * t498 + t438 * t499) * t558 - (-t433 * t555 + t435 * t498 + t437 * t499) * t557 + ((-t459 * t518 + t460 * t515 + t489 * t525 + t490 * t523) * t516 - t456 * t555 + t498 * t457 + t499 * t458 + t519 * t488) * t519) * t519 + ((t460 * t558 - t462 * t502 + t464 * t503) * t558 - (t459 * t558 - t461 * t502 + t463 * t503) * t557 + ((t434 * t502 + t436 * t484 + t438 * t485) * t515 - (t433 * t502 + t435 * t484 + t437 * t485) * t518) * t516 + (t457 * t484 + t458 * t485 + t488 * t558 + t490 * t503 + t502 * t569) * t519) * t558) * t568 / 0.2e1;
T  = t1;
