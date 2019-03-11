% Calculate kinetic energy for
% S6RRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:52:32
% EndTime: 2019-03-09 11:52:34
% DurationCPUTime: 2.69s
% Computational Cost: add. (1990->280), mult. (2201->430), div. (0->0), fcn. (2168->10), ass. (0->154)
t573 = Icges(6,1) + Icges(7,1);
t572 = -Icges(6,4) + Icges(7,5);
t571 = Icges(7,4) + Icges(6,5);
t570 = Icges(6,2) + Icges(7,3);
t569 = -Icges(7,6) + Icges(6,6);
t568 = Icges(3,3) + Icges(4,3);
t567 = -Icges(6,3) - Icges(7,2);
t477 = qJ(2) + pkin(10);
t471 = sin(t477);
t472 = cos(t477);
t481 = sin(qJ(2));
t484 = cos(qJ(2));
t566 = Icges(3,5) * t484 + Icges(4,5) * t472 - Icges(3,6) * t481 - Icges(4,6) * t471;
t565 = rSges(7,1) + pkin(5);
t564 = rSges(7,3) + qJ(6);
t478 = qJ(4) + qJ(5);
t476 = cos(t478);
t485 = cos(qJ(1));
t528 = t476 * t485;
t475 = sin(t478);
t482 = sin(qJ(1));
t530 = t475 * t482;
t430 = t472 * t530 + t528;
t525 = t482 * t476;
t529 = t475 * t485;
t431 = t472 * t525 - t529;
t532 = t471 * t482;
t563 = t570 * t430 + t572 * t431 - t569 * t532;
t432 = t472 * t529 - t525;
t433 = t472 * t528 + t530;
t531 = t471 * t485;
t562 = t570 * t432 + t572 * t433 - t569 * t531;
t561 = -t569 * t430 + t571 * t431 - t567 * t532;
t560 = -t569 * t432 + t571 * t433 - t567 * t531;
t559 = t572 * t430 + t573 * t431 + t571 * t532;
t558 = t572 * t432 + t573 * t433 + t571 * t531;
t557 = t569 * t472 + (t570 * t475 + t572 * t476) * t471;
t556 = t567 * t472 + (-t569 * t475 + t571 * t476) * t471;
t555 = -t571 * t472 + (t572 * t475 + t573 * t476) * t471;
t554 = t566 * t482 - t568 * t485;
t553 = t568 * t482 + t566 * t485;
t552 = Icges(3,5) * t481 + Icges(4,5) * t471 + Icges(3,6) * t484 + Icges(4,6) * t472;
t534 = Icges(4,4) * t471;
t453 = Icges(4,2) * t472 + t534;
t533 = Icges(4,4) * t472;
t454 = Icges(4,1) * t471 + t533;
t536 = Icges(3,4) * t481;
t459 = Icges(3,2) * t484 + t536;
t535 = Icges(3,4) * t484;
t460 = Icges(3,1) * t481 + t535;
t551 = -t453 * t471 + t454 * t472 - t459 * t481 + t460 * t484;
t502 = -Icges(4,2) * t471 + t533;
t422 = Icges(4,6) * t482 + t485 * t502;
t504 = Icges(4,1) * t472 - t534;
t424 = Icges(4,5) * t482 + t485 * t504;
t503 = -Icges(3,2) * t481 + t535;
t438 = Icges(3,6) * t482 + t485 * t503;
t505 = Icges(3,1) * t484 - t536;
t440 = Icges(3,5) * t482 + t485 * t505;
t550 = -t422 * t471 + t424 * t472 - t438 * t481 + t440 * t484;
t421 = -Icges(4,6) * t485 + t482 * t502;
t423 = -Icges(4,5) * t485 + t482 * t504;
t437 = -Icges(3,6) * t485 + t482 * t503;
t439 = -Icges(3,5) * t485 + t482 * t505;
t549 = t421 * t471 - t423 * t472 + t437 * t481 - t439 * t484;
t542 = pkin(2) * t481;
t540 = pkin(2) * t484;
t483 = cos(qJ(4));
t539 = pkin(4) * t483;
t480 = sin(qJ(4));
t527 = t480 * t482;
t526 = t480 * t485;
t524 = t482 * t483;
t523 = t483 * t485;
t522 = rSges(7,2) * t532 + t564 * t430 + t565 * t431;
t521 = rSges(7,2) * t531 + t564 * t432 + t565 * t433;
t520 = -rSges(7,2) * t472 + (t564 * t475 + t565 * t476) * t471;
t417 = -qJ(3) * t485 + t482 * t540;
t418 = qJ(3) * t482 + t485 * t540;
t474 = qJD(2) * t482;
t517 = qJD(2) * t485;
t519 = t417 * t474 + t418 * t517;
t466 = pkin(1) * t482 - pkin(7) * t485;
t518 = -t417 - t466;
t516 = qJD(4) * t471;
t450 = t485 * t516 + t474;
t515 = qJD(5) * t471;
t451 = t482 * t516 - t517;
t511 = pkin(3) * t472 + pkin(8) * t471;
t443 = t511 * t482;
t444 = t511 * t485;
t512 = t443 * t474 + t444 * t517 + t519;
t457 = qJD(1) * (pkin(1) * t485 + pkin(7) * t482);
t510 = qJD(1) * t418 - qJD(3) * t485 + t457;
t509 = rSges(3,1) * t484 - rSges(3,2) * t481;
t508 = rSges(4,1) * t472 - rSges(4,2) * t471;
t507 = qJD(2) * (-rSges(4,1) * t471 - rSges(4,2) * t472 - t542);
t506 = qJD(2) * (-pkin(3) * t471 + pkin(8) * t472 - t542);
t492 = pkin(9) * t471 + t472 * t539;
t383 = -pkin(4) * t526 + t482 * t492;
t384 = pkin(4) * t527 + t485 * t492;
t493 = t450 * t383 - t384 * t451 + t512;
t491 = qJD(1) * t444 + t482 * t506 + t510;
t473 = qJD(3) * t482;
t490 = t473 + (-t443 + t518) * qJD(1) + t485 * t506;
t399 = -pkin(9) * t472 + t471 * t539;
t467 = -qJD(4) * t472 + qJD(1);
t489 = t467 * t384 - t399 * t450 + t491;
t488 = -t383 * t467 + t451 * t399 + t490;
t465 = rSges(2,1) * t485 - rSges(2,2) * t482;
t464 = rSges(2,1) * t482 + rSges(2,2) * t485;
t463 = rSges(3,1) * t481 + rSges(3,2) * t484;
t449 = qJD(1) + (-qJD(4) - qJD(5)) * t472;
t448 = t472 * t523 + t527;
t447 = -t472 * t526 + t524;
t446 = t472 * t524 - t526;
t445 = -t472 * t527 - t523;
t442 = rSges(3,3) * t482 + t485 * t509;
t441 = -rSges(3,3) * t485 + t482 * t509;
t426 = rSges(4,3) * t482 + t485 * t508;
t425 = -rSges(4,3) * t485 + t482 * t508;
t416 = t482 * t515 + t451;
t415 = t485 * t515 + t450;
t414 = -rSges(5,3) * t472 + (rSges(5,1) * t483 - rSges(5,2) * t480) * t471;
t412 = -Icges(5,5) * t472 + (Icges(5,1) * t483 - Icges(5,4) * t480) * t471;
t411 = -Icges(5,6) * t472 + (Icges(5,4) * t483 - Icges(5,2) * t480) * t471;
t410 = -Icges(5,3) * t472 + (Icges(5,5) * t483 - Icges(5,6) * t480) * t471;
t407 = -rSges(6,3) * t472 + (rSges(6,1) * t476 - rSges(6,2) * t475) * t471;
t396 = qJD(1) * t442 - t463 * t474 + t457;
t395 = -t463 * t517 + (-t441 - t466) * qJD(1);
t394 = (t441 * t482 + t442 * t485) * qJD(2);
t393 = rSges(5,1) * t448 + rSges(5,2) * t447 + rSges(5,3) * t531;
t392 = rSges(5,1) * t446 + rSges(5,2) * t445 + rSges(5,3) * t532;
t390 = Icges(5,1) * t448 + Icges(5,4) * t447 + Icges(5,5) * t531;
t389 = Icges(5,1) * t446 + Icges(5,4) * t445 + Icges(5,5) * t532;
t388 = Icges(5,4) * t448 + Icges(5,2) * t447 + Icges(5,6) * t531;
t387 = Icges(5,4) * t446 + Icges(5,2) * t445 + Icges(5,6) * t532;
t386 = Icges(5,5) * t448 + Icges(5,6) * t447 + Icges(5,3) * t531;
t385 = Icges(5,5) * t446 + Icges(5,6) * t445 + Icges(5,3) * t532;
t382 = rSges(6,1) * t433 - rSges(6,2) * t432 + rSges(6,3) * t531;
t380 = rSges(6,1) * t431 - rSges(6,2) * t430 + rSges(6,3) * t532;
t364 = qJD(1) * t426 + t482 * t507 + t510;
t363 = t473 + t485 * t507 + (-t425 + t518) * qJD(1);
t362 = (t425 * t482 + t426 * t485) * qJD(2) + t519;
t361 = t393 * t467 - t414 * t450 + t491;
t360 = -t392 * t467 + t414 * t451 + t490;
t359 = t392 * t450 - t393 * t451 + t512;
t358 = t382 * t449 - t407 * t415 + t489;
t357 = -t380 * t449 + t407 * t416 + t488;
t356 = t380 * t415 - t382 * t416 + t493;
t355 = qJD(6) * t430 - t415 * t520 + t449 * t521 + t489;
t354 = qJD(6) * t432 + t416 * t520 - t449 * t522 + t488;
t353 = qJD(6) * t471 * t475 + t415 * t522 - t416 * t521 + t493;
t1 = m(3) * (t394 ^ 2 + t395 ^ 2 + t396 ^ 2) / 0.2e1 + m(4) * (t362 ^ 2 + t363 ^ 2 + t364 ^ 2) / 0.2e1 + m(5) * (t359 ^ 2 + t360 ^ 2 + t361 ^ 2) / 0.2e1 + m(6) * (t356 ^ 2 + t357 ^ 2 + t358 ^ 2) / 0.2e1 + m(7) * (t353 ^ 2 + t354 ^ 2 + t355 ^ 2) / 0.2e1 + t450 * ((t386 * t531 + t447 * t388 + t448 * t390) * t450 + (t385 * t531 + t387 * t447 + t389 * t448) * t451 + (t410 * t531 + t411 * t447 + t412 * t448) * t467) / 0.2e1 + t451 * ((t386 * t532 + t388 * t445 + t390 * t446) * t450 + (t385 * t532 + t445 * t387 + t446 * t389) * t451 + (t410 * t532 + t411 * t445 + t412 * t446) * t467) / 0.2e1 + t467 * ((-t385 * t451 - t386 * t450 - t410 * t467) * t472 + ((-t388 * t480 + t390 * t483) * t450 + (-t387 * t480 + t389 * t483) * t451 + (-t411 * t480 + t412 * t483) * t467) * t471) / 0.2e1 + ((t557 * t432 + t555 * t433 + t556 * t531) * t449 + (t563 * t432 + t559 * t433 + t561 * t531) * t416 + (t562 * t432 + t558 * t433 + t560 * t531) * t415) * t415 / 0.2e1 + ((t557 * t430 + t555 * t431 + t556 * t532) * t449 + (t563 * t430 + t559 * t431 + t561 * t532) * t416 + (t562 * t430 + t558 * t431 + t560 * t532) * t415) * t416 / 0.2e1 + ((-t560 * t415 - t561 * t416 - t556 * t449) * t472 + ((t557 * t475 + t555 * t476) * t449 + (t563 * t475 + t559 * t476) * t416 + (t562 * t475 + t558 * t476) * t415) * t471) * t449 / 0.2e1 + (Icges(2,3) + m(2) * (t464 ^ 2 + t465 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((-t421 * t472 - t423 * t471 - t437 * t484 - t439 * t481) * t485 + (t422 * t472 + t424 * t471 + t438 * t484 + t440 * t481) * t482) * qJD(2) + (t472 * t453 + t471 * t454 + t484 * t459 + t481 * t460) * qJD(1)) * qJD(1) / 0.2e1 + ((t553 * t482 ^ 2 + (t549 * t485 + (t550 - t554) * t482) * t485) * qJD(2) + (t552 * t482 + t551 * t485) * qJD(1)) * t474 / 0.2e1 - ((t554 * t485 ^ 2 + (t550 * t482 + (t549 - t553) * t485) * t482) * qJD(2) + (t551 * t482 - t552 * t485) * qJD(1)) * t517 / 0.2e1;
T  = t1;
