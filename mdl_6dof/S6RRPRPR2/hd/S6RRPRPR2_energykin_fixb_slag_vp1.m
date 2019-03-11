% Calculate kinetic energy for
% S6RRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:12:02
% EndTime: 2019-03-09 10:12:05
% DurationCPUTime: 3.00s
% Computational Cost: add. (1611->262), mult. (1637->396), div. (0->0), fcn. (1482->10), ass. (0->153)
t563 = Icges(5,4) + Icges(6,6);
t562 = Icges(5,1) + Icges(6,2);
t561 = -Icges(5,2) - Icges(6,3);
t457 = qJ(2) + pkin(10);
t451 = qJ(4) + t457;
t447 = cos(t451);
t560 = t563 * t447;
t446 = sin(t451);
t559 = t563 * t446;
t558 = Icges(6,4) - Icges(5,5);
t557 = Icges(6,5) - Icges(5,6);
t556 = t561 * t446 + t560;
t555 = t562 * t447 - t559;
t554 = Icges(6,1) + Icges(5,3);
t553 = Icges(3,3) + Icges(4,3);
t461 = sin(qJ(1));
t464 = cos(qJ(1));
t552 = t556 * t461 + t557 * t464;
t551 = -t557 * t461 + t556 * t464;
t550 = t555 * t461 + t558 * t464;
t549 = -t558 * t461 + t555 * t464;
t548 = t561 * t447 - t559;
t547 = t562 * t446 + t560;
t546 = t557 * t446 - t558 * t447;
t449 = sin(t457);
t450 = cos(t457);
t460 = sin(qJ(2));
t463 = cos(qJ(2));
t545 = Icges(3,5) * t463 + Icges(4,5) * t450 - Icges(3,6) * t460 - Icges(4,6) * t449;
t544 = t545 * t461 - t553 * t464;
t543 = t553 * t461 + t545 * t464;
t542 = Icges(3,5) * t460 + Icges(4,5) * t449 + Icges(3,6) * t463 + Icges(4,6) * t450;
t454 = qJD(2) * t461;
t437 = qJD(4) * t461 + t454;
t438 = (-qJD(2) - qJD(4)) * t464;
t541 = (-t552 * t446 + t550 * t447) * t438 + (-t551 * t446 + t549 * t447) * t437 + (t548 * t446 + t547 * t447) * qJD(1);
t540 = (t546 * t461 - t554 * t464) * t438 + (t554 * t461 + t546 * t464) * t437 + (-t558 * t446 - t557 * t447) * qJD(1);
t524 = Icges(4,4) * t449;
t429 = Icges(4,2) * t450 + t524;
t523 = Icges(4,4) * t450;
t430 = Icges(4,1) * t449 + t523;
t526 = Icges(3,4) * t460;
t440 = Icges(3,2) * t463 + t526;
t525 = Icges(3,4) * t463;
t441 = Icges(3,1) * t460 + t525;
t539 = -t429 * t449 + t430 * t450 - t440 * t460 + t441 * t463;
t489 = -Icges(4,2) * t449 + t523;
t395 = Icges(4,6) * t461 + t464 * t489;
t492 = Icges(4,1) * t450 - t524;
t397 = Icges(4,5) * t461 + t464 * t492;
t490 = -Icges(3,2) * t460 + t525;
t408 = Icges(3,6) * t461 + t464 * t490;
t493 = Icges(3,1) * t463 - t526;
t410 = Icges(3,5) * t461 + t464 * t493;
t538 = -t395 * t449 + t397 * t450 - t408 * t460 + t410 * t463;
t394 = -Icges(4,6) * t464 + t461 * t489;
t396 = -Icges(4,5) * t464 + t461 * t492;
t407 = -Icges(3,6) * t464 + t461 * t490;
t409 = -Icges(3,5) * t464 + t461 * t493;
t537 = t394 * t449 - t396 * t450 + t407 * t460 - t409 * t463;
t531 = pkin(2) * t460;
t530 = pkin(9) * t446;
t528 = t463 * pkin(2);
t518 = t447 * t461;
t517 = t447 * t464;
t459 = sin(qJ(6));
t516 = t459 * t461;
t515 = t459 * t464;
t462 = cos(qJ(6));
t514 = t461 * t462;
t513 = t462 * t464;
t390 = -qJ(3) * t464 + t461 * t528;
t391 = qJ(3) * t461 + t464 * t528;
t508 = qJD(2) * t464;
t512 = t390 * t454 + t391 * t508;
t445 = pkin(1) * t461 - pkin(7) * t464;
t511 = -t390 - t445;
t510 = pkin(3) * t450;
t507 = qJD(5) * t446;
t506 = qJD(6) * t447;
t363 = -pkin(8) * t464 + t461 * t510;
t505 = -t363 + t511;
t364 = pkin(8) * t461 + t464 * t510;
t502 = t363 * t454 + t364 * t508 + t512;
t495 = pkin(4) * t447 + qJ(5) * t446;
t401 = t495 * t461;
t501 = -t401 + t505;
t433 = qJD(1) * (pkin(1) * t464 + pkin(7) * t461);
t500 = qJD(1) * t391 - qJD(3) * t464 + t433;
t499 = rSges(3,1) * t463 - rSges(3,2) * t460;
t498 = rSges(4,1) * t450 - rSges(4,2) * t449;
t497 = rSges(5,1) * t447 - rSges(5,2) * t446;
t496 = -rSges(6,2) * t447 + rSges(6,3) * t446;
t494 = qJD(2) * (-rSges(4,1) * t449 - rSges(4,2) * t450 - t531);
t475 = qJD(2) * (-pkin(3) * t449 - t531);
t474 = -qJD(5) * t447 + t437 * t401 + t502;
t453 = qJD(3) * t461;
t471 = t464 * t475 + t453;
t423 = pkin(4) * t446 - qJ(5) * t447;
t470 = t438 * t423 + t464 * t507 + t471;
t469 = qJD(1) * t364 + t461 * t475 + t500;
t402 = t495 * t464;
t468 = qJD(1) * t402 + t461 * t507 + t469;
t444 = rSges(2,1) * t464 - rSges(2,2) * t461;
t443 = rSges(2,1) * t461 + rSges(2,2) * t464;
t442 = rSges(3,1) * t460 + rSges(3,2) * t463;
t436 = qJD(6) * t446 + qJD(1);
t427 = -pkin(5) * t464 + pkin(9) * t518;
t426 = pkin(5) * t461 + pkin(9) * t517;
t425 = rSges(5,1) * t446 + rSges(5,2) * t447;
t424 = -rSges(6,2) * t446 - rSges(6,3) * t447;
t416 = t446 * t516 - t513;
t415 = t446 * t514 + t515;
t414 = t446 * t515 + t514;
t413 = t446 * t513 - t516;
t412 = rSges(3,3) * t461 + t464 * t499;
t411 = -rSges(3,3) * t464 + t461 * t499;
t404 = t461 * t506 + t438;
t403 = t464 * t506 + t437;
t400 = rSges(4,3) * t461 + t464 * t498;
t399 = -rSges(4,3) * t464 + t461 * t498;
t389 = -rSges(6,1) * t464 + t461 * t496;
t388 = rSges(6,1) * t461 + t464 * t496;
t387 = rSges(5,3) * t461 + t464 * t497;
t386 = -rSges(5,3) * t464 + t461 * t497;
t369 = rSges(7,3) * t446 + (-rSges(7,1) * t459 - rSges(7,2) * t462) * t447;
t368 = Icges(7,5) * t446 + (-Icges(7,1) * t459 - Icges(7,4) * t462) * t447;
t367 = Icges(7,6) * t446 + (-Icges(7,4) * t459 - Icges(7,2) * t462) * t447;
t366 = Icges(7,3) * t446 + (-Icges(7,5) * t459 - Icges(7,6) * t462) * t447;
t359 = qJD(1) * t412 - t442 * t454 + t433;
t358 = -t442 * t508 + (-t411 - t445) * qJD(1);
t357 = (t411 * t461 + t412 * t464) * qJD(2);
t356 = rSges(7,1) * t416 + rSges(7,2) * t415 + rSges(7,3) * t518;
t355 = rSges(7,1) * t414 + rSges(7,2) * t413 + rSges(7,3) * t517;
t354 = Icges(7,1) * t416 + Icges(7,4) * t415 + Icges(7,5) * t518;
t353 = Icges(7,1) * t414 + Icges(7,4) * t413 + Icges(7,5) * t517;
t352 = Icges(7,4) * t416 + Icges(7,2) * t415 + Icges(7,6) * t518;
t351 = Icges(7,4) * t414 + Icges(7,2) * t413 + Icges(7,6) * t517;
t350 = Icges(7,5) * t416 + Icges(7,6) * t415 + Icges(7,3) * t518;
t349 = Icges(7,5) * t414 + Icges(7,6) * t413 + Icges(7,3) * t517;
t348 = qJD(1) * t400 + t461 * t494 + t500;
t347 = t453 + t464 * t494 + (-t399 + t511) * qJD(1);
t346 = (t399 * t461 + t400 * t464) * qJD(2) + t512;
t345 = qJD(1) * t387 - t425 * t437 + t469;
t344 = t425 * t438 + (-t386 + t505) * qJD(1) + t471;
t343 = t386 * t437 - t387 * t438 + t502;
t342 = qJD(1) * t388 + (-t423 - t424) * t437 + t468;
t341 = t424 * t438 + (-t389 + t501) * qJD(1) + t470;
t340 = t389 * t437 + (-t388 - t402) * t438 + t474;
t339 = qJD(1) * t426 + t355 * t436 - t369 * t403 + (-t423 - t530) * t437 + t468;
t338 = t438 * t530 - t356 * t436 + t369 * t404 + (-t427 + t501) * qJD(1) + t470;
t337 = -t355 * t404 + t356 * t403 + t427 * t437 + (-t402 - t426) * t438 + t474;
t1 = m(3) * (t357 ^ 2 + t358 ^ 2 + t359 ^ 2) / 0.2e1 + m(4) * (t346 ^ 2 + t347 ^ 2 + t348 ^ 2) / 0.2e1 + m(7) * (t337 ^ 2 + t338 ^ 2 + t339 ^ 2) / 0.2e1 + m(5) * (t343 ^ 2 + t344 ^ 2 + t345 ^ 2) / 0.2e1 + m(6) * (t340 ^ 2 + t341 ^ 2 + t342 ^ 2) / 0.2e1 + t404 * ((t349 * t518 + t351 * t415 + t353 * t416) * t403 + (t350 * t518 + t415 * t352 + t416 * t354) * t404 + (t366 * t518 + t367 * t415 + t368 * t416) * t436) / 0.2e1 + t436 * ((t349 * t403 + t350 * t404 + t366 * t436) * t446 + ((-t351 * t462 - t353 * t459) * t403 + (-t352 * t462 - t354 * t459) * t404 + (-t367 * t462 - t368 * t459) * t436) * t447) / 0.2e1 + t403 * ((t349 * t517 + t413 * t351 + t414 * t353) * t403 + (t350 * t517 + t352 * t413 + t354 * t414) * t404 + (t366 * t517 + t367 * t413 + t368 * t414) * t436) / 0.2e1 + (t540 * t461 + t541 * t464) * t437 / 0.2e1 + (t541 * t461 - t540 * t464) * t438 / 0.2e1 + (m(2) * (t443 ^ 2 + t444 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t543 * t461 ^ 2 + (t537 * t464 + (t538 - t544) * t461) * t464) * qJD(2) + (t461 * t542 + t464 * t539) * qJD(1)) * t454 / 0.2e1 - ((t544 * t464 ^ 2 + (t538 * t461 + (t537 - t543) * t464) * t461) * qJD(2) + (t461 * t539 - t464 * t542) * qJD(1)) * t508 / 0.2e1 + ((t550 * t446 + t552 * t447) * t438 + (t549 * t446 + t551 * t447) * t437 + ((-t394 * t450 - t396 * t449 - t407 * t463 - t409 * t460) * t464 + (t395 * t450 + t397 * t449 + t408 * t463 + t410 * t460) * t461) * qJD(2) + (t450 * t429 + t449 * t430 + t463 * t440 + t460 * t441 + t547 * t446 - t548 * t447) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
