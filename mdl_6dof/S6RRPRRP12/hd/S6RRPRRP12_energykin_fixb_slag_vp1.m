% Calculate kinetic energy for
% S6RRPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP12_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP12_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP12_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP12_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP12_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:50:21
% EndTime: 2019-03-09 12:50:24
% DurationCPUTime: 3.59s
% Computational Cost: add. (1255->263), mult. (2175->404), div. (0->0), fcn. (2142->8), ass. (0->145)
t581 = Icges(3,4) + Icges(4,6);
t580 = Icges(3,1) + Icges(4,2);
t579 = -Icges(3,2) - Icges(4,3);
t480 = cos(qJ(2));
t578 = t581 * t480;
t477 = sin(qJ(2));
t577 = t581 * t477;
t576 = -Icges(4,4) + Icges(3,5);
t575 = Icges(4,5) - Icges(3,6);
t574 = t477 * t579 + t578;
t573 = -t480 * t580 + t577;
t572 = Icges(4,1) + Icges(3,3);
t571 = Icges(6,1) + Icges(7,1);
t570 = Icges(6,4) - Icges(7,5);
t569 = Icges(7,4) + Icges(6,5);
t568 = Icges(6,2) + Icges(7,3);
t567 = Icges(7,6) - Icges(6,6);
t566 = Icges(6,3) + Icges(7,2);
t478 = sin(qJ(1));
t481 = cos(qJ(1));
t565 = t478 * t574 + t481 * t575;
t564 = -t478 * t575 + t481 * t574;
t563 = t573 * t478 + t481 * t576;
t562 = t478 * t576 - t573 * t481;
t561 = t480 * t579 - t577;
t560 = t477 * t580 + t578;
t559 = t575 * t477 + t480 * t576;
t558 = rSges(7,1) + pkin(5);
t557 = rSges(7,3) + qJ(6);
t475 = qJ(4) + qJ(5);
t473 = sin(t475);
t474 = cos(t475);
t525 = t477 * t481;
t430 = t473 * t478 - t474 * t525;
t431 = t473 * t525 + t474 * t478;
t521 = t480 * t481;
t556 = t568 * t430 - t570 * t431 + t567 * t521;
t526 = t477 * t478;
t432 = t473 * t481 + t474 * t526;
t433 = t473 * t526 - t474 * t481;
t523 = t478 * t480;
t555 = -t568 * t432 - t570 * t433 + t567 * t523;
t554 = t567 * t430 + t569 * t431 + t566 * t521;
t553 = -t567 * t432 + t569 * t433 + t566 * t523;
t552 = -t570 * t430 + t571 * t431 + t569 * t521;
t551 = t570 * t432 + t571 * t433 + t569 * t523;
t550 = (t570 * t473 + t568 * t474) * t480 + t567 * t477;
t549 = (-t569 * t473 + t567 * t474) * t480 + t566 * t477;
t548 = (-t571 * t473 - t570 * t474) * t480 + t569 * t477;
t547 = t572 * t478 + t559 * t481;
t546 = t559 * t478 - t572 * t481;
t448 = pkin(3) * t478 + pkin(8) * t521;
t449 = -pkin(3) * t481 + pkin(8) * t523;
t472 = qJD(2) * t478;
t515 = qJD(2) * t481;
t545 = t448 * t515 + t449 * t472;
t544 = t477 * t576 - t575 * t480;
t543 = t561 * t477 + t560 * t480;
t542 = t565 * t477 + t563 * t480;
t541 = -t564 * t477 + t562 * t480;
t476 = sin(qJ(4));
t534 = pkin(4) * t476;
t479 = cos(qJ(4));
t532 = pkin(4) * t479;
t524 = t478 * t479;
t522 = t479 * t481;
t520 = rSges(7,2) * t521 + t557 * t430 + t558 * t431;
t519 = rSges(7,2) * t523 - t557 * t432 + t558 * t433;
t518 = rSges(7,2) * t477 + (-t558 * t473 + t557 * t474) * t480;
t503 = pkin(2) * t480 + qJ(3) * t477;
t443 = t503 * t478;
t444 = t503 * t481;
t517 = t443 * t472 + t444 * t515;
t464 = pkin(1) * t478 - pkin(7) * t481;
t516 = -t443 - t464;
t512 = qJD(4) * t480;
t446 = t481 * t512 + t472;
t514 = qJD(3) * t477;
t513 = qJD(3) * t480;
t511 = qJD(5) * t480;
t469 = qJD(4) * t477 + qJD(1);
t451 = qJD(1) * (pkin(1) * t481 + pkin(7) * t478);
t510 = qJD(1) * t444 + t478 * t514 + t451;
t459 = pkin(2) * t477 - qJ(3) * t480;
t507 = qJD(2) * (rSges(4,2) * t477 + rSges(4,3) * t480 - t459);
t447 = t478 * t512 - t515;
t506 = -t513 + t517;
t505 = rSges(3,1) * t480 - rSges(3,2) * t477;
t504 = -rSges(4,2) * t480 + rSges(4,3) * t477;
t502 = qJD(2) * (-pkin(8) * t477 - t459);
t489 = pkin(9) * t480 + t477 * t534;
t395 = t478 * t532 + t481 * t489;
t396 = t478 * t489 - t481 * t532;
t488 = -t395 * t447 + t446 * t396 + t517 + t545;
t487 = qJD(1) * t448 + t478 * t502 + t510;
t468 = t481 * t514;
t486 = t468 + (-t449 + t516) * qJD(1) + t481 * t502;
t435 = pkin(9) * t477 - t480 * t534;
t485 = t469 * t395 - t435 * t446 + t487;
t484 = -t396 * t469 + t447 * t435 + t486;
t463 = rSges(2,1) * t481 - rSges(2,2) * t478;
t462 = rSges(2,1) * t478 + rSges(2,2) * t481;
t461 = rSges(3,1) * t477 + rSges(3,2) * t480;
t450 = qJD(5) * t477 + t469;
t442 = t476 * t526 - t522;
t441 = t476 * t481 + t477 * t524;
t440 = t476 * t525 + t524;
t439 = -t476 * t478 + t477 * t522;
t429 = -rSges(4,1) * t481 + t478 * t504;
t428 = rSges(4,1) * t478 + t481 * t504;
t427 = rSges(3,3) * t478 + t481 * t505;
t426 = rSges(5,3) * t477 + (-rSges(5,1) * t476 - rSges(5,2) * t479) * t480;
t425 = -rSges(3,3) * t481 + t478 * t505;
t414 = Icges(5,5) * t477 + (-Icges(5,1) * t476 - Icges(5,4) * t479) * t480;
t411 = Icges(5,6) * t477 + (-Icges(5,4) * t476 - Icges(5,2) * t479) * t480;
t408 = Icges(5,3) * t477 + (-Icges(5,5) * t476 - Icges(5,6) * t479) * t480;
t407 = t478 * t511 + t447;
t406 = t481 * t511 + t446;
t405 = rSges(6,3) * t477 + (-rSges(6,1) * t473 - rSges(6,2) * t474) * t480;
t394 = rSges(5,1) * t442 + rSges(5,2) * t441 + rSges(5,3) * t523;
t393 = rSges(5,1) * t440 + rSges(5,2) * t439 + rSges(5,3) * t521;
t392 = Icges(5,1) * t442 + Icges(5,4) * t441 + Icges(5,5) * t523;
t391 = Icges(5,1) * t440 + Icges(5,4) * t439 + Icges(5,5) * t521;
t390 = Icges(5,4) * t442 + Icges(5,2) * t441 + Icges(5,6) * t523;
t389 = Icges(5,4) * t440 + Icges(5,2) * t439 + Icges(5,6) * t521;
t388 = Icges(5,5) * t442 + Icges(5,6) * t441 + Icges(5,3) * t523;
t387 = Icges(5,5) * t440 + Icges(5,6) * t439 + Icges(5,3) * t521;
t383 = qJD(1) * t427 - t461 * t472 + t451;
t382 = -t461 * t515 + (-t425 - t464) * qJD(1);
t381 = (t425 * t478 + t427 * t481) * qJD(2);
t380 = rSges(6,1) * t433 + rSges(6,2) * t432 + rSges(6,3) * t523;
t378 = rSges(6,1) * t431 - rSges(6,2) * t430 + rSges(6,3) * t521;
t363 = qJD(1) * t428 + t478 * t507 + t510;
t362 = t468 + t481 * t507 + (-t429 + t516) * qJD(1);
t361 = (t428 * t481 + t429 * t478) * qJD(2) + t506;
t360 = t393 * t469 - t426 * t446 + t487;
t359 = -t394 * t469 + t426 * t447 + t486;
t358 = -t393 * t447 + t394 * t446 + t506 + t545;
t357 = t378 * t450 - t405 * t406 + t485;
t356 = -t380 * t450 + t405 * t407 + t484;
t355 = -t378 * t407 + t380 * t406 + t488 - t513;
t354 = -qJD(6) * t432 - t406 * t518 + t450 * t520 + t485;
t353 = qJD(6) * t430 + t407 * t518 - t450 * t519 + t484;
t352 = (qJD(6) * t474 - qJD(3)) * t480 - t520 * t407 + t519 * t406 + t488;
t1 = t446 * ((t387 * t521 + t439 * t389 + t440 * t391) * t446 + (t388 * t521 + t439 * t390 + t440 * t392) * t447 + (t408 * t521 + t439 * t411 + t440 * t414) * t469) / 0.2e1 + t447 * ((t387 * t523 + t441 * t389 + t442 * t391) * t446 + (t388 * t523 + t441 * t390 + t442 * t392) * t447 + (t408 * t523 + t441 * t411 + t442 * t414) * t469) / 0.2e1 + t469 * ((t387 * t446 + t388 * t447 + t408 * t469) * t477 + ((-t389 * t479 - t391 * t476) * t446 + (-t390 * t479 - t392 * t476) * t447 + (-t411 * t479 - t414 * t476) * t469) * t480) / 0.2e1 + m(7) * (t352 ^ 2 + t353 ^ 2 + t354 ^ 2) / 0.2e1 + m(6) * (t355 ^ 2 + t356 ^ 2 + t357 ^ 2) / 0.2e1 + m(5) * (t358 ^ 2 + t359 ^ 2 + t360 ^ 2) / 0.2e1 + m(4) * (t361 ^ 2 + t362 ^ 2 + t363 ^ 2) / 0.2e1 + m(3) * (t381 ^ 2 + t382 ^ 2 + t383 ^ 2) / 0.2e1 + ((t430 * t550 + t431 * t548 + t521 * t549) * t450 + (t430 * t555 + t431 * t551 + t521 * t553) * t407 + (t556 * t430 + t552 * t431 + t554 * t521) * t406) * t406 / 0.2e1 + ((-t432 * t550 + t433 * t548 + t523 * t549) * t450 + (-t555 * t432 + t551 * t433 + t553 * t523) * t407 + (-t432 * t556 + t552 * t433 + t554 * t523) * t406) * t407 / 0.2e1 + (((-t473 * t548 + t474 * t550) * t450 + (-t473 * t551 + t474 * t555) * t407 + (-t552 * t473 + t474 * t556) * t406) * t480 + (t406 * t554 + t407 * t553 + t450 * t549) * t477) * t450 / 0.2e1 + (m(2) * (t462 ^ 2 + t463 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t563 * t477 - t565 * t480) * t481 + (t562 * t477 + t564 * t480) * t478) * qJD(2) + (t560 * t477 - t561 * t480) * qJD(1)) * qJD(1) / 0.2e1 + ((t547 * t478 ^ 2 + (t542 * t481 + (t541 - t546) * t478) * t481) * qJD(2) + (t478 * t544 + t481 * t543) * qJD(1)) * t472 / 0.2e1 - ((t546 * t481 ^ 2 + (t541 * t478 + (t542 - t547) * t481) * t478) * qJD(2) + (t478 * t543 - t481 * t544) * qJD(1)) * t515 / 0.2e1;
T  = t1;
