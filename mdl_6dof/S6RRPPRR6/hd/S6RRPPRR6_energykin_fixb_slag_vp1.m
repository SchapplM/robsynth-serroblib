% Calculate kinetic energy for
% S6RRPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:13:02
% EndTime: 2019-03-09 09:13:05
% DurationCPUTime: 3.62s
% Computational Cost: add. (1505->288), mult. (2467->437), div. (0->0), fcn. (2624->10), ass. (0->150)
t563 = Icges(3,4) - Icges(4,5);
t562 = Icges(3,1) + Icges(4,1);
t561 = Icges(3,2) + Icges(4,3);
t479 = sin(qJ(2));
t560 = t563 * t479;
t482 = cos(qJ(2));
t559 = t563 * t482;
t558 = Icges(4,4) + Icges(3,5);
t557 = Icges(3,6) - Icges(4,6);
t556 = t479 * t561 - t559;
t555 = t482 * t562 - t560;
t554 = Icges(4,2) + Icges(3,3) + Icges(5,3);
t480 = sin(qJ(1));
t483 = cos(qJ(1));
t552 = t556 * t480 + t557 * t483;
t551 = -t557 * t480 + t556 * t483;
t550 = -t555 * t480 + t558 * t483;
t549 = t558 * t480 + t555 * t483;
t548 = -t482 * t561 - t560;
t547 = t479 * t562 + t559;
t546 = -t557 * t479 + t558 * t482;
t476 = cos(pkin(10));
t475 = sin(pkin(10));
t526 = t482 * t475;
t448 = t479 * t476 - t526;
t435 = t448 * t480;
t528 = t479 * t475;
t492 = t482 * t476 + t528;
t436 = t492 * t480;
t545 = -Icges(5,5) * t436 - Icges(5,6) * t435 + t546 * t480 - t483 * t554;
t437 = t448 * t483;
t438 = t492 * t483;
t544 = -Icges(5,5) * t438 - Icges(5,6) * t437 + t480 * t554 + t546 * t483;
t541 = t548 * t479 + t547 * t482;
t540 = t551 * t479 + t549 * t482;
t539 = -t552 * t479 + t550 * t482;
t520 = pkin(10) + qJ(5);
t471 = sin(t520);
t514 = cos(t520);
t506 = t479 * t514;
t440 = -t482 * t471 + t506;
t538 = Icges(5,5) * t448 - Icges(5,6) * t492 - t558 * t479 - t557 * t482;
t533 = t476 * pkin(4);
t525 = t482 * t483;
t507 = pkin(2) * t482 + qJ(3) * t479;
t444 = t507 * t480;
t467 = t480 * pkin(1) - t483 * pkin(7);
t523 = -t444 - t467;
t474 = qJD(2) * t480;
t522 = qJD(2) * t483;
t521 = qJD(3) * t479;
t445 = t507 * t483;
t451 = qJD(1) * (t483 * pkin(1) + t480 * pkin(7));
t519 = qJD(1) * t445 + t480 * t521 + t451;
t449 = t480 * t482 * pkin(3) + t483 * qJ(4);
t518 = -t449 + t523;
t462 = t479 * pkin(2) - t482 * qJ(3);
t515 = -t479 * pkin(3) - t462;
t513 = qJD(2) * (-t479 * rSges(4,1) + t482 * rSges(4,3) - t462);
t455 = qJD(5) * t483 - t522;
t469 = t483 * t521;
t512 = -qJD(4) * t480 + t469;
t454 = -qJD(5) * t480 + t474;
t488 = pkin(4) * t528 + t482 * t533;
t394 = pkin(8) * t483 + t480 * t488;
t511 = -t394 + t518;
t510 = -qJD(3) * t482 + t444 * t474 + t445 * t522;
t509 = rSges(3,1) * t482 - rSges(3,2) * t479;
t508 = rSges(4,1) * t482 + rSges(4,3) * t479;
t450 = pkin(3) * t525 - t480 * qJ(4);
t505 = qJD(1) * t450 + qJD(4) * t483 + t519;
t491 = qJD(2) * (-t448 * rSges(5,1) + rSges(5,2) * t492 + t515);
t490 = qJD(2) * (pkin(4) * t526 - t479 * t533 + t515);
t489 = t449 * t474 + t450 * t522 + t510;
t439 = t479 * t471 + t482 * t514;
t395 = -pkin(8) * t480 + t483 * t488;
t487 = t394 * t474 + t395 * t522 + t489;
t486 = t483 * t490 + t512;
t485 = qJD(1) * t395 + t480 * t490 + t505;
t481 = cos(qJ(6));
t478 = sin(qJ(6));
t466 = t483 * rSges(2,1) - t480 * rSges(2,2);
t465 = t480 * rSges(2,1) + t483 * rSges(2,2);
t464 = t479 * rSges(3,1) + t482 * rSges(3,2);
t434 = t480 * rSges(3,3) + t483 * t509;
t433 = t480 * rSges(4,2) + t483 * t508;
t432 = -t483 * rSges(3,3) + t480 * t509;
t431 = -t483 * rSges(4,2) + t480 * t508;
t416 = qJD(6) * t439 + qJD(1);
t414 = t439 * t483;
t413 = t471 * t525 - t483 * t506;
t412 = t439 * t480;
t411 = t440 * t480;
t409 = Icges(5,1) * t448 - Icges(5,4) * t492;
t408 = Icges(5,4) * t448 - Icges(5,2) * t492;
t406 = t414 * t481 - t480 * t478;
t405 = -t414 * t478 - t480 * t481;
t404 = t412 * t481 + t483 * t478;
t403 = -t412 * t478 + t483 * t481;
t402 = -qJD(6) * t411 + t455;
t401 = qJD(6) * t413 + t454;
t400 = t440 * pkin(5) + t439 * pkin(9);
t399 = t440 * rSges(6,1) - t439 * rSges(6,2);
t398 = Icges(6,1) * t440 - Icges(6,4) * t439;
t397 = Icges(6,4) * t440 - Icges(6,2) * t439;
t396 = Icges(6,5) * t440 - Icges(6,6) * t439;
t390 = t438 * rSges(5,1) + t437 * rSges(5,2) - t480 * rSges(5,3);
t389 = t436 * rSges(5,1) + t435 * rSges(5,2) + t483 * rSges(5,3);
t388 = Icges(5,1) * t438 + Icges(5,4) * t437 - Icges(5,5) * t480;
t387 = Icges(5,1) * t436 + Icges(5,4) * t435 + Icges(5,5) * t483;
t386 = Icges(5,4) * t438 + Icges(5,2) * t437 - Icges(5,6) * t480;
t385 = Icges(5,4) * t436 + Icges(5,2) * t435 + Icges(5,6) * t483;
t382 = qJD(1) * t434 - t464 * t474 + t451;
t381 = -t464 * t522 + (-t432 - t467) * qJD(1);
t380 = (t432 * t480 + t434 * t483) * qJD(2);
t379 = pkin(5) * t414 + pkin(9) * t413;
t378 = pkin(5) * t412 - pkin(9) * t411;
t377 = t414 * rSges(6,1) - t413 * rSges(6,2) - t480 * rSges(6,3);
t376 = t412 * rSges(6,1) + t411 * rSges(6,2) + t483 * rSges(6,3);
t375 = Icges(6,1) * t414 - Icges(6,4) * t413 - Icges(6,5) * t480;
t374 = Icges(6,1) * t412 + Icges(6,4) * t411 + Icges(6,5) * t483;
t373 = Icges(6,4) * t414 - Icges(6,2) * t413 - Icges(6,6) * t480;
t372 = Icges(6,4) * t412 + Icges(6,2) * t411 + Icges(6,6) * t483;
t371 = Icges(6,5) * t414 - Icges(6,6) * t413 - Icges(6,3) * t480;
t370 = Icges(6,5) * t412 + Icges(6,6) * t411 + Icges(6,3) * t483;
t369 = t439 * rSges(7,3) + (rSges(7,1) * t481 - rSges(7,2) * t478) * t440;
t368 = Icges(7,5) * t439 + (Icges(7,1) * t481 - Icges(7,4) * t478) * t440;
t367 = Icges(7,6) * t439 + (Icges(7,4) * t481 - Icges(7,2) * t478) * t440;
t366 = Icges(7,3) * t439 + (Icges(7,5) * t481 - Icges(7,6) * t478) * t440;
t365 = qJD(1) * t433 + t480 * t513 + t519;
t364 = t469 + t483 * t513 + (-t431 + t523) * qJD(1);
t363 = (t431 * t480 + t433 * t483) * qJD(2) + t510;
t362 = rSges(7,1) * t406 + rSges(7,2) * t405 + rSges(7,3) * t413;
t361 = rSges(7,1) * t404 + rSges(7,2) * t403 - rSges(7,3) * t411;
t360 = Icges(7,1) * t406 + Icges(7,4) * t405 + Icges(7,5) * t413;
t359 = Icges(7,1) * t404 + Icges(7,4) * t403 - Icges(7,5) * t411;
t358 = Icges(7,4) * t406 + Icges(7,2) * t405 + Icges(7,6) * t413;
t357 = Icges(7,4) * t404 + Icges(7,2) * t403 - Icges(7,6) * t411;
t356 = Icges(7,5) * t406 + Icges(7,6) * t405 + Icges(7,3) * t413;
t355 = Icges(7,5) * t404 + Icges(7,6) * t403 - Icges(7,3) * t411;
t354 = qJD(1) * t390 + t480 * t491 + t505;
t353 = t483 * t491 + (-t389 + t518) * qJD(1) + t512;
t352 = (t389 * t480 + t390 * t483) * qJD(2) + t489;
t351 = qJD(1) * t377 - t454 * t399 + t485;
t350 = t455 * t399 + (-t376 + t511) * qJD(1) + t486;
t349 = t454 * t376 - t455 * t377 + t487;
t348 = qJD(1) * t379 + t416 * t362 - t401 * t369 - t454 * t400 + t485;
t347 = -t416 * t361 + t402 * t369 + t455 * t400 + (-t378 + t511) * qJD(1) + t486;
t346 = t401 * t361 - t402 * t362 + t454 * t378 - t455 * t379 + t487;
t1 = m(7) * (t346 ^ 2 + t347 ^ 2 + t348 ^ 2) / 0.2e1 + m(6) * (t349 ^ 2 + t350 ^ 2 + t351 ^ 2) / 0.2e1 + m(5) * (t352 ^ 2 + t353 ^ 2 + t354 ^ 2) / 0.2e1 + m(3) * (t380 ^ 2 + t381 ^ 2 + t382 ^ 2) / 0.2e1 + m(4) * (t363 ^ 2 + t364 ^ 2 + t365 ^ 2) / 0.2e1 + t401 * ((t356 * t413 + t358 * t405 + t360 * t406) * t401 + (t355 * t413 + t357 * t405 + t359 * t406) * t402 + (t366 * t413 + t367 * t405 + t368 * t406) * t416) / 0.2e1 + t416 * ((t355 * t402 + t356 * t401 + t366 * t416) * t439 + ((-t358 * t478 + t360 * t481) * t401 + (-t357 * t478 + t359 * t481) * t402 + (-t367 * t478 + t368 * t481) * t416) * t440) / 0.2e1 + t402 * ((-t356 * t411 + t358 * t403 + t360 * t404) * t401 + (-t355 * t411 + t357 * t403 + t359 * t404) * t402 + (-t366 * t411 + t367 * t403 + t368 * t404) * t416) / 0.2e1 + t454 * ((-t480 * t371 - t413 * t373 + t414 * t375) * t454 + (-t480 * t370 - t413 * t372 + t414 * t374) * t455 + (-t480 * t396 - t413 * t397 + t414 * t398) * qJD(1)) / 0.2e1 + t455 * ((t483 * t371 + t411 * t373 + t412 * t375) * t454 + (t483 * t370 + t411 * t372 + t412 * t374) * t455 + (t483 * t396 + t411 * t397 + t412 * t398) * qJD(1)) / 0.2e1 + (Icges(2,3) + m(2) * (t465 ^ 2 + t466 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((-t437 * t385 - t438 * t387 + t539 * t483) * t483 + (t437 * t386 + t438 * t388 + (t540 - t545) * t483 + t544 * t480) * t480) * qJD(2) + (t437 * t408 + t438 * t409 - t538 * t480 + t541 * t483) * qJD(1)) * t474 / 0.2e1 - (((t435 * t386 + t436 * t388 + t540 * t480) * t480 + (-t435 * t385 - t436 * t387 + (t539 - t544) * t480 + t545 * t483) * t483) * qJD(2) + (t435 * t408 + t436 * t409 + t541 * t480 + t538 * t483) * qJD(1)) * t522 / 0.2e1 + ((-t439 * t373 + t440 * t375) * t454 + (-t439 * t372 + t440 * t374) * t455 + ((t492 * t385 - t448 * t387 + t550 * t479 + t552 * t482) * t483 + (-t492 * t386 + t448 * t388 + t549 * t479 - t551 * t482) * t480) * qJD(2) + (-t439 * t397 + t440 * t398 - t492 * t408 + t448 * t409 + t547 * t479 - t548 * t482) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
