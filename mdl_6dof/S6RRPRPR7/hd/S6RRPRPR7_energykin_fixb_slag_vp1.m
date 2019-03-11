% Calculate kinetic energy for
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 10:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR7_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR7_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR7_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR7_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:45:16
% EndTime: 2019-03-09 10:45:20
% DurationCPUTime: 3.60s
% Computational Cost: add. (1523->291), mult. (2485->435), div. (0->0), fcn. (2642->10), ass. (0->151)
t568 = Icges(3,4) - Icges(4,5);
t567 = Icges(3,1) + Icges(4,1);
t566 = Icges(3,2) + Icges(4,3);
t486 = cos(qJ(2));
t565 = t568 * t486;
t482 = sin(qJ(2));
t564 = t568 * t482;
t563 = Icges(4,4) + Icges(3,5);
t562 = Icges(3,6) - Icges(4,6);
t561 = t566 * t482 - t565;
t560 = t567 * t486 - t564;
t559 = Icges(4,2) + Icges(3,3);
t558 = Icges(5,3) + Icges(6,3);
t483 = sin(qJ(1));
t487 = cos(qJ(1));
t557 = t561 * t483 + t562 * t487;
t556 = -t562 * t483 + t561 * t487;
t555 = -t560 * t483 + t563 * t487;
t554 = t563 * t483 + t560 * t487;
t553 = -t566 * t486 - t564;
t552 = t567 * t482 + t565;
t551 = -t562 * t482 + t563 * t486;
t522 = qJ(4) + pkin(10);
t475 = sin(t522);
t517 = cos(t522);
t509 = t482 * t517;
t444 = -t486 * t475 + t509;
t415 = t444 * t483;
t443 = t482 * t475 + t486 * t517;
t416 = t443 * t483;
t485 = cos(qJ(4));
t481 = sin(qJ(4));
t528 = t481 * t486;
t452 = t482 * t485 - t528;
t439 = t452 * t483;
t529 = t481 * t482;
t496 = t485 * t486 + t529;
t440 = t496 * t483;
t550 = Icges(5,5) * t440 + Icges(6,5) * t416 + Icges(5,6) * t439 + Icges(6,6) * t415 + t558 * t487;
t526 = t486 * t487;
t417 = t475 * t526 - t487 * t509;
t418 = t443 * t487;
t441 = t452 * t487;
t442 = t496 * t487;
t549 = Icges(5,5) * t442 + Icges(6,5) * t418 + Icges(5,6) * t441 - Icges(6,6) * t417 - t558 * t483;
t548 = Icges(5,5) * t452 + Icges(6,5) * t444 - Icges(5,6) * t496 - Icges(6,6) * t443;
t547 = t551 * t483 - t559 * t487;
t546 = t559 * t483 + t551 * t487;
t545 = t563 * t482 + t562 * t486;
t544 = t553 * t482 + t552 * t486;
t543 = t556 * t482 + t554 * t486;
t542 = -t557 * t482 + t555 * t486;
t535 = pkin(4) * t485;
t511 = pkin(2) * t486 + qJ(3) * t482;
t448 = t511 * t483;
t471 = pkin(1) * t483 - pkin(7) * t487;
t525 = -t448 - t471;
t478 = qJD(2) * t483;
t524 = qJD(2) * t487;
t523 = qJD(3) * t482;
t449 = t511 * t487;
t455 = qJD(1) * (pkin(1) * t487 + pkin(7) * t483);
t521 = qJD(1) * t449 + t483 * t523 + t455;
t453 = pkin(3) * t483 * t486 + pkin(8) * t487;
t520 = -t453 + t525;
t466 = pkin(2) * t482 - qJ(3) * t486;
t516 = qJD(2) * (-rSges(4,1) * t482 + rSges(4,3) * t486 - t466);
t459 = qJD(4) * t487 - t524;
t458 = -qJD(4) * t483 + t478;
t494 = pkin(4) * t529 + t486 * t535;
t397 = qJ(5) * t487 + t483 * t494;
t515 = -t397 + t520;
t514 = -qJD(3) * t486 + t448 * t478 + t449 * t524;
t513 = rSges(3,1) * t486 - rSges(3,2) * t482;
t512 = rSges(4,1) * t486 + rSges(4,3) * t482;
t510 = qJD(2) * (-pkin(3) * t482 - t466);
t454 = pkin(3) * t526 - pkin(8) * t483;
t495 = t453 * t478 + t454 * t524 + t514;
t473 = t487 * t523;
t493 = t487 * t510 + t473;
t492 = t458 * t397 + t495;
t491 = qJD(1) * t454 + t483 * t510 + t521;
t420 = -pkin(4) * t528 + t482 * t535;
t490 = -qJD(5) * t483 + t459 * t420 + t493;
t398 = -qJ(5) * t483 + t487 * t494;
t489 = qJD(1) * t398 + qJD(5) * t487 + t491;
t484 = cos(qJ(6));
t480 = sin(qJ(6));
t470 = rSges(2,1) * t487 - rSges(2,2) * t483;
t469 = rSges(2,1) * t483 + rSges(2,2) * t487;
t468 = rSges(3,1) * t482 + rSges(3,2) * t486;
t438 = rSges(3,3) * t483 + t487 * t513;
t437 = rSges(4,2) * t483 + t487 * t512;
t436 = -rSges(3,3) * t487 + t483 * t513;
t435 = -rSges(4,2) * t487 + t483 * t512;
t419 = qJD(6) * t443 + qJD(1);
t414 = rSges(5,1) * t452 - rSges(5,2) * t496;
t413 = Icges(5,1) * t452 - Icges(5,4) * t496;
t412 = Icges(5,4) * t452 - Icges(5,2) * t496;
t409 = t418 * t484 - t480 * t483;
t408 = -t418 * t480 - t483 * t484;
t407 = t416 * t484 + t480 * t487;
t406 = -t416 * t480 + t484 * t487;
t405 = -qJD(6) * t415 + t459;
t404 = qJD(6) * t417 + t458;
t403 = pkin(5) * t444 + pkin(9) * t443;
t402 = rSges(6,1) * t444 - rSges(6,2) * t443;
t401 = Icges(6,1) * t444 - Icges(6,4) * t443;
t400 = Icges(6,4) * t444 - Icges(6,2) * t443;
t395 = rSges(5,1) * t442 + rSges(5,2) * t441 - rSges(5,3) * t483;
t394 = rSges(5,1) * t440 + rSges(5,2) * t439 + rSges(5,3) * t487;
t393 = Icges(5,1) * t442 + Icges(5,4) * t441 - Icges(5,5) * t483;
t392 = Icges(5,1) * t440 + Icges(5,4) * t439 + Icges(5,5) * t487;
t391 = Icges(5,4) * t442 + Icges(5,2) * t441 - Icges(5,6) * t483;
t390 = Icges(5,4) * t440 + Icges(5,2) * t439 + Icges(5,6) * t487;
t387 = qJD(1) * t438 - t468 * t478 + t455;
t386 = -t468 * t524 + (-t436 - t471) * qJD(1);
t385 = (t436 * t483 + t438 * t487) * qJD(2);
t384 = pkin(5) * t418 + pkin(9) * t417;
t383 = pkin(5) * t416 - pkin(9) * t415;
t381 = rSges(6,1) * t418 - rSges(6,2) * t417 - rSges(6,3) * t483;
t380 = rSges(6,1) * t416 + rSges(6,2) * t415 + rSges(6,3) * t487;
t379 = Icges(6,1) * t418 - Icges(6,4) * t417 - Icges(6,5) * t483;
t378 = Icges(6,1) * t416 + Icges(6,4) * t415 + Icges(6,5) * t487;
t377 = Icges(6,4) * t418 - Icges(6,2) * t417 - Icges(6,6) * t483;
t376 = Icges(6,4) * t416 + Icges(6,2) * t415 + Icges(6,6) * t487;
t373 = rSges(7,3) * t443 + (rSges(7,1) * t484 - rSges(7,2) * t480) * t444;
t372 = Icges(7,5) * t443 + (Icges(7,1) * t484 - Icges(7,4) * t480) * t444;
t371 = Icges(7,6) * t443 + (Icges(7,4) * t484 - Icges(7,2) * t480) * t444;
t370 = Icges(7,3) * t443 + (Icges(7,5) * t484 - Icges(7,6) * t480) * t444;
t369 = qJD(1) * t437 + t483 * t516 + t521;
t368 = t473 + t487 * t516 + (-t435 + t525) * qJD(1);
t367 = (t435 * t483 + t437 * t487) * qJD(2) + t514;
t366 = rSges(7,1) * t409 + rSges(7,2) * t408 + rSges(7,3) * t417;
t365 = rSges(7,1) * t407 + rSges(7,2) * t406 - rSges(7,3) * t415;
t364 = Icges(7,1) * t409 + Icges(7,4) * t408 + Icges(7,5) * t417;
t363 = Icges(7,1) * t407 + Icges(7,4) * t406 - Icges(7,5) * t415;
t362 = Icges(7,4) * t409 + Icges(7,2) * t408 + Icges(7,6) * t417;
t361 = Icges(7,4) * t407 + Icges(7,2) * t406 - Icges(7,6) * t415;
t360 = Icges(7,5) * t409 + Icges(7,6) * t408 + Icges(7,3) * t417;
t359 = Icges(7,5) * t407 + Icges(7,6) * t406 - Icges(7,3) * t415;
t358 = qJD(1) * t395 - t414 * t458 + t491;
t357 = t414 * t459 + (-t394 + t520) * qJD(1) + t493;
t356 = t394 * t458 - t395 * t459 + t495;
t355 = qJD(1) * t381 + (-t402 - t420) * t458 + t489;
t354 = t402 * t459 + (-t380 + t515) * qJD(1) + t490;
t353 = t380 * t458 + (-t381 - t398) * t459 + t492;
t352 = qJD(1) * t384 + t366 * t419 - t373 * t404 + (-t403 - t420) * t458 + t489;
t351 = -t365 * t419 + t373 * t405 + t403 * t459 + (-t383 + t515) * qJD(1) + t490;
t350 = t365 * t404 - t366 * t405 + t383 * t458 + (-t384 - t398) * t459 + t492;
t1 = m(7) * (t350 ^ 2 + t351 ^ 2 + t352 ^ 2) / 0.2e1 + m(6) * (t353 ^ 2 + t354 ^ 2 + t355 ^ 2) / 0.2e1 + m(5) * (t356 ^ 2 + t357 ^ 2 + t358 ^ 2) / 0.2e1 + m(4) * (t367 ^ 2 + t368 ^ 2 + t369 ^ 2) / 0.2e1 + m(3) * (t385 ^ 2 + t386 ^ 2 + t387 ^ 2) / 0.2e1 + t419 * ((t359 * t405 + t360 * t404 + t370 * t419) * t443 + ((-t362 * t480 + t364 * t484) * t404 + (-t361 * t480 + t363 * t484) * t405 + (-t371 * t480 + t372 * t484) * t419) * t444) / 0.2e1 + t405 * ((-t360 * t415 + t362 * t406 + t364 * t407) * t404 + (-t359 * t415 + t361 * t406 + t407 * t363) * t405 + (-t370 * t415 + t371 * t406 + t372 * t407) * t419) / 0.2e1 + t404 * ((t360 * t417 + t362 * t408 + t364 * t409) * t404 + (t359 * t417 + t361 * t408 + t363 * t409) * t405 + (t370 * t417 + t371 * t408 + t372 * t409) * t419) / 0.2e1 + ((-t376 * t417 + t378 * t418 + t390 * t441 + t392 * t442 - t550 * t483) * t459 + (-t377 * t417 + t379 * t418 + t391 * t441 + t393 * t442 - t549 * t483) * t458 + (-t400 * t417 + t401 * t418 + t412 * t441 + t413 * t442 - t548 * t483) * qJD(1)) * t458 / 0.2e1 + ((t376 * t415 + t378 * t416 + t390 * t439 + t392 * t440 + t550 * t487) * t459 + (t377 * t415 + t379 * t416 + t391 * t439 + t393 * t440 + t549 * t487) * t458 + (t400 * t415 + t401 * t416 + t412 * t439 + t413 * t440 + t548 * t487) * qJD(1)) * t459 / 0.2e1 + (Icges(2,3) + m(2) * (t469 ^ 2 + t470 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((t546 * t483 ^ 2 + (t542 * t487 + (t543 - t547) * t483) * t487) * qJD(2) + (t545 * t483 + t544 * t487) * qJD(1)) * t478 / 0.2e1 - ((t547 * t487 ^ 2 + (t543 * t483 + (t542 - t546) * t487) * t483) * qJD(2) + (t544 * t483 - t545 * t487) * qJD(1)) * t524 / 0.2e1 + ((-t376 * t443 + t378 * t444 - t390 * t496 + t392 * t452) * t459 + (-t377 * t443 + t379 * t444 - t391 * t496 + t393 * t452) * t458 + ((t555 * t482 + t557 * t486) * t487 + (t554 * t482 - t556 * t486) * t483) * qJD(2) + (-t443 * t400 + t444 * t401 - t496 * t412 + t452 * t413 + t552 * t482 - t553 * t486) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
