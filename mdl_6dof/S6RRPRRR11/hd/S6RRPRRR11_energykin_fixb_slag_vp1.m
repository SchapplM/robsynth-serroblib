% Calculate kinetic energy for
% S6RRPRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR11_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR11_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR11_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR11_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR11_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:29:12
% EndTime: 2019-03-09 14:29:16
% DurationCPUTime: 4.12s
% Computational Cost: add. (1441->312), mult. (2233->502), div. (0->0), fcn. (2190->10), ass. (0->161)
t568 = Icges(3,4) + Icges(4,6);
t567 = Icges(3,1) + Icges(4,2);
t566 = -Icges(3,2) - Icges(4,3);
t487 = cos(qJ(2));
t565 = t568 * t487;
t484 = sin(qJ(2));
t564 = t568 * t484;
t563 = -Icges(4,4) + Icges(3,5);
t562 = Icges(4,5) - Icges(3,6);
t561 = t566 * t484 + t565;
t560 = -t567 * t487 + t564;
t559 = Icges(4,1) + Icges(3,3);
t485 = sin(qJ(1));
t488 = cos(qJ(1));
t558 = t561 * t485 + t562 * t488;
t557 = -t562 * t485 + t561 * t488;
t556 = t560 * t485 + t563 * t488;
t555 = t563 * t485 - t560 * t488;
t554 = t566 * t487 - t564;
t553 = t567 * t484 + t565;
t552 = t562 * t484 + t563 * t487;
t551 = t552 * t485 - t488 * t559;
t550 = t485 * t559 + t552 * t488;
t549 = t563 * t484 - t562 * t487;
t548 = t484 * t554 + t487 * t553;
t547 = -t484 * t557 + t487 * t555;
t546 = t484 * t558 + t487 * t556;
t483 = sin(qJ(4));
t542 = pkin(4) * t483;
t486 = cos(qJ(4));
t540 = t486 * pkin(4);
t534 = t484 * t485;
t533 = t484 * t488;
t532 = t485 * t486;
t531 = t485 * t487;
t530 = t486 * t488;
t529 = t487 * t488;
t482 = qJ(4) + qJ(5);
t512 = pkin(2) * t487 + qJ(3) * t484;
t441 = t512 * t485;
t463 = pkin(1) * t485 - pkin(7) * t488;
t528 = -t441 - t463;
t478 = cos(t482);
t527 = pkin(5) * t478;
t476 = qJD(2) * t485;
t523 = qJD(4) * t487;
t444 = t488 * t523 + t476;
t525 = qJD(2) * t488;
t524 = qJD(3) * t484;
t522 = qJD(5) * t487;
t521 = qJD(6) * t487;
t470 = qJD(4) * t484 + qJD(1);
t442 = t512 * t488;
t450 = qJD(1) * (pkin(1) * t488 + pkin(7) * t485);
t520 = qJD(1) * t442 + t485 * t524 + t450;
t404 = t488 * t522 + t444;
t449 = qJD(5) * t484 + t470;
t477 = sin(t482);
t517 = pkin(5) * t477;
t458 = pkin(2) * t484 - qJ(3) * t487;
t516 = qJD(2) * (rSges(4,2) * t484 + rSges(4,3) * t487 - t458);
t445 = t485 * t523 - t525;
t405 = t485 * t522 + t445;
t515 = -qJD(3) * t487 + t441 * t476 + t442 * t525;
t514 = rSges(3,1) * t487 - rSges(3,2) * t484;
t513 = -rSges(4,2) * t487 + rSges(4,3) * t484;
t511 = qJD(2) * (-pkin(8) * t484 - t458);
t447 = pkin(3) * t485 + pkin(8) * t529;
t448 = -pkin(3) * t488 + pkin(8) * t531;
t498 = t447 * t525 + t448 * t476 + t515;
t497 = pkin(9) * t487 + t484 * t542;
t496 = pkin(10) * t487 + t484 * t517;
t386 = t485 * t540 + t488 * t497;
t387 = t485 * t497 - t488 * t540;
t495 = -t386 * t445 + t444 * t387 + t498;
t494 = qJD(1) * t447 + t485 * t511 + t520;
t469 = t488 * t524;
t493 = t469 + (-t448 + t528) * qJD(1) + t488 * t511;
t432 = pkin(9) * t484 - t487 * t542;
t492 = t470 * t386 - t432 * t444 + t494;
t491 = -t387 * t470 + t445 * t432 + t493;
t479 = qJ(6) + t482;
t472 = cos(t479);
t471 = sin(t479);
t462 = rSges(2,1) * t488 - rSges(2,2) * t485;
t461 = rSges(2,1) * t485 + rSges(2,2) * t488;
t460 = rSges(3,1) * t484 + rSges(3,2) * t487;
t440 = t483 * t534 - t530;
t439 = t483 * t488 + t484 * t532;
t438 = t483 * t533 + t532;
t437 = -t483 * t485 + t484 * t530;
t436 = qJD(6) * t484 + t449;
t431 = t477 * t534 - t478 * t488;
t430 = t477 * t488 + t478 * t534;
t429 = t477 * t533 + t478 * t485;
t428 = -t477 * t485 + t478 * t533;
t427 = -rSges(4,1) * t488 + t485 * t513;
t426 = rSges(4,1) * t485 + t488 * t513;
t425 = rSges(3,3) * t485 + t488 * t514;
t424 = rSges(5,3) * t484 + (-rSges(5,1) * t483 - rSges(5,2) * t486) * t487;
t423 = -rSges(3,3) * t488 + t485 * t514;
t412 = Icges(5,5) * t484 + (-Icges(5,1) * t483 - Icges(5,4) * t486) * t487;
t409 = Icges(5,6) * t484 + (-Icges(5,4) * t483 - Icges(5,2) * t486) * t487;
t406 = Icges(5,3) * t484 + (-Icges(5,5) * t483 - Icges(5,6) * t486) * t487;
t403 = t471 * t534 - t472 * t488;
t402 = t471 * t488 + t472 * t534;
t401 = t471 * t533 + t472 * t485;
t400 = -t471 * t485 + t472 * t533;
t399 = rSges(6,3) * t484 + (-rSges(6,1) * t477 - rSges(6,2) * t478) * t487;
t398 = Icges(6,5) * t484 + (-Icges(6,1) * t477 - Icges(6,4) * t478) * t487;
t397 = Icges(6,6) * t484 + (-Icges(6,4) * t477 - Icges(6,2) * t478) * t487;
t396 = Icges(6,3) * t484 + (-Icges(6,5) * t477 - Icges(6,6) * t478) * t487;
t395 = rSges(7,3) * t484 + (-rSges(7,1) * t471 - rSges(7,2) * t472) * t487;
t394 = Icges(7,5) * t484 + (-Icges(7,1) * t471 - Icges(7,4) * t472) * t487;
t393 = Icges(7,6) * t484 + (-Icges(7,4) * t471 - Icges(7,2) * t472) * t487;
t392 = Icges(7,3) * t484 + (-Icges(7,5) * t471 - Icges(7,6) * t472) * t487;
t391 = t485 * t521 + t405;
t390 = t488 * t521 + t404;
t388 = pkin(10) * t484 - t487 * t517;
t385 = rSges(5,1) * t440 + rSges(5,2) * t439 + rSges(5,3) * t531;
t384 = rSges(5,1) * t438 + rSges(5,2) * t437 + rSges(5,3) * t529;
t383 = Icges(5,1) * t440 + Icges(5,4) * t439 + Icges(5,5) * t531;
t382 = Icges(5,1) * t438 + Icges(5,4) * t437 + Icges(5,5) * t529;
t381 = Icges(5,4) * t440 + Icges(5,2) * t439 + Icges(5,6) * t531;
t380 = Icges(5,4) * t438 + Icges(5,2) * t437 + Icges(5,6) * t529;
t379 = Icges(5,5) * t440 + Icges(5,6) * t439 + Icges(5,3) * t531;
t378 = Icges(5,5) * t438 + Icges(5,6) * t437 + Icges(5,3) * t529;
t376 = qJD(1) * t425 - t460 * t476 + t450;
t375 = -t460 * t525 + (-t423 - t463) * qJD(1);
t374 = (t423 * t485 + t425 * t488) * qJD(2);
t373 = rSges(6,1) * t431 + rSges(6,2) * t430 + rSges(6,3) * t531;
t372 = rSges(6,1) * t429 + rSges(6,2) * t428 + rSges(6,3) * t529;
t371 = Icges(6,1) * t431 + Icges(6,4) * t430 + Icges(6,5) * t531;
t370 = Icges(6,1) * t429 + Icges(6,4) * t428 + Icges(6,5) * t529;
t369 = Icges(6,4) * t431 + Icges(6,2) * t430 + Icges(6,6) * t531;
t368 = Icges(6,4) * t429 + Icges(6,2) * t428 + Icges(6,6) * t529;
t367 = Icges(6,5) * t431 + Icges(6,6) * t430 + Icges(6,3) * t531;
t366 = Icges(6,5) * t429 + Icges(6,6) * t428 + Icges(6,3) * t529;
t364 = rSges(7,1) * t403 + rSges(7,2) * t402 + rSges(7,3) * t531;
t363 = rSges(7,1) * t401 + rSges(7,2) * t400 + rSges(7,3) * t529;
t362 = Icges(7,1) * t403 + Icges(7,4) * t402 + Icges(7,5) * t531;
t361 = Icges(7,1) * t401 + Icges(7,4) * t400 + Icges(7,5) * t529;
t360 = Icges(7,4) * t403 + Icges(7,2) * t402 + Icges(7,6) * t531;
t359 = Icges(7,4) * t401 + Icges(7,2) * t400 + Icges(7,6) * t529;
t358 = Icges(7,5) * t403 + Icges(7,6) * t402 + Icges(7,3) * t531;
t357 = Icges(7,5) * t401 + Icges(7,6) * t400 + Icges(7,3) * t529;
t356 = t485 * t496 - t488 * t527;
t355 = t485 * t527 + t488 * t496;
t354 = qJD(1) * t426 + t485 * t516 + t520;
t353 = t469 + t488 * t516 + (-t427 + t528) * qJD(1);
t352 = (t426 * t488 + t427 * t485) * qJD(2) + t515;
t351 = t384 * t470 - t424 * t444 + t494;
t350 = -t385 * t470 + t424 * t445 + t493;
t349 = -t384 * t445 + t385 * t444 + t498;
t348 = t372 * t449 - t399 * t404 + t492;
t347 = -t373 * t449 + t399 * t405 + t491;
t346 = -t372 * t405 + t373 * t404 + t495;
t345 = t355 * t449 + t363 * t436 - t388 * t404 - t390 * t395 + t492;
t344 = -t356 * t449 - t364 * t436 + t388 * t405 + t391 * t395 + t491;
t343 = -t355 * t405 + t356 * t404 - t363 * t391 + t364 * t390 + t495;
t1 = t436 * ((t357 * t390 + t358 * t391 + t392 * t436) * t484 + ((-t359 * t472 - t361 * t471) * t390 + (-t360 * t472 - t362 * t471) * t391 + (-t393 * t472 - t394 * t471) * t436) * t487) / 0.2e1 + t390 * ((t357 * t529 + t400 * t359 + t401 * t361) * t390 + (t358 * t529 + t360 * t400 + t362 * t401) * t391 + (t392 * t529 + t393 * t400 + t394 * t401) * t436) / 0.2e1 + t391 * ((t357 * t531 + t359 * t402 + t361 * t403) * t390 + (t358 * t531 + t402 * t360 + t403 * t362) * t391 + (t392 * t531 + t393 * t402 + t394 * t403) * t436) / 0.2e1 + t404 * ((t366 * t529 + t428 * t368 + t429 * t370) * t404 + (t367 * t529 + t369 * t428 + t371 * t429) * t405 + (t396 * t529 + t397 * t428 + t398 * t429) * t449) / 0.2e1 + t405 * ((t366 * t531 + t368 * t430 + t370 * t431) * t404 + (t367 * t531 + t430 * t369 + t431 * t371) * t405 + (t396 * t531 + t397 * t430 + t398 * t431) * t449) / 0.2e1 + t449 * ((t366 * t404 + t367 * t405 + t396 * t449) * t484 + ((-t368 * t478 - t370 * t477) * t404 + (-t369 * t478 - t371 * t477) * t405 + (-t397 * t478 - t398 * t477) * t449) * t487) / 0.2e1 + t470 * ((t378 * t444 + t379 * t445 + t406 * t470) * t484 + ((-t380 * t486 - t382 * t483) * t444 + (-t381 * t486 - t383 * t483) * t445 + (-t409 * t486 - t412 * t483) * t470) * t487) / 0.2e1 + t444 * ((t378 * t529 + t437 * t380 + t438 * t382) * t444 + (t379 * t529 + t381 * t437 + t383 * t438) * t445 + (t406 * t529 + t409 * t437 + t412 * t438) * t470) / 0.2e1 + t445 * ((t378 * t531 + t380 * t439 + t382 * t440) * t444 + (t379 * t531 + t439 * t381 + t440 * t383) * t445 + (t406 * t531 + t409 * t439 + t412 * t440) * t470) / 0.2e1 + m(7) * (t343 ^ 2 + t344 ^ 2 + t345 ^ 2) / 0.2e1 + m(6) * (t346 ^ 2 + t347 ^ 2 + t348 ^ 2) / 0.2e1 + m(5) * (t349 ^ 2 + t350 ^ 2 + t351 ^ 2) / 0.2e1 + m(4) * (t352 ^ 2 + t353 ^ 2 + t354 ^ 2) / 0.2e1 + m(3) * (t374 ^ 2 + t375 ^ 2 + t376 ^ 2) / 0.2e1 + (Icges(2,3) + m(2) * (t461 ^ 2 + t462 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((t484 * t556 - t487 * t558) * t488 + (t484 * t555 + t487 * t557) * t485) * qJD(2) + (t553 * t484 - t554 * t487) * qJD(1)) * qJD(1) / 0.2e1 + ((t550 * t485 ^ 2 + (t546 * t488 + (t547 - t551) * t485) * t488) * qJD(2) + (t485 * t549 + t488 * t548) * qJD(1)) * t476 / 0.2e1 - ((t551 * t488 ^ 2 + (t547 * t485 + (t546 - t550) * t488) * t485) * qJD(2) + (t485 * t548 - t488 * t549) * qJD(1)) * t525 / 0.2e1;
T  = t1;
