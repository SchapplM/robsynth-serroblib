% Calculate kinetic energy for
% S6RPRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:51:11
% EndTime: 2019-03-09 03:51:13
% DurationCPUTime: 2.09s
% Computational Cost: add. (1957->303), mult. (1907->480), div. (0->0), fcn. (1897->12), ass. (0->149)
t458 = cos(pkin(10));
t508 = pkin(2) * t458;
t457 = cos(pkin(11));
t507 = t457 * pkin(4);
t454 = pkin(10) + qJ(3);
t445 = sin(t454);
t506 = Icges(4,4) * t445;
t447 = cos(t454);
t505 = Icges(4,4) * t447;
t461 = sin(qJ(1));
t504 = t445 * t461;
t462 = cos(qJ(1));
t503 = t445 * t462;
t502 = t447 * t461;
t501 = t447 * t462;
t455 = sin(pkin(11));
t500 = t455 * t461;
t499 = t455 * t462;
t498 = t457 * t461;
t497 = t457 * t462;
t436 = pkin(1) * t461 - qJ(2) * t462;
t494 = pkin(7) * t462 - t461 * t508 - t436;
t453 = pkin(11) + qJ(5);
t446 = cos(t453);
t493 = pkin(5) * t446;
t450 = qJD(2) * t461;
t489 = qJD(4) * t445;
t492 = t462 * t489 + t450;
t449 = qJD(3) * t461;
t488 = qJD(5) * t445;
t422 = t462 * t488 + t449;
t490 = qJD(3) * t462;
t487 = qJD(6) * t445;
t476 = pkin(3) * t447 + qJ(4) * t445;
t415 = t476 * t461;
t486 = -t415 + t494;
t444 = sin(t453);
t483 = pkin(5) * t444;
t427 = pkin(3) * t445 - qJ(4) * t447;
t482 = qJD(3) * (pkin(8) * t447 - t445 * t507 - t427);
t481 = qJD(3) * (rSges(5,3) * t447 - (rSges(5,1) * t457 - rSges(5,2) * t455) * t445 - t427);
t423 = t461 * t488 - t490;
t431 = qJD(1) * (pkin(1) * t462 + qJ(2) * t461);
t480 = -qJD(2) * t462 + qJD(1) * (pkin(7) * t461 + t462 * t508) + t431;
t416 = t476 * t462;
t479 = -qJD(4) * t447 + t415 * t449 + t416 * t490;
t456 = sin(pkin(10));
t478 = rSges(3,1) * t458 - rSges(3,2) * t456;
t477 = rSges(4,1) * t447 - rSges(4,2) * t445;
t475 = Icges(4,1) * t447 - t506;
t474 = -Icges(4,2) * t445 + t505;
t473 = Icges(4,5) * t447 - Icges(4,6) * t445;
t402 = -Icges(4,6) * t462 + t461 * t474;
t404 = -Icges(4,5) * t462 + t461 * t475;
t472 = t402 * t445 - t404 * t447;
t403 = Icges(4,6) * t461 + t462 * t474;
t405 = Icges(4,5) * t461 + t462 * t475;
t471 = -t403 * t445 + t405 * t447;
t425 = Icges(4,2) * t447 + t506;
t426 = Icges(4,1) * t445 + t505;
t470 = -t425 * t445 + t426 * t447;
t469 = qJD(1) * t416 + t461 * t489 + t480;
t467 = pkin(8) * t445 + t447 * t507;
t366 = -pkin(4) * t499 + t461 * t467;
t367 = pkin(4) * t500 + t462 * t467;
t468 = t366 * t449 + t367 * t490 + t479;
t466 = pkin(9) * t445 + t447 * t493;
t465 = qJD(1) * t367 + t461 * t482 + t469;
t464 = (-t366 + t486) * qJD(1) + t462 * t482 + t492;
t448 = qJ(6) + t453;
t441 = cos(t448);
t440 = sin(t448);
t439 = -qJD(5) * t447 + qJD(1);
t438 = rSges(2,1) * t462 - rSges(2,2) * t461;
t437 = rSges(2,1) * t461 + rSges(2,2) * t462;
t428 = rSges(4,1) * t445 + rSges(4,2) * t447;
t424 = Icges(4,5) * t445 + Icges(4,6) * t447;
t421 = qJD(1) + (-qJD(5) - qJD(6)) * t447;
t420 = t447 * t497 + t500;
t419 = -t447 * t499 + t498;
t418 = t447 * t498 - t499;
t417 = -t447 * t500 - t497;
t413 = t444 * t461 + t446 * t501;
t412 = -t444 * t501 + t446 * t461;
t411 = -t444 * t462 + t446 * t502;
t410 = -t444 * t502 - t446 * t462;
t409 = rSges(4,3) * t461 + t462 * t477;
t408 = -rSges(4,3) * t462 + t461 * t477;
t401 = Icges(4,3) * t461 + t462 * t473;
t400 = -Icges(4,3) * t462 + t461 * t473;
t399 = t440 * t461 + t441 * t501;
t398 = -t440 * t501 + t441 * t461;
t397 = -t440 * t462 + t441 * t502;
t396 = -t440 * t502 - t441 * t462;
t395 = t461 * t487 + t423;
t394 = t462 * t487 + t422;
t390 = -Icges(5,5) * t447 + (Icges(5,1) * t457 - Icges(5,4) * t455) * t445;
t389 = -Icges(5,6) * t447 + (Icges(5,4) * t457 - Icges(5,2) * t455) * t445;
t388 = -Icges(5,3) * t447 + (Icges(5,5) * t457 - Icges(5,6) * t455) * t445;
t387 = -rSges(6,3) * t447 + (rSges(6,1) * t446 - rSges(6,2) * t444) * t445;
t386 = -Icges(6,5) * t447 + (Icges(6,1) * t446 - Icges(6,4) * t444) * t445;
t385 = -Icges(6,6) * t447 + (Icges(6,4) * t446 - Icges(6,2) * t444) * t445;
t384 = -Icges(6,3) * t447 + (Icges(6,5) * t446 - Icges(6,6) * t444) * t445;
t383 = -rSges(7,3) * t447 + (rSges(7,1) * t441 - rSges(7,2) * t440) * t445;
t382 = -Icges(7,5) * t447 + (Icges(7,1) * t441 - Icges(7,4) * t440) * t445;
t381 = -Icges(7,6) * t447 + (Icges(7,4) * t441 - Icges(7,2) * t440) * t445;
t380 = -Icges(7,3) * t447 + (Icges(7,5) * t441 - Icges(7,6) * t440) * t445;
t378 = qJD(1) * t461 * rSges(3,3) + t431 + (qJD(1) * t478 - qJD(2)) * t462;
t377 = t450 + (t462 * rSges(3,3) - t461 * t478 - t436) * qJD(1);
t376 = -pkin(9) * t447 + t445 * t493;
t375 = rSges(5,1) * t420 + rSges(5,2) * t419 + rSges(5,3) * t503;
t374 = rSges(5,1) * t418 + rSges(5,2) * t417 + rSges(5,3) * t504;
t373 = Icges(5,1) * t420 + Icges(5,4) * t419 + Icges(5,5) * t503;
t372 = Icges(5,1) * t418 + Icges(5,4) * t417 + Icges(5,5) * t504;
t371 = Icges(5,4) * t420 + Icges(5,2) * t419 + Icges(5,6) * t503;
t370 = Icges(5,4) * t418 + Icges(5,2) * t417 + Icges(5,6) * t504;
t369 = Icges(5,5) * t420 + Icges(5,6) * t419 + Icges(5,3) * t503;
t368 = Icges(5,5) * t418 + Icges(5,6) * t417 + Icges(5,3) * t504;
t362 = (t408 * t461 + t409 * t462) * qJD(3);
t361 = rSges(6,1) * t413 + rSges(6,2) * t412 + rSges(6,3) * t503;
t360 = rSges(6,1) * t411 + rSges(6,2) * t410 + rSges(6,3) * t504;
t359 = Icges(6,1) * t413 + Icges(6,4) * t412 + Icges(6,5) * t503;
t358 = Icges(6,1) * t411 + Icges(6,4) * t410 + Icges(6,5) * t504;
t357 = Icges(6,4) * t413 + Icges(6,2) * t412 + Icges(6,6) * t503;
t356 = Icges(6,4) * t411 + Icges(6,2) * t410 + Icges(6,6) * t504;
t355 = Icges(6,5) * t413 + Icges(6,6) * t412 + Icges(6,3) * t503;
t354 = Icges(6,5) * t411 + Icges(6,6) * t410 + Icges(6,3) * t504;
t353 = rSges(7,1) * t399 + rSges(7,2) * t398 + rSges(7,3) * t503;
t352 = rSges(7,1) * t397 + rSges(7,2) * t396 + rSges(7,3) * t504;
t351 = Icges(7,1) * t399 + Icges(7,4) * t398 + Icges(7,5) * t503;
t350 = Icges(7,1) * t397 + Icges(7,4) * t396 + Icges(7,5) * t504;
t349 = Icges(7,4) * t399 + Icges(7,2) * t398 + Icges(7,6) * t503;
t348 = Icges(7,4) * t397 + Icges(7,2) * t396 + Icges(7,6) * t504;
t347 = Icges(7,5) * t399 + Icges(7,6) * t398 + Icges(7,3) * t503;
t346 = Icges(7,5) * t397 + Icges(7,6) * t396 + Icges(7,3) * t504;
t345 = t461 * t483 + t462 * t466;
t344 = t461 * t466 - t462 * t483;
t343 = qJD(1) * t409 - t428 * t449 + t480;
t342 = -t428 * t490 + t450 + (-t408 + t494) * qJD(1);
t341 = (t374 * t461 + t375 * t462) * qJD(3) + t479;
t340 = qJD(1) * t375 + t461 * t481 + t469;
t339 = t462 * t481 + (-t374 + t486) * qJD(1) + t492;
t338 = t361 * t439 - t387 * t422 + t465;
t337 = -t360 * t439 + t387 * t423 + t464;
t336 = t360 * t422 - t361 * t423 + t468;
t335 = t345 * t439 + t353 * t421 - t376 * t422 - t383 * t394 + t465;
t334 = -t344 * t439 - t352 * t421 + t376 * t423 + t383 * t395 + t464;
t333 = t344 * t422 - t345 * t423 + t352 * t394 - t353 * t395 + t468;
t1 = m(7) * (t333 ^ 2 + t334 ^ 2 + t335 ^ 2) / 0.2e1 + m(6) * (t336 ^ 2 + t337 ^ 2 + t338 ^ 2) / 0.2e1 + m(5) * (t339 ^ 2 + t340 ^ 2 + t341 ^ 2) / 0.2e1 + m(4) * (t342 ^ 2 + t343 ^ 2 + t362 ^ 2) / 0.2e1 + m(3) * (t377 ^ 2 + t378 ^ 2) / 0.2e1 + t394 * ((t347 * t503 + t398 * t349 + t399 * t351) * t394 + (t346 * t503 + t348 * t398 + t350 * t399) * t395 + (t380 * t503 + t381 * t398 + t382 * t399) * t421) / 0.2e1 + t395 * ((t347 * t504 + t349 * t396 + t351 * t397) * t394 + (t346 * t504 + t396 * t348 + t397 * t350) * t395 + (t380 * t504 + t381 * t396 + t382 * t397) * t421) / 0.2e1 + t439 * ((-t354 * t423 - t355 * t422 - t384 * t439) * t447 + ((-t357 * t444 + t359 * t446) * t422 + (-t356 * t444 + t358 * t446) * t423 + (-t385 * t444 + t386 * t446) * t439) * t445) / 0.2e1 + t421 * ((-t346 * t395 - t347 * t394 - t380 * t421) * t447 + ((-t349 * t440 + t351 * t441) * t394 + (-t348 * t440 + t350 * t441) * t395 + (-t381 * t440 + t382 * t441) * t421) * t445) / 0.2e1 + t422 * ((t355 * t503 + t412 * t357 + t413 * t359) * t422 + (t354 * t503 + t356 * t412 + t358 * t413) * t423 + (t384 * t503 + t385 * t412 + t386 * t413) * t439) / 0.2e1 + t423 * ((t355 * t504 + t357 * t410 + t359 * t411) * t422 + (t354 * t504 + t410 * t356 + t411 * t358) * t423 + (t384 * t504 + t385 * t410 + t386 * t411) * t439) / 0.2e1 + (((t403 * t447 + t405 * t445) * t461 - (t402 * t447 + t404 * t445) * t462 + (t368 * t462 - t369 * t461) * t447 + ((-t371 * t455 + t373 * t457) * t461 - (-t370 * t455 + t372 * t457) * t462) * t445) * qJD(3) + ((t425 - t388) * t447 + (-t389 * t455 + t390 * t457 + t426) * t445) * qJD(1)) * qJD(1) / 0.2e1 + (((-t368 * t503 - t370 * t419 - t372 * t420 + t472 * t462) * t462 + ((-t400 + t471) * t462 + t369 * t503 + t371 * t419 + t373 * t420 + t401 * t461) * t461) * qJD(3) + (t388 * t503 + t389 * t419 + t390 * t420 + t461 * t424 + t462 * t470) * qJD(1)) * t449 / 0.2e1 - (((-t368 * t504 - t370 * t417 - t372 * t418 + t400 * t462) * t462 + ((-t401 + t472) * t462 + t369 * t504 + t371 * t417 + t373 * t418 + t471 * t461) * t461) * qJD(3) + (t388 * t504 + t389 * t417 + t390 * t418 - t462 * t424 + t461 * t470) * qJD(1)) * t490 / 0.2e1 + (Icges(2,3) + Icges(3,2) * t458 ^ 2 + (Icges(3,1) * t456 + 0.2e1 * Icges(3,4) * t458) * t456 + m(2) * (t437 ^ 2 + t438 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
