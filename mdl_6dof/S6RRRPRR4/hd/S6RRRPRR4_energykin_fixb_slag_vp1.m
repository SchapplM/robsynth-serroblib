% Calculate kinetic energy for
% S6RRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:15:22
% EndTime: 2019-03-09 18:15:24
% DurationCPUTime: 2.62s
% Computational Cost: add. (2145->339), mult. (2218->538), div. (0->0), fcn. (2175->12), ass. (0->172)
t470 = cos(qJ(2));
t524 = pkin(2) * t470;
t466 = cos(pkin(11));
t523 = t466 * pkin(4);
t468 = sin(qJ(2));
t521 = Icges(3,4) * t468;
t520 = Icges(3,4) * t470;
t464 = qJ(2) + qJ(3);
t459 = sin(t464);
t519 = Icges(4,4) * t459;
t460 = cos(t464);
t518 = Icges(4,4) * t460;
t469 = sin(qJ(1));
t517 = t459 * t469;
t471 = cos(qJ(1));
t516 = t459 * t471;
t515 = t460 * t469;
t514 = t460 * t471;
t465 = sin(pkin(11));
t513 = t465 * t469;
t512 = t465 * t471;
t511 = t466 * t469;
t510 = t466 * t471;
t396 = -pkin(8) * t471 + t469 * t524;
t397 = pkin(8) * t469 + t471 * t524;
t458 = qJD(2) * t469;
t504 = qJD(2) * t471;
t508 = t396 * t458 + t397 * t504;
t448 = pkin(1) * t469 - pkin(7) * t471;
t507 = -t396 - t448;
t463 = pkin(11) + qJ(5);
t455 = cos(t463);
t506 = pkin(5) * t455;
t436 = qJD(3) * t469 + t458;
t503 = qJD(4) * t459;
t502 = qJD(5) * t459;
t501 = qJD(6) * t459;
t500 = pkin(2) * qJD(2) * t468;
t494 = pkin(3) * t460 + qJ(4) * t459;
t421 = t494 * t469;
t499 = -t421 + t507;
t417 = t471 * t502 + t436;
t454 = sin(t463);
t498 = pkin(5) * t454;
t437 = (-qJD(2) - qJD(3)) * t471;
t497 = t471 * t500;
t496 = rSges(3,1) * t470 - rSges(3,2) * t468;
t495 = rSges(4,1) * t460 - rSges(4,2) * t459;
t493 = Icges(3,1) * t470 - t521;
t492 = Icges(4,1) * t460 - t519;
t491 = -Icges(3,2) * t468 + t520;
t490 = -Icges(4,2) * t459 + t518;
t489 = Icges(3,5) * t470 - Icges(3,6) * t468;
t488 = Icges(4,5) * t460 - Icges(4,6) * t459;
t413 = -Icges(3,6) * t471 + t469 * t491;
t415 = -Icges(3,5) * t471 + t469 * t493;
t487 = t413 * t468 - t415 * t470;
t414 = Icges(3,6) * t469 + t471 * t491;
t416 = Icges(3,5) * t469 + t471 * t493;
t486 = -t414 * t468 + t416 * t470;
t439 = Icges(3,2) * t470 + t521;
t440 = Icges(3,1) * t468 + t520;
t485 = -t439 * t468 + t440 * t470;
t418 = t469 * t502 + t437;
t484 = -qJD(4) * t460 + t436 * t421 + t508;
t435 = qJD(1) * (pkin(1) * t471 + pkin(7) * t469);
t483 = qJD(1) * t397 - t469 * t500 + t435;
t431 = pkin(3) * t459 - qJ(4) * t460;
t482 = t437 * t431 + t471 * t503 - t497;
t481 = qJD(1) * (Icges(4,5) * t459 + Icges(4,6) * t460) + (-Icges(4,3) * t471 + t469 * t488) * t437 + (Icges(4,3) * t469 + t471 * t488) * t436;
t480 = pkin(9) * t459 + t460 * t523;
t479 = pkin(10) * t459 + t460 * t506;
t422 = t494 * t471;
t478 = qJD(1) * t422 + t469 * t503 + t483;
t357 = -pkin(4) * t512 + t469 * t480;
t358 = pkin(4) * t513 + t471 * t480;
t477 = t436 * t357 + (-t358 - t422) * t437 + t484;
t372 = -pkin(9) * t460 + t459 * t523;
t476 = qJD(1) * t358 + (-t372 - t431) * t436 + t478;
t475 = t437 * t372 + (-t357 + t499) * qJD(1) + t482;
t400 = -Icges(4,6) * t471 + t469 * t490;
t401 = Icges(4,6) * t469 + t471 * t490;
t402 = -Icges(4,5) * t471 + t469 * t492;
t403 = Icges(4,5) * t469 + t471 * t492;
t429 = Icges(4,2) * t460 + t519;
t430 = Icges(4,1) * t459 + t518;
t474 = (-t401 * t459 + t403 * t460) * t436 + (-t400 * t459 + t402 * t460) * t437 + (-t429 * t459 + t430 * t460) * qJD(1);
t456 = qJ(6) + t463;
t451 = cos(t456);
t450 = sin(t456);
t449 = -qJD(5) * t460 + qJD(1);
t443 = rSges(2,1) * t471 - rSges(2,2) * t469;
t442 = rSges(2,1) * t469 + rSges(2,2) * t471;
t441 = rSges(3,1) * t468 + rSges(3,2) * t470;
t438 = Icges(3,5) * t468 + Icges(3,6) * t470;
t432 = rSges(4,1) * t459 + rSges(4,2) * t460;
t427 = qJD(1) + (-qJD(5) - qJD(6)) * t460;
t426 = t460 * t510 + t513;
t425 = -t460 * t512 + t511;
t424 = t460 * t511 - t512;
t423 = -t460 * t513 - t510;
t420 = rSges(3,3) * t469 + t471 * t496;
t419 = -rSges(3,3) * t471 + t469 * t496;
t412 = Icges(3,3) * t469 + t471 * t489;
t411 = -Icges(3,3) * t471 + t469 * t489;
t409 = t454 * t469 + t455 * t514;
t408 = -t454 * t514 + t455 * t469;
t407 = -t454 * t471 + t455 * t515;
t406 = -t454 * t515 - t455 * t471;
t405 = rSges(4,3) * t469 + t471 * t495;
t404 = -rSges(4,3) * t471 + t469 * t495;
t394 = t450 * t469 + t451 * t514;
t393 = -t450 * t514 + t451 * t469;
t392 = -t450 * t471 + t451 * t515;
t391 = -t450 * t515 - t451 * t471;
t390 = -rSges(5,3) * t460 + (rSges(5,1) * t466 - rSges(5,2) * t465) * t459;
t388 = -Icges(5,5) * t460 + (Icges(5,1) * t466 - Icges(5,4) * t465) * t459;
t387 = -Icges(5,6) * t460 + (Icges(5,4) * t466 - Icges(5,2) * t465) * t459;
t386 = -Icges(5,3) * t460 + (Icges(5,5) * t466 - Icges(5,6) * t465) * t459;
t383 = t469 * t501 + t418;
t382 = t471 * t501 + t417;
t380 = -rSges(6,3) * t460 + (rSges(6,1) * t455 - rSges(6,2) * t454) * t459;
t379 = -Icges(6,5) * t460 + (Icges(6,1) * t455 - Icges(6,4) * t454) * t459;
t378 = -Icges(6,6) * t460 + (Icges(6,4) * t455 - Icges(6,2) * t454) * t459;
t377 = -Icges(6,3) * t460 + (Icges(6,5) * t455 - Icges(6,6) * t454) * t459;
t376 = -rSges(7,3) * t460 + (rSges(7,1) * t451 - rSges(7,2) * t450) * t459;
t375 = -Icges(7,5) * t460 + (Icges(7,1) * t451 - Icges(7,4) * t450) * t459;
t374 = -Icges(7,6) * t460 + (Icges(7,4) * t451 - Icges(7,2) * t450) * t459;
t373 = -Icges(7,3) * t460 + (Icges(7,5) * t451 - Icges(7,6) * t450) * t459;
t370 = -pkin(10) * t460 + t459 * t506;
t369 = qJD(1) * t420 - t441 * t458 + t435;
t368 = -t441 * t504 + (-t419 - t448) * qJD(1);
t367 = (t419 * t469 + t420 * t471) * qJD(2);
t366 = rSges(5,1) * t426 + rSges(5,2) * t425 + rSges(5,3) * t516;
t365 = rSges(5,1) * t424 + rSges(5,2) * t423 + rSges(5,3) * t517;
t364 = Icges(5,1) * t426 + Icges(5,4) * t425 + Icges(5,5) * t516;
t363 = Icges(5,1) * t424 + Icges(5,4) * t423 + Icges(5,5) * t517;
t362 = Icges(5,4) * t426 + Icges(5,2) * t425 + Icges(5,6) * t516;
t361 = Icges(5,4) * t424 + Icges(5,2) * t423 + Icges(5,6) * t517;
t360 = Icges(5,5) * t426 + Icges(5,6) * t425 + Icges(5,3) * t516;
t359 = Icges(5,5) * t424 + Icges(5,6) * t423 + Icges(5,3) * t517;
t355 = rSges(6,1) * t409 + rSges(6,2) * t408 + rSges(6,3) * t516;
t354 = rSges(6,1) * t407 + rSges(6,2) * t406 + rSges(6,3) * t517;
t353 = Icges(6,1) * t409 + Icges(6,4) * t408 + Icges(6,5) * t516;
t352 = Icges(6,1) * t407 + Icges(6,4) * t406 + Icges(6,5) * t517;
t351 = Icges(6,4) * t409 + Icges(6,2) * t408 + Icges(6,6) * t516;
t350 = Icges(6,4) * t407 + Icges(6,2) * t406 + Icges(6,6) * t517;
t349 = Icges(6,5) * t409 + Icges(6,6) * t408 + Icges(6,3) * t516;
t348 = Icges(6,5) * t407 + Icges(6,6) * t406 + Icges(6,3) * t517;
t346 = rSges(7,1) * t394 + rSges(7,2) * t393 + rSges(7,3) * t516;
t345 = rSges(7,1) * t392 + rSges(7,2) * t391 + rSges(7,3) * t517;
t344 = Icges(7,1) * t394 + Icges(7,4) * t393 + Icges(7,5) * t516;
t343 = Icges(7,1) * t392 + Icges(7,4) * t391 + Icges(7,5) * t517;
t342 = Icges(7,4) * t394 + Icges(7,2) * t393 + Icges(7,6) * t516;
t341 = Icges(7,4) * t392 + Icges(7,2) * t391 + Icges(7,6) * t517;
t340 = Icges(7,5) * t394 + Icges(7,6) * t393 + Icges(7,3) * t516;
t339 = Icges(7,5) * t392 + Icges(7,6) * t391 + Icges(7,3) * t517;
t338 = t469 * t498 + t471 * t479;
t337 = t469 * t479 - t471 * t498;
t336 = qJD(1) * t405 - t432 * t436 + t483;
t335 = -t497 + t432 * t437 + (-t404 + t507) * qJD(1);
t334 = t404 * t436 - t405 * t437 + t508;
t333 = qJD(1) * t366 + (-t390 - t431) * t436 + t478;
t332 = t390 * t437 + (-t365 + t499) * qJD(1) + t482;
t331 = t365 * t436 + (-t366 - t422) * t437 + t484;
t330 = t355 * t449 - t380 * t417 + t476;
t329 = -t354 * t449 + t380 * t418 + t475;
t328 = t354 * t417 - t355 * t418 + t477;
t327 = t338 * t449 + t346 * t427 - t370 * t417 - t376 * t382 + t476;
t326 = -t337 * t449 - t345 * t427 + t370 * t418 + t376 * t383 + t475;
t325 = t337 * t417 - t338 * t418 + t345 * t382 - t346 * t383 + t477;
t1 = m(3) * (t367 ^ 2 + t368 ^ 2 + t369 ^ 2) / 0.2e1 + m(4) * (t334 ^ 2 + t335 ^ 2 + t336 ^ 2) / 0.2e1 + m(5) * (t331 ^ 2 + t332 ^ 2 + t333 ^ 2) / 0.2e1 + m(6) * (t328 ^ 2 + t329 ^ 2 + t330 ^ 2) / 0.2e1 + m(7) * (t325 ^ 2 + t326 ^ 2 + t327 ^ 2) / 0.2e1 + t417 * ((t349 * t516 + t408 * t351 + t409 * t353) * t417 + (t348 * t516 + t350 * t408 + t352 * t409) * t418 + (t377 * t516 + t378 * t408 + t379 * t409) * t449) / 0.2e1 + t418 * ((t349 * t517 + t351 * t406 + t353 * t407) * t417 + (t348 * t517 + t406 * t350 + t407 * t352) * t418 + (t377 * t517 + t378 * t406 + t379 * t407) * t449) / 0.2e1 + t449 * ((-t348 * t418 - t349 * t417 - t377 * t449) * t460 + ((-t351 * t454 + t353 * t455) * t417 + (-t350 * t454 + t352 * t455) * t418 + (-t378 * t454 + t379 * t455) * t449) * t459) / 0.2e1 + t382 * ((t340 * t516 + t393 * t342 + t394 * t344) * t382 + (t339 * t516 + t341 * t393 + t343 * t394) * t383 + (t373 * t516 + t374 * t393 + t375 * t394) * t427) / 0.2e1 + t383 * ((t340 * t517 + t342 * t391 + t344 * t392) * t382 + (t339 * t517 + t391 * t341 + t392 * t343) * t383 + (t373 * t517 + t374 * t391 + t375 * t392) * t427) / 0.2e1 + t427 * ((-t339 * t383 - t340 * t382 - t373 * t427) * t460 + ((-t342 * t450 + t344 * t451) * t382 + (-t341 * t450 + t343 * t451) * t383 + (-t374 * t450 + t375 * t451) * t427) * t459) / 0.2e1 - ((-t471 * t438 + t469 * t485) * qJD(1) + (t471 ^ 2 * t411 + (t486 * t469 + (-t412 + t487) * t471) * t469) * qJD(2)) * t504 / 0.2e1 + ((t469 * t438 + t471 * t485) * qJD(1) + (t469 ^ 2 * t412 + (t487 * t471 + (-t411 + t486) * t469) * t471) * qJD(2)) * t458 / 0.2e1 + (t469 * t481 + t471 * t474 + (t360 * t516 + t425 * t362 + t426 * t364) * t436 + (t359 * t516 + t361 * t425 + t363 * t426) * t437 + (t386 * t516 + t387 * t425 + t388 * t426) * qJD(1)) * t436 / 0.2e1 + (t469 * t474 - t481 * t471 + (t360 * t517 + t362 * t423 + t364 * t424) * t436 + (t359 * t517 + t361 * t423 + t424 * t363) * t437 + (t386 * t517 + t387 * t423 + t388 * t424) * qJD(1)) * t437 / 0.2e1 + (m(2) * (t442 ^ 2 + t443 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t414 * t470 + t416 * t468) * t469 - (t413 * t470 + t415 * t468) * t471) * qJD(2) + (t401 * t460 + t403 * t459) * t436 + (t400 * t460 + t402 * t459) * t437 + (-t359 * t437 - t360 * t436) * t460 + ((-t362 * t465 + t364 * t466) * t436 + (-t361 * t465 + t363 * t466) * t437) * t459 + (t470 * t439 + t468 * t440 + (t429 - t386) * t460 + (-t387 * t465 + t388 * t466 + t430) * t459) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
