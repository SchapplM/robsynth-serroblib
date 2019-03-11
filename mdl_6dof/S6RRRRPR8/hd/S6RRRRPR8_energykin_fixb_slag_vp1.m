% Calculate kinetic energy for
% S6RRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR8_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR8_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR8_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR8_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:38:49
% EndTime: 2019-03-09 22:38:52
% DurationCPUTime: 3.14s
% Computational Cost: add. (1916->319), mult. (2957->499), div. (0->0), fcn. (3134->10), ass. (0->156)
t524 = Icges(5,1) + Icges(6,1);
t523 = -Icges(5,4) + Icges(6,5);
t522 = Icges(6,4) + Icges(5,5);
t521 = Icges(5,2) + Icges(6,3);
t520 = -Icges(6,6) + Icges(5,6);
t519 = -Icges(5,3) - Icges(6,2);
t462 = qJ(3) + qJ(4);
t460 = sin(t462);
t461 = cos(t462);
t470 = cos(qJ(1));
t466 = sin(qJ(1));
t469 = cos(qJ(2));
t497 = t466 * t469;
t423 = t460 * t497 + t461 * t470;
t424 = -t460 * t470 + t461 * t497;
t465 = sin(qJ(2));
t499 = t465 * t466;
t518 = t521 * t423 + t523 * t424 - t520 * t499;
t496 = t469 * t470;
t425 = t460 * t496 - t466 * t461;
t426 = t460 * t466 + t461 * t496;
t498 = t465 * t470;
t517 = t521 * t425 + t523 * t426 - t520 * t498;
t516 = -t520 * t423 + t522 * t424 - t519 * t499;
t515 = -t520 * t425 + t522 * t426 - t519 * t498;
t514 = t523 * t423 + t524 * t424 + t522 * t499;
t513 = t523 * t425 + t524 * t426 + t522 * t498;
t512 = t520 * t469 + (t521 * t460 + t523 * t461) * t465;
t511 = t519 * t469 + (-t520 * t460 + t522 * t461) * t465;
t510 = -t522 * t469 + (t523 * t460 + t524 * t461) * t465;
t468 = cos(qJ(3));
t505 = pkin(3) * t468;
t503 = Icges(3,4) * t465;
t502 = Icges(3,4) * t469;
t464 = sin(qJ(3));
t501 = t464 * t466;
t500 = t464 * t470;
t489 = pkin(2) * t469 + pkin(8) * t465;
t436 = t489 * t466;
t437 = t489 * t470;
t459 = qJD(2) * t466;
t494 = qJD(2) * t470;
t495 = t436 * t459 + t437 * t494;
t493 = qJD(3) * t465;
t439 = t470 * t493 + t459;
t492 = qJD(4) * t465;
t491 = qJD(6) * t465;
t490 = -qJD(3) - qJD(4);
t409 = t470 * t492 + t439;
t440 = t466 * t493 - t494;
t410 = t466 * t492 + t440;
t488 = rSges(3,1) * t469 - rSges(3,2) * t465;
t487 = Icges(3,1) * t469 - t503;
t486 = -Icges(3,2) * t465 + t502;
t485 = Icges(3,5) * t469 - Icges(3,6) * t465;
t415 = -Icges(3,6) * t470 + t466 * t486;
t418 = -Icges(3,5) * t470 + t466 * t487;
t484 = t415 * t465 - t418 * t469;
t416 = Icges(3,6) * t466 + t470 * t486;
t419 = Icges(3,5) * t466 + t470 * t487;
t483 = -t416 * t465 + t419 * t469;
t444 = Icges(3,2) * t469 + t503;
t445 = Icges(3,1) * t465 + t502;
t482 = -t444 * t465 + t445 * t469;
t478 = pkin(9) * t465 + t469 * t505;
t388 = -pkin(3) * t500 + t466 * t478;
t389 = pkin(3) * t501 + t470 * t478;
t481 = t439 * t388 - t389 * t440 + t495;
t442 = qJD(1) * (pkin(1) * t470 + pkin(7) * t466);
t450 = pkin(2) * t465 - pkin(8) * t469;
t480 = qJD(1) * t437 - t450 * t459 + t442;
t378 = pkin(4) * t424 + qJ(5) * t423;
t479 = qJD(5) * t465 * t460 + t409 * t378 + t481;
t451 = pkin(1) * t466 - pkin(7) * t470;
t477 = (-t436 - t451) * qJD(1) - t450 * t494;
t396 = -pkin(9) * t469 + t465 * t505;
t457 = -qJD(3) * t469 + qJD(1);
t476 = t457 * t389 - t396 * t439 + t480;
t379 = pkin(4) * t426 + qJ(5) * t425;
t441 = t469 * t490 + qJD(1);
t475 = qJD(5) * t423 + t441 * t379 + t476;
t474 = -t388 * t457 + t440 * t396 + t477;
t429 = (pkin(4) * t461 + qJ(5) * t460) * t465;
t473 = qJD(5) * t425 + t410 * t429 + t474;
t467 = cos(qJ(6));
t463 = sin(qJ(6));
t448 = rSges(2,1) * t470 - rSges(2,2) * t466;
t447 = rSges(2,1) * t466 + rSges(2,2) * t470;
t446 = rSges(3,1) * t465 + rSges(3,2) * t469;
t443 = Icges(3,5) * t465 + Icges(3,6) * t469;
t438 = pkin(5) * t461 * t465 + pkin(10) * t469;
t435 = t468 * t496 + t501;
t434 = -t464 * t496 + t466 * t468;
t433 = t468 * t497 - t500;
t432 = -t464 * t497 - t468 * t470;
t431 = qJD(1) + (qJD(6) + t490) * t469;
t422 = rSges(3,3) * t466 + t470 * t488;
t421 = -rSges(3,3) * t470 + t466 * t488;
t420 = -rSges(4,3) * t469 + (rSges(4,1) * t468 - rSges(4,2) * t464) * t465;
t417 = -Icges(4,5) * t469 + (Icges(4,1) * t468 - Icges(4,4) * t464) * t465;
t414 = -Icges(4,6) * t469 + (Icges(4,4) * t468 - Icges(4,2) * t464) * t465;
t413 = Icges(3,3) * t466 + t470 * t485;
t412 = -Icges(3,3) * t470 + t466 * t485;
t411 = -Icges(4,3) * t469 + (Icges(4,5) * t468 - Icges(4,6) * t464) * t465;
t406 = (t460 * t463 + t461 * t467) * t465;
t405 = (t460 * t467 - t461 * t463) * t465;
t404 = -rSges(5,3) * t469 + (rSges(5,1) * t461 - rSges(5,2) * t460) * t465;
t403 = -rSges(6,2) * t469 + (rSges(6,1) * t461 + rSges(6,3) * t460) * t465;
t395 = -t466 * t491 + t410;
t394 = -t470 * t491 + t409;
t393 = pkin(5) * t426 - pkin(10) * t498;
t392 = pkin(5) * t424 - pkin(10) * t499;
t387 = rSges(4,1) * t435 + rSges(4,2) * t434 + rSges(4,3) * t498;
t386 = rSges(4,1) * t433 + rSges(4,2) * t432 + rSges(4,3) * t499;
t385 = Icges(4,1) * t435 + Icges(4,4) * t434 + Icges(4,5) * t498;
t384 = Icges(4,1) * t433 + Icges(4,4) * t432 + Icges(4,5) * t499;
t383 = Icges(4,4) * t435 + Icges(4,2) * t434 + Icges(4,6) * t498;
t382 = Icges(4,4) * t433 + Icges(4,2) * t432 + Icges(4,6) * t499;
t381 = Icges(4,5) * t435 + Icges(4,6) * t434 + Icges(4,3) * t498;
t380 = Icges(4,5) * t433 + Icges(4,6) * t432 + Icges(4,3) * t499;
t377 = t425 * t463 + t426 * t467;
t376 = t425 * t467 - t426 * t463;
t375 = t423 * t463 + t424 * t467;
t374 = t423 * t467 - t424 * t463;
t373 = qJD(1) * t422 - t446 * t459 + t442;
t372 = -t446 * t494 + (-t421 - t451) * qJD(1);
t371 = (t421 * t466 + t422 * t470) * qJD(2);
t369 = rSges(5,1) * t426 - rSges(5,2) * t425 + rSges(5,3) * t498;
t368 = rSges(6,1) * t426 + rSges(6,2) * t498 + rSges(6,3) * t425;
t367 = rSges(5,1) * t424 - rSges(5,2) * t423 + rSges(5,3) * t499;
t366 = rSges(6,1) * t424 + rSges(6,2) * t499 + rSges(6,3) * t423;
t353 = rSges(7,1) * t406 + rSges(7,2) * t405 + rSges(7,3) * t469;
t352 = Icges(7,1) * t406 + Icges(7,4) * t405 + Icges(7,5) * t469;
t351 = Icges(7,4) * t406 + Icges(7,2) * t405 + Icges(7,6) * t469;
t350 = Icges(7,5) * t406 + Icges(7,6) * t405 + Icges(7,3) * t469;
t346 = rSges(7,1) * t377 + rSges(7,2) * t376 - rSges(7,3) * t498;
t345 = rSges(7,1) * t375 + rSges(7,2) * t374 - rSges(7,3) * t499;
t344 = Icges(7,1) * t377 + Icges(7,4) * t376 - Icges(7,5) * t498;
t343 = Icges(7,1) * t375 + Icges(7,4) * t374 - Icges(7,5) * t499;
t342 = Icges(7,4) * t377 + Icges(7,2) * t376 - Icges(7,6) * t498;
t341 = Icges(7,4) * t375 + Icges(7,2) * t374 - Icges(7,6) * t499;
t340 = Icges(7,5) * t377 + Icges(7,6) * t376 - Icges(7,3) * t498;
t339 = Icges(7,5) * t375 + Icges(7,6) * t374 - Icges(7,3) * t499;
t338 = t387 * t457 - t420 * t439 + t480;
t337 = -t386 * t457 + t420 * t440 + t477;
t336 = t386 * t439 - t387 * t440 + t495;
t335 = t369 * t441 - t404 * t409 + t476;
t334 = -t367 * t441 + t404 * t410 + t474;
t333 = t367 * t409 - t369 * t410 + t481;
t332 = t368 * t441 + (-t403 - t429) * t409 + t475;
t331 = t403 * t410 + (-t366 - t378) * t441 + t473;
t330 = t366 * t409 + (-t368 - t379) * t410 + t479;
t329 = t346 * t431 - t353 * t394 + t393 * t441 + (-t429 - t438) * t409 + t475;
t328 = -t345 * t431 + t353 * t395 + t410 * t438 + (-t378 - t392) * t441 + t473;
t327 = t345 * t394 - t346 * t395 + t392 * t409 + (-t379 - t393) * t410 + t479;
t1 = qJD(1) * ((t469 * t444 + t465 * t445) * qJD(1) + ((t416 * t469 + t419 * t465) * t466 - (t415 * t469 + t465 * t418) * t470) * qJD(2)) / 0.2e1 + t439 * ((t381 * t498 + t434 * t383 + t435 * t385) * t439 + (t380 * t498 + t382 * t434 + t384 * t435) * t440 + (t411 * t498 + t414 * t434 + t417 * t435) * t457) / 0.2e1 + t440 * ((t381 * t499 + t383 * t432 + t385 * t433) * t439 + (t380 * t499 + t432 * t382 + t433 * t384) * t440 + (t411 * t499 + t414 * t432 + t417 * t433) * t457) / 0.2e1 + t457 * ((-t380 * t440 - t381 * t439 - t411 * t457) * t469 + ((-t383 * t464 + t385 * t468) * t439 + (-t382 * t464 + t384 * t468) * t440 + (-t414 * t464 + t417 * t468) * t457) * t465) / 0.2e1 + m(7) * (t327 ^ 2 + t328 ^ 2 + t329 ^ 2) / 0.2e1 + m(5) * (t333 ^ 2 + t334 ^ 2 + t335 ^ 2) / 0.2e1 + m(6) * (t330 ^ 2 + t331 ^ 2 + t332 ^ 2) / 0.2e1 + m(3) * (t371 ^ 2 + t372 ^ 2 + t373 ^ 2) / 0.2e1 + m(4) * (t336 ^ 2 + t337 ^ 2 + t338 ^ 2) / 0.2e1 + t394 * ((-t340 * t498 + t376 * t342 + t377 * t344) * t394 + (-t339 * t498 + t341 * t376 + t343 * t377) * t395 + (-t350 * t498 + t351 * t376 + t352 * t377) * t431) / 0.2e1 + t395 * ((-t340 * t499 + t342 * t374 + t344 * t375) * t394 + (-t339 * t499 + t374 * t341 + t375 * t343) * t395 + (-t350 * t499 + t351 * t374 + t352 * t375) * t431) / 0.2e1 + t431 * ((t340 * t469 + t342 * t405 + t344 * t406) * t394 + (t339 * t469 + t341 * t405 + t343 * t406) * t395 + (t469 * t350 + t405 * t351 + t406 * t352) * t431) / 0.2e1 - ((-t470 * t443 + t466 * t482) * qJD(1) + (t470 ^ 2 * t412 + (t483 * t466 + (-t413 + t484) * t470) * t466) * qJD(2)) * t494 / 0.2e1 + ((t466 * t443 + t470 * t482) * qJD(1) + (t466 ^ 2 * t413 + (t484 * t470 + (-t412 + t483) * t466) * t470) * qJD(2)) * t459 / 0.2e1 + ((t425 * t512 + t426 * t510 + t498 * t511) * t441 + (t425 * t518 + t514 * t426 + t516 * t498) * t410 + (t425 * t517 + t426 * t513 + t498 * t515) * t409) * t409 / 0.2e1 + ((t423 * t512 + t424 * t510 + t499 * t511) * t441 + (t423 * t518 + t514 * t424 + t516 * t499) * t410 + (t423 * t517 + t424 * t513 + t499 * t515) * t409) * t410 / 0.2e1 + ((-t409 * t515 - t410 * t516 - t441 * t511) * t469 + ((t460 * t512 + t461 * t510) * t441 + (t460 * t518 + t514 * t461) * t410 + (t460 * t517 + t461 * t513) * t409) * t465) * t441 / 0.2e1 + (m(2) * (t447 ^ 2 + t448 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
