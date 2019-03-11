% Calculate kinetic energy for
% S6RRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:31:44
% EndTime: 2019-03-09 15:31:46
% DurationCPUTime: 2.90s
% Computational Cost: add. (1853->312), mult. (2852->466), div. (0->0), fcn. (3029->10), ass. (0->152)
t522 = Icges(5,1) + Icges(6,1);
t521 = -Icges(5,4) + Icges(6,5);
t520 = Icges(6,4) + Icges(5,5);
t519 = Icges(5,2) + Icges(6,3);
t518 = -Icges(6,6) + Icges(5,6);
t517 = -Icges(5,3) - Icges(6,2) - Icges(4,3);
t458 = qJ(3) + pkin(10);
t455 = sin(t458);
t456 = cos(t458);
t467 = cos(qJ(1));
t463 = sin(qJ(1));
t466 = cos(qJ(2));
t495 = t463 * t466;
t408 = t455 * t495 + t456 * t467;
t409 = -t455 * t467 + t456 * t495;
t462 = sin(qJ(2));
t497 = t462 * t463;
t516 = t519 * t408 + t521 * t409 - t518 * t497;
t494 = t466 * t467;
t410 = t455 * t494 - t463 * t456;
t411 = t455 * t463 + t456 * t494;
t496 = t462 * t467;
t515 = t519 * t410 + t521 * t411 - t518 * t496;
t514 = t521 * t408 + t522 * t409 + t520 * t497;
t513 = t521 * t410 + t522 * t411 + t520 * t496;
t512 = t518 * t466 + (t519 * t455 + t521 * t456) * t462;
t511 = -t520 * t466 + (t521 * t455 + t522 * t456) * t462;
t461 = sin(qJ(3));
t465 = cos(qJ(3));
t428 = -t461 * t495 - t465 * t467;
t498 = t461 * t467;
t429 = t465 * t495 - t498;
t510 = Icges(4,5) * t429 + Icges(4,6) * t428 - t518 * t408 + t520 * t409 - t517 * t497;
t430 = -t461 * t494 + t463 * t465;
t499 = t461 * t463;
t431 = t465 * t494 + t499;
t509 = Icges(4,5) * t431 + Icges(4,6) * t430 - t518 * t410 + t520 * t411 - t517 * t496;
t508 = t517 * t466 + (Icges(4,5) * t465 - Icges(4,6) * t461 - t518 * t455 + t520 * t456) * t462;
t503 = pkin(3) * t465;
t501 = Icges(3,4) * t462;
t500 = Icges(3,4) * t466;
t377 = pkin(4) * t409 + qJ(5) * t408;
t475 = qJ(4) * t462 + t466 * t503;
t385 = -pkin(3) * t498 + t463 * t475;
t493 = -t377 - t385;
t378 = pkin(4) * t411 + qJ(5) * t410;
t386 = pkin(3) * t499 + t467 * t475;
t492 = -t378 - t386;
t393 = -qJ(4) * t466 + t462 * t503;
t421 = (pkin(4) * t456 + qJ(5) * t455) * t462;
t491 = -t393 - t421;
t485 = pkin(2) * t466 + pkin(8) * t462;
t432 = t485 * t463;
t433 = t485 * t467;
t457 = qJD(2) * t463;
t489 = qJD(2) * t467;
t490 = t432 * t457 + t433 * t489;
t488 = qJD(3) * t462;
t435 = t467 * t488 + t457;
t487 = qJD(4) * t462;
t486 = qJD(6) * t462;
t436 = t463 * t488 - t489;
t484 = rSges(3,1) * t466 - rSges(3,2) * t462;
t483 = Icges(3,1) * t466 - t501;
t482 = -Icges(3,2) * t462 + t500;
t481 = Icges(3,5) * t466 - Icges(3,6) * t462;
t416 = -Icges(3,6) * t467 + t463 * t482;
t419 = -Icges(3,5) * t467 + t463 * t483;
t480 = t416 * t462 - t419 * t466;
t417 = Icges(3,6) * t463 + t467 * t482;
t420 = Icges(3,5) * t463 + t467 * t483;
t479 = -t417 * t462 + t420 * t466;
t440 = Icges(3,2) * t466 + t501;
t441 = Icges(3,1) * t462 + t500;
t478 = -t440 * t462 + t441 * t466;
t438 = qJD(1) * (pkin(1) * t467 + pkin(7) * t463);
t446 = pkin(2) * t462 - pkin(8) * t466;
t477 = qJD(1) * t433 - t446 * t457 + t438;
t476 = -qJD(4) * t466 + t435 * t385 + t490;
t453 = -qJD(3) * t466 + qJD(1);
t474 = t453 * t386 + t463 * t487 + t477;
t473 = qJD(5) * t462 * t455 + t435 * t377 + t476;
t447 = pkin(1) * t463 - pkin(7) * t467;
t472 = (-t432 - t447) * qJD(1) - t446 * t489;
t471 = qJD(5) * t408 + t453 * t378 + t474;
t470 = t436 * t393 + t467 * t487 + t472;
t469 = qJD(5) * t410 + t436 * t421 + t470;
t464 = cos(qJ(6));
t460 = sin(qJ(6));
t445 = rSges(2,1) * t467 - rSges(2,2) * t463;
t444 = rSges(2,1) * t463 + rSges(2,2) * t467;
t443 = rSges(3,1) * t462 + rSges(3,2) * t466;
t439 = Icges(3,5) * t462 + Icges(3,6) * t466;
t437 = qJD(1) + (-qJD(3) + qJD(6)) * t466;
t434 = pkin(5) * t456 * t462 + pkin(9) * t466;
t424 = rSges(3,3) * t463 + t467 * t484;
t423 = -rSges(3,3) * t467 + t463 * t484;
t422 = -rSges(4,3) * t466 + (rSges(4,1) * t465 - rSges(4,2) * t461) * t462;
t418 = -Icges(4,5) * t466 + (Icges(4,1) * t465 - Icges(4,4) * t461) * t462;
t415 = -Icges(4,6) * t466 + (Icges(4,4) * t465 - Icges(4,2) * t461) * t462;
t414 = Icges(3,3) * t463 + t467 * t481;
t413 = -Icges(3,3) * t467 + t463 * t481;
t407 = -t463 * t486 + t436;
t406 = -t467 * t486 + t435;
t403 = (t455 * t460 + t456 * t464) * t462;
t402 = (t455 * t464 - t456 * t460) * t462;
t401 = -rSges(5,3) * t466 + (rSges(5,1) * t456 - rSges(5,2) * t455) * t462;
t400 = -rSges(6,2) * t466 + (rSges(6,1) * t456 + rSges(6,3) * t455) * t462;
t392 = pkin(5) * t411 - pkin(9) * t496;
t391 = pkin(5) * t409 - pkin(9) * t497;
t388 = rSges(4,1) * t431 + rSges(4,2) * t430 + rSges(4,3) * t496;
t387 = rSges(4,1) * t429 + rSges(4,2) * t428 + rSges(4,3) * t497;
t384 = Icges(4,1) * t431 + Icges(4,4) * t430 + Icges(4,5) * t496;
t383 = Icges(4,1) * t429 + Icges(4,4) * t428 + Icges(4,5) * t497;
t382 = Icges(4,4) * t431 + Icges(4,2) * t430 + Icges(4,6) * t496;
t381 = Icges(4,4) * t429 + Icges(4,2) * t428 + Icges(4,6) * t497;
t376 = t410 * t460 + t411 * t464;
t375 = t410 * t464 - t411 * t460;
t374 = t408 * t460 + t409 * t464;
t373 = t408 * t464 - t409 * t460;
t372 = qJD(1) * t424 - t443 * t457 + t438;
t371 = -t443 * t489 + (-t423 - t447) * qJD(1);
t370 = (t423 * t463 + t424 * t467) * qJD(2);
t368 = rSges(5,1) * t411 - rSges(5,2) * t410 + rSges(5,3) * t496;
t367 = rSges(6,1) * t411 + rSges(6,2) * t496 + rSges(6,3) * t410;
t366 = rSges(5,1) * t409 - rSges(5,2) * t408 + rSges(5,3) * t497;
t365 = rSges(6,1) * t409 + rSges(6,2) * t497 + rSges(6,3) * t408;
t351 = rSges(7,1) * t403 + rSges(7,2) * t402 + rSges(7,3) * t466;
t350 = Icges(7,1) * t403 + Icges(7,4) * t402 + Icges(7,5) * t466;
t349 = Icges(7,4) * t403 + Icges(7,2) * t402 + Icges(7,6) * t466;
t348 = Icges(7,5) * t403 + Icges(7,6) * t402 + Icges(7,3) * t466;
t345 = rSges(7,1) * t376 + rSges(7,2) * t375 - rSges(7,3) * t496;
t344 = rSges(7,1) * t374 + rSges(7,2) * t373 - rSges(7,3) * t497;
t343 = Icges(7,1) * t376 + Icges(7,4) * t375 - Icges(7,5) * t496;
t342 = Icges(7,1) * t374 + Icges(7,4) * t373 - Icges(7,5) * t497;
t341 = Icges(7,4) * t376 + Icges(7,2) * t375 - Icges(7,6) * t496;
t340 = Icges(7,4) * t374 + Icges(7,2) * t373 - Icges(7,6) * t497;
t339 = Icges(7,5) * t376 + Icges(7,6) * t375 - Icges(7,3) * t496;
t338 = Icges(7,5) * t374 + Icges(7,6) * t373 - Icges(7,3) * t497;
t337 = t388 * t453 - t422 * t435 + t477;
t336 = -t387 * t453 + t422 * t436 + t472;
t335 = t387 * t435 - t388 * t436 + t490;
t334 = t368 * t453 + (-t393 - t401) * t435 + t474;
t333 = t401 * t436 + (-t366 - t385) * t453 + t470;
t332 = t366 * t435 + (-t368 - t386) * t436 + t476;
t331 = t367 * t453 + (-t400 + t491) * t435 + t471;
t330 = t400 * t436 + (-t365 + t493) * t453 + t469;
t329 = t365 * t435 + (-t367 + t492) * t436 + t473;
t328 = t345 * t437 - t351 * t406 + t392 * t453 + (-t434 + t491) * t435 + t471;
t327 = -t344 * t437 + t351 * t407 + t434 * t436 + (-t391 + t493) * t453 + t469;
t326 = t344 * t406 - t345 * t407 + t391 * t435 + (-t392 + t492) * t436 + t473;
t1 = t437 * ((t339 * t466 + t341 * t402 + t343 * t403) * t406 + (t338 * t466 + t340 * t402 + t342 * t403) * t407 + (t466 * t348 + t402 * t349 + t403 * t350) * t437) / 0.2e1 + qJD(1) * ((t466 * t440 + t462 * t441) * qJD(1) + ((t417 * t466 + t420 * t462) * t463 - (t416 * t466 + t419 * t462) * t467) * qJD(2)) / 0.2e1 + m(6) * (t329 ^ 2 + t330 ^ 2 + t331 ^ 2) / 0.2e1 + m(7) * (t326 ^ 2 + t327 ^ 2 + t328 ^ 2) / 0.2e1 + m(5) * (t332 ^ 2 + t333 ^ 2 + t334 ^ 2) / 0.2e1 + m(4) * (t335 ^ 2 + t336 ^ 2 + t337 ^ 2) / 0.2e1 + m(3) * (t370 ^ 2 + t371 ^ 2 + t372 ^ 2) / 0.2e1 + ((t463 * t439 + t467 * t478) * qJD(1) + (t463 ^ 2 * t414 + (t480 * t467 + (-t413 + t479) * t463) * t467) * qJD(2)) * t457 / 0.2e1 + t406 * ((-t339 * t496 + t375 * t341 + t376 * t343) * t406 + (-t338 * t496 + t340 * t375 + t342 * t376) * t407 + (-t348 * t496 + t349 * t375 + t350 * t376) * t437) / 0.2e1 + t407 * ((-t339 * t497 + t341 * t373 + t343 * t374) * t406 + (-t338 * t497 + t373 * t340 + t374 * t342) * t407 + (-t348 * t497 + t349 * t373 + t350 * t374) * t437) / 0.2e1 - ((-t467 * t439 + t463 * t478) * qJD(1) + (t467 ^ 2 * t413 + (t479 * t463 + (-t414 + t480) * t467) * t463) * qJD(2)) * t489 / 0.2e1 + (m(2) * (t444 ^ 2 + t445 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t512 * t410 + t511 * t411 + t415 * t430 + t418 * t431 + t508 * t496) * t453 + (t381 * t430 + t383 * t431 + t516 * t410 + t514 * t411 + t510 * t496) * t436 + (t430 * t382 + t431 * t384 + t515 * t410 + t513 * t411 + t509 * t496) * t435) * t435 / 0.2e1 + ((t512 * t408 + t511 * t409 + t415 * t428 + t418 * t429 + t508 * t497) * t453 + (t428 * t381 + t429 * t383 + t516 * t408 + t514 * t409 + t510 * t497) * t436 + (t382 * t428 + t384 * t429 + t515 * t408 + t513 * t409 + t509 * t497) * t435) * t436 / 0.2e1 + ((-t509 * t435 - t510 * t436 - t508 * t453) * t466 + ((-t415 * t461 + t418 * t465 + t512 * t455 + t511 * t456) * t453 + (-t381 * t461 + t383 * t465 + t516 * t455 + t514 * t456) * t436 + (-t382 * t461 + t384 * t465 + t515 * t455 + t513 * t456) * t435) * t462) * t453 / 0.2e1;
T  = t1;
