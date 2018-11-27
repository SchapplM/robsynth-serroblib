% Calculate kinetic energy for
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% m_mdh [8x1]
%   mass of all robot links (including the base)
% rSges [8x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [8x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S7RRRRRRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(4,1),zeros(8,1),zeros(8,3),zeros(8,6)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_energykin_fixb_slag_vp1: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_energykin_fixb_slag_vp1: qJD has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_energykin_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [8 1]), ...
  'S7RRRRRRR1_energykin_fixb_slag_vp1: m has to be [8x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [8,3]), ...
  'S7RRRRRRR1_energykin_fixb_slag_vp1: rSges has to be [8x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [8 6]), ...
  'S7RRRRRRR1_energykin_fixb_slag_vp1: Icges has to be [8x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-26 19:14:21
% EndTime: 2018-11-26 19:14:25
% DurationCPUTime: 3.98s
% Computational Cost: add. (3779->405), mult. (9507->652), div. (0->0), fcn. (12092->14), ass. (0->179)
t492 = sin(qJ(1));
t484 = t492 ^ 2;
t497 = cos(qJ(1));
t485 = t497 ^ 2;
t526 = cos(qJ(4));
t525 = cos(qJ(5));
t489 = sin(qJ(4));
t491 = sin(qJ(2));
t495 = cos(qJ(3));
t496 = cos(qJ(2));
t463 = t491 * t495 * t489 + t496 * t526;
t524 = pkin(3) * t463;
t481 = qJD(3) * t496 + qJD(1);
t490 = sin(qJ(3));
t520 = t490 * t491;
t469 = -qJD(4) * t520 + t481;
t523 = pkin(3) * t469;
t522 = Icges(3,4) * t491;
t521 = Icges(3,4) * t496;
t519 = t491 * t492;
t518 = t491 * t497;
t517 = t496 * t492;
t516 = t496 * t497;
t515 = qJD(1) * t491;
t483 = qJD(2) * t492;
t514 = qJD(2) * t497;
t513 = qJD(3) * t491;
t512 = t491 * t526;
t470 = -t497 * t513 + t483;
t511 = rSges(3,1) * t496 - rSges(3,2) * t491;
t510 = Icges(3,1) * t496 - t522;
t509 = -Icges(3,2) * t491 + t521;
t508 = Icges(3,5) * t496 - Icges(3,6) * t491;
t450 = -Icges(3,6) * t497 + t492 * t509;
t453 = -Icges(3,5) * t497 + t492 * t510;
t507 = t450 * t491 - t453 * t496;
t451 = Icges(3,6) * t492 + t497 * t509;
t454 = Icges(3,5) * t492 + t497 * t510;
t506 = -t451 * t491 + t454 * t496;
t473 = Icges(3,2) * t496 + t522;
t474 = Icges(3,1) * t491 + t521;
t505 = -t473 * t491 + t474 * t496;
t504 = (-t484 - t485) * t491 * qJD(2) * pkin(2);
t503 = (t492 * t515 - t496 * t514) * pkin(2);
t467 = -t490 * t516 - t492 * t495;
t435 = qJD(4) * t467 + t470;
t471 = -t492 * t513 - t514;
t468 = -t490 * t492 + t495 * t516;
t444 = t468 * t489 - t497 * t512;
t406 = qJD(5) * t444 + t435;
t434 = qJD(5) * t463 + t469;
t465 = -t490 * t517 + t495 * t497;
t436 = qJD(4) * t465 + t471;
t445 = t468 * t526 + t489 * t518;
t488 = sin(qJ(5));
t415 = t445 * t488 - t467 * t525;
t383 = qJD(6) * t415 + t406;
t464 = -t496 * t489 + t495 * t512;
t440 = t464 * t488 + t520 * t525;
t405 = qJD(6) * t440 + t434;
t466 = t490 * t497 + t495 * t517;
t442 = t466 * t489 - t492 * t512;
t407 = qJD(5) * t442 + t436;
t502 = (-t483 * t496 - t497 * t515) * pkin(2);
t443 = t466 * t526 + t489 * t519;
t413 = t443 * t488 - t465 * t525;
t384 = qJD(6) * t413 + t407;
t501 = t436 * t524 - t442 * t523 + t503;
t500 = t504 + (t435 * t442 - t436 * t444) * pkin(3);
t499 = -t435 * t524 + t444 * t523 + t502;
t494 = cos(qJ(6));
t493 = cos(qJ(7));
t487 = sin(qJ(6));
t486 = sin(qJ(7));
t477 = rSges(2,1) * t497 - rSges(2,2) * t492;
t476 = rSges(2,1) * t492 + rSges(2,2) * t497;
t475 = rSges(3,1) * t491 + rSges(3,2) * t496;
t472 = Icges(3,5) * t491 + Icges(3,6) * t496;
t457 = rSges(3,3) * t492 + t497 * t511;
t456 = -rSges(3,3) * t497 + t492 * t511;
t455 = rSges(4,3) * t496 + (rSges(4,1) * t495 - rSges(4,2) * t490) * t491;
t452 = Icges(4,5) * t496 + (Icges(4,1) * t495 - Icges(4,4) * t490) * t491;
t449 = Icges(4,6) * t496 + (Icges(4,4) * t495 - Icges(4,2) * t490) * t491;
t448 = Icges(3,3) * t492 + t497 * t508;
t447 = -Icges(3,3) * t497 + t492 * t508;
t446 = Icges(4,3) * t496 + (Icges(4,5) * t495 - Icges(4,6) * t490) * t491;
t441 = t464 * t525 - t488 * t520;
t433 = -qJD(1) * t456 - t475 * t514;
t432 = qJD(1) * t457 - t475 * t483;
t431 = rSges(4,1) * t468 + rSges(4,2) * t467 - rSges(4,3) * t518;
t430 = rSges(4,1) * t466 + rSges(4,2) * t465 - rSges(4,3) * t519;
t429 = rSges(5,1) * t464 - rSges(5,2) * t463 - rSges(5,3) * t520;
t428 = Icges(4,1) * t468 + Icges(4,4) * t467 - Icges(4,5) * t518;
t427 = Icges(4,1) * t466 + Icges(4,4) * t465 - Icges(4,5) * t519;
t426 = Icges(5,1) * t464 - Icges(5,4) * t463 - Icges(5,5) * t520;
t425 = Icges(4,4) * t468 + Icges(4,2) * t467 - Icges(4,6) * t518;
t424 = Icges(4,4) * t466 + Icges(4,2) * t465 - Icges(4,6) * t519;
t423 = Icges(5,4) * t464 - Icges(5,2) * t463 - Icges(5,6) * t520;
t422 = Icges(4,5) * t468 + Icges(4,6) * t467 - Icges(4,3) * t518;
t421 = Icges(4,5) * t466 + Icges(4,6) * t465 - Icges(4,3) * t519;
t420 = Icges(5,5) * t464 - Icges(5,6) * t463 - Icges(5,3) * t520;
t417 = (t456 * t492 + t457 * t497) * qJD(2);
t416 = t445 * t525 + t467 * t488;
t414 = t443 * t525 + t465 * t488;
t412 = t441 * t494 + t463 * t487;
t411 = -t441 * t487 + t463 * t494;
t404 = rSges(5,1) * t445 - rSges(5,2) * t444 + rSges(5,3) * t467;
t403 = rSges(5,1) * t443 - rSges(5,2) * t442 + rSges(5,3) * t465;
t402 = rSges(6,1) * t441 - rSges(6,2) * t440 + rSges(6,3) * t463;
t401 = Icges(5,1) * t445 - Icges(5,4) * t444 + Icges(5,5) * t467;
t400 = Icges(5,1) * t443 - Icges(5,4) * t442 + Icges(5,5) * t465;
t399 = Icges(6,1) * t441 - Icges(6,4) * t440 + Icges(6,5) * t463;
t398 = Icges(5,4) * t445 - Icges(5,2) * t444 + Icges(5,6) * t467;
t397 = Icges(5,4) * t443 - Icges(5,2) * t442 + Icges(5,6) * t465;
t396 = Icges(6,4) * t441 - Icges(6,2) * t440 + Icges(6,6) * t463;
t395 = Icges(5,5) * t445 - Icges(5,6) * t444 + Icges(5,3) * t467;
t394 = Icges(5,5) * t443 - Icges(5,6) * t442 + Icges(5,3) * t465;
t393 = Icges(6,5) * t441 - Icges(6,6) * t440 + Icges(6,3) * t463;
t392 = t416 * t494 + t444 * t487;
t391 = -t416 * t487 + t444 * t494;
t390 = t414 * t494 + t442 * t487;
t389 = -t414 * t487 + t442 * t494;
t388 = t412 * t493 - t440 * t486;
t387 = -t412 * t486 - t440 * t493;
t386 = -t430 * t481 + t455 * t471 + t503;
t385 = t431 * t481 - t455 * t470 + t502;
t382 = qJD(7) * t411 + t405;
t381 = t430 * t470 - t431 * t471 + t504;
t380 = rSges(6,1) * t416 - rSges(6,2) * t415 + rSges(6,3) * t444;
t379 = rSges(6,1) * t414 - rSges(6,2) * t413 + rSges(6,3) * t442;
t378 = rSges(7,1) * t412 + rSges(7,2) * t411 + rSges(7,3) * t440;
t377 = Icges(6,1) * t416 - Icges(6,4) * t415 + Icges(6,5) * t444;
t376 = Icges(6,1) * t414 - Icges(6,4) * t413 + Icges(6,5) * t442;
t375 = Icges(7,1) * t412 + Icges(7,4) * t411 + Icges(7,5) * t440;
t374 = Icges(6,4) * t416 - Icges(6,2) * t415 + Icges(6,6) * t444;
t373 = Icges(6,4) * t414 - Icges(6,2) * t413 + Icges(6,6) * t442;
t372 = Icges(7,4) * t412 + Icges(7,2) * t411 + Icges(7,6) * t440;
t371 = Icges(6,5) * t416 - Icges(6,6) * t415 + Icges(6,3) * t444;
t370 = Icges(6,5) * t414 - Icges(6,6) * t413 + Icges(6,3) * t442;
t369 = Icges(7,5) * t412 + Icges(7,6) * t411 + Icges(7,3) * t440;
t368 = t392 * t493 - t415 * t486;
t367 = -t392 * t486 - t415 * t493;
t366 = t390 * t493 - t413 * t486;
t365 = -t390 * t486 - t413 * t493;
t364 = -t403 * t469 + t429 * t436 + t503;
t363 = t404 * t469 - t429 * t435 + t502;
t362 = qJD(7) * t389 + t384;
t361 = qJD(7) * t391 + t383;
t360 = rSges(7,1) * t392 + rSges(7,2) * t391 + rSges(7,3) * t415;
t359 = rSges(7,1) * t390 + rSges(7,2) * t389 + rSges(7,3) * t413;
t358 = rSges(8,1) * t388 + rSges(8,2) * t387 + rSges(8,3) * t411;
t357 = Icges(7,1) * t392 + Icges(7,4) * t391 + Icges(7,5) * t415;
t356 = Icges(7,1) * t390 + Icges(7,4) * t389 + Icges(7,5) * t413;
t355 = Icges(8,1) * t388 + Icges(8,4) * t387 + Icges(8,5) * t411;
t354 = Icges(7,4) * t392 + Icges(7,2) * t391 + Icges(7,6) * t415;
t353 = Icges(7,4) * t390 + Icges(7,2) * t389 + Icges(7,6) * t413;
t352 = Icges(8,4) * t388 + Icges(8,2) * t387 + Icges(8,6) * t411;
t351 = Icges(7,5) * t392 + Icges(7,6) * t391 + Icges(7,3) * t415;
t350 = Icges(7,5) * t390 + Icges(7,6) * t389 + Icges(7,3) * t413;
t349 = Icges(8,5) * t388 + Icges(8,6) * t387 + Icges(8,3) * t411;
t348 = t403 * t435 - t404 * t436 + t504;
t347 = rSges(8,1) * t368 + rSges(8,2) * t367 + rSges(8,3) * t391;
t346 = rSges(8,1) * t366 + rSges(8,2) * t365 + rSges(8,3) * t389;
t345 = Icges(8,1) * t368 + Icges(8,4) * t367 + Icges(8,5) * t391;
t344 = Icges(8,1) * t366 + Icges(8,4) * t365 + Icges(8,5) * t389;
t343 = Icges(8,4) * t368 + Icges(8,2) * t367 + Icges(8,6) * t391;
t342 = Icges(8,4) * t366 + Icges(8,2) * t365 + Icges(8,6) * t389;
t341 = Icges(8,5) * t368 + Icges(8,6) * t367 + Icges(8,3) * t391;
t340 = Icges(8,5) * t366 + Icges(8,6) * t365 + Icges(8,3) * t389;
t339 = -t379 * t434 + t402 * t407 + t501;
t338 = t380 * t434 - t402 * t406 + t499;
t337 = t379 * t406 - t380 * t407 + t500;
t336 = -t359 * t405 + t378 * t384 + t501;
t335 = t360 * t405 - t378 * t383 + t499;
t334 = t359 * t383 - t360 * t384 + t500;
t333 = -t346 * t382 + t358 * t362 + (t384 * t411 - t389 * t405) * pkin(4) + t501;
t332 = t347 * t382 - t358 * t361 + (-t383 * t411 + t391 * t405) * pkin(4) + t499;
t331 = t346 * t361 - t347 * t362 + (t383 * t389 - t384 * t391) * pkin(4) + t500;
t1 = t383 * ((t351 * t415 + t354 * t391 + t357 * t392) * t383 + (t350 * t415 + t353 * t391 + t356 * t392) * t384 + (t369 * t415 + t372 * t391 + t375 * t392) * t405) / 0.2e1 + t382 * ((t341 * t411 + t343 * t387 + t345 * t388) * t361 + (t340 * t411 + t342 * t387 + t344 * t388) * t362 + (t349 * t411 + t352 * t387 + t355 * t388) * t382) / 0.2e1 + t384 * ((t351 * t413 + t354 * t389 + t357 * t390) * t383 + (t350 * t413 + t353 * t389 + t356 * t390) * t384 + (t369 * t413 + t372 * t389 + t375 * t390) * t405) / 0.2e1 + t434 * ((t371 * t463 - t374 * t440 + t377 * t441) * t406 + (t370 * t463 - t373 * t440 + t376 * t441) * t407 + (t393 * t463 - t396 * t440 + t399 * t441) * t434) / 0.2e1 + t407 * ((t371 * t442 - t374 * t413 + t377 * t414) * t406 + (t370 * t442 - t373 * t413 + t376 * t414) * t407 + (t393 * t442 - t396 * t413 + t399 * t414) * t434) / 0.2e1 + t406 * ((t371 * t444 - t374 * t415 + t377 * t416) * t406 + (t370 * t444 - t373 * t415 + t376 * t416) * t407 + (t393 * t444 - t396 * t415 + t399 * t416) * t434) / 0.2e1 + t405 * ((t351 * t440 + t354 * t411 + t357 * t412) * t383 + (t350 * t440 + t353 * t411 + t356 * t412) * t384 + (t369 * t440 + t372 * t411 + t375 * t412) * t405) / 0.2e1 + m(3) * (t417 ^ 2 + t432 ^ 2 + t433 ^ 2) / 0.2e1 + t471 * ((-t422 * t519 + t425 * t465 + t428 * t466) * t470 + (-t421 * t519 + t424 * t465 + t427 * t466) * t471 + (-t446 * t519 + t449 * t465 + t452 * t466) * t481) / 0.2e1 + t469 * ((-t395 * t520 - t398 * t463 + t401 * t464) * t435 + (-t394 * t520 - t397 * t463 + t400 * t464) * t436 + (-t420 * t520 - t423 * t463 + t426 * t464) * t469) / 0.2e1 - ((-t497 * t472 + t492 * t505) * qJD(1) + (t485 * t447 + (t506 * t492 + (-t448 + t507) * t497) * t492) * qJD(2)) * t514 / 0.2e1 + ((t492 * t472 + t497 * t505) * qJD(1) + (t484 * t448 + (t507 * t497 + (-t447 + t506) * t492) * t497) * qJD(2)) * t483 / 0.2e1 + m(5) * (t348 ^ 2 + t363 ^ 2 + t364 ^ 2) / 0.2e1 + t470 * ((-t422 * t518 + t425 * t467 + t428 * t468) * t470 + (-t421 * t518 + t424 * t467 + t427 * t468) * t471 + (-t446 * t518 + t449 * t467 + t452 * t468) * t481) / 0.2e1 + m(4) * (t381 ^ 2 + t385 ^ 2 + t386 ^ 2) / 0.2e1 + t362 * ((t341 * t389 + t343 * t365 + t345 * t366) * t361 + (t340 * t389 + t342 * t365 + t344 * t366) * t362 + (t349 * t389 + t352 * t365 + t355 * t366) * t382) / 0.2e1 + t361 * ((t341 * t391 + t343 * t367 + t345 * t368) * t361 + (t340 * t391 + t342 * t367 + t344 * t368) * t362 + (t349 * t391 + t352 * t367 + t355 * t368) * t382) / 0.2e1 + m(8) * (t331 ^ 2 + t332 ^ 2 + t333 ^ 2) / 0.2e1 + t436 * ((t395 * t465 - t398 * t442 + t401 * t443) * t435 + (t394 * t465 - t397 * t442 + t400 * t443) * t436 + (t420 * t465 - t423 * t442 + t426 * t443) * t469) / 0.2e1 + t435 * ((t395 * t467 - t398 * t444 + t401 * t445) * t435 + (t394 * t467 - t397 * t444 + t400 * t445) * t436 + (t420 * t467 - t423 * t444 + t426 * t445) * t469) / 0.2e1 + qJD(1) * ((t473 * t496 + t474 * t491) * qJD(1) + ((t451 * t496 + t454 * t491) * t492 - (t450 * t496 + t453 * t491) * t497) * qJD(2)) / 0.2e1 + m(6) * (t337 ^ 2 + t338 ^ 2 + t339 ^ 2) / 0.2e1 + t481 * ((t421 * t471 + t422 * t470 + t446 * t481) * t496 + ((-t425 * t490 + t428 * t495) * t470 + (-t424 * t490 + t427 * t495) * t471 + (-t449 * t490 + t452 * t495) * t481) * t491) / 0.2e1 + (m(2) * (t476 ^ 2 + t477 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + m(7) * (t334 ^ 2 + t335 ^ 2 + t336 ^ 2) / 0.2e1;
T  = t1;
