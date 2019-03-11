% Calculate kinetic energy for
% S6RPRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR8_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR8_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR8_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:57:43
% EndTime: 2019-03-09 03:57:45
% DurationCPUTime: 2.08s
% Computational Cost: add. (1279->272), mult. (1562->420), div. (0->0), fcn. (1496->10), ass. (0->141)
t490 = Icges(4,3) + Icges(5,3);
t422 = qJ(3) + pkin(10);
t413 = sin(t422);
t414 = cos(t422);
t426 = sin(qJ(3));
t429 = cos(qJ(3));
t489 = Icges(4,5) * t426 + Icges(5,5) * t413 + Icges(4,6) * t429 + Icges(5,6) * t414;
t427 = sin(qJ(1));
t430 = cos(qJ(1));
t488 = t490 * t427 - t489 * t430;
t487 = t489 * t427 + t490 * t430;
t486 = Icges(4,5) * t429 + Icges(5,5) * t414 - Icges(4,6) * t426 - Icges(5,6) * t413;
t469 = Icges(5,4) * t414;
t395 = -Icges(5,2) * t413 + t469;
t470 = Icges(5,4) * t413;
t396 = Icges(5,1) * t414 - t470;
t471 = Icges(4,4) * t429;
t401 = -Icges(4,2) * t426 + t471;
t472 = Icges(4,4) * t426;
t402 = Icges(4,1) * t429 - t472;
t485 = t395 * t414 + t396 * t413 + t401 * t429 + t402 * t426;
t444 = Icges(5,2) * t414 + t470;
t361 = Icges(5,6) * t427 - t430 * t444;
t446 = Icges(5,1) * t413 + t469;
t363 = Icges(5,5) * t427 - t430 * t446;
t445 = Icges(4,2) * t429 + t472;
t377 = Icges(4,6) * t427 - t430 * t445;
t447 = Icges(4,1) * t426 + t471;
t379 = Icges(4,5) * t427 - t430 * t447;
t484 = t361 * t414 + t363 * t413 + t377 * t429 + t379 * t426;
t360 = Icges(5,6) * t430 + t427 * t444;
t362 = Icges(5,5) * t430 + t427 * t446;
t376 = Icges(4,6) * t430 + t427 * t445;
t378 = Icges(4,5) * t430 + t427 * t447;
t483 = -t360 * t414 - t362 * t413 - t376 * t429 - t378 * t426;
t428 = cos(qJ(5));
t475 = pkin(5) * t428;
t482 = -pkin(9) * t414 + t413 * t475;
t478 = pkin(3) * t426;
t477 = pkin(3) * t429;
t468 = t414 * t427;
t467 = t414 * t430;
t423 = qJ(5) + qJ(6);
t420 = sin(t423);
t466 = t420 * t427;
t465 = t420 * t430;
t421 = cos(t423);
t464 = t421 * t427;
t463 = t421 * t430;
t425 = sin(qJ(5));
t462 = t425 * t427;
t461 = t425 * t430;
t460 = t427 * t428;
t459 = t428 * t430;
t399 = qJD(1) * (pkin(1) * t430 + qJ(2) * t427);
t458 = qJD(1) * t430 * pkin(7) + t399;
t417 = qJD(3) * t427;
t457 = qJD(5) * t414;
t392 = t430 * t457 + t417;
t418 = qJD(3) * t430;
t408 = qJD(5) * t413 + qJD(1);
t419 = qJD(2) * t427;
t456 = qJD(4) * t430 + t417 * t477 + t419;
t404 = pkin(1) * t427 - qJ(2) * t430;
t453 = -pkin(7) * t427 - t404;
t386 = qJ(4) * t430 + t427 * t478;
t452 = qJD(1) * t386 + qJD(4) * t427 + t458;
t385 = qJ(4) * t427 - t430 * t478;
t451 = -t385 + t453;
t450 = pkin(4) * t413 - pkin(8) * t414;
t449 = rSges(4,1) * t426 + rSges(4,2) * t429;
t448 = rSges(5,1) * t413 + rSges(5,2) * t414;
t367 = t385 * t418;
t382 = t450 * t427;
t383 = t450 * t430;
t435 = -t383 * t418 + t367 + (-t382 - t386) * t417;
t398 = pkin(4) * t414 + pkin(8) * t413;
t434 = t398 * t417 + (t383 + t451) * qJD(1) + t456;
t433 = qJD(1) * t382 + (-qJD(2) + (-t398 - t477) * qJD(3)) * t430 + t452;
t407 = rSges(2,1) * t430 - rSges(2,2) * t427;
t406 = rSges(4,1) * t429 - rSges(4,2) * t426;
t405 = rSges(2,1) * t427 + rSges(2,2) * t430;
t397 = rSges(5,1) * t414 - rSges(5,2) * t413;
t393 = -t427 * t457 + t418;
t391 = qJD(6) * t413 + t408;
t390 = -t413 * t459 + t462;
t389 = t413 * t461 + t460;
t388 = t413 * t460 + t461;
t387 = -t413 * t462 + t459;
t381 = rSges(4,3) * t427 - t430 * t449;
t380 = rSges(4,3) * t430 + t427 * t449;
t371 = -t413 * t463 + t466;
t370 = t413 * t465 + t464;
t369 = t413 * t464 + t465;
t368 = -t413 * t466 + t463;
t365 = rSges(5,3) * t427 - t430 * t448;
t364 = rSges(5,3) * t430 + t427 * t448;
t357 = t418 + (-qJD(5) - qJD(6)) * t468;
t356 = qJD(6) * t467 + t392;
t355 = rSges(6,3) * t413 + (rSges(6,1) * t428 - rSges(6,2) * t425) * t414;
t354 = Icges(6,5) * t413 + (Icges(6,1) * t428 - Icges(6,4) * t425) * t414;
t353 = Icges(6,6) * t413 + (Icges(6,4) * t428 - Icges(6,2) * t425) * t414;
t352 = Icges(6,3) * t413 + (Icges(6,5) * t428 - Icges(6,6) * t425) * t414;
t351 = t399 - qJD(2) * t430 + qJD(1) * (-rSges(3,2) * t430 + rSges(3,3) * t427);
t350 = t419 + (rSges(3,2) * t427 + rSges(3,3) * t430 - t404) * qJD(1);
t349 = rSges(7,3) * t413 + (rSges(7,1) * t421 - rSges(7,2) * t420) * t414;
t348 = Icges(7,5) * t413 + (Icges(7,1) * t421 - Icges(7,4) * t420) * t414;
t347 = Icges(7,6) * t413 + (Icges(7,4) * t421 - Icges(7,2) * t420) * t414;
t346 = Icges(7,3) * t413 + (Icges(7,5) * t421 - Icges(7,6) * t420) * t414;
t345 = pkin(9) * t413 + t414 * t475;
t344 = (-t380 * t427 + t381 * t430) * qJD(3);
t343 = rSges(6,1) * t390 + rSges(6,2) * t389 + rSges(6,3) * t467;
t342 = rSges(6,1) * t388 + rSges(6,2) * t387 - rSges(6,3) * t468;
t341 = Icges(6,1) * t390 + Icges(6,4) * t389 + Icges(6,5) * t467;
t340 = Icges(6,1) * t388 + Icges(6,4) * t387 - Icges(6,5) * t468;
t339 = Icges(6,4) * t390 + Icges(6,2) * t389 + Icges(6,6) * t467;
t338 = Icges(6,4) * t388 + Icges(6,2) * t387 - Icges(6,6) * t468;
t337 = Icges(6,5) * t390 + Icges(6,6) * t389 + Icges(6,3) * t467;
t336 = Icges(6,5) * t388 + Icges(6,6) * t387 - Icges(6,3) * t468;
t335 = pkin(5) * t462 - t482 * t430;
t334 = pkin(5) * t461 + t482 * t427;
t333 = rSges(7,1) * t371 + rSges(7,2) * t370 + rSges(7,3) * t467;
t332 = rSges(7,1) * t369 + rSges(7,2) * t368 - rSges(7,3) * t468;
t331 = Icges(7,1) * t371 + Icges(7,4) * t370 + Icges(7,5) * t467;
t330 = Icges(7,1) * t369 + Icges(7,4) * t368 - Icges(7,5) * t468;
t329 = Icges(7,4) * t371 + Icges(7,2) * t370 + Icges(7,6) * t467;
t328 = Icges(7,4) * t369 + Icges(7,2) * t368 - Icges(7,6) * t468;
t327 = Icges(7,5) * t371 + Icges(7,6) * t370 + Icges(7,3) * t467;
t326 = Icges(7,5) * t369 + Icges(7,6) * t368 - Icges(7,3) * t468;
t325 = qJD(1) * t380 + (-qJD(3) * t406 - qJD(2)) * t430 + t458;
t324 = t406 * t417 + t419 + (-t381 + t453) * qJD(1);
t323 = qJD(1) * t364 + (-qJD(2) + (-t397 - t477) * qJD(3)) * t430 + t452;
t322 = t397 * t417 + (-t365 + t451) * qJD(1) + t456;
t321 = t367 + (t365 * t430 + (-t364 - t386) * t427) * qJD(3);
t320 = t342 * t408 - t355 * t393 + t433;
t319 = -t343 * t408 + t355 * t392 + t434;
t318 = -t342 * t392 + t343 * t393 + t435;
t317 = t332 * t391 + t334 * t408 - t345 * t393 - t349 * t357 + t433;
t316 = -t333 * t391 - t335 * t408 + t345 * t392 + t349 * t356 + t434;
t315 = -t332 * t356 + t333 * t357 - t334 * t392 + t335 * t393 + t435;
t1 = m(7) * (t315 ^ 2 + t316 ^ 2 + t317 ^ 2) / 0.2e1 + m(6) * (t318 ^ 2 + t319 ^ 2 + t320 ^ 2) / 0.2e1 + m(5) * (t321 ^ 2 + t322 ^ 2 + t323 ^ 2) / 0.2e1 + m(4) * (t324 ^ 2 + t325 ^ 2 + t344 ^ 2) / 0.2e1 + m(3) * (t350 ^ 2 + t351 ^ 2) / 0.2e1 + t391 * ((t326 * t357 + t327 * t356 + t346 * t391) * t413 + ((-t328 * t420 + t330 * t421) * t357 + (-t329 * t420 + t331 * t421) * t356 + (-t347 * t420 + t348 * t421) * t391) * t414) / 0.2e1 + t357 * ((-t326 * t468 + t368 * t328 + t369 * t330) * t357 + (-t327 * t468 + t329 * t368 + t331 * t369) * t356 + (-t346 * t468 + t347 * t368 + t348 * t369) * t391) / 0.2e1 + t356 * ((t326 * t467 + t328 * t370 + t330 * t371) * t357 + (t327 * t467 + t370 * t329 + t371 * t331) * t356 + (t346 * t467 + t347 * t370 + t348 * t371) * t391) / 0.2e1 + t393 * ((-t336 * t468 + t387 * t338 + t388 * t340) * t393 + (-t337 * t468 + t339 * t387 + t341 * t388) * t392 + (-t352 * t468 + t353 * t387 + t354 * t388) * t408) / 0.2e1 + t392 * ((t336 * t467 + t338 * t389 + t340 * t390) * t393 + (t337 * t467 + t389 * t339 + t390 * t341) * t392 + (t352 * t467 + t353 * t389 + t354 * t390) * t408) / 0.2e1 + t408 * ((t336 * t393 + t337 * t392 + t352 * t408) * t413 + ((-t338 * t425 + t340 * t428) * t393 + (-t339 * t425 + t341 * t428) * t392 + (-t353 * t425 + t354 * t428) * t408) * t414) / 0.2e1 + (((-t360 * t413 + t362 * t414 - t376 * t426 + t378 * t429) * t430 + (-t413 * t361 + t414 * t363 - t377 * t426 + t379 * t429) * t427) * qJD(3) + (-t413 * t395 + t414 * t396 - t426 * t401 + t429 * t402) * qJD(1)) * qJD(1) / 0.2e1 + ((t488 * t427 ^ 2 + (t483 * t430 + (-t484 + t487) * t427) * t430) * qJD(3) + (t486 * t427 - t485 * t430) * qJD(1)) * t417 / 0.2e1 + ((t487 * t430 ^ 2 + (t484 * t427 + (-t483 + t488) * t430) * t427) * qJD(3) + (t485 * t427 + t486 * t430) * qJD(1)) * t418 / 0.2e1 + (Icges(2,3) + Icges(3,1) + m(2) * (t405 ^ 2 + t407 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
