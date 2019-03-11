% Calculate kinetic energy for
% S6RPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:00:16
% EndTime: 2019-03-09 07:00:17
% DurationCPUTime: 1.99s
% Computational Cost: add. (2101->296), mult. (1957->480), div. (0->0), fcn. (1946->12), ass. (0->147)
t442 = sin(qJ(1));
t485 = pkin(1) * t442;
t443 = cos(qJ(4));
t483 = t443 * pkin(4);
t441 = sin(qJ(3));
t481 = Icges(4,4) * t441;
t444 = cos(qJ(3));
t480 = Icges(4,4) * t444;
t437 = qJ(1) + pkin(11);
t431 = sin(t437);
t440 = sin(qJ(4));
t479 = t431 * t440;
t478 = t431 * t441;
t432 = cos(t437);
t477 = t432 * t440;
t476 = t432 * t441;
t475 = t444 * t431;
t474 = t444 * t432;
t439 = qJ(4) + qJ(5);
t445 = cos(qJ(1));
t430 = qJD(1) * t445 * pkin(1);
t473 = qJD(1) * (pkin(2) * t432 + pkin(7) * t431) + t430;
t434 = cos(t439);
t472 = pkin(5) * t434;
t425 = qJD(3) * t431;
t469 = qJD(4) * t441;
t406 = t432 * t469 + t425;
t470 = qJD(3) * t432;
t468 = qJD(5) * t441;
t467 = qJD(6) * t441;
t466 = -qJD(4) - qJD(5);
t364 = t432 * t468 + t406;
t462 = pkin(3) * t444 + pkin(8) * t441;
t404 = t462 * t431;
t405 = t462 * t432;
t465 = t404 * t425 + t405 * t470 + qJD(2);
t464 = -pkin(2) * t431 + pkin(7) * t432 - t485;
t433 = sin(t439);
t463 = pkin(5) * t433;
t407 = t431 * t469 - t470;
t365 = t431 * t468 + t407;
t461 = rSges(4,1) * t444 - rSges(4,2) * t441;
t460 = Icges(4,1) * t444 - t481;
t459 = -Icges(4,2) * t441 + t480;
t458 = Icges(4,5) * t444 - Icges(4,6) * t441;
t368 = -Icges(4,6) * t432 + t431 * t459;
t370 = -Icges(4,5) * t432 + t431 * t460;
t457 = t368 * t441 - t370 * t444;
t369 = Icges(4,6) * t431 + t432 * t459;
t371 = Icges(4,5) * t431 + t432 * t460;
t456 = -t369 * t441 + t371 * t444;
t415 = Icges(4,2) * t444 + t481;
t416 = Icges(4,1) * t441 + t480;
t455 = -t415 * t441 + t416 * t444;
t424 = pkin(3) * t441 - pkin(8) * t444;
t454 = qJD(1) * t405 - t424 * t425 + t473;
t452 = pkin(9) * t441 + t444 * t483;
t358 = -pkin(4) * t477 + t431 * t452;
t359 = pkin(4) * t479 + t432 * t452;
t453 = t406 * t358 - t359 * t407 + t465;
t451 = pkin(10) * t441 + t444 * t472;
t450 = (-t404 + t464) * qJD(1) - t424 * t470;
t384 = -pkin(9) * t444 + t441 * t483;
t426 = -qJD(4) * t444 + qJD(1);
t449 = t426 * t359 - t384 * t406 + t454;
t448 = -t358 * t426 + t407 * t384 + t450;
t435 = qJ(6) + t439;
t428 = cos(t435);
t427 = sin(t435);
t423 = rSges(2,1) * t445 - rSges(2,2) * t442;
t422 = rSges(2,1) * t442 + rSges(2,2) * t445;
t421 = rSges(4,1) * t441 + rSges(4,2) * t444;
t414 = Icges(4,5) * t441 + Icges(4,6) * t444;
t412 = t444 * t466 + qJD(1);
t408 = qJD(1) + (-qJD(6) + t466) * t444;
t403 = -rSges(5,3) * t444 + (rSges(5,1) * t443 - rSges(5,2) * t440) * t441;
t402 = -Icges(5,5) * t444 + (Icges(5,1) * t443 - Icges(5,4) * t440) * t441;
t401 = -Icges(5,6) * t444 + (Icges(5,4) * t443 - Icges(5,2) * t440) * t441;
t400 = -Icges(5,3) * t444 + (Icges(5,5) * t443 - Icges(5,6) * t440) * t441;
t399 = t443 * t474 + t479;
t398 = t431 * t443 - t440 * t474;
t397 = t443 * t475 - t477;
t396 = -t432 * t443 - t440 * t475;
t394 = t430 + qJD(1) * (rSges(3,1) * t432 - rSges(3,2) * t431);
t393 = (-rSges(3,1) * t431 - rSges(3,2) * t432 - t485) * qJD(1);
t392 = -rSges(6,3) * t444 + (rSges(6,1) * t434 - rSges(6,2) * t433) * t441;
t391 = -Icges(6,5) * t444 + (Icges(6,1) * t434 - Icges(6,4) * t433) * t441;
t390 = -Icges(6,6) * t444 + (Icges(6,4) * t434 - Icges(6,2) * t433) * t441;
t389 = -Icges(6,3) * t444 + (Icges(6,5) * t434 - Icges(6,6) * t433) * t441;
t388 = t431 * t433 + t434 * t474;
t387 = t431 * t434 - t433 * t474;
t386 = -t432 * t433 + t434 * t475;
t385 = -t432 * t434 - t433 * t475;
t383 = -rSges(7,3) * t444 + (rSges(7,1) * t428 - rSges(7,2) * t427) * t441;
t380 = -Icges(7,5) * t444 + (Icges(7,1) * t428 - Icges(7,4) * t427) * t441;
t379 = -Icges(7,6) * t444 + (Icges(7,4) * t428 - Icges(7,2) * t427) * t441;
t378 = -Icges(7,3) * t444 + (Icges(7,5) * t428 - Icges(7,6) * t427) * t441;
t377 = rSges(4,3) * t431 + t432 * t461;
t376 = -rSges(4,3) * t432 + t431 * t461;
t375 = t427 * t431 + t428 * t474;
t374 = -t427 * t474 + t428 * t431;
t373 = -t427 * t432 + t428 * t475;
t372 = -t427 * t475 - t428 * t432;
t367 = Icges(4,3) * t431 + t432 * t458;
t366 = -Icges(4,3) * t432 + t431 * t458;
t363 = t431 * t467 + t365;
t362 = t432 * t467 + t364;
t361 = -pkin(10) * t444 + t441 * t472;
t357 = rSges(5,1) * t399 + rSges(5,2) * t398 + rSges(5,3) * t476;
t356 = rSges(5,1) * t397 + rSges(5,2) * t396 + rSges(5,3) * t478;
t355 = Icges(5,1) * t399 + Icges(5,4) * t398 + Icges(5,5) * t476;
t354 = Icges(5,1) * t397 + Icges(5,4) * t396 + Icges(5,5) * t478;
t353 = Icges(5,4) * t399 + Icges(5,2) * t398 + Icges(5,6) * t476;
t352 = Icges(5,4) * t397 + Icges(5,2) * t396 + Icges(5,6) * t478;
t351 = Icges(5,5) * t399 + Icges(5,6) * t398 + Icges(5,3) * t476;
t350 = Icges(5,5) * t397 + Icges(5,6) * t396 + Icges(5,3) * t478;
t348 = rSges(6,1) * t388 + rSges(6,2) * t387 + rSges(6,3) * t476;
t347 = rSges(6,1) * t386 + rSges(6,2) * t385 + rSges(6,3) * t478;
t346 = Icges(6,1) * t388 + Icges(6,4) * t387 + Icges(6,5) * t476;
t345 = Icges(6,1) * t386 + Icges(6,4) * t385 + Icges(6,5) * t478;
t344 = Icges(6,4) * t388 + Icges(6,2) * t387 + Icges(6,6) * t476;
t343 = Icges(6,4) * t386 + Icges(6,2) * t385 + Icges(6,6) * t478;
t342 = Icges(6,5) * t388 + Icges(6,6) * t387 + Icges(6,3) * t476;
t341 = Icges(6,5) * t386 + Icges(6,6) * t385 + Icges(6,3) * t478;
t340 = rSges(7,1) * t375 + rSges(7,2) * t374 + rSges(7,3) * t476;
t339 = rSges(7,1) * t373 + rSges(7,2) * t372 + rSges(7,3) * t478;
t338 = Icges(7,1) * t375 + Icges(7,4) * t374 + Icges(7,5) * t476;
t337 = Icges(7,1) * t373 + Icges(7,4) * t372 + Icges(7,5) * t478;
t336 = Icges(7,4) * t375 + Icges(7,2) * t374 + Icges(7,6) * t476;
t335 = Icges(7,4) * t373 + Icges(7,2) * t372 + Icges(7,6) * t478;
t334 = Icges(7,5) * t375 + Icges(7,6) * t374 + Icges(7,3) * t476;
t333 = Icges(7,5) * t373 + Icges(7,6) * t372 + Icges(7,3) * t478;
t332 = qJD(1) * t377 - t421 * t425 + t473;
t331 = -t421 * t470 + (-t376 + t464) * qJD(1);
t330 = qJD(2) + (t376 * t431 + t377 * t432) * qJD(3);
t328 = t431 * t463 + t432 * t451;
t327 = t431 * t451 - t432 * t463;
t326 = t357 * t426 - t403 * t406 + t454;
t325 = -t356 * t426 + t403 * t407 + t450;
t324 = t356 * t406 - t357 * t407 + t465;
t323 = t348 * t412 - t364 * t392 + t449;
t322 = -t347 * t412 + t365 * t392 + t448;
t321 = t347 * t364 - t348 * t365 + t453;
t320 = t328 * t412 + t340 * t408 - t361 * t364 - t362 * t383 + t449;
t319 = -t327 * t412 - t339 * t408 + t361 * t365 + t363 * t383 + t448;
t318 = t327 * t364 - t328 * t365 + t339 * t362 - t340 * t363 + t453;
t1 = t365 * ((t342 * t478 + t344 * t385 + t346 * t386) * t364 + (t341 * t478 + t385 * t343 + t386 * t345) * t365 + (t385 * t390 + t386 * t391 + t389 * t478) * t412) / 0.2e1 + t412 * ((-t341 * t365 - t342 * t364 - t389 * t412) * t444 + ((-t344 * t433 + t346 * t434) * t364 + (-t343 * t433 + t345 * t434) * t365 + (-t390 * t433 + t391 * t434) * t412) * t441) / 0.2e1 + t364 * ((t342 * t476 + t387 * t344 + t388 * t346) * t364 + (t341 * t476 + t343 * t387 + t345 * t388) * t365 + (t387 * t390 + t388 * t391 + t389 * t476) * t412) / 0.2e1 + t406 * ((t351 * t476 + t398 * t353 + t399 * t355) * t406 + (t350 * t476 + t352 * t398 + t354 * t399) * t407 + (t398 * t401 + t399 * t402 + t400 * t476) * t426) / 0.2e1 + t407 * ((t351 * t478 + t353 * t396 + t355 * t397) * t406 + (t350 * t478 + t396 * t352 + t397 * t354) * t407 + (t396 * t401 + t397 * t402 + t400 * t478) * t426) / 0.2e1 + t426 * ((-t350 * t407 - t351 * t406 - t400 * t426) * t444 + ((-t353 * t440 + t355 * t443) * t406 + (-t352 * t440 + t354 * t443) * t407 + (-t401 * t440 + t402 * t443) * t426) * t441) / 0.2e1 - ((-t432 * t414 + t431 * t455) * qJD(1) + (t432 ^ 2 * t366 + (t456 * t431 + (-t367 + t457) * t432) * t431) * qJD(3)) * t470 / 0.2e1 + qJD(1) * ((t444 * t415 + t441 * t416) * qJD(1) + ((t369 * t444 + t371 * t441) * t431 - (t368 * t444 + t370 * t441) * t432) * qJD(3)) / 0.2e1 + m(7) * (t318 ^ 2 + t319 ^ 2 + t320 ^ 2) / 0.2e1 + m(6) * (t321 ^ 2 + t322 ^ 2 + t323 ^ 2) / 0.2e1 + m(5) * (t324 ^ 2 + t325 ^ 2 + t326 ^ 2) / 0.2e1 + m(4) * (t330 ^ 2 + t331 ^ 2 + t332 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 + t393 ^ 2 + t394 ^ 2) / 0.2e1 + t362 * ((t334 * t476 + t374 * t336 + t375 * t338) * t362 + (t333 * t476 + t335 * t374 + t337 * t375) * t363 + (t374 * t379 + t375 * t380 + t378 * t476) * t408) / 0.2e1 + t363 * ((t334 * t478 + t336 * t372 + t338 * t373) * t362 + (t333 * t478 + t372 * t335 + t373 * t337) * t363 + (t372 * t379 + t373 * t380 + t378 * t478) * t408) / 0.2e1 + t408 * ((-t333 * t363 - t334 * t362 - t378 * t408) * t444 + ((-t336 * t427 + t338 * t428) * t362 + (-t335 * t427 + t337 * t428) * t363 + (-t379 * t427 + t380 * t428) * t408) * t441) / 0.2e1 + ((t431 * t414 + t432 * t455) * qJD(1) + (t431 ^ 2 * t367 + (t457 * t432 + (-t366 + t456) * t431) * t432) * qJD(3)) * t425 / 0.2e1 + (Icges(2,3) + Icges(3,3) + m(2) * (t422 ^ 2 + t423 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
