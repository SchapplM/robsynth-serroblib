% Calculate kinetic energy for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR10_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR10_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR10_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR10_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:23:45
% EndTime: 2019-12-31 20:23:48
% DurationCPUTime: 2.39s
% Computational Cost: add. (2237->300), mult. (5620->467), div. (0->0), fcn. (7045->12), ass. (0->138)
t447 = sin(qJ(1));
t450 = cos(qJ(1));
t446 = sin(qJ(2));
t449 = cos(qJ(2));
t480 = sin(pkin(10));
t481 = cos(pkin(10));
t428 = -t446 * t480 + t449 * t481;
t443 = cos(pkin(5));
t454 = t443 * t428;
t455 = t446 * t481 + t449 * t480;
t405 = -t447 * t455 + t450 * t454;
t420 = t455 * t443;
t406 = t420 * t450 + t428 * t447;
t442 = sin(pkin(5));
t474 = t450 * t442;
t358 = Icges(4,5) * t406 + Icges(4,6) * t405 - Icges(4,3) * t474;
t407 = -t447 * t454 - t450 * t455;
t408 = -t420 * t447 + t428 * t450;
t477 = t447 * t442;
t359 = Icges(4,5) * t408 + Icges(4,6) * t407 + Icges(4,3) * t477;
t475 = t449 * t450;
t479 = t446 * t447;
t422 = t443 * t475 - t479;
t476 = t447 * t449;
t478 = t446 * t450;
t423 = t443 * t478 + t476;
t424 = -t443 * t476 - t478;
t425 = -t443 * t479 + t475;
t456 = (Icges(3,5) * t423 + Icges(3,6) * t422 - Icges(3,3) * t474) * t450 - (Icges(3,5) * t425 + Icges(3,6) * t424 + Icges(3,3) * t477) * t447;
t489 = (t358 * t450 - t359 * t447 + t456) * t442;
t418 = t428 * t442;
t419 = t455 * t442;
t488 = Icges(4,5) * t419 + Icges(4,6) * t418 + (Icges(3,5) * t446 + Icges(3,6) * t449) * t442 + (Icges(4,3) + Icges(3,3)) * t443;
t484 = cos(qJ(4));
t483 = pkin(2) * t446;
t482 = pkin(2) * t449;
t466 = -qJ(3) * t442 + t443 * t483;
t404 = -t447 * t466 + t450 * t482;
t426 = qJD(1) * (pkin(1) * t450 + pkin(7) * t477);
t437 = qJD(2) * t443 + qJD(1);
t473 = t437 * t404 + t426;
t471 = qJD(2) * t442;
t436 = t447 * t471;
t383 = -qJD(4) * t407 + t436;
t472 = qJD(1) * (pkin(1) * t447 - pkin(7) * t474);
t470 = qJD(3) * t450;
t403 = t447 * t482 + t450 * t466;
t467 = t450 * t471;
t469 = qJD(3) * t443 + t403 * t436 + t404 * t467;
t468 = t442 * t484;
t411 = -qJD(4) * t418 + t437;
t429 = qJ(3) * t443 + t442 * t483;
t463 = qJD(2) * (-rSges(4,1) * t419 - rSges(4,2) * t418 - rSges(4,3) * t443 - t429);
t462 = qJD(2) * (-pkin(3) * t419 + pkin(8) * t418 - t429);
t461 = qJD(3) * t477 - t472;
t384 = -qJD(4) * t405 - t467;
t373 = pkin(3) * t406 - pkin(8) * t405;
t374 = pkin(3) * t408 - pkin(8) * t407;
t458 = t373 * t436 + t374 * t467 + t469;
t453 = t437 * t374 + (t447 * t462 - t470) * t442 + t473;
t452 = (-t373 - t403) * t437 + t462 * t474 + t461;
t448 = cos(qJ(5));
t445 = sin(qJ(4));
t444 = sin(qJ(5));
t432 = rSges(2,1) * t450 - rSges(2,2) * t447;
t431 = rSges(2,1) * t447 + rSges(2,2) * t450;
t417 = rSges(3,3) * t443 + (rSges(3,1) * t446 + rSges(3,2) * t449) * t442;
t416 = Icges(3,5) * t443 + (Icges(3,1) * t446 + Icges(3,4) * t449) * t442;
t415 = Icges(3,6) * t443 + (Icges(3,4) * t446 + Icges(3,2) * t449) * t442;
t410 = t419 * t484 + t443 * t445;
t409 = t419 * t445 - t443 * t484;
t399 = rSges(3,1) * t425 + rSges(3,2) * t424 + rSges(3,3) * t477;
t398 = rSges(3,1) * t423 + rSges(3,2) * t422 - rSges(3,3) * t474;
t397 = Icges(3,1) * t425 + Icges(3,4) * t424 + Icges(3,5) * t477;
t396 = Icges(3,1) * t423 + Icges(3,4) * t422 - Icges(3,5) * t474;
t395 = Icges(3,4) * t425 + Icges(3,2) * t424 + Icges(3,6) * t477;
t394 = Icges(3,4) * t423 + Icges(3,2) * t422 - Icges(3,6) * t474;
t390 = Icges(4,1) * t419 + Icges(4,4) * t418 + Icges(4,5) * t443;
t389 = Icges(4,4) * t419 + Icges(4,2) * t418 + Icges(4,6) * t443;
t382 = t408 * t484 + t445 * t477;
t381 = t408 * t445 - t447 * t468;
t380 = t406 * t484 - t445 * t474;
t379 = t406 * t445 + t450 * t468;
t378 = t410 * t448 - t418 * t444;
t377 = -t410 * t444 - t418 * t448;
t376 = qJD(5) * t409 + t411;
t375 = pkin(4) * t410 + pkin(9) * t409;
t372 = rSges(5,1) * t410 - rSges(5,2) * t409 - rSges(5,3) * t418;
t371 = Icges(5,1) * t410 - Icges(5,4) * t409 - Icges(5,5) * t418;
t370 = Icges(5,4) * t410 - Icges(5,2) * t409 - Icges(5,6) * t418;
t369 = Icges(5,5) * t410 - Icges(5,6) * t409 - Icges(5,3) * t418;
t366 = rSges(4,1) * t408 + rSges(4,2) * t407 + rSges(4,3) * t477;
t365 = rSges(4,1) * t406 + rSges(4,2) * t405 - rSges(4,3) * t474;
t363 = Icges(4,1) * t408 + Icges(4,4) * t407 + Icges(4,5) * t477;
t362 = Icges(4,1) * t406 + Icges(4,4) * t405 - Icges(4,5) * t474;
t361 = Icges(4,4) * t408 + Icges(4,2) * t407 + Icges(4,6) * t477;
t360 = Icges(4,4) * t406 + Icges(4,2) * t405 - Icges(4,6) * t474;
t357 = t399 * t437 - t417 * t436 + t426;
t356 = -t398 * t437 - t417 * t467 - t472;
t355 = t382 * t448 - t407 * t444;
t354 = -t382 * t444 - t407 * t448;
t353 = t380 * t448 - t405 * t444;
t352 = -t380 * t444 - t405 * t448;
t351 = (t398 * t447 + t399 * t450) * t471;
t350 = qJD(5) * t379 + t384;
t349 = qJD(5) * t381 + t383;
t348 = pkin(4) * t382 + pkin(9) * t381;
t347 = pkin(4) * t380 + pkin(9) * t379;
t346 = rSges(6,1) * t378 + rSges(6,2) * t377 + rSges(6,3) * t409;
t345 = Icges(6,1) * t378 + Icges(6,4) * t377 + Icges(6,5) * t409;
t344 = Icges(6,4) * t378 + Icges(6,2) * t377 + Icges(6,6) * t409;
t343 = Icges(6,5) * t378 + Icges(6,6) * t377 + Icges(6,3) * t409;
t342 = rSges(5,1) * t382 - rSges(5,2) * t381 - rSges(5,3) * t407;
t341 = rSges(5,1) * t380 - rSges(5,2) * t379 - rSges(5,3) * t405;
t340 = Icges(5,1) * t382 - Icges(5,4) * t381 - Icges(5,5) * t407;
t339 = Icges(5,1) * t380 - Icges(5,4) * t379 - Icges(5,5) * t405;
t338 = Icges(5,4) * t382 - Icges(5,2) * t381 - Icges(5,6) * t407;
t337 = Icges(5,4) * t380 - Icges(5,2) * t379 - Icges(5,6) * t405;
t336 = Icges(5,5) * t382 - Icges(5,6) * t381 - Icges(5,3) * t407;
t335 = Icges(5,5) * t380 - Icges(5,6) * t379 - Icges(5,3) * t405;
t334 = t366 * t437 + (t447 * t463 - t470) * t442 + t473;
t333 = (-t365 - t403) * t437 + t463 * t474 + t461;
t332 = rSges(6,1) * t355 + rSges(6,2) * t354 + rSges(6,3) * t381;
t331 = rSges(6,1) * t353 + rSges(6,2) * t352 + rSges(6,3) * t379;
t330 = Icges(6,1) * t355 + Icges(6,4) * t354 + Icges(6,5) * t381;
t329 = Icges(6,1) * t353 + Icges(6,4) * t352 + Icges(6,5) * t379;
t328 = Icges(6,4) * t355 + Icges(6,2) * t354 + Icges(6,6) * t381;
t327 = Icges(6,4) * t353 + Icges(6,2) * t352 + Icges(6,6) * t379;
t326 = Icges(6,5) * t355 + Icges(6,6) * t354 + Icges(6,3) * t381;
t325 = Icges(6,5) * t353 + Icges(6,6) * t352 + Icges(6,3) * t379;
t324 = (t365 * t447 + t366 * t450) * t471 + t469;
t323 = t342 * t411 - t372 * t383 + t453;
t322 = -t341 * t411 + t372 * t384 + t452;
t321 = t341 * t383 - t342 * t384 + t458;
t320 = t332 * t376 - t346 * t349 + t348 * t411 - t375 * t383 + t453;
t319 = -t331 * t376 + t346 * t350 - t347 * t411 + t375 * t384 + t452;
t318 = t331 * t349 - t332 * t350 + t347 * t383 - t348 * t384 + t458;
t1 = m(3) * (t351 ^ 2 + t356 ^ 2 + t357 ^ 2) / 0.2e1 + m(4) * (t324 ^ 2 + t333 ^ 2 + t334 ^ 2) / 0.2e1 + m(5) * (t321 ^ 2 + t322 ^ 2 + t323 ^ 2) / 0.2e1 + t383 * ((-t407 * t336 - t381 * t338 + t382 * t340) * t383 + (-t335 * t407 - t337 * t381 + t339 * t382) * t384 + (-t369 * t407 - t370 * t381 + t371 * t382) * t411) / 0.2e1 + t384 * ((-t336 * t405 - t338 * t379 + t340 * t380) * t383 + (-t405 * t335 - t379 * t337 + t380 * t339) * t384 + (-t369 * t405 - t370 * t379 + t371 * t380) * t411) / 0.2e1 + t411 * ((-t336 * t418 - t338 * t409 + t340 * t410) * t383 + (-t335 * t418 - t337 * t409 + t339 * t410) * t384 + (-t418 * t369 - t409 * t370 + t410 * t371) * t411) / 0.2e1 + m(6) * (t318 ^ 2 + t319 ^ 2 + t320 ^ 2) / 0.2e1 + t349 * ((t381 * t326 + t354 * t328 + t355 * t330) * t349 + (t325 * t381 + t327 * t354 + t329 * t355) * t350 + (t343 * t381 + t344 * t354 + t345 * t355) * t376) / 0.2e1 + t350 * ((t326 * t379 + t328 * t352 + t330 * t353) * t349 + (t379 * t325 + t352 * t327 + t353 * t329) * t350 + (t343 * t379 + t344 * t352 + t345 * t353) * t376) / 0.2e1 + t376 * ((t326 * t409 + t328 * t377 + t330 * t378) * t349 + (t325 * t409 + t327 * t377 + t329 * t378) * t350 + (t409 * t343 + t377 * t344 + t378 * t345) * t376) / 0.2e1 + ((((t395 * t449 + t397 * t446) * t447 - (t394 * t449 + t396 * t446) * t450) * t442 - t456 * t443 + (t359 * t443 + t361 * t418 + t363 * t419) * t447 - (t358 * t443 + t360 * t418 + t362 * t419) * t450) * t471 + ((t415 * t449 + t416 * t446) * t442 + t418 * t389 + t419 * t390 + t488 * t443) * t437) * t437 / 0.2e1 + (m(2) * (t431 ^ 2 + t432 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t360 * t407 - t362 * t408 - t394 * t424 - t396 * t425) * t450 + (t407 * t361 + t408 * t363 + t424 * t395 + t425 * t397 - t489) * t447) * t471 + (t389 * t407 + t390 * t408 + t415 * t424 + t416 * t425 + t488 * t477) * t437) * t436 / 0.2e1 - (((-t405 * t360 - t406 * t362 - t422 * t394 - t423 * t396 + t489) * t450 + (t361 * t405 + t363 * t406 + t395 * t422 + t397 * t423) * t447) * t471 + (t389 * t405 + t390 * t406 + t415 * t422 + t416 * t423 - t488 * t474) * t437) * t467 / 0.2e1;
T = t1;
