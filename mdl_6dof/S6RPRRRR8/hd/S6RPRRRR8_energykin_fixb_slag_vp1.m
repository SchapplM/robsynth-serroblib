% Calculate kinetic energy for
% S6RPRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-03-09 07:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR8_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR8_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR8_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR8_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:19:44
% EndTime: 2019-03-09 07:19:46
% DurationCPUTime: 1.76s
% Computational Cost: add. (1321->273), mult. (1604->439), div. (0->0), fcn. (1538->10), ass. (0->145)
t420 = qJ(3) + qJ(4);
t416 = sin(t420);
t418 = cos(t420);
t424 = cos(qJ(5));
t473 = pkin(5) * t424;
t478 = -pkin(10) * t418 + t416 * t473;
t423 = sin(qJ(1));
t426 = cos(qJ(1));
t468 = Icges(5,4) * t416;
t443 = Icges(5,2) * t418 + t468;
t357 = Icges(5,6) * t426 + t423 * t443;
t358 = Icges(5,6) * t423 - t426 * t443;
t467 = Icges(5,4) * t418;
t445 = Icges(5,1) * t416 + t467;
t359 = Icges(5,5) * t426 + t423 * t445;
t360 = Icges(5,5) * t423 - t426 * t445;
t390 = -Icges(5,2) * t416 + t467;
t391 = Icges(5,1) * t418 - t468;
t412 = qJD(3) * t423;
t395 = qJD(4) * t423 + t412;
t413 = qJD(3) * t426;
t396 = qJD(4) * t426 + t413;
t477 = (t357 * t418 + t359 * t416) * t396 + (t358 * t418 + t360 * t416) * t395 + (t390 * t418 + t391 * t416) * qJD(1);
t422 = sin(qJ(3));
t475 = pkin(3) * t422;
t470 = Icges(4,4) * t422;
t425 = cos(qJ(3));
t469 = Icges(4,4) * t425;
t419 = qJ(5) + qJ(6);
t415 = sin(t419);
t466 = t415 * t423;
t465 = t415 * t426;
t417 = cos(t419);
t464 = t417 * t423;
t463 = t417 * t426;
t462 = t418 * t423;
t461 = t418 * t426;
t421 = sin(qJ(5));
t460 = t421 * t423;
t459 = t421 * t426;
t458 = t423 * t424;
t457 = t424 * t426;
t394 = qJD(1) * (pkin(1) * t426 + qJ(2) * t423);
t456 = qJD(1) * t426 * pkin(7) + t394;
t414 = qJD(2) * t423;
t453 = pkin(3) * qJD(3) * t425;
t455 = t423 * t453 + t414;
t454 = qJD(5) * t418;
t405 = qJD(5) * t416 + qJD(1);
t376 = t426 * t454 + t395;
t400 = pkin(1) * t423 - qJ(2) * t426;
t452 = -pkin(7) * t423 - t400;
t382 = pkin(8) * t423 - t426 * t475;
t451 = -t382 + t452;
t450 = pkin(4) * t416 - pkin(9) * t418;
t383 = pkin(8) * t426 + t423 * t475;
t449 = t382 * t413 - t383 * t412;
t448 = rSges(4,1) * t422 + rSges(4,2) * t425;
t447 = rSges(5,1) * t416 + rSges(5,2) * t418;
t446 = Icges(4,1) * t422 + t469;
t444 = Icges(4,2) * t425 + t470;
t442 = Icges(4,5) * t422 + Icges(4,6) * t425;
t441 = Icges(5,5) * t416 + Icges(5,6) * t418;
t371 = Icges(4,6) * t426 + t423 * t444;
t373 = Icges(4,5) * t426 + t423 * t446;
t438 = -t371 * t425 - t373 * t422;
t372 = Icges(4,6) * t423 - t426 * t444;
t374 = Icges(4,5) * t423 - t426 * t446;
t437 = t372 * t425 + t374 * t422;
t398 = -Icges(4,2) * t422 + t469;
t399 = Icges(4,1) * t425 - t470;
t435 = t398 * t425 + t399 * t422;
t434 = (Icges(5,5) * t418 - Icges(5,6) * t416) * qJD(1) + (Icges(5,3) * t426 + t423 * t441) * t396 + (Icges(5,3) * t423 - t426 * t441) * t395;
t380 = t450 * t423;
t381 = t450 * t426;
t433 = -t380 * t395 - t396 * t381 + t449;
t432 = qJD(1) * t383 + (-qJD(2) - t453) * t426 + t456;
t393 = pkin(4) * t418 + pkin(9) * t416;
t431 = t395 * t393 + (t381 + t451) * qJD(1) + t455;
t430 = qJD(1) * t380 - t393 * t396 + t432;
t403 = rSges(2,1) * t426 - rSges(2,2) * t423;
t402 = rSges(4,1) * t425 - rSges(4,2) * t422;
t401 = rSges(2,1) * t423 + rSges(2,2) * t426;
t397 = Icges(4,5) * t425 - Icges(4,6) * t422;
t392 = rSges(5,1) * t418 - rSges(5,2) * t416;
t388 = qJD(6) * t416 + t405;
t387 = -t416 * t457 + t460;
t386 = t416 * t459 + t458;
t385 = t416 * t458 + t459;
t384 = -t416 * t460 + t457;
t379 = rSges(4,3) * t423 - t426 * t448;
t378 = rSges(4,3) * t426 + t423 * t448;
t377 = -t423 * t454 + t396;
t370 = Icges(4,3) * t423 - t426 * t442;
t369 = Icges(4,3) * t426 + t423 * t442;
t367 = -t416 * t463 + t466;
t366 = t416 * t465 + t464;
t365 = t416 * t464 + t465;
t364 = -t416 * t466 + t463;
t362 = rSges(5,3) * t423 - t426 * t447;
t361 = rSges(5,3) * t426 + t423 * t447;
t353 = rSges(6,3) * t416 + (rSges(6,1) * t424 - rSges(6,2) * t421) * t418;
t352 = Icges(6,5) * t416 + (Icges(6,1) * t424 - Icges(6,4) * t421) * t418;
t351 = Icges(6,6) * t416 + (Icges(6,4) * t424 - Icges(6,2) * t421) * t418;
t350 = Icges(6,3) * t416 + (Icges(6,5) * t424 - Icges(6,6) * t421) * t418;
t349 = t394 - qJD(2) * t426 + qJD(1) * (-rSges(3,2) * t426 + rSges(3,3) * t423);
t348 = t414 + (rSges(3,2) * t423 + rSges(3,3) * t426 - t400) * qJD(1);
t347 = rSges(7,3) * t416 + (rSges(7,1) * t417 - rSges(7,2) * t415) * t418;
t346 = (-qJD(5) - qJD(6)) * t462 + t396;
t345 = qJD(6) * t461 + t376;
t343 = Icges(7,5) * t416 + (Icges(7,1) * t417 - Icges(7,4) * t415) * t418;
t342 = Icges(7,6) * t416 + (Icges(7,4) * t417 - Icges(7,2) * t415) * t418;
t341 = Icges(7,3) * t416 + (Icges(7,5) * t417 - Icges(7,6) * t415) * t418;
t340 = pkin(10) * t416 + t418 * t473;
t339 = rSges(6,1) * t387 + rSges(6,2) * t386 + rSges(6,3) * t461;
t338 = rSges(6,1) * t385 + rSges(6,2) * t384 - rSges(6,3) * t462;
t337 = Icges(6,1) * t387 + Icges(6,4) * t386 + Icges(6,5) * t461;
t336 = Icges(6,1) * t385 + Icges(6,4) * t384 - Icges(6,5) * t462;
t335 = Icges(6,4) * t387 + Icges(6,2) * t386 + Icges(6,6) * t461;
t334 = Icges(6,4) * t385 + Icges(6,2) * t384 - Icges(6,6) * t462;
t333 = Icges(6,5) * t387 + Icges(6,6) * t386 + Icges(6,3) * t461;
t332 = Icges(6,5) * t385 + Icges(6,6) * t384 - Icges(6,3) * t462;
t331 = (-t378 * t423 + t379 * t426) * qJD(3);
t330 = pkin(5) * t460 - t478 * t426;
t329 = pkin(5) * t459 + t478 * t423;
t328 = rSges(7,1) * t367 + rSges(7,2) * t366 + rSges(7,3) * t461;
t327 = rSges(7,1) * t365 + rSges(7,2) * t364 - rSges(7,3) * t462;
t326 = Icges(7,1) * t367 + Icges(7,4) * t366 + Icges(7,5) * t461;
t325 = Icges(7,1) * t365 + Icges(7,4) * t364 - Icges(7,5) * t462;
t324 = Icges(7,4) * t367 + Icges(7,2) * t366 + Icges(7,6) * t461;
t323 = Icges(7,4) * t365 + Icges(7,2) * t364 - Icges(7,6) * t462;
t322 = Icges(7,5) * t367 + Icges(7,6) * t366 + Icges(7,3) * t461;
t321 = Icges(7,5) * t365 + Icges(7,6) * t364 - Icges(7,3) * t462;
t320 = qJD(1) * t378 + (-qJD(3) * t402 - qJD(2)) * t426 + t456;
t319 = t402 * t412 + t414 + (-t379 + t452) * qJD(1);
t318 = qJD(1) * t361 - t392 * t396 + t432;
t317 = t392 * t395 + (-t362 + t451) * qJD(1) + t455;
t316 = -t361 * t395 + t362 * t396 + t449;
t315 = t338 * t405 - t353 * t377 + t430;
t314 = -t339 * t405 + t353 * t376 + t431;
t313 = -t338 * t376 + t339 * t377 + t433;
t312 = t327 * t388 + t329 * t405 - t340 * t377 - t346 * t347 + t430;
t311 = -t328 * t388 - t330 * t405 + t340 * t376 + t345 * t347 + t431;
t310 = -t327 * t345 + t328 * t346 - t329 * t376 + t330 * t377 + t433;
t1 = t395 * (t434 * t423 - t477 * t426) / 0.2e1 + m(7) * (t310 ^ 2 + t311 ^ 2 + t312 ^ 2) / 0.2e1 + m(6) * (t313 ^ 2 + t314 ^ 2 + t315 ^ 2) / 0.2e1 + m(5) * (t316 ^ 2 + t317 ^ 2 + t318 ^ 2) / 0.2e1 + m(4) * (t319 ^ 2 + t320 ^ 2 + t331 ^ 2) / 0.2e1 + m(3) * (t348 ^ 2 + t349 ^ 2) / 0.2e1 + ((t423 * t397 - t426 * t435) * qJD(1) + (t423 ^ 2 * t370 + (t438 * t426 + (t369 - t437) * t423) * t426) * qJD(3)) * t412 / 0.2e1 + t388 * ((t321 * t346 + t322 * t345 + t341 * t388) * t416 + ((-t323 * t415 + t325 * t417) * t346 + (-t324 * t415 + t326 * t417) * t345 + (-t342 * t415 + t343 * t417) * t388) * t418) / 0.2e1 + t346 * ((-t321 * t462 + t364 * t323 + t365 * t325) * t346 + (-t322 * t462 + t324 * t364 + t326 * t365) * t345 + (-t341 * t462 + t342 * t364 + t343 * t365) * t388) / 0.2e1 + t345 * ((t321 * t461 + t323 * t366 + t325 * t367) * t346 + (t322 * t461 + t366 * t324 + t367 * t326) * t345 + (t341 * t461 + t342 * t366 + t343 * t367) * t388) / 0.2e1 + t377 * ((-t332 * t462 + t384 * t334 + t385 * t336) * t377 + (-t333 * t462 + t335 * t384 + t337 * t385) * t376 + (-t350 * t462 + t384 * t351 + t385 * t352) * t405) / 0.2e1 + t376 * ((t332 * t461 + t334 * t386 + t336 * t387) * t377 + (t333 * t461 + t386 * t335 + t387 * t337) * t376 + (t350 * t461 + t351 * t386 + t352 * t387) * t405) / 0.2e1 + t405 * ((t332 * t377 + t333 * t376 + t350 * t405) * t416 + ((-t334 * t421 + t336 * t424) * t377 + (-t335 * t421 + t337 * t424) * t376 + (-t351 * t421 + t352 * t424) * t405) * t418) / 0.2e1 + t396 * (t477 * t423 + t434 * t426) / 0.2e1 + ((t426 * t397 + t423 * t435) * qJD(1) + (t426 ^ 2 * t369 + (t437 * t423 + (t370 - t438) * t426) * t423) * qJD(3)) * t413 / 0.2e1 + ((-t357 * t416 + t359 * t418) * t396 + (-t358 * t416 + t360 * t418) * t395 + ((-t371 * t422 + t373 * t425) * t426 + (-t372 * t422 + t374 * t425) * t423) * qJD(3) + (-t416 * t390 + t418 * t391 - t422 * t398 + t425 * t399) * qJD(1)) * qJD(1) / 0.2e1 + (m(2) * (t401 ^ 2 + t403 ^ 2) + Icges(2,3) + Icges(3,1)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
