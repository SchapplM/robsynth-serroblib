% Calculate kinetic energy for
% S6RPRPRR10
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
% Datum: 2019-03-09 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR10_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR10_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR10_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR10_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:07:37
% EndTime: 2019-03-09 04:07:39
% DurationCPUTime: 2.04s
% Computational Cost: add. (1279->302), mult. (1863->469), div. (0->0), fcn. (1853->10), ass. (0->141)
t433 = sin(qJ(3));
t435 = cos(qJ(3));
t429 = pkin(10) + qJ(5);
t420 = cos(t429);
t459 = pkin(5) * t420;
t476 = -pkin(9) * t435 + t433 * t459;
t431 = cos(pkin(10));
t471 = t431 * pkin(4);
t475 = -pkin(8) * t435 + t433 * t471;
t470 = Icges(4,4) * t433;
t469 = Icges(4,4) * t435;
t430 = sin(pkin(10));
t434 = sin(qJ(1));
t468 = t430 * t434;
t436 = cos(qJ(1));
t467 = t430 * t436;
t466 = t433 * t434;
t465 = t433 * t436;
t464 = t434 * t435;
t463 = t435 * t436;
t447 = pkin(3) * t433 - qJ(4) * t435;
t396 = t447 * t436;
t425 = qJD(3) * t436;
t461 = qJD(4) * t433 - t396 * t425;
t409 = pkin(3) * t435 + qJ(4) * t433;
t424 = qJD(3) * t434;
t426 = qJD(2) * t434;
t460 = t409 * t424 + t426;
t403 = qJD(1) * (pkin(1) * t436 + qJ(2) * t434);
t458 = qJD(1) * t436 * pkin(7) + t403;
t455 = qJD(5) * t435;
t398 = t436 * t455 + t424;
t456 = qJD(4) * t435;
t414 = qJD(5) * t433 + qJD(1);
t419 = sin(t429);
t452 = pkin(5) * t419;
t407 = pkin(1) * t434 - qJ(2) * t436;
t451 = -pkin(7) * t434 - t407;
t395 = t447 * t434;
t450 = qJD(1) * t395 + t436 * t456 + t458;
t449 = t396 + t451;
t448 = rSges(4,1) * t433 + rSges(4,2) * t435;
t446 = Icges(4,1) * t433 + t469;
t445 = Icges(4,2) * t435 + t470;
t444 = Icges(4,5) * t433 + Icges(4,6) * t435;
t383 = Icges(4,6) * t436 + t434 * t445;
t385 = Icges(4,5) * t436 + t434 * t446;
t443 = -t383 * t435 - t385 * t433;
t384 = Icges(4,6) * t434 - t436 * t445;
t386 = Icges(4,5) * t434 - t436 * t446;
t442 = t384 * t435 + t386 * t433;
t405 = -Icges(4,2) * t433 + t469;
t406 = Icges(4,1) * t435 - t470;
t441 = t405 * t435 + t406 * t433;
t352 = pkin(4) * t467 + t475 * t434;
t353 = pkin(4) * t468 - t475 * t436;
t440 = t353 * t425 + (-t352 - t395) * t424 + t461;
t362 = pkin(8) * t433 + t435 * t471;
t439 = qJD(1) * t352 + (-qJD(2) + (-t362 - t409) * qJD(3)) * t436 + t450;
t438 = t362 * t424 + (-t353 + t449) * qJD(1) - t434 * t456 + t460;
t421 = qJ(6) + t429;
t416 = cos(t421);
t415 = sin(t421);
t411 = rSges(2,1) * t436 - rSges(2,2) * t434;
t410 = rSges(4,1) * t435 - rSges(4,2) * t433;
t408 = rSges(2,1) * t434 + rSges(2,2) * t436;
t404 = Icges(4,5) * t435 - Icges(4,6) * t433;
t402 = qJD(6) * t433 + t414;
t399 = -t434 * t455 + t425;
t394 = -t431 * t465 + t468;
t393 = t430 * t465 + t431 * t434;
t392 = t431 * t466 + t467;
t391 = -t430 * t466 + t431 * t436;
t389 = rSges(4,3) * t434 - t436 * t448;
t388 = rSges(4,3) * t436 + t434 * t448;
t382 = Icges(4,3) * t434 - t436 * t444;
t381 = Icges(4,3) * t436 + t434 * t444;
t380 = t419 * t434 - t420 * t465;
t379 = t419 * t465 + t420 * t434;
t378 = t419 * t436 + t420 * t466;
t377 = -t419 * t466 + t420 * t436;
t376 = rSges(5,3) * t433 + (rSges(5,1) * t431 - rSges(5,2) * t430) * t435;
t375 = t425 + (-qJD(5) - qJD(6)) * t464;
t374 = qJD(6) * t463 + t398;
t373 = Icges(5,5) * t433 + (Icges(5,1) * t431 - Icges(5,4) * t430) * t435;
t372 = Icges(5,6) * t433 + (Icges(5,4) * t431 - Icges(5,2) * t430) * t435;
t371 = Icges(5,3) * t433 + (Icges(5,5) * t431 - Icges(5,6) * t430) * t435;
t370 = t415 * t434 - t416 * t465;
t369 = t415 * t465 + t416 * t434;
t368 = t415 * t436 + t416 * t466;
t367 = -t415 * t466 + t416 * t436;
t366 = rSges(6,3) * t433 + (rSges(6,1) * t420 - rSges(6,2) * t419) * t435;
t365 = Icges(6,5) * t433 + (Icges(6,1) * t420 - Icges(6,4) * t419) * t435;
t364 = Icges(6,6) * t433 + (Icges(6,4) * t420 - Icges(6,2) * t419) * t435;
t363 = Icges(6,3) * t433 + (Icges(6,5) * t420 - Icges(6,6) * t419) * t435;
t361 = rSges(7,3) * t433 + (rSges(7,1) * t416 - rSges(7,2) * t415) * t435;
t360 = Icges(7,5) * t433 + (Icges(7,1) * t416 - Icges(7,4) * t415) * t435;
t359 = Icges(7,6) * t433 + (Icges(7,4) * t416 - Icges(7,2) * t415) * t435;
t358 = Icges(7,3) * t433 + (Icges(7,5) * t416 - Icges(7,6) * t415) * t435;
t357 = t403 - qJD(2) * t436 + qJD(1) * (-rSges(3,2) * t436 + rSges(3,3) * t434);
t356 = t426 + (rSges(3,2) * t434 + rSges(3,3) * t436 - t407) * qJD(1);
t354 = pkin(9) * t433 + t435 * t459;
t351 = rSges(5,1) * t394 + rSges(5,2) * t393 + rSges(5,3) * t463;
t350 = rSges(5,1) * t392 + rSges(5,2) * t391 - rSges(5,3) * t464;
t349 = Icges(5,1) * t394 + Icges(5,4) * t393 + Icges(5,5) * t463;
t348 = Icges(5,1) * t392 + Icges(5,4) * t391 - Icges(5,5) * t464;
t347 = Icges(5,4) * t394 + Icges(5,2) * t393 + Icges(5,6) * t463;
t346 = Icges(5,4) * t392 + Icges(5,2) * t391 - Icges(5,6) * t464;
t345 = Icges(5,5) * t394 + Icges(5,6) * t393 + Icges(5,3) * t463;
t344 = Icges(5,5) * t392 + Icges(5,6) * t391 - Icges(5,3) * t464;
t341 = (-t388 * t434 + t389 * t436) * qJD(3);
t340 = rSges(6,1) * t380 + rSges(6,2) * t379 + rSges(6,3) * t463;
t339 = rSges(6,1) * t378 + rSges(6,2) * t377 - rSges(6,3) * t464;
t338 = Icges(6,1) * t380 + Icges(6,4) * t379 + Icges(6,5) * t463;
t337 = Icges(6,1) * t378 + Icges(6,4) * t377 - Icges(6,5) * t464;
t336 = Icges(6,4) * t380 + Icges(6,2) * t379 + Icges(6,6) * t463;
t335 = Icges(6,4) * t378 + Icges(6,2) * t377 - Icges(6,6) * t464;
t334 = Icges(6,5) * t380 + Icges(6,6) * t379 + Icges(6,3) * t463;
t333 = Icges(6,5) * t378 + Icges(6,6) * t377 - Icges(6,3) * t464;
t332 = rSges(7,1) * t370 + rSges(7,2) * t369 + rSges(7,3) * t463;
t331 = rSges(7,1) * t368 + rSges(7,2) * t367 - rSges(7,3) * t464;
t330 = Icges(7,1) * t370 + Icges(7,4) * t369 + Icges(7,5) * t463;
t329 = Icges(7,1) * t368 + Icges(7,4) * t367 - Icges(7,5) * t464;
t328 = Icges(7,4) * t370 + Icges(7,2) * t369 + Icges(7,6) * t463;
t327 = Icges(7,4) * t368 + Icges(7,2) * t367 - Icges(7,6) * t464;
t326 = Icges(7,5) * t370 + Icges(7,6) * t369 + Icges(7,3) * t463;
t325 = Icges(7,5) * t368 + Icges(7,6) * t367 - Icges(7,3) * t464;
t324 = qJD(1) * t388 + (-qJD(3) * t410 - qJD(2)) * t436 + t458;
t323 = t410 * t424 + t426 + (-t389 + t451) * qJD(1);
t322 = t452 * t434 - t476 * t436;
t321 = t476 * t434 + t452 * t436;
t320 = qJD(1) * t350 + (-qJD(2) + (-t376 - t409) * qJD(3)) * t436 + t450;
t319 = (qJD(3) * t376 - t456) * t434 + (-t351 + t449) * qJD(1) + t460;
t318 = (t351 * t436 + (-t350 - t395) * t434) * qJD(3) + t461;
t317 = t339 * t414 - t366 * t399 + t439;
t316 = -t340 * t414 + t366 * t398 + t438;
t315 = -t339 * t398 + t340 * t399 + t440;
t314 = t321 * t414 + t331 * t402 - t354 * t399 - t361 * t375 + t439;
t313 = -t322 * t414 - t332 * t402 + t354 * t398 + t361 * t374 + t438;
t312 = -t321 * t398 + t322 * t399 - t331 * t374 + t332 * t375 + t440;
t1 = t402 * ((t325 * t375 + t326 * t374 + t358 * t402) * t433 + ((-t327 * t415 + t329 * t416) * t375 + (-t328 * t415 + t330 * t416) * t374 + (-t359 * t415 + t360 * t416) * t402) * t435) / 0.2e1 + t375 * ((-t325 * t464 + t367 * t327 + t368 * t329) * t375 + (-t326 * t464 + t328 * t367 + t330 * t368) * t374 + (-t358 * t464 + t359 * t367 + t360 * t368) * t402) / 0.2e1 + t399 * ((-t333 * t464 + t377 * t335 + t378 * t337) * t399 + (-t334 * t464 + t336 * t377 + t338 * t378) * t398 + (-t363 * t464 + t364 * t377 + t365 * t378) * t414) / 0.2e1 + t398 * ((t333 * t463 + t335 * t379 + t337 * t380) * t399 + (t334 * t463 + t379 * t336 + t380 * t338) * t398 + (t363 * t463 + t364 * t379 + t365 * t380) * t414) / 0.2e1 + t414 * ((t333 * t399 + t334 * t398 + t363 * t414) * t433 + ((-t335 * t419 + t337 * t420) * t399 + (-t336 * t419 + t338 * t420) * t398 + (-t364 * t419 + t365 * t420) * t414) * t435) / 0.2e1 + t374 * ((t325 * t463 + t327 * t369 + t329 * t370) * t375 + (t326 * t463 + t369 * t328 + t370 * t330) * t374 + (t358 * t463 + t359 * t369 + t360 * t370) * t402) / 0.2e1 + m(6) * (t315 ^ 2 + t316 ^ 2 + t317 ^ 2) / 0.2e1 + m(5) * (t318 ^ 2 + t319 ^ 2 + t320 ^ 2) / 0.2e1 + m(4) * (t323 ^ 2 + t324 ^ 2 + t341 ^ 2) / 0.2e1 + m(7) * (t312 ^ 2 + t313 ^ 2 + t314 ^ 2) / 0.2e1 + m(3) * (t356 ^ 2 + t357 ^ 2) / 0.2e1 + (((t344 * t436 + t345 * t434) * t433 + ((-t346 * t430 + t348 * t431) * t436 + (-t347 * t430 + t349 * t431) * t434) * t435 + (-t383 * t433 + t385 * t435) * t436 + (-t384 * t433 + t386 * t435) * t434) * qJD(3) + ((-t372 * t430 + t373 * t431 + t406) * t435 + (t371 - t405) * t433) * qJD(1)) * qJD(1) / 0.2e1 + (((t344 * t463 + t346 * t393 + t348 * t394 + t443 * t436) * t436 + (t345 * t463 + t347 * t393 + t349 * t394 + (t381 - t442) * t436 + t382 * t434) * t434) * qJD(3) + (t371 * t463 + t372 * t393 + t373 * t394 + t434 * t404 - t436 * t441) * qJD(1)) * t424 / 0.2e1 + (((-t344 * t464 + t346 * t391 + t348 * t392 + t381 * t436) * t436 + (-t345 * t464 + t347 * t391 + t349 * t392 + (t382 - t443) * t436 + t442 * t434) * t434) * qJD(3) + (-t371 * t464 + t372 * t391 + t373 * t392 + t436 * t404 + t434 * t441) * qJD(1)) * t425 / 0.2e1 + (m(2) * (t408 ^ 2 + t411 ^ 2) + Icges(3,1) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
