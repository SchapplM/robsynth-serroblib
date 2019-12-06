% Calculate kinetic energy for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:49:36
% EndTime: 2019-12-05 15:49:39
% DurationCPUTime: 2.87s
% Computational Cost: add. (2185->291), mult. (5576->457), div. (0->0), fcn. (7017->12), ass. (0->133)
t435 = sin(pkin(9));
t437 = cos(pkin(9));
t441 = sin(qJ(2));
t443 = cos(qJ(2));
t467 = sin(pkin(10));
t468 = cos(pkin(10));
t425 = -t441 * t467 + t443 * t468;
t438 = cos(pkin(5));
t446 = t438 * t425;
t447 = t441 * t468 + t443 * t467;
t403 = -t435 * t447 + t437 * t446;
t418 = t447 * t438;
t404 = t418 * t437 + t425 * t435;
t436 = sin(pkin(5));
t465 = t437 * t436;
t354 = Icges(4,5) * t404 + Icges(4,6) * t403 - Icges(4,3) * t465;
t405 = -t435 * t446 - t437 * t447;
t406 = -t418 * t435 + t425 * t437;
t466 = t435 * t436;
t355 = Icges(4,5) * t406 + Icges(4,6) * t405 + Icges(4,3) * t466;
t463 = t438 * t443;
t420 = -t435 * t441 + t437 * t463;
t464 = t438 * t441;
t421 = t435 * t443 + t437 * t464;
t390 = Icges(3,5) * t421 + Icges(3,6) * t420 - Icges(3,3) * t465;
t422 = -t435 * t463 - t437 * t441;
t423 = -t435 * t464 + t437 * t443;
t391 = Icges(3,5) * t423 + Icges(3,6) * t422 + Icges(3,3) * t466;
t416 = t425 * t436;
t417 = t447 * t436;
t472 = (Icges(4,5) * t417 + Icges(4,6) * t416 + (Icges(3,5) * t441 + Icges(3,6) * t443) * t436 + (Icges(4,3) + Icges(3,3)) * t438) * t438;
t477 = -t472 + (t390 + t354) * t465 - (t391 + t355) * t466;
t471 = qJD(2) ^ 2;
t470 = cos(qJ(4));
t469 = pkin(2) * t443;
t426 = pkin(2) * t436 * t441 + qJ(3) * t438;
t462 = -pkin(3) * t417 + pkin(7) * t416 - t426;
t461 = qJD(2) * t436;
t430 = t435 * t461;
t381 = -qJD(4) * t405 + t430;
t434 = qJD(2) * t438;
t409 = -qJD(4) * t416 + t434;
t460 = qJD(3) * t437;
t458 = t436 * t470;
t457 = t437 * t461;
t455 = pkin(2) * t464 - qJ(3) * t436;
t452 = (-rSges(4,1) * t417 - rSges(4,2) * t416 - rSges(4,3) * t438 - t426) * t436;
t398 = t435 * t469 + t455 * t437;
t399 = -t455 * t435 + t437 * t469;
t451 = qJD(3) * t438 + t398 * t430 + t399 * t457 + qJD(1);
t382 = -qJD(4) * t403 - t457;
t365 = pkin(3) * t404 - pkin(7) * t403;
t366 = pkin(3) * t406 - pkin(7) * t405;
t448 = t365 * t430 + t366 * t457 + t451;
t385 = t399 * t434;
t445 = t366 * t434 + t385 + (qJD(2) * t435 * t462 - t460) * t436;
t429 = qJD(3) * t466;
t444 = t429 + ((-t365 - t398) * t438 + t462 * t465) * qJD(2);
t442 = cos(qJ(5));
t440 = sin(qJ(4));
t439 = sin(qJ(5));
t415 = t438 * rSges(3,3) + (rSges(3,1) * t441 + rSges(3,2) * t443) * t436;
t414 = Icges(3,5) * t438 + (Icges(3,1) * t441 + Icges(3,4) * t443) * t436;
t413 = Icges(3,6) * t438 + (Icges(3,4) * t441 + Icges(3,2) * t443) * t436;
t408 = t417 * t470 + t438 * t440;
t407 = t417 * t440 - t438 * t470;
t397 = rSges(3,1) * t423 + rSges(3,2) * t422 + rSges(3,3) * t466;
t396 = rSges(3,1) * t421 + rSges(3,2) * t420 - rSges(3,3) * t465;
t395 = Icges(3,1) * t423 + Icges(3,4) * t422 + Icges(3,5) * t466;
t394 = Icges(3,1) * t421 + Icges(3,4) * t420 - Icges(3,5) * t465;
t393 = Icges(3,4) * t423 + Icges(3,2) * t422 + Icges(3,6) * t466;
t392 = Icges(3,4) * t421 + Icges(3,2) * t420 - Icges(3,6) * t465;
t388 = Icges(4,1) * t417 + Icges(4,4) * t416 + Icges(4,5) * t438;
t387 = Icges(4,4) * t417 + Icges(4,2) * t416 + Icges(4,6) * t438;
t380 = t406 * t470 + t440 * t466;
t379 = t406 * t440 - t435 * t458;
t378 = t404 * t470 - t440 * t465;
t377 = t404 * t440 + t437 * t458;
t376 = t408 * t442 - t416 * t439;
t375 = -t408 * t439 - t416 * t442;
t374 = qJD(5) * t407 + t409;
t373 = pkin(4) * t408 + pkin(8) * t407;
t372 = (-t396 * t438 - t415 * t465) * qJD(2);
t371 = (t397 * t438 - t415 * t466) * qJD(2);
t370 = rSges(5,1) * t408 - rSges(5,2) * t407 - rSges(5,3) * t416;
t369 = Icges(5,1) * t408 - Icges(5,4) * t407 - Icges(5,5) * t416;
t368 = Icges(5,4) * t408 - Icges(5,2) * t407 - Icges(5,6) * t416;
t367 = Icges(5,5) * t408 - Icges(5,6) * t407 - Icges(5,3) * t416;
t361 = rSges(4,1) * t406 + rSges(4,2) * t405 + rSges(4,3) * t466;
t360 = rSges(4,1) * t404 + rSges(4,2) * t403 - rSges(4,3) * t465;
t359 = Icges(4,1) * t406 + Icges(4,4) * t405 + Icges(4,5) * t466;
t358 = Icges(4,1) * t404 + Icges(4,4) * t403 - Icges(4,5) * t465;
t357 = Icges(4,4) * t406 + Icges(4,2) * t405 + Icges(4,6) * t466;
t356 = Icges(4,4) * t404 + Icges(4,2) * t403 - Icges(4,6) * t465;
t353 = t380 * t442 - t405 * t439;
t352 = -t380 * t439 - t405 * t442;
t351 = t378 * t442 - t403 * t439;
t350 = -t378 * t439 - t403 * t442;
t349 = qJD(1) + (t396 * t435 + t397 * t437) * t461;
t348 = qJD(5) * t377 + t382;
t347 = qJD(5) * t379 + t381;
t346 = pkin(4) * t380 + pkin(8) * t379;
t345 = pkin(4) * t378 + pkin(8) * t377;
t344 = rSges(6,1) * t376 + rSges(6,2) * t375 + rSges(6,3) * t407;
t343 = Icges(6,1) * t376 + Icges(6,4) * t375 + Icges(6,5) * t407;
t342 = Icges(6,4) * t376 + Icges(6,2) * t375 + Icges(6,6) * t407;
t341 = Icges(6,5) * t376 + Icges(6,6) * t375 + Icges(6,3) * t407;
t340 = rSges(5,1) * t380 - rSges(5,2) * t379 - rSges(5,3) * t405;
t339 = rSges(5,1) * t378 - rSges(5,2) * t377 - rSges(5,3) * t403;
t338 = Icges(5,1) * t380 - Icges(5,4) * t379 - Icges(5,5) * t405;
t337 = Icges(5,1) * t378 - Icges(5,4) * t377 - Icges(5,5) * t403;
t336 = Icges(5,4) * t380 - Icges(5,2) * t379 - Icges(5,6) * t405;
t335 = Icges(5,4) * t378 - Icges(5,2) * t377 - Icges(5,6) * t403;
t334 = Icges(5,5) * t380 - Icges(5,6) * t379 - Icges(5,3) * t405;
t333 = Icges(5,5) * t378 - Icges(5,6) * t377 - Icges(5,3) * t403;
t332 = t429 + ((-t360 - t398) * t438 + t437 * t452) * qJD(2);
t331 = -t436 * t460 + t385 + (t361 * t438 + t435 * t452) * qJD(2);
t330 = rSges(6,1) * t353 + rSges(6,2) * t352 + rSges(6,3) * t379;
t329 = rSges(6,1) * t351 + rSges(6,2) * t350 + rSges(6,3) * t377;
t328 = Icges(6,1) * t353 + Icges(6,4) * t352 + Icges(6,5) * t379;
t327 = Icges(6,1) * t351 + Icges(6,4) * t350 + Icges(6,5) * t377;
t326 = Icges(6,4) * t353 + Icges(6,2) * t352 + Icges(6,6) * t379;
t325 = Icges(6,4) * t351 + Icges(6,2) * t350 + Icges(6,6) * t377;
t324 = Icges(6,5) * t353 + Icges(6,6) * t352 + Icges(6,3) * t379;
t323 = Icges(6,5) * t351 + Icges(6,6) * t350 + Icges(6,3) * t377;
t322 = (t360 * t435 + t361 * t437) * t461 + t451;
t321 = -t339 * t409 + t370 * t382 + t444;
t320 = t340 * t409 - t370 * t381 + t445;
t319 = t339 * t381 - t340 * t382 + t448;
t318 = -t329 * t374 + t344 * t348 - t345 * t409 + t373 * t382 + t444;
t317 = t330 * t374 - t344 * t347 + t346 * t409 - t373 * t381 + t445;
t316 = t329 * t347 - t330 * t348 + t345 * t381 - t346 * t382 + t448;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t349 ^ 2 + t371 ^ 2 + t372 ^ 2) / 0.2e1 + m(4) * (t322 ^ 2 + t331 ^ 2 + t332 ^ 2) / 0.2e1 + m(5) * (t319 ^ 2 + t320 ^ 2 + t321 ^ 2) / 0.2e1 + t381 * ((-t405 * t334 - t379 * t336 + t380 * t338) * t381 + (-t333 * t405 - t335 * t379 + t337 * t380) * t382 + (-t367 * t405 - t368 * t379 + t369 * t380) * t409) / 0.2e1 + t382 * ((-t334 * t403 - t336 * t377 + t338 * t378) * t381 + (-t403 * t333 - t377 * t335 + t378 * t337) * t382 + (-t367 * t403 - t368 * t377 + t369 * t378) * t409) / 0.2e1 + t409 * ((-t334 * t416 - t336 * t407 + t338 * t408) * t381 + (-t333 * t416 - t335 * t407 + t337 * t408) * t382 + (-t416 * t367 - t407 * t368 + t408 * t369) * t409) / 0.2e1 + m(6) * (t316 ^ 2 + t317 ^ 2 + t318 ^ 2) / 0.2e1 + t347 * ((t379 * t324 + t352 * t326 + t353 * t328) * t347 + (t323 * t379 + t325 * t352 + t327 * t353) * t348 + (t341 * t379 + t342 * t352 + t343 * t353) * t374) / 0.2e1 + t348 * ((t324 * t377 + t326 * t350 + t328 * t351) * t347 + (t377 * t323 + t350 * t325 + t351 * t327) * t348 + (t341 * t377 + t342 * t350 + t343 * t351) * t374) / 0.2e1 + t374 * ((t324 * t407 + t326 * t375 + t328 * t376) * t347 + (t323 * t407 + t325 * t375 + t327 * t376) * t348 + (t407 * t341 + t375 * t342 + t376 * t343) * t374) / 0.2e1 - ((t357 * t403 + t359 * t404 + t393 * t420 + t395 * t421) * t466 + (t387 * t403 + t388 * t404 + t413 * t420 + t414 * t421) * t438 + (-t356 * t403 - t358 * t404 - t392 * t420 - t394 * t421 + t477) * t465) * t471 * t465 / 0.2e1 + (((t387 * t416 + t388 * t417 + t472) * t438 + ((-t390 * t437 + t391 * t435 + t413 * t443 + t414 * t441) * t438 + (t355 * t438 + t357 * t416 + t359 * t417) * t435 - (t354 * t438 + t356 * t416 + t358 * t417) * t437 + ((t393 * t443 + t395 * t441) * t435 - (t392 * t443 + t394 * t441) * t437) * t436) * t436) * t438 + ((-t356 * t405 - t358 * t406 - t392 * t422 - t394 * t423) * t465 + (t387 * t405 + t388 * t406 + t413 * t422 + t414 * t423) * t438 + (t357 * t405 + t359 * t406 + t393 * t422 + t395 * t423 - t477) * t466) * t466) * t471 / 0.2e1;
T = t1;
