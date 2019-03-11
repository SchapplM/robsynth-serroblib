% Calculate kinetic energy for
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
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
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:37:57
% EndTime: 2019-03-09 02:37:59
% DurationCPUTime: 2.04s
% Computational Cost: add. (1858->267), mult. (1503->399), div. (0->0), fcn. (1436->12), ass. (0->141)
t488 = Icges(4,3) + Icges(5,3);
t412 = qJ(3) + pkin(10);
t406 = sin(t412);
t409 = cos(t412);
t418 = sin(qJ(3));
t420 = cos(qJ(3));
t487 = Icges(4,5) * t420 + Icges(5,5) * t409 - Icges(4,6) * t418 - Icges(5,6) * t406;
t413 = qJ(1) + pkin(9);
t407 = sin(t413);
t410 = cos(t413);
t486 = t487 * t407 - t410 * t488;
t485 = t407 * t488 + t487 * t410;
t484 = Icges(4,5) * t418 + Icges(5,5) * t406 + Icges(4,6) * t420 + Icges(5,6) * t409;
t469 = Icges(5,4) * t406;
t387 = Icges(5,2) * t409 + t469;
t468 = Icges(5,4) * t409;
t388 = Icges(5,1) * t406 + t468;
t471 = Icges(4,4) * t418;
t395 = Icges(4,2) * t420 + t471;
t470 = Icges(4,4) * t420;
t396 = Icges(4,1) * t418 + t470;
t483 = -t387 * t406 + t388 * t409 - t395 * t418 + t396 * t420;
t438 = -Icges(5,2) * t406 + t468;
t349 = Icges(5,6) * t407 + t438 * t410;
t440 = Icges(5,1) * t409 - t469;
t352 = Icges(5,5) * t407 + t440 * t410;
t439 = -Icges(4,2) * t418 + t470;
t369 = Icges(4,6) * t407 + t439 * t410;
t441 = Icges(4,1) * t420 - t471;
t371 = Icges(4,5) * t407 + t441 * t410;
t482 = -t349 * t406 + t352 * t409 - t369 * t418 + t371 * t420;
t348 = -Icges(5,6) * t410 + t438 * t407;
t351 = -Icges(5,5) * t410 + t440 * t407;
t368 = -Icges(4,6) * t410 + t439 * t407;
t370 = -Icges(4,5) * t410 + t441 * t407;
t481 = t348 * t406 - t351 * t409 + t368 * t418 - t370 * t420;
t477 = pkin(3) * t418;
t419 = sin(qJ(1));
t476 = t419 * pkin(1);
t474 = t420 * pkin(3);
t415 = cos(pkin(11));
t473 = t415 * pkin(5);
t467 = t406 * t407;
t466 = t406 * t410;
t465 = t407 * t409;
t414 = sin(pkin(11));
t464 = t407 * t414;
t463 = t407 * t415;
t411 = pkin(11) + qJ(6);
t405 = sin(t411);
t462 = t410 * t405;
t408 = cos(t411);
t461 = t410 * t408;
t460 = t410 * t414;
t459 = t410 * t415;
t421 = cos(qJ(1));
t404 = qJD(1) * t421 * pkin(1);
t457 = qJD(1) * (t410 * pkin(2) + t407 * pkin(7)) + t404;
t401 = qJD(4) * t407;
t453 = qJD(5) * t406;
t456 = t410 * t453 + t401;
t455 = qJD(3) * t407;
t454 = qJD(3) * t410;
t452 = qJD(6) * t406;
t342 = -qJ(4) * t410 + t474 * t407;
t343 = qJ(4) * t407 + t474 * t410;
t451 = t342 * t455 + t343 * t454 + qJD(2);
t448 = -t407 * pkin(2) + t410 * pkin(7) - t476;
t447 = -t406 * pkin(4) + t409 * qJ(5) - t477;
t446 = -t342 + t448;
t445 = rSges(4,1) * t420 - rSges(4,2) * t418;
t444 = rSges(5,1) * t409 - rSges(5,2) * t406;
t443 = pkin(4) * t409 + qJ(5) * t406;
t442 = qJD(3) * (-t406 * rSges(5,1) - t409 * rSges(5,2) - t477);
t375 = t443 * t407;
t429 = -t375 + t446;
t428 = qJD(1) * t343 - qJD(4) * t410 + t457;
t427 = qJD(3) * (pkin(8) * t409 - t473 * t406 + t447);
t426 = qJD(3) * (t409 * rSges(6,3) - (rSges(6,1) * t415 - rSges(6,2) * t414) * t406 + t447);
t376 = t443 * t410;
t425 = qJD(1) * t376 + t407 * t453 + t428;
t424 = -qJD(5) * t409 + t375 * t455 + t376 * t454 + t451;
t423 = pkin(8) * t406 + t473 * t409;
t400 = -qJD(6) * t409 + qJD(1);
t399 = t421 * rSges(2,1) - t419 * rSges(2,2);
t398 = t419 * rSges(2,1) + t421 * rSges(2,2);
t397 = t418 * rSges(4,1) + t420 * rSges(4,2);
t384 = t407 * t452 - t454;
t383 = t410 * t452 + t455;
t382 = t404 + qJD(1) * (t410 * rSges(3,1) - t407 * rSges(3,2));
t381 = (-t407 * rSges(3,1) - t410 * rSges(3,2) - t476) * qJD(1);
t380 = t409 * t459 + t464;
t379 = -t409 * t460 + t463;
t378 = t409 * t463 - t460;
t377 = -t409 * t464 - t459;
t373 = t407 * rSges(4,3) + t445 * t410;
t372 = -t410 * rSges(4,3) + t445 * t407;
t365 = t407 * t405 + t409 * t461;
t364 = t407 * t408 - t409 * t462;
t363 = t408 * t465 - t462;
t362 = -t405 * t465 - t461;
t360 = -Icges(6,5) * t409 + (Icges(6,1) * t415 - Icges(6,4) * t414) * t406;
t359 = -Icges(6,6) * t409 + (Icges(6,4) * t415 - Icges(6,2) * t414) * t406;
t358 = -Icges(6,3) * t409 + (Icges(6,5) * t415 - Icges(6,6) * t414) * t406;
t357 = t407 * rSges(5,3) + t444 * t410;
t356 = -t410 * rSges(5,3) + t444 * t407;
t355 = -t409 * rSges(7,3) + (rSges(7,1) * t408 - rSges(7,2) * t405) * t406;
t350 = -Icges(7,5) * t409 + (Icges(7,1) * t408 - Icges(7,4) * t405) * t406;
t347 = -Icges(7,6) * t409 + (Icges(7,4) * t408 - Icges(7,2) * t405) * t406;
t344 = -Icges(7,3) * t409 + (Icges(7,5) * t408 - Icges(7,6) * t405) * t406;
t337 = t380 * rSges(6,1) + t379 * rSges(6,2) + rSges(6,3) * t466;
t336 = t378 * rSges(6,1) + t377 * rSges(6,2) + rSges(6,3) * t467;
t335 = Icges(6,1) * t380 + Icges(6,4) * t379 + Icges(6,5) * t466;
t334 = Icges(6,1) * t378 + Icges(6,4) * t377 + Icges(6,5) * t467;
t333 = Icges(6,4) * t380 + Icges(6,2) * t379 + Icges(6,6) * t466;
t332 = Icges(6,4) * t378 + Icges(6,2) * t377 + Icges(6,6) * t467;
t331 = Icges(6,5) * t380 + Icges(6,6) * t379 + Icges(6,3) * t466;
t330 = Icges(6,5) * t378 + Icges(6,6) * t377 + Icges(6,3) * t467;
t329 = pkin(5) * t464 + t423 * t410;
t328 = -pkin(5) * t460 + t423 * t407;
t327 = qJD(1) * t373 - t397 * t455 + t457;
t326 = -t397 * t454 + (-t372 + t448) * qJD(1);
t325 = qJD(2) + (t372 * t407 + t373 * t410) * qJD(3);
t324 = t365 * rSges(7,1) + t364 * rSges(7,2) + rSges(7,3) * t466;
t323 = t363 * rSges(7,1) + t362 * rSges(7,2) + rSges(7,3) * t467;
t322 = Icges(7,1) * t365 + Icges(7,4) * t364 + Icges(7,5) * t466;
t321 = Icges(7,1) * t363 + Icges(7,4) * t362 + Icges(7,5) * t467;
t320 = Icges(7,4) * t365 + Icges(7,2) * t364 + Icges(7,6) * t466;
t319 = Icges(7,4) * t363 + Icges(7,2) * t362 + Icges(7,6) * t467;
t318 = Icges(7,5) * t365 + Icges(7,6) * t364 + Icges(7,3) * t466;
t317 = Icges(7,5) * t363 + Icges(7,6) * t362 + Icges(7,3) * t467;
t316 = qJD(1) * t357 + t407 * t442 + t428;
t315 = t401 + t410 * t442 + (-t356 + t446) * qJD(1);
t314 = (t356 * t407 + t357 * t410) * qJD(3) + t451;
t313 = qJD(1) * t337 + t407 * t426 + t425;
t312 = t410 * t426 + (-t336 + t429) * qJD(1) + t456;
t311 = (t336 * t407 + t337 * t410) * qJD(3) + t424;
t310 = qJD(1) * t329 + t400 * t324 - t383 * t355 + t407 * t427 + t425;
t309 = -t400 * t323 + t384 * t355 + t410 * t427 + (-t328 + t429) * qJD(1) + t456;
t308 = t383 * t323 - t384 * t324 + (t328 * t407 + t329 * t410) * qJD(3) + t424;
t1 = t383 * ((t318 * t466 + t364 * t320 + t365 * t322) * t383 + (t317 * t466 + t364 * t319 + t365 * t321) * t384 + (t344 * t466 + t364 * t347 + t365 * t350) * t400) / 0.2e1 + t384 * ((t318 * t467 + t362 * t320 + t363 * t322) * t383 + (t317 * t467 + t362 * t319 + t363 * t321) * t384 + (t344 * t467 + t362 * t347 + t363 * t350) * t400) / 0.2e1 + t400 * ((-t317 * t384 - t318 * t383 - t344 * t400) * t409 + ((-t320 * t405 + t322 * t408) * t383 + (-t319 * t405 + t321 * t408) * t384 + (-t347 * t405 + t350 * t408) * t400) * t406) / 0.2e1 + m(3) * (qJD(2) ^ 2 + t381 ^ 2 + t382 ^ 2) / 0.2e1 + m(4) * (t325 ^ 2 + t326 ^ 2 + t327 ^ 2) / 0.2e1 + m(5) * (t314 ^ 2 + t315 ^ 2 + t316 ^ 2) / 0.2e1 + m(6) * (t311 ^ 2 + t312 ^ 2 + t313 ^ 2) / 0.2e1 + m(7) * (t308 ^ 2 + t309 ^ 2 + t310 ^ 2) / 0.2e1 + (Icges(2,3) + Icges(3,3) + m(2) * (t398 ^ 2 + t399 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((-t420 * t368 - t418 * t370 + (-t348 + t330) * t409 + (t332 * t414 - t334 * t415 - t351) * t406) * t410 + (t420 * t369 + t418 * t371 + (t349 - t331) * t409 + (-t333 * t414 + t335 * t415 + t352) * t406) * t407) * qJD(3) + (t420 * t395 + t418 * t396 + (t387 - t358) * t409 + (-t359 * t414 + t360 * t415 + t388) * t406) * qJD(1)) * qJD(1) / 0.2e1 + (((-t330 * t466 - t379 * t332 - t380 * t334 + t481 * t410) * t410 + (t331 * t466 + t379 * t333 + t380 * t335 + (t482 - t486) * t410 + t485 * t407) * t407) * qJD(3) + (t358 * t466 + t379 * t359 + t380 * t360 + t484 * t407 + t483 * t410) * qJD(1)) * t455 / 0.2e1 - (((t331 * t467 + t377 * t333 + t378 * t335 + t482 * t407) * t407 + (-t330 * t467 - t377 * t332 - t378 * t334 + (t481 - t485) * t407 + t486 * t410) * t410) * qJD(3) + (t358 * t467 + t377 * t359 + t378 * t360 + t483 * t407 - t484 * t410) * qJD(1)) * t454 / 0.2e1;
T  = t1;
