% Calculate kinetic energy for
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
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
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR10V2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10V2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10V2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR10V2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR10V2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:41:35
% EndTime: 2019-04-11 14:41:37
% DurationCPUTime: 2.32s
% Computational Cost: add. (2331->319), mult. (3295->517), div. (0->0), fcn. (3738->12), ass. (0->160)
t487 = cos(qJ(5));
t442 = sin(qJ(1));
t486 = pkin(1) * t442;
t445 = cos(qJ(2));
t484 = pkin(2) * t445;
t441 = sin(qJ(2));
t483 = Icges(3,4) * t441;
t482 = Icges(3,4) * t445;
t437 = qJ(2) + qJ(3);
t435 = sin(t437);
t481 = Icges(4,4) * t435;
t436 = cos(t437);
t480 = Icges(4,4) * t436;
t440 = sin(qJ(4));
t479 = t435 * t440;
t478 = t435 * t442;
t446 = cos(qJ(1));
t477 = t435 * t446;
t476 = t440 * t442;
t475 = t440 * t446;
t444 = cos(qJ(4));
t474 = t442 * t444;
t473 = t444 * t446;
t409 = t484 * t442;
t410 = t484 * t446;
t434 = qJD(2) * t442;
t471 = qJD(2) * t446;
t472 = t409 * t434 + t410 * t471;
t418 = qJD(3) * t442 + t434;
t470 = qJD(4) * t435;
t469 = pkin(2) * qJD(2) * t441;
t391 = t446 * t470 + t418;
t468 = t435 * t487;
t467 = -t409 - t486;
t419 = (-qJD(2) - qJD(3)) * t446;
t466 = t446 * t469;
t429 = -qJD(4) * t436 + qJD(1);
t406 = t436 * t475 - t474;
t362 = qJD(5) * t406 + t391;
t465 = pkin(3) * t436 + pkin(5) * t435;
t464 = rSges(3,1) * t445 - rSges(3,2) * t441;
t463 = rSges(4,1) * t436 - rSges(4,2) * t435;
t403 = qJD(5) * t479 + t429;
t462 = Icges(3,1) * t445 - t483;
t461 = Icges(4,1) * t436 - t481;
t460 = -Icges(3,2) * t441 + t482;
t459 = -Icges(4,2) * t435 + t480;
t458 = Icges(3,5) * t445 - Icges(3,6) * t441;
t457 = Icges(4,5) * t436 - Icges(4,6) * t435;
t387 = -Icges(3,6) * t446 + t442 * t460;
t389 = -Icges(3,5) * t446 + t442 * t462;
t456 = t387 * t441 - t389 * t445;
t388 = Icges(3,6) * t442 + t446 * t460;
t390 = Icges(3,5) * t442 + t446 * t462;
t455 = -t388 * t441 + t390 * t445;
t421 = Icges(3,2) * t445 + t483;
t422 = Icges(3,1) * t441 + t482;
t454 = -t421 * t441 + t422 * t445;
t392 = t442 * t470 + t419;
t401 = t465 * t442;
t402 = t465 * t446;
t453 = t418 * t401 - t402 * t419 + t472;
t404 = t436 * t476 + t473;
t363 = qJD(5) * t404 + t392;
t432 = qJD(1) * t446 * pkin(1);
t452 = qJD(1) * t410 - t442 * t469 + t432;
t451 = (Icges(4,5) * t435 + Icges(4,6) * t436) * qJD(1) + (-Icges(4,3) * t446 + t442 * t457) * t419 + (Icges(4,3) * t442 + t446 * t457) * t418;
t415 = pkin(3) * t435 - pkin(5) * t436;
t450 = qJD(1) * t402 - t415 * t418 + t452;
t449 = t419 * t415 + (-t401 + t467) * qJD(1) - t466;
t378 = -Icges(4,6) * t446 + t442 * t459;
t379 = Icges(4,6) * t442 + t446 * t459;
t380 = -Icges(4,5) * t446 + t442 * t461;
t381 = Icges(4,5) * t442 + t446 * t461;
t412 = Icges(4,2) * t436 + t481;
t413 = Icges(4,1) * t435 + t480;
t448 = (-t379 * t435 + t381 * t436) * t418 + (-t378 * t435 + t380 * t436) * t419 + (-t412 * t435 + t413 * t436) * qJD(1);
t443 = cos(qJ(6));
t439 = sin(qJ(5));
t438 = sin(qJ(6));
t425 = rSges(2,1) * t446 - rSges(2,2) * t442;
t424 = rSges(2,1) * t442 + rSges(2,2) * t446;
t423 = rSges(3,1) * t441 + rSges(3,2) * t445;
t420 = Icges(3,5) * t441 + Icges(3,6) * t445;
t414 = rSges(4,1) * t435 + rSges(4,2) * t436;
t407 = t436 * t473 + t476;
t405 = t436 * t474 - t475;
t398 = -t436 * t439 + t444 * t468;
t397 = t435 * t444 * t439 + t436 * t487;
t396 = rSges(3,3) * t442 + t446 * t464;
t395 = -rSges(3,3) * t446 + t442 * t464;
t386 = Icges(3,3) * t442 + t446 * t458;
t385 = -Icges(3,3) * t446 + t442 * t458;
t383 = rSges(4,3) * t442 + t446 * t463;
t382 = -rSges(4,3) * t446 + t442 * t463;
t374 = -rSges(5,3) * t436 + (rSges(5,1) * t444 - rSges(5,2) * t440) * t435;
t373 = -Icges(5,5) * t436 + (Icges(5,1) * t444 - Icges(5,4) * t440) * t435;
t372 = -Icges(5,6) * t436 + (Icges(5,4) * t444 - Icges(5,2) * t440) * t435;
t371 = -Icges(5,3) * t436 + (Icges(5,5) * t444 - Icges(5,6) * t440) * t435;
t369 = t407 * t487 + t439 * t477;
t368 = t407 * t439 - t446 * t468;
t367 = t405 * t487 + t439 * t478;
t366 = t405 * t439 - t442 * t468;
t365 = t398 * t443 + t438 * t479;
t364 = -t398 * t438 + t443 * t479;
t361 = qJD(6) * t397 + t403;
t360 = qJD(1) * t396 - t423 * t434 + t432;
t359 = -t423 * t471 + (-t395 - t486) * qJD(1);
t358 = rSges(5,1) * t407 - rSges(5,2) * t406 + rSges(5,3) * t477;
t357 = rSges(5,1) * t405 - rSges(5,2) * t404 + rSges(5,3) * t478;
t356 = Icges(5,1) * t407 - Icges(5,4) * t406 + Icges(5,5) * t477;
t355 = Icges(5,1) * t405 - Icges(5,4) * t404 + Icges(5,5) * t478;
t354 = Icges(5,4) * t407 - Icges(5,2) * t406 + Icges(5,6) * t477;
t353 = Icges(5,4) * t405 - Icges(5,2) * t404 + Icges(5,6) * t478;
t352 = Icges(5,5) * t407 - Icges(5,6) * t406 + Icges(5,3) * t477;
t351 = Icges(5,5) * t405 - Icges(5,6) * t404 + Icges(5,3) * t478;
t350 = (t395 * t442 + t396 * t446) * qJD(2);
t349 = rSges(6,1) * t398 - rSges(6,2) * t397 + rSges(6,3) * t479;
t348 = Icges(6,1) * t398 - Icges(6,4) * t397 + Icges(6,5) * t479;
t347 = Icges(6,4) * t398 - Icges(6,2) * t397 + Icges(6,6) * t479;
t346 = Icges(6,5) * t398 - Icges(6,6) * t397 + Icges(6,3) * t479;
t345 = t369 * t443 + t406 * t438;
t344 = -t369 * t438 + t406 * t443;
t343 = t367 * t443 + t404 * t438;
t342 = -t367 * t438 + t404 * t443;
t341 = qJD(6) * t366 + t363;
t340 = qJD(6) * t368 + t362;
t339 = qJD(1) * t383 - t414 * t418 + t452;
t338 = -t466 + t414 * t419 + (-t382 + t467) * qJD(1);
t337 = rSges(6,1) * t369 - rSges(6,2) * t368 + rSges(6,3) * t406;
t336 = rSges(6,1) * t367 - rSges(6,2) * t366 + rSges(6,3) * t404;
t335 = Icges(6,1) * t369 - Icges(6,4) * t368 + Icges(6,5) * t406;
t334 = Icges(6,1) * t367 - Icges(6,4) * t366 + Icges(6,5) * t404;
t333 = Icges(6,4) * t369 - Icges(6,2) * t368 + Icges(6,6) * t406;
t332 = Icges(6,4) * t367 - Icges(6,2) * t366 + Icges(6,6) * t404;
t331 = Icges(6,5) * t369 - Icges(6,6) * t368 + Icges(6,3) * t406;
t330 = Icges(6,5) * t367 - Icges(6,6) * t366 + Icges(6,3) * t404;
t329 = rSges(7,1) * t365 + rSges(7,2) * t364 + rSges(7,3) * t397;
t328 = Icges(7,1) * t365 + Icges(7,4) * t364 + Icges(7,5) * t397;
t327 = Icges(7,4) * t365 + Icges(7,2) * t364 + Icges(7,6) * t397;
t326 = Icges(7,5) * t365 + Icges(7,6) * t364 + Icges(7,3) * t397;
t325 = t382 * t418 - t383 * t419 + t472;
t324 = rSges(7,1) * t345 + rSges(7,2) * t344 + rSges(7,3) * t368;
t323 = rSges(7,1) * t343 + rSges(7,2) * t342 + rSges(7,3) * t366;
t322 = Icges(7,1) * t345 + Icges(7,4) * t344 + Icges(7,5) * t368;
t321 = Icges(7,1) * t343 + Icges(7,4) * t342 + Icges(7,5) * t366;
t320 = Icges(7,4) * t345 + Icges(7,2) * t344 + Icges(7,6) * t368;
t319 = Icges(7,4) * t343 + Icges(7,2) * t342 + Icges(7,6) * t366;
t318 = Icges(7,5) * t345 + Icges(7,6) * t344 + Icges(7,3) * t368;
t317 = Icges(7,5) * t343 + Icges(7,6) * t342 + Icges(7,3) * t366;
t316 = t358 * t429 - t374 * t391 + t450;
t315 = -t357 * t429 + t374 * t392 + t449;
t314 = t357 * t391 - t358 * t392 + t453;
t313 = t337 * t403 - t349 * t362 + t450;
t312 = -t336 * t403 + t349 * t363 + t449;
t311 = t336 * t362 - t337 * t363 + t453;
t310 = t324 * t361 - t329 * t340 + (-t362 * t397 + t368 * t403) * pkin(6) + t450;
t309 = -t323 * t361 + t329 * t341 + (t363 * t397 - t366 * t403) * pkin(6) + t449;
t308 = t323 * t340 - t324 * t341 + (t362 * t366 - t363 * t368) * pkin(6) + t453;
t1 = m(5) * (t314 ^ 2 + t315 ^ 2 + t316 ^ 2) / 0.2e1 + m(3) * (t350 ^ 2 + t359 ^ 2 + t360 ^ 2) / 0.2e1 + t340 * ((t318 * t368 + t320 * t344 + t322 * t345) * t340 + (t317 * t368 + t319 * t344 + t321 * t345) * t341 + (t326 * t368 + t327 * t344 + t328 * t345) * t361) / 0.2e1 + t403 * ((t331 * t479 - t333 * t397 + t335 * t398) * t362 + (t330 * t479 - t332 * t397 + t334 * t398) * t363 + (t346 * t479 - t347 * t397 + t348 * t398) * t403) / 0.2e1 + t341 * ((t318 * t366 + t320 * t342 + t322 * t343) * t340 + (t317 * t366 + t319 * t342 + t321 * t343) * t341 + (t326 * t366 + t327 * t342 + t328 * t343) * t361) / 0.2e1 + t361 * ((t318 * t397 + t320 * t364 + t322 * t365) * t340 + (t317 * t397 + t319 * t364 + t321 * t365) * t341 + (t326 * t397 + t327 * t364 + t328 * t365) * t361) / 0.2e1 + t362 * ((t331 * t406 - t333 * t368 + t335 * t369) * t362 + (t330 * t406 - t332 * t368 + t334 * t369) * t363 + (t346 * t406 - t347 * t368 + t348 * t369) * t403) / 0.2e1 + t363 * ((t331 * t404 - t333 * t366 + t335 * t367) * t362 + (t330 * t404 - t332 * t366 + t334 * t367) * t363 + (t346 * t404 - t347 * t366 + t348 * t367) * t403) / 0.2e1 + t419 * (t448 * t442 - t451 * t446) / 0.2e1 + t392 * ((t352 * t478 - t354 * t404 + t356 * t405) * t391 + (t351 * t478 - t404 * t353 + t405 * t355) * t392 + (t371 * t478 - t372 * t404 + t373 * t405) * t429) / 0.2e1 + t429 * ((-t351 * t392 - t352 * t391 - t371 * t429) * t436 + ((-t354 * t440 + t356 * t444) * t391 + (-t353 * t440 + t355 * t444) * t392 + (-t372 * t440 + t373 * t444) * t429) * t435) / 0.2e1 + t391 * ((t352 * t477 - t354 * t406 + t356 * t407) * t391 + (t351 * t477 - t353 * t406 + t355 * t407) * t392 + (t371 * t477 - t372 * t406 + t373 * t407) * t429) / 0.2e1 + t418 * (t451 * t442 + t448 * t446) / 0.2e1 + m(7) * (t308 ^ 2 + t309 ^ 2 + t310 ^ 2) / 0.2e1 + m(6) * (t311 ^ 2 + t312 ^ 2 + t313 ^ 2) / 0.2e1 + m(4) * (t325 ^ 2 + t338 ^ 2 + t339 ^ 2) / 0.2e1 - ((-t446 * t420 + t442 * t454) * qJD(1) + (t446 ^ 2 * t385 + (t455 * t442 + (-t386 + t456) * t446) * t442) * qJD(2)) * t471 / 0.2e1 + ((t442 * t420 + t454 * t446) * qJD(1) + (t442 ^ 2 * t386 + (t456 * t446 + (-t385 + t455) * t442) * t446) * qJD(2)) * t434 / 0.2e1 + (m(2) * (t424 ^ 2 + t425 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t388 * t445 + t390 * t441) * t442 - (t387 * t445 + t389 * t441) * t446) * qJD(2) + (t379 * t436 + t381 * t435) * t418 + (t378 * t436 + t380 * t435) * t419 + (t436 * t412 + t435 * t413 + t445 * t421 + t441 * t422) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
