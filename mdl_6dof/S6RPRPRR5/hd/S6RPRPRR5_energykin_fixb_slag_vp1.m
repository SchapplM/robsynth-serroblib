% Calculate kinetic energy for
% S6RPRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:47:56
% EndTime: 2019-03-09 03:47:59
% DurationCPUTime: 2.74s
% Computational Cost: add. (1525->245), mult. (1963->382), div. (0->0), fcn. (2075->10), ass. (0->129)
t518 = Icges(4,4) - Icges(5,5);
t517 = Icges(4,1) + Icges(5,1);
t516 = Icges(4,2) + Icges(5,3);
t440 = pkin(10) + qJ(3);
t436 = cos(t440);
t515 = t518 * t436;
t435 = sin(t440);
t514 = t518 * t435;
t513 = Icges(5,4) + Icges(4,5);
t512 = Icges(4,6) - Icges(5,6);
t511 = t516 * t435 - t515;
t510 = t517 * t436 - t514;
t509 = Icges(5,2) + Icges(4,3);
t446 = sin(qJ(1));
t448 = cos(qJ(1));
t508 = t511 * t446 + t512 * t448;
t507 = -t512 * t446 + t511 * t448;
t506 = -t510 * t446 + t513 * t448;
t505 = t513 * t446 + t510 * t448;
t504 = -t516 * t436 - t514;
t503 = t517 * t435 + t515;
t502 = -t512 * t435 + t513 * t436;
t501 = t502 * t446 - t448 * t509;
t500 = t446 * t509 + t502 * t448;
t499 = t513 * t435 + t512 * t436;
t498 = t435 * t504 + t436 * t503;
t497 = t435 * t507 + t436 * t505;
t496 = -t435 * t508 + t436 * t506;
t492 = cos(qJ(5));
t442 = cos(pkin(10));
t490 = pkin(2) * t442;
t485 = t436 * t446;
t484 = t436 * t448;
t431 = pkin(1) * t446 - qJ(2) * t448;
t482 = pkin(7) * t448 - t446 * t490 - t431;
t439 = qJD(2) * t446;
t479 = qJD(4) * t435;
t481 = t448 * t479 + t439;
t438 = qJD(3) * t446;
t480 = qJD(3) * t448;
t467 = pkin(3) * t436 + qJ(4) * t435;
t406 = t467 * t446;
t478 = -t406 + t482;
t477 = t435 * t492;
t421 = pkin(3) * t435 - qJ(4) * t436;
t474 = qJD(3) * (-rSges(5,1) * t435 + rSges(5,3) * t436 - t421);
t428 = qJD(5) * t448 - t480;
t427 = -qJD(5) * t446 + t438;
t419 = pkin(4) * t485 + pkin(8) * t448;
t473 = -t419 + t478;
t424 = qJD(1) * (pkin(1) * t448 + qJ(2) * t446);
t472 = -qJD(2) * t448 + qJD(1) * (pkin(7) * t446 + t448 * t490) + t424;
t407 = t467 * t448;
t471 = -qJD(4) * t436 + t406 * t438 + t407 * t480;
t441 = sin(pkin(10));
t470 = rSges(3,1) * t442 - rSges(3,2) * t441;
t469 = rSges(4,1) * t436 - rSges(4,2) * t435;
t468 = rSges(5,1) * t436 + rSges(5,3) * t435;
t466 = qJD(3) * (-pkin(4) * t435 - t421);
t445 = sin(qJ(5));
t410 = t435 * t445 + t436 * t492;
t453 = qJD(1) * t407 + t446 * t479 + t472;
t420 = pkin(4) * t484 - pkin(8) * t446;
t452 = t419 * t438 + t420 * t480 + t471;
t451 = t448 * t466 + t481;
t450 = qJD(1) * t420 + t446 * t466 + t453;
t447 = cos(qJ(6));
t444 = sin(qJ(6));
t433 = rSges(2,1) * t448 - rSges(2,2) * t446;
t432 = rSges(2,1) * t446 + rSges(2,2) * t448;
t423 = rSges(4,1) * t435 + rSges(4,2) * t436;
t411 = -t436 * t445 + t477;
t404 = qJD(6) * t410 + qJD(1);
t403 = t410 * t448;
t402 = t445 * t484 - t448 * t477;
t401 = t410 * t446;
t400 = t445 * t485 - t446 * t477;
t399 = rSges(4,3) * t446 + t448 * t469;
t398 = rSges(5,2) * t446 + t448 * t468;
t397 = -rSges(4,3) * t448 + t446 * t469;
t396 = -rSges(5,2) * t448 + t446 * t468;
t379 = t403 * t447 - t444 * t446;
t378 = -t403 * t444 - t446 * t447;
t377 = t401 * t447 + t444 * t448;
t376 = -t401 * t444 + t447 * t448;
t375 = qJD(6) * t400 + t428;
t374 = qJD(6) * t402 + t427;
t373 = pkin(5) * t411 + pkin(9) * t410;
t372 = rSges(6,1) * t411 - rSges(6,2) * t410;
t371 = Icges(6,1) * t411 - Icges(6,4) * t410;
t370 = Icges(6,4) * t411 - Icges(6,2) * t410;
t369 = Icges(6,5) * t411 - Icges(6,6) * t410;
t368 = qJD(1) * t446 * rSges(3,3) + t424 + (qJD(1) * t470 - qJD(2)) * t448;
t367 = t439 + (t448 * rSges(3,3) - t446 * t470 - t431) * qJD(1);
t366 = pkin(5) * t403 + pkin(9) * t402;
t365 = pkin(5) * t401 + pkin(9) * t400;
t364 = rSges(6,1) * t403 - rSges(6,2) * t402 - rSges(6,3) * t446;
t363 = rSges(6,1) * t401 - rSges(6,2) * t400 + rSges(6,3) * t448;
t362 = Icges(6,1) * t403 - Icges(6,4) * t402 - Icges(6,5) * t446;
t361 = Icges(6,1) * t401 - Icges(6,4) * t400 + Icges(6,5) * t448;
t360 = Icges(6,4) * t403 - Icges(6,2) * t402 - Icges(6,6) * t446;
t359 = Icges(6,4) * t401 - Icges(6,2) * t400 + Icges(6,6) * t448;
t358 = Icges(6,5) * t403 - Icges(6,6) * t402 - Icges(6,3) * t446;
t357 = Icges(6,5) * t401 - Icges(6,6) * t400 + Icges(6,3) * t448;
t356 = (t397 * t446 + t399 * t448) * qJD(3);
t355 = rSges(7,3) * t410 + (rSges(7,1) * t447 - rSges(7,2) * t444) * t411;
t354 = Icges(7,5) * t410 + (Icges(7,1) * t447 - Icges(7,4) * t444) * t411;
t353 = Icges(7,6) * t410 + (Icges(7,4) * t447 - Icges(7,2) * t444) * t411;
t352 = Icges(7,3) * t410 + (Icges(7,5) * t447 - Icges(7,6) * t444) * t411;
t351 = qJD(1) * t399 - t423 * t438 + t472;
t350 = -t423 * t480 + t439 + (-t397 + t482) * qJD(1);
t349 = rSges(7,1) * t379 + rSges(7,2) * t378 + rSges(7,3) * t402;
t348 = rSges(7,1) * t377 + rSges(7,2) * t376 + rSges(7,3) * t400;
t347 = Icges(7,1) * t379 + Icges(7,4) * t378 + Icges(7,5) * t402;
t346 = Icges(7,1) * t377 + Icges(7,4) * t376 + Icges(7,5) * t400;
t345 = Icges(7,4) * t379 + Icges(7,2) * t378 + Icges(7,6) * t402;
t344 = Icges(7,4) * t377 + Icges(7,2) * t376 + Icges(7,6) * t400;
t343 = Icges(7,5) * t379 + Icges(7,6) * t378 + Icges(7,3) * t402;
t342 = Icges(7,5) * t377 + Icges(7,6) * t376 + Icges(7,3) * t400;
t341 = (t396 * t446 + t398 * t448) * qJD(3) + t471;
t340 = qJD(1) * t398 + t446 * t474 + t453;
t339 = t448 * t474 + (-t396 + t478) * qJD(1) + t481;
t338 = t363 * t427 - t364 * t428 + t452;
t337 = qJD(1) * t364 - t372 * t427 + t450;
t336 = t372 * t428 + (-t363 + t473) * qJD(1) + t451;
t335 = qJD(1) * t366 + t349 * t404 - t355 * t374 - t373 * t427 + t450;
t334 = -t348 * t404 + t355 * t375 + t373 * t428 + (-t365 + t473) * qJD(1) + t451;
t333 = t348 * t374 - t349 * t375 + t365 * t427 - t366 * t428 + t452;
t1 = t404 * ((t342 * t375 + t343 * t374 + t352 * t404) * t410 + ((-t345 * t444 + t347 * t447) * t374 + (-t344 * t444 + t346 * t447) * t375 + (-t353 * t444 + t354 * t447) * t404) * t411) / 0.2e1 + t374 * ((t402 * t343 + t378 * t345 + t379 * t347) * t374 + (t342 * t402 + t344 * t378 + t346 * t379) * t375 + (t352 * t402 + t353 * t378 + t354 * t379) * t404) / 0.2e1 + t427 * ((-t446 * t358 - t402 * t360 + t403 * t362) * t427 + (-t357 * t446 - t359 * t402 + t361 * t403) * t428 + (-t369 * t446 - t370 * t402 + t371 * t403) * qJD(1)) / 0.2e1 + t428 * ((t358 * t448 - t360 * t400 + t362 * t401) * t427 + (t448 * t357 - t400 * t359 + t401 * t361) * t428 + (t369 * t448 - t370 * t400 + t371 * t401) * qJD(1)) / 0.2e1 + m(7) * (t333 ^ 2 + t334 ^ 2 + t335 ^ 2) / 0.2e1 + m(6) * (t336 ^ 2 + t337 ^ 2 + t338 ^ 2) / 0.2e1 + m(4) * (t350 ^ 2 + t351 ^ 2 + t356 ^ 2) / 0.2e1 + m(5) * (t339 ^ 2 + t340 ^ 2 + t341 ^ 2) / 0.2e1 + m(3) * (t367 ^ 2 + t368 ^ 2) / 0.2e1 + t375 * ((t343 * t400 + t345 * t376 + t347 * t377) * t374 + (t400 * t342 + t376 * t344 + t377 * t346) * t375 + (t352 * t400 + t353 * t376 + t354 * t377) * t404) / 0.2e1 + ((t500 * t446 ^ 2 + (t496 * t448 + (t497 - t501) * t446) * t448) * qJD(3) + (t446 * t499 + t448 * t498) * qJD(1)) * t438 / 0.2e1 - ((t501 * t448 ^ 2 + (t497 * t446 + (t496 - t500) * t448) * t446) * qJD(3) + (t446 * t498 - t448 * t499) * qJD(1)) * t480 / 0.2e1 + (Icges(2,3) + Icges(3,2) * t442 ^ 2 + (Icges(3,1) * t441 + 0.2e1 * Icges(3,4) * t442) * t441 + m(2) * (t432 ^ 2 + t433 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((-t360 * t410 + t362 * t411) * t427 + (-t359 * t410 + t361 * t411) * t428 + ((t435 * t506 + t436 * t508) * t448 + (t435 * t505 - t436 * t507) * t446) * qJD(3) + (-t410 * t370 + t411 * t371 + t503 * t435 - t504 * t436) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
