% Calculate kinetic energy for
% S6RPRRRR5
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
% Datum: 2019-03-09 07:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:08:10
% EndTime: 2019-03-09 07:08:12
% DurationCPUTime: 1.86s
% Computational Cost: add. (1972->278), mult. (1684->451), div. (0->0), fcn. (1618->12), ass. (0->151)
t453 = cos(pkin(11));
t509 = t453 * pkin(2);
t457 = cos(qJ(5));
t508 = pkin(5) * t457;
t450 = pkin(11) + qJ(3);
t440 = sin(t450);
t506 = Icges(4,4) * t440;
t441 = cos(t450);
t505 = Icges(4,4) * t441;
t442 = qJ(4) + t450;
t436 = sin(t442);
t504 = Icges(5,4) * t436;
t437 = cos(t442);
t503 = Icges(5,4) * t437;
t456 = sin(qJ(1));
t502 = t436 * t456;
t458 = cos(qJ(1));
t501 = t436 * t458;
t451 = qJ(5) + qJ(6);
t446 = sin(t451);
t500 = t446 * t456;
t499 = t446 * t458;
t447 = cos(t451);
t498 = t447 * t456;
t497 = t447 * t458;
t455 = sin(qJ(5));
t496 = t455 * t456;
t495 = t455 * t458;
t494 = t456 * t457;
t493 = t457 * t458;
t489 = pkin(3) * t441;
t368 = -pkin(8) * t458 + t456 * t489;
t369 = pkin(8) * t456 + t458 * t489;
t444 = qJD(3) * t456;
t487 = qJD(3) * t458;
t491 = t368 * t444 + t369 * t487;
t433 = pkin(1) * t456 - qJ(2) * t458;
t490 = pkin(7) * t458 - t456 * t509 - t433;
t431 = qJD(4) * t456 + t444;
t486 = qJD(5) * t436;
t485 = qJD(6) * t436;
t484 = pkin(3) * qJD(3) * t440;
t483 = -t368 + t490;
t410 = t458 * t486 + t431;
t432 = (-qJD(3) - qJD(4)) * t458;
t482 = pkin(4) * t437 + pkin(9) * t436;
t427 = qJD(1) * (pkin(1) * t458 + qJ(2) * t456);
t481 = -qJD(2) * t458 + qJD(1) * (pkin(7) * t456 + t458 * t509) + t427;
t452 = sin(pkin(11));
t480 = rSges(3,1) * t453 - rSges(3,2) * t452;
t479 = rSges(4,1) * t441 - rSges(4,2) * t440;
t478 = rSges(5,1) * t437 - rSges(5,2) * t436;
t477 = Icges(4,1) * t441 - t506;
t476 = Icges(5,1) * t437 - t504;
t475 = -Icges(4,2) * t440 + t505;
t474 = -Icges(5,2) * t436 + t503;
t473 = Icges(4,5) * t441 - Icges(4,6) * t440;
t472 = Icges(5,5) * t437 - Icges(5,6) * t436;
t397 = -Icges(4,6) * t458 + t456 * t475;
t399 = -Icges(4,5) * t458 + t456 * t477;
t471 = t397 * t440 - t399 * t441;
t398 = Icges(4,6) * t456 + t458 * t475;
t400 = Icges(4,5) * t456 + t458 * t477;
t470 = -t398 * t440 + t400 * t441;
t423 = Icges(4,2) * t441 + t506;
t424 = Icges(4,1) * t440 + t505;
t469 = -t423 * t440 + t424 * t441;
t411 = t456 * t486 + t432;
t445 = qJD(2) * t456;
t468 = -t458 * t484 + t445;
t408 = t482 * t456;
t409 = t482 * t458;
t467 = t431 * t408 - t409 * t432 + t491;
t466 = pkin(10) * t436 + t437 * t508;
t465 = (Icges(5,5) * t436 + Icges(5,6) * t437) * qJD(1) + (-Icges(5,3) * t458 + t456 * t472) * t432 + (Icges(5,3) * t456 + t458 * t472) * t431;
t464 = qJD(1) * t369 - t456 * t484 + t481;
t421 = pkin(4) * t436 - pkin(9) * t437;
t463 = t432 * t421 + (-t408 + t483) * qJD(1) + t468;
t462 = qJD(1) * t409 - t421 * t431 + t464;
t388 = -Icges(5,6) * t458 + t456 * t474;
t389 = Icges(5,6) * t456 + t458 * t474;
t390 = -Icges(5,5) * t458 + t456 * t476;
t391 = Icges(5,5) * t456 + t458 * t476;
t418 = Icges(5,2) * t437 + t504;
t419 = Icges(5,1) * t436 + t503;
t461 = (-t389 * t436 + t391 * t437) * t431 + (-t388 * t436 + t390 * t437) * t432 + (-t418 * t436 + t419 * t437) * qJD(1);
t435 = rSges(2,1) * t458 - rSges(2,2) * t456;
t434 = rSges(2,1) * t456 + rSges(2,2) * t458;
t430 = -qJD(5) * t437 + qJD(1);
t425 = rSges(4,1) * t440 + rSges(4,2) * t441;
t422 = Icges(4,5) * t440 + Icges(4,6) * t441;
t420 = rSges(5,1) * t436 + rSges(5,2) * t437;
t416 = qJD(1) + (-qJD(5) - qJD(6)) * t437;
t415 = t437 * t493 + t496;
t414 = -t437 * t495 + t494;
t413 = t437 * t494 - t495;
t412 = -t437 * t496 - t493;
t407 = t437 * t497 + t500;
t406 = -t437 * t499 + t498;
t405 = t437 * t498 - t499;
t404 = -t437 * t500 - t497;
t402 = rSges(4,3) * t456 + t458 * t479;
t401 = -rSges(4,3) * t458 + t456 * t479;
t396 = Icges(4,3) * t456 + t458 * t473;
t395 = -Icges(4,3) * t458 + t456 * t473;
t393 = rSges(5,3) * t456 + t458 * t478;
t392 = -rSges(5,3) * t458 + t456 * t478;
t383 = -rSges(6,3) * t437 + (rSges(6,1) * t457 - rSges(6,2) * t455) * t436;
t382 = -Icges(6,5) * t437 + (Icges(6,1) * t457 - Icges(6,4) * t455) * t436;
t381 = -Icges(6,6) * t437 + (Icges(6,4) * t457 - Icges(6,2) * t455) * t436;
t380 = -Icges(6,3) * t437 + (Icges(6,5) * t457 - Icges(6,6) * t455) * t436;
t379 = t456 * t485 + t411;
t378 = t458 * t485 + t410;
t376 = -rSges(7,3) * t437 + (rSges(7,1) * t447 - rSges(7,2) * t446) * t436;
t375 = -Icges(7,5) * t437 + (Icges(7,1) * t447 - Icges(7,4) * t446) * t436;
t374 = -Icges(7,6) * t437 + (Icges(7,4) * t447 - Icges(7,2) * t446) * t436;
t373 = -Icges(7,3) * t437 + (Icges(7,5) * t447 - Icges(7,6) * t446) * t436;
t372 = -pkin(10) * t437 + t436 * t508;
t371 = qJD(1) * t456 * rSges(3,3) + t427 + (qJD(1) * t480 - qJD(2)) * t458;
t370 = t445 + (t458 * rSges(3,3) - t456 * t480 - t433) * qJD(1);
t364 = rSges(6,1) * t415 + rSges(6,2) * t414 + rSges(6,3) * t501;
t363 = rSges(6,1) * t413 + rSges(6,2) * t412 + rSges(6,3) * t502;
t362 = Icges(6,1) * t415 + Icges(6,4) * t414 + Icges(6,5) * t501;
t361 = Icges(6,1) * t413 + Icges(6,4) * t412 + Icges(6,5) * t502;
t360 = Icges(6,4) * t415 + Icges(6,2) * t414 + Icges(6,6) * t501;
t359 = Icges(6,4) * t413 + Icges(6,2) * t412 + Icges(6,6) * t502;
t358 = Icges(6,5) * t415 + Icges(6,6) * t414 + Icges(6,3) * t501;
t357 = Icges(6,5) * t413 + Icges(6,6) * t412 + Icges(6,3) * t502;
t356 = pkin(5) * t496 + t458 * t466;
t355 = -pkin(5) * t495 + t456 * t466;
t354 = (t401 * t456 + t402 * t458) * qJD(3);
t353 = rSges(7,1) * t407 + rSges(7,2) * t406 + rSges(7,3) * t501;
t352 = rSges(7,1) * t405 + rSges(7,2) * t404 + rSges(7,3) * t502;
t351 = Icges(7,1) * t407 + Icges(7,4) * t406 + Icges(7,5) * t501;
t350 = Icges(7,1) * t405 + Icges(7,4) * t404 + Icges(7,5) * t502;
t349 = Icges(7,4) * t407 + Icges(7,2) * t406 + Icges(7,6) * t501;
t348 = Icges(7,4) * t405 + Icges(7,2) * t404 + Icges(7,6) * t502;
t347 = Icges(7,5) * t407 + Icges(7,6) * t406 + Icges(7,3) * t501;
t346 = Icges(7,5) * t405 + Icges(7,6) * t404 + Icges(7,3) * t502;
t345 = qJD(1) * t402 - t425 * t444 + t481;
t344 = -t425 * t487 + t445 + (-t401 + t490) * qJD(1);
t343 = qJD(1) * t393 - t420 * t431 + t464;
t342 = t420 * t432 + (-t392 + t483) * qJD(1) + t468;
t341 = t392 * t431 - t393 * t432 + t491;
t340 = t364 * t430 - t383 * t410 + t462;
t339 = -t363 * t430 + t383 * t411 + t463;
t338 = t363 * t410 - t364 * t411 + t467;
t337 = t353 * t416 + t356 * t430 - t372 * t410 - t376 * t378 + t462;
t336 = -t352 * t416 - t355 * t430 + t372 * t411 + t376 * t379 + t463;
t335 = t352 * t378 - t353 * t379 + t355 * t410 - t356 * t411 + t467;
t1 = ((t456 * t422 + t458 * t469) * qJD(1) + (t456 ^ 2 * t396 + (t471 * t458 + (-t395 + t470) * t456) * t458) * qJD(3)) * t444 / 0.2e1 + m(6) * (t338 ^ 2 + t339 ^ 2 + t340 ^ 2) / 0.2e1 + m(7) * (t335 ^ 2 + t336 ^ 2 + t337 ^ 2) / 0.2e1 - ((-t458 * t422 + t456 * t469) * qJD(1) + (t458 ^ 2 * t395 + (t470 * t456 + (-t396 + t471) * t458) * t456) * qJD(3)) * t487 / 0.2e1 + t431 * (t456 * t465 + t458 * t461) / 0.2e1 + t432 * (t456 * t461 - t458 * t465) / 0.2e1 + t410 * ((t358 * t501 + t414 * t360 + t415 * t362) * t410 + (t357 * t501 + t359 * t414 + t361 * t415) * t411 + (t380 * t501 + t381 * t414 + t382 * t415) * t430) / 0.2e1 + t411 * ((t358 * t502 + t360 * t412 + t362 * t413) * t410 + (t357 * t502 + t412 * t359 + t413 * t361) * t411 + (t380 * t502 + t381 * t412 + t382 * t413) * t430) / 0.2e1 + t430 * ((-t357 * t411 - t358 * t410 - t380 * t430) * t437 + ((-t360 * t455 + t362 * t457) * t410 + (-t359 * t455 + t361 * t457) * t411 + (-t381 * t455 + t382 * t457) * t430) * t436) / 0.2e1 + t378 * ((t347 * t501 + t406 * t349 + t407 * t351) * t378 + (t346 * t501 + t348 * t406 + t350 * t407) * t379 + (t373 * t501 + t374 * t406 + t375 * t407) * t416) / 0.2e1 + t379 * ((t347 * t502 + t349 * t404 + t351 * t405) * t378 + (t346 * t502 + t404 * t348 + t405 * t350) * t379 + (t373 * t502 + t404 * t374 + t405 * t375) * t416) / 0.2e1 + t416 * ((-t346 * t379 - t347 * t378 - t373 * t416) * t437 + ((-t349 * t446 + t351 * t447) * t378 + (-t348 * t446 + t350 * t447) * t379 + (-t374 * t446 + t375 * t447) * t416) * t436) / 0.2e1 + m(3) * (t370 ^ 2 + t371 ^ 2) / 0.2e1 + m(4) * (t344 ^ 2 + t345 ^ 2 + t354 ^ 2) / 0.2e1 + m(5) * (t341 ^ 2 + t342 ^ 2 + t343 ^ 2) / 0.2e1 + (((t398 * t441 + t400 * t440) * t456 - (t397 * t441 + t440 * t399) * t458) * qJD(3) + (t389 * t437 + t391 * t436) * t431 + (t388 * t437 + t390 * t436) * t432 + (t437 * t418 + t436 * t419 + t441 * t423 + t440 * t424) * qJD(1)) * qJD(1) / 0.2e1 + (Icges(2,3) + Icges(3,2) * t453 ^ 2 + (Icges(3,1) * t452 + 0.2e1 * Icges(3,4) * t453) * t452 + m(2) * (t434 ^ 2 + t435 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
