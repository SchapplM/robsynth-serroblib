% Calculate kinetic energy for
% S6RRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:57:22
% EndTime: 2019-03-09 21:57:24
% DurationCPUTime: 2.23s
% Computational Cost: add. (2052->314), mult. (1935->495), div. (0->0), fcn. (1836->12), ass. (0->173)
t468 = qJ(2) + qJ(3);
t462 = sin(t468);
t539 = pkin(3) * t462;
t474 = cos(qJ(2));
t537 = t474 * pkin(2);
t470 = cos(pkin(11));
t536 = pkin(5) * t470;
t472 = sin(qJ(2));
t534 = Icges(3,4) * t472;
t533 = Icges(3,4) * t474;
t532 = Icges(4,4) * t462;
t463 = cos(t468);
t531 = Icges(4,4) * t463;
t464 = qJ(4) + t468;
t454 = sin(t464);
t530 = Icges(5,4) * t454;
t455 = cos(t464);
t529 = Icges(5,4) * t455;
t473 = sin(qJ(1));
t528 = t454 * t473;
t475 = cos(qJ(1));
t527 = t454 * t475;
t466 = pkin(11) + qJ(6);
t457 = sin(t466);
t526 = t457 * t473;
t525 = t457 * t475;
t458 = cos(t466);
t524 = t458 * t473;
t523 = t458 * t475;
t469 = sin(pkin(11));
t522 = t469 * t473;
t521 = t469 * t475;
t520 = t470 * t473;
t519 = t470 * t475;
t398 = -pkin(8) * t475 + t473 * t537;
t399 = pkin(8) * t473 + t475 * t537;
t461 = qJD(2) * t473;
t513 = qJD(2) * t475;
t517 = t398 * t461 + t399 * t513;
t452 = pkin(1) * t473 - pkin(7) * t475;
t516 = -t398 - t452;
t515 = pkin(3) * t463;
t443 = qJD(3) * t473 + t461;
t512 = qJD(5) * t454;
t511 = qJD(6) * t454;
t510 = -qJD(2) - qJD(3);
t509 = pkin(2) * qJD(2) * t472;
t373 = -pkin(9) * t475 + t473 * t515;
t508 = -t373 + t516;
t433 = qJD(4) * t473 + t443;
t507 = t475 * t509;
t502 = pkin(4) * t455 + qJ(5) * t454;
t413 = t502 * t473;
t506 = -t413 + t508;
t505 = rSges(3,1) * t474 - rSges(3,2) * t472;
t504 = rSges(4,1) * t463 - rSges(4,2) * t462;
t503 = rSges(5,1) * t455 - rSges(5,2) * t454;
t434 = (-qJD(4) + t510) * t475;
t501 = Icges(3,1) * t474 - t534;
t500 = Icges(4,1) * t463 - t532;
t499 = Icges(5,1) * t455 - t530;
t498 = -Icges(3,2) * t472 + t533;
t497 = -Icges(4,2) * t462 + t531;
t496 = -Icges(5,2) * t454 + t529;
t495 = Icges(3,5) * t474 - Icges(3,6) * t472;
t494 = Icges(4,5) * t463 - Icges(4,6) * t462;
t493 = Icges(5,5) * t455 - Icges(5,6) * t454;
t417 = -Icges(3,6) * t475 + t473 * t498;
t419 = -Icges(3,5) * t475 + t473 * t501;
t492 = t417 * t472 - t419 * t474;
t418 = Icges(3,6) * t473 + t475 * t498;
t420 = Icges(3,5) * t473 + t475 * t501;
t491 = -t418 * t472 + t420 * t474;
t447 = Icges(3,2) * t474 + t534;
t448 = Icges(3,1) * t472 + t533;
t490 = -t447 * t472 + t448 * t474;
t444 = t510 * t475;
t489 = t444 * t539 - t507;
t374 = pkin(9) * t473 + t475 * t515;
t488 = t443 * t373 - t374 * t444 + t517;
t440 = qJD(1) * (pkin(1) * t475 + pkin(7) * t473);
t487 = qJD(1) * t399 - t473 * t509 + t440;
t431 = pkin(4) * t454 - qJ(5) * t455;
t486 = t434 * t431 + t475 * t512 + t489;
t485 = (Icges(5,5) * t454 + Icges(5,6) * t455) * qJD(1) + (-Icges(5,3) * t475 + t473 * t493) * t434 + (Icges(5,3) * t473 + t475 * t493) * t433;
t484 = (Icges(4,5) * t462 + Icges(4,6) * t463) * qJD(1) + (-Icges(4,3) * t475 + t473 * t494) * t444 + (Icges(4,3) * t473 + t475 * t494) * t443;
t483 = pkin(10) * t454 + t455 * t536;
t482 = -qJD(5) * t455 + t433 * t413 + t488;
t481 = qJD(1) * t374 - t443 * t539 + t487;
t414 = t502 * t475;
t480 = qJD(1) * t414 + t473 * t512 + t481;
t390 = -Icges(5,6) * t475 + t473 * t496;
t391 = Icges(5,6) * t473 + t475 * t496;
t392 = -Icges(5,5) * t475 + t473 * t499;
t393 = Icges(5,5) * t473 + t475 * t499;
t429 = Icges(5,2) * t455 + t530;
t430 = Icges(5,1) * t454 + t529;
t479 = (-t391 * t454 + t393 * t455) * t433 + (-t390 * t454 + t392 * t455) * t434 + (-t429 * t454 + t430 * t455) * qJD(1);
t402 = -Icges(4,6) * t475 + t473 * t497;
t403 = Icges(4,6) * t473 + t475 * t497;
t404 = -Icges(4,5) * t475 + t473 * t500;
t405 = Icges(4,5) * t473 + t475 * t500;
t436 = Icges(4,2) * t463 + t532;
t437 = Icges(4,1) * t462 + t531;
t478 = (-t403 * t462 + t405 * t463) * t443 + (-t402 * t462 + t404 * t463) * t444 + (-t436 * t462 + t437 * t463) * qJD(1);
t451 = rSges(2,1) * t475 - rSges(2,2) * t473;
t450 = rSges(2,1) * t473 + rSges(2,2) * t475;
t449 = rSges(3,1) * t472 + rSges(3,2) * t474;
t446 = Icges(3,5) * t472 + Icges(3,6) * t474;
t445 = -qJD(6) * t455 + qJD(1);
t438 = rSges(4,1) * t462 + rSges(4,2) * t463;
t432 = rSges(5,1) * t454 + rSges(5,2) * t455;
t426 = t455 * t519 + t522;
t425 = -t455 * t521 + t520;
t424 = t455 * t520 - t521;
t423 = -t455 * t522 - t519;
t422 = rSges(3,3) * t473 + t475 * t505;
t421 = -rSges(3,3) * t475 + t473 * t505;
t416 = Icges(3,3) * t473 + t475 * t495;
t415 = -Icges(3,3) * t475 + t473 * t495;
t412 = t455 * t523 + t526;
t411 = -t455 * t525 + t524;
t410 = t455 * t524 - t525;
t409 = -t455 * t526 - t523;
t408 = rSges(4,3) * t473 + t475 * t504;
t407 = -rSges(4,3) * t475 + t473 * t504;
t397 = rSges(5,3) * t473 + t475 * t503;
t396 = -rSges(5,3) * t475 + t473 * t503;
t395 = t473 * t511 + t434;
t394 = t475 * t511 + t433;
t384 = -rSges(6,3) * t455 + (rSges(6,1) * t470 - rSges(6,2) * t469) * t454;
t383 = -Icges(6,5) * t455 + (Icges(6,1) * t470 - Icges(6,4) * t469) * t454;
t382 = -Icges(6,6) * t455 + (Icges(6,4) * t470 - Icges(6,2) * t469) * t454;
t381 = -Icges(6,3) * t455 + (Icges(6,5) * t470 - Icges(6,6) * t469) * t454;
t379 = -rSges(7,3) * t455 + (rSges(7,1) * t458 - rSges(7,2) * t457) * t454;
t378 = -Icges(7,5) * t455 + (Icges(7,1) * t458 - Icges(7,4) * t457) * t454;
t377 = -Icges(7,6) * t455 + (Icges(7,4) * t458 - Icges(7,2) * t457) * t454;
t376 = -Icges(7,3) * t455 + (Icges(7,5) * t458 - Icges(7,6) * t457) * t454;
t372 = -pkin(10) * t455 + t454 * t536;
t370 = qJD(1) * t422 - t449 * t461 + t440;
t369 = -t449 * t513 + (-t421 - t452) * qJD(1);
t367 = (t421 * t473 + t422 * t475) * qJD(2);
t366 = rSges(6,1) * t426 + rSges(6,2) * t425 + rSges(6,3) * t527;
t365 = rSges(6,1) * t424 + rSges(6,2) * t423 + rSges(6,3) * t528;
t364 = Icges(6,1) * t426 + Icges(6,4) * t425 + Icges(6,5) * t527;
t363 = Icges(6,1) * t424 + Icges(6,4) * t423 + Icges(6,5) * t528;
t362 = Icges(6,4) * t426 + Icges(6,2) * t425 + Icges(6,6) * t527;
t361 = Icges(6,4) * t424 + Icges(6,2) * t423 + Icges(6,6) * t528;
t360 = Icges(6,5) * t426 + Icges(6,6) * t425 + Icges(6,3) * t527;
t359 = Icges(6,5) * t424 + Icges(6,6) * t423 + Icges(6,3) * t528;
t358 = pkin(5) * t522 + t475 * t483;
t357 = -pkin(5) * t521 + t473 * t483;
t356 = rSges(7,1) * t412 + rSges(7,2) * t411 + rSges(7,3) * t527;
t355 = rSges(7,1) * t410 + rSges(7,2) * t409 + rSges(7,3) * t528;
t354 = Icges(7,1) * t412 + Icges(7,4) * t411 + Icges(7,5) * t527;
t353 = Icges(7,1) * t410 + Icges(7,4) * t409 + Icges(7,5) * t528;
t352 = Icges(7,4) * t412 + Icges(7,2) * t411 + Icges(7,6) * t527;
t351 = Icges(7,4) * t410 + Icges(7,2) * t409 + Icges(7,6) * t528;
t350 = Icges(7,5) * t412 + Icges(7,6) * t411 + Icges(7,3) * t527;
t349 = Icges(7,5) * t410 + Icges(7,6) * t409 + Icges(7,3) * t528;
t348 = qJD(1) * t408 - t438 * t443 + t487;
t347 = -t507 + t438 * t444 + (-t407 + t516) * qJD(1);
t346 = t407 * t443 - t408 * t444 + t517;
t345 = qJD(1) * t397 - t432 * t433 + t481;
t344 = t432 * t434 + (-t396 + t508) * qJD(1) + t489;
t343 = t396 * t433 - t397 * t434 + t488;
t342 = qJD(1) * t366 + (-t384 - t431) * t433 + t480;
t341 = t384 * t434 + (-t365 + t506) * qJD(1) + t486;
t340 = t365 * t433 + (-t366 - t414) * t434 + t482;
t339 = qJD(1) * t358 + t356 * t445 - t379 * t394 + (-t372 - t431) * t433 + t480;
t338 = -t355 * t445 + t372 * t434 + t379 * t395 + (-t357 + t506) * qJD(1) + t486;
t337 = t355 * t394 - t356 * t395 + t357 * t433 + (-t358 - t414) * t434 + t482;
t1 = m(6) * (t340 ^ 2 + t341 ^ 2 + t342 ^ 2) / 0.2e1 + ((t473 * t446 + t475 * t490) * qJD(1) + (t473 ^ 2 * t416 + (t492 * t475 + (-t415 + t491) * t473) * t475) * qJD(2)) * t461 / 0.2e1 - ((-t475 * t446 + t473 * t490) * qJD(1) + (t475 ^ 2 * t415 + (t491 * t473 + (-t416 + t492) * t475) * t473) * qJD(2)) * t513 / 0.2e1 + m(7) * (t337 ^ 2 + t338 ^ 2 + t339 ^ 2) / 0.2e1 + m(3) * (t367 ^ 2 + t369 ^ 2 + t370 ^ 2) / 0.2e1 + m(4) * (t346 ^ 2 + t347 ^ 2 + t348 ^ 2) / 0.2e1 + m(5) * (t343 ^ 2 + t344 ^ 2 + t345 ^ 2) / 0.2e1 + t443 * (t484 * t473 + t478 * t475) / 0.2e1 + t444 * (t478 * t473 - t484 * t475) / 0.2e1 + t394 * ((t350 * t527 + t411 * t352 + t412 * t354) * t394 + (t349 * t527 + t351 * t411 + t353 * t412) * t395 + (t376 * t527 + t377 * t411 + t378 * t412) * t445) / 0.2e1 + t395 * ((t350 * t528 + t352 * t409 + t354 * t410) * t394 + (t349 * t528 + t409 * t351 + t410 * t353) * t395 + (t376 * t528 + t377 * t409 + t378 * t410) * t445) / 0.2e1 + t445 * ((-t349 * t395 - t350 * t394 - t376 * t445) * t455 + ((-t352 * t457 + t354 * t458) * t394 + (-t351 * t457 + t353 * t458) * t395 + (-t377 * t457 + t378 * t458) * t445) * t454) / 0.2e1 + (t473 * t485 + t475 * t479 + (t360 * t527 + t425 * t362 + t426 * t364) * t433 + (t359 * t527 + t361 * t425 + t363 * t426) * t434 + (t381 * t527 + t382 * t425 + t383 * t426) * qJD(1)) * t433 / 0.2e1 + (t479 * t473 - t485 * t475 + (t360 * t528 + t362 * t423 + t364 * t424) * t433 + (t359 * t528 + t423 * t361 + t424 * t363) * t434 + (t381 * t528 + t382 * t423 + t383 * t424) * qJD(1)) * t434 / 0.2e1 + (Icges(2,3) + m(2) * (t450 ^ 2 + t451 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((t418 * t474 + t420 * t472) * t473 - (t417 * t474 + t419 * t472) * t475) * qJD(2) + (t403 * t463 + t405 * t462) * t443 + (t402 * t463 + t404 * t462) * t444 + (t391 * t455 + t393 * t454) * t433 + (t390 * t455 + t392 * t454) * t434 + (-t359 * t434 - t360 * t433) * t455 + ((-t362 * t469 + t364 * t470) * t433 + (-t361 * t469 + t363 * t470) * t434) * t454 + (t463 * t436 + t462 * t437 + t474 * t447 + t472 * t448 + (t429 - t381) * t455 + (-t382 * t469 + t383 * t470 + t430) * t454) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
