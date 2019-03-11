% Calculate kinetic energy for
% S6RRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:01:53
% EndTime: 2019-03-09 22:01:56
% DurationCPUTime: 3.24s
% Computational Cost: add. (1675->266), mult. (1701->415), div. (0->0), fcn. (1546->10), ass. (0->158)
t550 = Icges(5,4) + Icges(6,6);
t549 = Icges(5,1) + Icges(6,2);
t548 = -Icges(5,2) - Icges(6,3);
t456 = qJ(2) + qJ(3);
t453 = qJ(4) + t456;
t446 = cos(t453);
t547 = t550 * t446;
t445 = sin(t453);
t546 = t550 * t445;
t545 = Icges(6,4) - Icges(5,5);
t544 = Icges(6,5) - Icges(5,6);
t543 = t548 * t445 + t547;
t542 = t549 * t446 - t546;
t541 = Icges(6,1) + Icges(5,3);
t459 = sin(qJ(1));
t462 = cos(qJ(1));
t540 = t543 * t459 + t544 * t462;
t539 = -t544 * t459 + t543 * t462;
t538 = t542 * t459 + t545 * t462;
t537 = -t545 * t459 + t542 * t462;
t536 = t548 * t446 - t546;
t535 = t549 * t445 + t547;
t534 = t544 * t445 - t545 * t446;
t450 = qJD(2) * t459;
t435 = qJD(3) * t459 + t450;
t423 = qJD(4) * t459 + t435;
t502 = -qJD(2) - qJD(3);
t424 = (-qJD(4) + t502) * t462;
t533 = (-t540 * t445 + t538 * t446) * t424 + (-t539 * t445 + t537 * t446) * t423 + (t536 * t445 + t535 * t446) * qJD(1);
t532 = (t534 * t459 - t541 * t462) * t424 + (t541 * t459 + t534 * t462) * t423 + (-t545 * t445 - t544 * t446) * qJD(1);
t451 = sin(t456);
t528 = pkin(3) * t451;
t527 = pkin(10) * t445;
t461 = cos(qJ(2));
t525 = t461 * pkin(2);
t458 = sin(qJ(2));
t523 = Icges(3,4) * t458;
t522 = Icges(3,4) * t461;
t521 = Icges(4,4) * t451;
t452 = cos(t456);
t520 = Icges(4,4) * t452;
t515 = t446 * t459;
t514 = t446 * t462;
t457 = sin(qJ(6));
t513 = t457 * t459;
t512 = t457 * t462;
t460 = cos(qJ(6));
t511 = t459 * t460;
t510 = t460 * t462;
t388 = -pkin(8) * t462 + t459 * t525;
t389 = pkin(8) * t459 + t462 * t525;
t505 = qJD(2) * t462;
t509 = t388 * t450 + t389 * t505;
t444 = pkin(1) * t459 - pkin(7) * t462;
t508 = -t388 - t444;
t507 = pkin(3) * t452;
t504 = qJD(5) * t445;
t503 = qJD(6) * t446;
t501 = pkin(2) * qJD(2) * t458;
t359 = -pkin(9) * t462 + t459 * t507;
t500 = -t359 + t508;
t499 = t462 * t501;
t493 = pkin(4) * t446 + qJ(5) * t445;
t399 = t493 * t459;
t498 = -t399 + t500;
t497 = rSges(3,1) * t461 - rSges(3,2) * t458;
t496 = rSges(4,1) * t452 - rSges(4,2) * t451;
t495 = rSges(5,1) * t446 - rSges(5,2) * t445;
t494 = -rSges(6,2) * t446 + rSges(6,3) * t445;
t492 = Icges(3,1) * t461 - t523;
t491 = Icges(4,1) * t452 - t521;
t489 = -Icges(3,2) * t458 + t522;
t488 = -Icges(4,2) * t451 + t520;
t485 = Icges(3,5) * t461 - Icges(3,6) * t458;
t484 = Icges(4,5) * t452 - Icges(4,6) * t451;
t403 = -Icges(3,6) * t462 + t459 * t489;
t405 = -Icges(3,5) * t462 + t459 * t492;
t480 = t403 * t458 - t405 * t461;
t404 = Icges(3,6) * t459 + t462 * t489;
t406 = Icges(3,5) * t459 + t462 * t492;
t479 = -t404 * t458 + t406 * t461;
t439 = Icges(3,2) * t461 + t523;
t440 = Icges(3,1) * t458 + t522;
t478 = -t439 * t458 + t440 * t461;
t436 = t502 * t462;
t477 = t436 * t528 - t499;
t360 = pkin(9) * t459 + t462 * t507;
t476 = t435 * t359 - t360 * t436 + t509;
t432 = qJD(1) * (pkin(1) * t462 + pkin(7) * t459);
t475 = qJD(1) * t389 - t459 * t501 + t432;
t420 = pkin(4) * t445 - qJ(5) * t446;
t474 = t424 * t420 + t462 * t504 + t477;
t471 = (Icges(4,5) * t451 + Icges(4,6) * t452) * qJD(1) + (-Icges(4,3) * t462 + t459 * t484) * t436 + (Icges(4,3) * t459 + t462 * t484) * t435;
t470 = -qJD(5) * t446 + t423 * t399 + t476;
t469 = qJD(1) * t360 - t435 * t528 + t475;
t400 = t493 * t462;
t468 = qJD(1) * t400 + t459 * t504 + t469;
t392 = -Icges(4,6) * t462 + t459 * t488;
t393 = Icges(4,6) * t459 + t462 * t488;
t394 = -Icges(4,5) * t462 + t459 * t491;
t395 = Icges(4,5) * t459 + t462 * t491;
t428 = Icges(4,2) * t452 + t521;
t429 = Icges(4,1) * t451 + t520;
t465 = (-t393 * t451 + t395 * t452) * t435 + (-t392 * t451 + t394 * t452) * t436 + (-t428 * t451 + t429 * t452) * qJD(1);
t443 = rSges(2,1) * t462 - rSges(2,2) * t459;
t442 = rSges(2,1) * t459 + rSges(2,2) * t462;
t441 = rSges(3,1) * t458 + rSges(3,2) * t461;
t438 = Icges(3,5) * t458 + Icges(3,6) * t461;
t437 = qJD(6) * t445 + qJD(1);
t430 = rSges(4,1) * t451 + rSges(4,2) * t452;
t426 = -pkin(5) * t462 + pkin(10) * t515;
t425 = pkin(5) * t459 + pkin(10) * t514;
t422 = rSges(5,1) * t445 + rSges(5,2) * t446;
t421 = -rSges(6,2) * t445 - rSges(6,3) * t446;
t412 = t445 * t513 - t510;
t411 = t445 * t511 + t512;
t410 = t445 * t512 + t511;
t409 = t445 * t510 - t513;
t408 = rSges(3,3) * t459 + t462 * t497;
t407 = -rSges(3,3) * t462 + t459 * t497;
t402 = Icges(3,3) * t459 + t462 * t485;
t401 = -Icges(3,3) * t462 + t459 * t485;
t398 = rSges(4,3) * t459 + t462 * t496;
t397 = -rSges(4,3) * t462 + t459 * t496;
t387 = -rSges(6,1) * t462 + t459 * t494;
t386 = rSges(6,1) * t459 + t462 * t494;
t385 = rSges(5,3) * t459 + t462 * t495;
t384 = -rSges(5,3) * t462 + t459 * t495;
t383 = t459 * t503 + t424;
t382 = t462 * t503 + t423;
t368 = rSges(7,3) * t445 + (-rSges(7,1) * t457 - rSges(7,2) * t460) * t446;
t365 = Icges(7,5) * t445 + (-Icges(7,1) * t457 - Icges(7,4) * t460) * t446;
t364 = Icges(7,6) * t445 + (-Icges(7,4) * t457 - Icges(7,2) * t460) * t446;
t363 = Icges(7,3) * t445 + (-Icges(7,5) * t457 - Icges(7,6) * t460) * t446;
t357 = qJD(1) * t408 - t441 * t450 + t432;
t356 = -t441 * t505 + (-t407 - t444) * qJD(1);
t354 = (t407 * t459 + t408 * t462) * qJD(2);
t353 = rSges(7,1) * t412 + rSges(7,2) * t411 + rSges(7,3) * t515;
t352 = rSges(7,1) * t410 + rSges(7,2) * t409 + rSges(7,3) * t514;
t351 = Icges(7,1) * t412 + Icges(7,4) * t411 + Icges(7,5) * t515;
t350 = Icges(7,1) * t410 + Icges(7,4) * t409 + Icges(7,5) * t514;
t349 = Icges(7,4) * t412 + Icges(7,2) * t411 + Icges(7,6) * t515;
t348 = Icges(7,4) * t410 + Icges(7,2) * t409 + Icges(7,6) * t514;
t347 = Icges(7,5) * t412 + Icges(7,6) * t411 + Icges(7,3) * t515;
t346 = Icges(7,5) * t410 + Icges(7,6) * t409 + Icges(7,3) * t514;
t345 = qJD(1) * t398 - t430 * t435 + t475;
t344 = -t499 + t430 * t436 + (-t397 + t508) * qJD(1);
t343 = t397 * t435 - t398 * t436 + t509;
t342 = qJD(1) * t385 - t422 * t423 + t469;
t341 = t422 * t424 + (-t384 + t500) * qJD(1) + t477;
t340 = qJD(1) * t386 + (-t420 - t421) * t423 + t468;
t339 = t421 * t424 + (-t387 + t498) * qJD(1) + t474;
t338 = t384 * t423 - t385 * t424 + t476;
t337 = t387 * t423 + (-t386 - t400) * t424 + t470;
t336 = qJD(1) * t425 + t352 * t437 - t368 * t382 + (-t420 - t527) * t423 + t468;
t335 = t424 * t527 - t353 * t437 + t368 * t383 + (-t426 + t498) * qJD(1) + t474;
t334 = -t352 * t383 + t353 * t382 + t423 * t426 + (-t400 - t425) * t424 + t470;
t1 = ((t459 * t438 + t462 * t478) * qJD(1) + (t459 ^ 2 * t402 + (t480 * t462 + (-t401 + t479) * t459) * t462) * qJD(2)) * t450 / 0.2e1 - ((-t462 * t438 + t459 * t478) * qJD(1) + (t462 ^ 2 * t401 + (t479 * t459 + (-t402 + t480) * t462) * t459) * qJD(2)) * t505 / 0.2e1 + t436 * (t465 * t459 - t471 * t462) / 0.2e1 + t435 * (t471 * t459 + t465 * t462) / 0.2e1 + t383 * ((t346 * t515 + t348 * t411 + t350 * t412) * t382 + (t347 * t515 + t411 * t349 + t412 * t351) * t383 + (t363 * t515 + t364 * t411 + t365 * t412) * t437) / 0.2e1 + t437 * ((t346 * t382 + t347 * t383 + t363 * t437) * t445 + ((-t348 * t460 - t350 * t457) * t382 + (-t349 * t460 - t351 * t457) * t383 + (-t364 * t460 - t365 * t457) * t437) * t446) / 0.2e1 + t382 * ((t346 * t514 + t409 * t348 + t410 * t350) * t382 + (t347 * t514 + t349 * t409 + t351 * t410) * t383 + (t363 * t514 + t364 * t409 + t365 * t410) * t437) / 0.2e1 + m(4) * (t343 ^ 2 + t344 ^ 2 + t345 ^ 2) / 0.2e1 + m(3) * (t354 ^ 2 + t356 ^ 2 + t357 ^ 2) / 0.2e1 + m(7) * (t334 ^ 2 + t335 ^ 2 + t336 ^ 2) / 0.2e1 + m(5) * (t338 ^ 2 + t341 ^ 2 + t342 ^ 2) / 0.2e1 + m(6) * (t337 ^ 2 + t339 ^ 2 + t340 ^ 2) / 0.2e1 + (t532 * t459 + t533 * t462) * t423 / 0.2e1 + (t533 * t459 - t532 * t462) * t424 / 0.2e1 + (Icges(2,3) + m(2) * (t442 ^ 2 + t443 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((t393 * t452 + t395 * t451) * t435 + (t392 * t452 + t394 * t451) * t436 + ((t404 * t461 + t406 * t458) * t459 - (t403 * t461 + t405 * t458) * t462) * qJD(2) + (t538 * t445 + t540 * t446) * t424 + (t537 * t445 + t539 * t446) * t423 + (t452 * t428 + t451 * t429 + t461 * t439 + t458 * t440 + t535 * t445 - t536 * t446) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
