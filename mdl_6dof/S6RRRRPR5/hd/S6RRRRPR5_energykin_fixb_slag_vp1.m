% Calculate kinetic energy for
% S6RRRRPR5
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
% Datum: 2019-03-09 22:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:11:59
% EndTime: 2019-03-09 22:12:02
% DurationCPUTime: 2.76s
% Computational Cost: add. (1865->294), mult. (2561->457), div. (0->0), fcn. (2682->10), ass. (0->154)
t521 = Icges(5,1) + Icges(6,1);
t520 = -Icges(5,4) + Icges(6,5);
t519 = Icges(6,4) + Icges(5,5);
t518 = Icges(5,2) + Icges(6,3);
t517 = -Icges(6,6) + Icges(5,6);
t516 = -Icges(5,3) - Icges(6,2);
t451 = qJ(2) + qJ(3);
t450 = cos(t451);
t457 = cos(qJ(4));
t459 = cos(qJ(1));
t490 = t457 * t459;
t453 = sin(qJ(4));
t455 = sin(qJ(1));
t493 = t453 * t455;
t420 = t450 * t493 + t490;
t491 = t455 * t457;
t492 = t453 * t459;
t421 = t450 * t491 - t492;
t449 = sin(t451);
t495 = t449 * t455;
t515 = t420 * t518 + t421 * t520 - t495 * t517;
t422 = t450 * t492 - t491;
t423 = t450 * t490 + t493;
t494 = t449 * t459;
t514 = t422 * t518 + t423 * t520 - t494 * t517;
t513 = -t420 * t517 + t421 * t519 - t495 * t516;
t512 = -t422 * t517 + t423 * t519 - t494 * t516;
t511 = t420 * t520 + t421 * t521 + t495 * t519;
t510 = t422 * t520 + t423 * t521 + t494 * t519;
t509 = t517 * t450 + (t453 * t518 + t457 * t520) * t449;
t508 = t516 * t450 + (-t453 * t517 + t457 * t519) * t449;
t507 = -t519 * t450 + (t453 * t520 + t457 * t521) * t449;
t458 = cos(qJ(2));
t501 = pkin(2) * t458;
t454 = sin(qJ(2));
t499 = Icges(3,4) * t454;
t498 = Icges(3,4) * t458;
t497 = Icges(4,4) * t449;
t496 = Icges(4,4) * t450;
t391 = -pkin(8) * t459 + t455 * t501;
t392 = pkin(8) * t455 + t459 * t501;
t448 = qJD(2) * t455;
t487 = qJD(2) * t459;
t489 = t391 * t448 + t392 * t487;
t443 = pkin(1) * t455 - pkin(7) * t459;
t488 = -t391 - t443;
t432 = qJD(3) * t455 + t448;
t486 = qJD(4) * t449;
t485 = qJD(6) * t449;
t484 = pkin(2) * qJD(2) * t454;
t411 = t459 * t486 + t432;
t433 = (-qJD(2) - qJD(3)) * t459;
t483 = t459 * t484;
t482 = pkin(3) * t450 + pkin(9) * t449;
t481 = rSges(3,1) * t458 - rSges(3,2) * t454;
t480 = rSges(4,1) * t450 - rSges(4,2) * t449;
t479 = Icges(3,1) * t458 - t499;
t478 = Icges(4,1) * t450 - t497;
t477 = -Icges(3,2) * t454 + t498;
t476 = -Icges(4,2) * t449 + t496;
t475 = Icges(3,5) * t458 - Icges(3,6) * t454;
t474 = Icges(4,5) * t450 - Icges(4,6) * t449;
t407 = -Icges(3,6) * t459 + t455 * t477;
t409 = -Icges(3,5) * t459 + t455 * t479;
t473 = t407 * t454 - t409 * t458;
t408 = Icges(3,6) * t455 + t459 * t477;
t410 = Icges(3,5) * t455 + t459 * t479;
t472 = -t408 * t454 + t410 * t458;
t435 = Icges(3,2) * t458 + t499;
t436 = Icges(3,1) * t454 + t498;
t471 = -t435 * t454 + t436 * t458;
t412 = t455 * t486 + t433;
t418 = t482 * t455;
t419 = t482 * t459;
t470 = t418 * t432 - t419 * t433 + t489;
t431 = qJD(1) * (pkin(1) * t459 + pkin(7) * t455);
t469 = qJD(1) * t392 - t455 * t484 + t431;
t373 = pkin(4) * t421 + qJ(5) * t420;
t468 = qJD(5) * t449 * t453 + t373 * t411 + t470;
t467 = qJD(1) * (Icges(4,5) * t449 + Icges(4,6) * t450) + (-Icges(4,3) * t459 + t455 * t474) * t433 + (Icges(4,3) * t455 + t459 * t474) * t432;
t430 = pkin(3) * t449 - pkin(9) * t450;
t466 = qJD(1) * t419 - t430 * t432 + t469;
t465 = t433 * t430 + (-t418 + t488) * qJD(1) - t483;
t374 = pkin(4) * t423 + qJ(5) * t422;
t444 = -qJD(4) * t450 + qJD(1);
t464 = qJD(5) * t420 + t444 * t374 + t466;
t417 = (pkin(4) * t457 + qJ(5) * t453) * t449;
t463 = qJD(5) * t422 + t412 * t417 + t465;
t396 = -Icges(4,6) * t459 + t455 * t476;
t397 = Icges(4,6) * t455 + t459 * t476;
t398 = -Icges(4,5) * t459 + t455 * t478;
t399 = Icges(4,5) * t455 + t459 * t478;
t427 = Icges(4,2) * t450 + t497;
t428 = Icges(4,1) * t449 + t496;
t462 = (-t397 * t449 + t399 * t450) * t432 + (-t396 * t449 + t398 * t450) * t433 + (-t427 * t449 + t428 * t450) * qJD(1);
t456 = cos(qJ(6));
t452 = sin(qJ(6));
t439 = rSges(2,1) * t459 - rSges(2,2) * t455;
t438 = rSges(2,1) * t455 + rSges(2,2) * t459;
t437 = rSges(3,1) * t454 + rSges(3,2) * t458;
t434 = Icges(3,5) * t454 + Icges(3,6) * t458;
t429 = rSges(4,1) * t449 + rSges(4,2) * t450;
t425 = qJD(1) + (-qJD(4) + qJD(6)) * t450;
t424 = pkin(5) * t449 * t457 + pkin(10) * t450;
t416 = rSges(3,3) * t455 + t459 * t481;
t415 = -rSges(3,3) * t459 + t455 * t481;
t406 = Icges(3,3) * t455 + t459 * t475;
t405 = -Icges(3,3) * t459 + t455 * t475;
t403 = (t452 * t453 + t456 * t457) * t449;
t402 = (-t452 * t457 + t453 * t456) * t449;
t401 = rSges(4,3) * t455 + t459 * t480;
t400 = -rSges(4,3) * t459 + t455 * t480;
t390 = -rSges(5,3) * t450 + (rSges(5,1) * t457 - rSges(5,2) * t453) * t449;
t389 = -rSges(6,2) * t450 + (rSges(6,1) * t457 + rSges(6,3) * t453) * t449;
t381 = pkin(5) * t423 - pkin(10) * t494;
t380 = pkin(5) * t421 - pkin(10) * t495;
t377 = -t455 * t485 + t412;
t376 = -t459 * t485 + t411;
t371 = t422 * t452 + t423 * t456;
t370 = t422 * t456 - t423 * t452;
t369 = t420 * t452 + t421 * t456;
t368 = t420 * t456 - t421 * t452;
t367 = qJD(1) * t416 - t437 * t448 + t431;
t366 = -t437 * t487 + (-t415 - t443) * qJD(1);
t365 = rSges(5,1) * t423 - rSges(5,2) * t422 + rSges(5,3) * t494;
t364 = rSges(6,1) * t423 + rSges(6,2) * t494 + rSges(6,3) * t422;
t363 = rSges(5,1) * t421 - rSges(5,2) * t420 + rSges(5,3) * t495;
t362 = rSges(6,1) * t421 + rSges(6,2) * t495 + rSges(6,3) * t420;
t348 = (t415 * t455 + t416 * t459) * qJD(2);
t347 = rSges(7,1) * t403 + rSges(7,2) * t402 + rSges(7,3) * t450;
t346 = Icges(7,1) * t403 + Icges(7,4) * t402 + Icges(7,5) * t450;
t345 = Icges(7,4) * t403 + Icges(7,2) * t402 + Icges(7,6) * t450;
t344 = Icges(7,5) * t403 + Icges(7,6) * t402 + Icges(7,3) * t450;
t342 = qJD(1) * t401 - t429 * t432 + t469;
t341 = -t483 + t429 * t433 + (-t400 + t488) * qJD(1);
t340 = rSges(7,1) * t371 + rSges(7,2) * t370 - rSges(7,3) * t494;
t339 = rSges(7,1) * t369 + rSges(7,2) * t368 - rSges(7,3) * t495;
t338 = Icges(7,1) * t371 + Icges(7,4) * t370 - Icges(7,5) * t494;
t337 = Icges(7,1) * t369 + Icges(7,4) * t368 - Icges(7,5) * t495;
t336 = Icges(7,4) * t371 + Icges(7,2) * t370 - Icges(7,6) * t494;
t335 = Icges(7,4) * t369 + Icges(7,2) * t368 - Icges(7,6) * t495;
t334 = Icges(7,5) * t371 + Icges(7,6) * t370 - Icges(7,3) * t494;
t333 = Icges(7,5) * t369 + Icges(7,6) * t368 - Icges(7,3) * t495;
t332 = t400 * t432 - t401 * t433 + t489;
t331 = t365 * t444 - t390 * t411 + t466;
t330 = -t363 * t444 + t390 * t412 + t465;
t329 = t363 * t411 - t365 * t412 + t470;
t328 = t364 * t444 + (-t389 - t417) * t411 + t464;
t327 = t389 * t412 + (-t362 - t373) * t444 + t463;
t326 = t362 * t411 + (-t364 - t374) * t412 + t468;
t325 = t340 * t425 - t347 * t376 + t381 * t444 + (-t417 - t424) * t411 + t464;
t324 = -t339 * t425 + t347 * t377 + t412 * t424 + (-t373 - t380) * t444 + t463;
t323 = t339 * t376 - t340 * t377 + t380 * t411 + (-t374 - t381) * t412 + t468;
t1 = t376 * ((-t334 * t494 + t370 * t336 + t371 * t338) * t376 + (-t333 * t494 + t335 * t370 + t337 * t371) * t377 + (-t344 * t494 + t345 * t370 + t346 * t371) * t425) / 0.2e1 + t377 * ((-t334 * t495 + t336 * t368 + t338 * t369) * t376 + (-t333 * t495 + t368 * t335 + t369 * t337) * t377 + (-t344 * t495 + t345 * t368 + t346 * t369) * t425) / 0.2e1 + t425 * ((t334 * t450 + t336 * t402 + t338 * t403) * t376 + (t333 * t450 + t335 * t402 + t337 * t403) * t377 + (t450 * t344 + t402 * t345 + t403 * t346) * t425) / 0.2e1 + t432 * (t467 * t455 + t462 * t459) / 0.2e1 + t433 * (t462 * t455 - t467 * t459) / 0.2e1 + m(6) * (t326 ^ 2 + t327 ^ 2 + t328 ^ 2) / 0.2e1 + m(7) * (t323 ^ 2 + t324 ^ 2 + t325 ^ 2) / 0.2e1 + m(5) * (t329 ^ 2 + t330 ^ 2 + t331 ^ 2) / 0.2e1 + m(4) * (t332 ^ 2 + t341 ^ 2 + t342 ^ 2) / 0.2e1 + m(3) * (t348 ^ 2 + t366 ^ 2 + t367 ^ 2) / 0.2e1 - ((-t434 * t459 + t455 * t471) * qJD(1) + (t459 ^ 2 * t405 + (t472 * t455 + (-t406 + t473) * t459) * t455) * qJD(2)) * t487 / 0.2e1 + ((t434 * t455 + t459 * t471) * qJD(1) + (t455 ^ 2 * t406 + (t473 * t459 + (-t405 + t472) * t455) * t459) * qJD(2)) * t448 / 0.2e1 + ((t422 * t509 + t423 * t507 + t494 * t508) * t444 + (t422 * t515 + t423 * t511 + t494 * t513) * t412 + (t514 * t422 + t510 * t423 + t512 * t494) * t411) * t411 / 0.2e1 + ((t420 * t509 + t421 * t507 + t495 * t508) * t444 + (t515 * t420 + t511 * t421 + t513 * t495) * t412 + (t420 * t514 + t421 * t510 + t495 * t512) * t411) * t412 / 0.2e1 + ((-t411 * t512 - t412 * t513 - t444 * t508) * t450 + ((t453 * t509 + t457 * t507) * t444 + (t453 * t515 + t457 * t511) * t412 + (t453 * t514 + t457 * t510) * t411) * t449) * t444 / 0.2e1 + (m(2) * (t438 ^ 2 + t439 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t397 * t450 + t399 * t449) * t432 + (t396 * t450 + t398 * t449) * t433 + ((t408 * t458 + t410 * t454) * t455 - (t407 * t458 + t409 * t454) * t459) * qJD(2) + (t450 * t427 + t449 * t428 + t458 * t435 + t454 * t436) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
