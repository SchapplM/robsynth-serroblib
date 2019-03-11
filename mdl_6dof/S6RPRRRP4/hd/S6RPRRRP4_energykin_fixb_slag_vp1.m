% Calculate kinetic energy for
% S6RPRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:06:56
% EndTime: 2019-03-09 06:06:58
% DurationCPUTime: 1.96s
% Computational Cost: add. (1787->233), mult. (1659->361), div. (0->0), fcn. (1593->10), ass. (0->132)
t517 = Icges(6,1) + Icges(7,1);
t516 = Icges(6,4) + Icges(7,4);
t515 = -Icges(7,5) - Icges(6,5);
t514 = Icges(6,2) + Icges(7,2);
t513 = -Icges(7,6) - Icges(6,6);
t512 = -Icges(7,3) - Icges(6,3);
t438 = pkin(10) + qJ(3);
t432 = qJ(4) + t438;
t427 = cos(t432);
t445 = cos(qJ(5));
t446 = cos(qJ(1));
t483 = t445 * t446;
t443 = sin(qJ(5));
t444 = sin(qJ(1));
t486 = t443 * t444;
t405 = -t427 * t486 - t483;
t484 = t444 * t445;
t485 = t443 * t446;
t406 = t427 * t484 - t485;
t426 = sin(t432);
t488 = t426 * t444;
t511 = -t513 * t405 - t515 * t406 - t512 * t488;
t407 = -t427 * t485 + t484;
t408 = t427 * t483 + t486;
t487 = t426 * t446;
t510 = -t513 * t407 - t515 * t408 - t512 * t487;
t509 = t514 * t405 + t516 * t406 - t513 * t488;
t508 = t514 * t407 + t516 * t408 - t513 * t487;
t507 = t516 * t405 + t517 * t406 - t515 * t488;
t506 = t516 * t407 + t517 * t408 - t515 * t487;
t505 = t512 * t427 + (t513 * t443 - t515 * t445) * t426;
t504 = t513 * t427 + (-t514 * t443 + t516 * t445) * t426;
t503 = t515 * t427 + (-t516 * t443 + t517 * t445) * t426;
t469 = pkin(4) * t427 + pkin(9) * t426;
t402 = t469 * t446;
t413 = t426 * pkin(4) - t427 * pkin(9);
t434 = qJD(3) * t444;
t421 = qJD(4) * t444 + t434;
t502 = qJD(1) * t402 - t421 * t413;
t440 = cos(pkin(10));
t496 = t440 * pkin(2);
t495 = pkin(5) * t445;
t430 = sin(t438);
t493 = Icges(4,4) * t430;
t431 = cos(t438);
t492 = Icges(4,4) * t431;
t491 = Icges(5,4) * t426;
t490 = Icges(5,4) * t427;
t452 = qJ(6) * t426 + t427 * t495;
t481 = rSges(7,1) * t406 + rSges(7,2) * t405 + rSges(7,3) * t488 - pkin(5) * t485 + t444 * t452;
t480 = rSges(7,1) * t408 + rSges(7,2) * t407 + rSges(7,3) * t487 + pkin(5) * t486 + t446 * t452;
t476 = pkin(3) * t431;
t367 = -pkin(8) * t446 + t444 * t476;
t368 = pkin(8) * t444 + t446 * t476;
t474 = qJD(3) * t446;
t479 = t367 * t434 + t368 * t474;
t478 = (-qJ(6) - rSges(7,3)) * t427 + (rSges(7,1) * t445 - rSges(7,2) * t443 + t495) * t426;
t423 = t444 * pkin(1) - t446 * qJ(2);
t477 = pkin(7) * t446 - t444 * t496 - t423;
t473 = qJD(5) * t426;
t472 = pkin(3) * qJD(3) * t430;
t471 = -t367 + t477;
t422 = (-qJD(3) - qJD(4)) * t446;
t470 = t446 * t472;
t419 = qJD(1) * (t446 * pkin(1) + t444 * qJ(2));
t468 = -qJD(2) * t446 + qJD(1) * (pkin(7) * t444 + t446 * t496) + t419;
t439 = sin(pkin(10));
t467 = rSges(3,1) * t440 - rSges(3,2) * t439;
t466 = rSges(4,1) * t431 - rSges(4,2) * t430;
t465 = rSges(5,1) * t427 - rSges(5,2) * t426;
t464 = Icges(4,1) * t431 - t493;
t463 = Icges(5,1) * t427 - t491;
t462 = -Icges(4,2) * t430 + t492;
t461 = -Icges(5,2) * t426 + t490;
t460 = Icges(4,5) * t431 - Icges(4,6) * t430;
t459 = Icges(5,5) * t427 - Icges(5,6) * t426;
t394 = -Icges(4,6) * t446 + t444 * t462;
t396 = -Icges(4,5) * t446 + t444 * t464;
t458 = t394 * t430 - t396 * t431;
t395 = Icges(4,6) * t444 + t446 * t462;
t397 = Icges(4,5) * t444 + t446 * t464;
t457 = -t395 * t430 + t397 * t431;
t415 = Icges(4,2) * t431 + t493;
t416 = Icges(4,1) * t430 + t492;
t456 = -t415 * t430 + t416 * t431;
t401 = t469 * t444;
t455 = t421 * t401 - t422 * t402 + t479;
t454 = qJD(1) * t368 + t468;
t453 = qJD(6) * t426 - t472;
t451 = qJD(1) * (Icges(5,5) * t426 + Icges(5,6) * t427) + (-Icges(5,3) * t446 + t444 * t459) * t422 + (Icges(5,3) * t444 + t446 * t459) * t421;
t435 = qJD(2) * t444;
t450 = t422 * t413 + t435 + (-t401 + t471) * qJD(1);
t449 = -t444 * t472 + t454;
t385 = -Icges(5,6) * t446 + t444 * t461;
t386 = Icges(5,6) * t444 + t446 * t461;
t387 = -Icges(5,5) * t446 + t444 * t463;
t388 = Icges(5,5) * t444 + t446 * t463;
t410 = Icges(5,2) * t427 + t491;
t411 = Icges(5,1) * t426 + t490;
t448 = (-t386 * t426 + t388 * t427) * t421 + (-t385 * t426 + t387 * t427) * t422 + (-t410 * t426 + t411 * t427) * qJD(1);
t425 = rSges(2,1) * t446 - rSges(2,2) * t444;
t424 = rSges(2,1) * t444 + rSges(2,2) * t446;
t420 = -qJD(5) * t427 + qJD(1);
t417 = rSges(4,1) * t430 + rSges(4,2) * t431;
t414 = Icges(4,5) * t430 + Icges(4,6) * t431;
t412 = rSges(5,1) * t426 + rSges(5,2) * t427;
t404 = t444 * t473 + t422;
t403 = t446 * t473 + t421;
t399 = t444 * rSges(4,3) + t446 * t466;
t398 = -t446 * rSges(4,3) + t444 * t466;
t393 = Icges(4,3) * t444 + t446 * t460;
t392 = -Icges(4,3) * t446 + t444 * t460;
t390 = t444 * rSges(5,3) + t446 * t465;
t389 = -t446 * rSges(5,3) + t444 * t465;
t380 = -t427 * rSges(6,3) + (rSges(6,1) * t445 - rSges(6,2) * t443) * t426;
t370 = qJD(1) * t444 * rSges(3,3) + t419 + (qJD(1) * t467 - qJD(2)) * t446;
t369 = t435 + (t446 * rSges(3,3) - t444 * t467 - t423) * qJD(1);
t363 = rSges(6,1) * t408 + rSges(6,2) * t407 + rSges(6,3) * t487;
t361 = rSges(6,1) * t406 + rSges(6,2) * t405 + rSges(6,3) * t488;
t345 = (t398 * t444 + t399 * t446) * qJD(3);
t344 = qJD(1) * t399 - t417 * t434 + t468;
t343 = -t417 * t474 + t435 + (-t398 + t477) * qJD(1);
t342 = qJD(1) * t390 - t412 * t421 + t449;
t341 = -t470 + t422 * t412 + t435 + (-t389 + t471) * qJD(1);
t340 = t389 * t421 - t390 * t422 + t479;
t339 = t363 * t420 - t380 * t403 + t449 + t502;
t338 = -t420 * t361 + t404 * t380 + t450 - t470;
t337 = t361 * t403 - t363 * t404 + t455;
t336 = -t403 * t478 + t420 * t480 + t444 * t453 + t454 + t502;
t335 = t404 * t478 - t420 * t481 + t446 * t453 + t450;
t334 = -qJD(6) * t427 + t403 * t481 - t404 * t480 + t455;
t1 = m(5) * (t340 ^ 2 + t341 ^ 2 + t342 ^ 2) / 0.2e1 + m(6) * (t337 ^ 2 + t338 ^ 2 + t339 ^ 2) / 0.2e1 + m(7) * (t334 ^ 2 + t335 ^ 2 + t336 ^ 2) / 0.2e1 - ((-t446 * t414 + t444 * t456) * qJD(1) + (t446 ^ 2 * t392 + (t457 * t444 + (-t393 + t458) * t446) * t444) * qJD(3)) * t474 / 0.2e1 + t421 * (t451 * t444 + t448 * t446) / 0.2e1 + t422 * (t448 * t444 - t451 * t446) / 0.2e1 + m(3) * (t369 ^ 2 + t370 ^ 2) / 0.2e1 + m(4) * (t343 ^ 2 + t344 ^ 2 + t345 ^ 2) / 0.2e1 + ((t444 * t414 + t446 * t456) * qJD(1) + (t444 ^ 2 * t393 + (t458 * t446 + (-t392 + t457) * t444) * t446) * qJD(3)) * t434 / 0.2e1 + ((t504 * t407 + t503 * t408 + t505 * t487) * t420 + (t509 * t407 + t507 * t408 + t511 * t487) * t404 + (t508 * t407 + t506 * t408 + t510 * t487) * t403) * t403 / 0.2e1 + ((t504 * t405 + t503 * t406 + t505 * t488) * t420 + (t509 * t405 + t507 * t406 + t511 * t488) * t404 + (t508 * t405 + t506 * t406 + t510 * t488) * t403) * t404 / 0.2e1 + ((-t510 * t403 - t511 * t404 - t505 * t420) * t427 + ((-t504 * t443 + t503 * t445) * t420 + (-t509 * t443 + t507 * t445) * t404 + (-t508 * t443 + t506 * t445) * t403) * t426) * t420 / 0.2e1 + (((t431 * t395 + t430 * t397) * t444 - (t394 * t431 + t396 * t430) * t446) * qJD(3) + (t386 * t427 + t388 * t426) * t421 + (t385 * t427 + t387 * t426) * t422 + (t427 * t410 + t426 * t411 + t415 * t431 + t430 * t416) * qJD(1)) * qJD(1) / 0.2e1 + (Icges(2,3) + Icges(3,2) * t440 ^ 2 + (Icges(3,1) * t439 + 0.2e1 * Icges(3,4) * t440) * t439 + m(2) * (t424 ^ 2 + t425 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
