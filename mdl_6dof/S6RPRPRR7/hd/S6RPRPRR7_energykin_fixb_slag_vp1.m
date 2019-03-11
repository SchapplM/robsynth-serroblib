% Calculate kinetic energy for
% S6RPRPRR7
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
% Datum: 2019-03-09 03:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR7_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR7_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR7_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR7_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:54:55
% EndTime: 2019-03-09 03:54:57
% DurationCPUTime: 1.77s
% Computational Cost: add. (1186->247), mult. (1312->381), div. (0->0), fcn. (1190->10), ass. (0->137)
t494 = Icges(4,3) + Icges(5,3);
t420 = qJ(3) + pkin(10);
t409 = sin(t420);
t410 = cos(t420);
t423 = sin(qJ(3));
t426 = cos(qJ(3));
t493 = Icges(4,5) * t423 + Icges(5,5) * t409 + Icges(4,6) * t426 + Icges(5,6) * t410;
t424 = sin(qJ(1));
t427 = cos(qJ(1));
t492 = t494 * t424 - t493 * t427;
t491 = t493 * t424 + t494 * t427;
t490 = Icges(4,5) * t426 + Icges(5,5) * t410 - Icges(4,6) * t423 - Icges(5,6) * t409;
t474 = Icges(5,4) * t410;
t389 = -Icges(5,2) * t409 + t474;
t475 = Icges(5,4) * t409;
t390 = Icges(5,1) * t410 - t475;
t476 = Icges(4,4) * t426;
t399 = -Icges(4,2) * t423 + t476;
t477 = Icges(4,4) * t423;
t400 = Icges(4,1) * t426 - t477;
t489 = t389 * t410 + t390 * t409 + t399 * t426 + t400 * t423;
t446 = Icges(5,2) * t410 + t475;
t358 = Icges(5,6) * t424 - t427 * t446;
t449 = Icges(5,1) * t409 + t474;
t360 = Icges(5,5) * t424 - t427 * t449;
t447 = Icges(4,2) * t426 + t477;
t372 = Icges(4,6) * t424 - t427 * t447;
t450 = Icges(4,1) * t423 + t476;
t374 = Icges(4,5) * t424 - t427 * t450;
t488 = t358 * t410 + t360 * t409 + t372 * t426 + t374 * t423;
t357 = Icges(5,6) * t427 + t424 * t446;
t359 = Icges(5,5) * t427 + t424 * t449;
t371 = Icges(4,6) * t427 + t424 * t447;
t373 = Icges(4,5) * t427 + t424 * t450;
t487 = -t357 * t410 - t359 * t409 - t371 * t426 - t373 * t423;
t411 = qJ(5) + t420;
t407 = cos(t411);
t406 = sin(t411);
t473 = Icges(6,4) * t406;
t445 = Icges(6,2) * t407 + t473;
t349 = Icges(6,6) * t427 + t424 * t445;
t350 = Icges(6,6) * t424 - t427 * t445;
t472 = Icges(6,4) * t407;
t448 = Icges(6,1) * t406 + t472;
t351 = Icges(6,5) * t427 + t424 * t448;
t352 = Icges(6,5) * t424 - t427 * t448;
t384 = -Icges(6,2) * t406 + t472;
t385 = Icges(6,1) * t407 - t473;
t416 = qJD(3) * t424;
t396 = qJD(5) * t424 + t416;
t417 = qJD(3) * t427;
t397 = qJD(5) * t427 + t417;
t486 = (t349 * t407 + t351 * t406) * t397 + (t350 * t407 + t352 * t406) * t396 + (t384 * t407 + t385 * t406) * qJD(1);
t482 = pkin(3) * t423;
t481 = pkin(3) * t426;
t480 = pkin(4) * t410;
t471 = t407 * t424;
t470 = t407 * t427;
t422 = sin(qJ(6));
t469 = t422 * t424;
t468 = t422 * t427;
t425 = cos(qJ(6));
t467 = t424 * t425;
t466 = t425 * t427;
t392 = qJD(1) * (pkin(1) * t427 + qJ(2) * t424);
t465 = qJD(1) * t427 * pkin(7) + t392;
t463 = qJD(6) * t407;
t418 = qJD(2) * t424;
t462 = qJD(4) * t427 + t416 * t481 + t418;
t459 = pkin(4) * t409;
t401 = pkin(1) * t424 - qJ(2) * t427;
t458 = -pkin(7) * t424 - t401;
t382 = qJ(4) * t427 + t424 * t482;
t457 = qJD(1) * t382 + qJD(4) * t424 + t465;
t456 = t416 * t480 + t462;
t381 = qJ(4) * t424 - t427 * t482;
t455 = -t381 + t458;
t454 = pkin(5) * t406 - pkin(9) * t407;
t453 = rSges(4,1) * t423 + rSges(4,2) * t426;
t452 = rSges(5,1) * t409 + rSges(5,2) * t410;
t451 = rSges(6,1) * t406 + rSges(6,2) * t407;
t442 = Icges(6,5) * t406 + Icges(6,6) * t407;
t339 = pkin(8) * t424 - t427 * t459;
t432 = -t339 + t455;
t340 = pkin(8) * t427 + t424 * t459;
t363 = t381 * t417;
t431 = t339 * t417 + t363 + (-t340 - t382) * t416;
t430 = qJD(1) * (Icges(6,5) * t407 - Icges(6,6) * t406) + (Icges(6,3) * t427 + t424 * t442) * t397 + (Icges(6,3) * t424 - t427 * t442) * t396;
t429 = qJD(1) * t340 + (-qJD(2) + (-t480 - t481) * qJD(3)) * t427 + t457;
t404 = rSges(2,1) * t427 - rSges(2,2) * t424;
t403 = rSges(4,1) * t426 - rSges(4,2) * t423;
t402 = rSges(2,1) * t424 + rSges(2,2) * t427;
t395 = qJD(6) * t406 + qJD(1);
t391 = rSges(5,1) * t410 - rSges(5,2) * t409;
t387 = pkin(5) * t407 + pkin(9) * t406;
t386 = rSges(6,1) * t407 - rSges(6,2) * t406;
t380 = -t406 * t466 + t469;
t379 = t406 * t468 + t467;
t378 = t406 * t467 + t468;
t377 = -t406 * t469 + t466;
t376 = rSges(4,3) * t424 - t427 * t453;
t375 = rSges(4,3) * t427 + t424 * t453;
t367 = -t424 * t463 + t397;
t366 = t427 * t463 + t396;
t365 = t454 * t427;
t364 = t454 * t424;
t362 = rSges(5,3) * t424 - t427 * t452;
t361 = rSges(5,3) * t427 + t424 * t452;
t354 = rSges(6,3) * t424 - t427 * t451;
t353 = rSges(6,3) * t427 + t424 * t451;
t346 = t392 - qJD(2) * t427 + qJD(1) * (-rSges(3,2) * t427 + rSges(3,3) * t424);
t345 = t418 + (rSges(3,2) * t424 + rSges(3,3) * t427 - t401) * qJD(1);
t344 = rSges(7,3) * t406 + (rSges(7,1) * t425 - rSges(7,2) * t422) * t407;
t343 = Icges(7,5) * t406 + (Icges(7,1) * t425 - Icges(7,4) * t422) * t407;
t342 = Icges(7,6) * t406 + (Icges(7,4) * t425 - Icges(7,2) * t422) * t407;
t341 = Icges(7,3) * t406 + (Icges(7,5) * t425 - Icges(7,6) * t422) * t407;
t336 = (-t375 * t424 + t376 * t427) * qJD(3);
t335 = rSges(7,1) * t380 + rSges(7,2) * t379 + rSges(7,3) * t470;
t334 = rSges(7,1) * t378 + rSges(7,2) * t377 - rSges(7,3) * t471;
t333 = Icges(7,1) * t380 + Icges(7,4) * t379 + Icges(7,5) * t470;
t332 = Icges(7,1) * t378 + Icges(7,4) * t377 - Icges(7,5) * t471;
t331 = Icges(7,4) * t380 + Icges(7,2) * t379 + Icges(7,6) * t470;
t330 = Icges(7,4) * t378 + Icges(7,2) * t377 - Icges(7,6) * t471;
t329 = Icges(7,5) * t380 + Icges(7,6) * t379 + Icges(7,3) * t470;
t328 = Icges(7,5) * t378 + Icges(7,6) * t377 - Icges(7,3) * t471;
t327 = qJD(1) * t375 + (-qJD(3) * t403 - qJD(2)) * t427 + t465;
t326 = t403 * t416 + t418 + (-t376 + t458) * qJD(1);
t325 = qJD(1) * t361 + (-qJD(2) + (-t391 - t481) * qJD(3)) * t427 + t457;
t324 = t391 * t416 + (-t362 + t455) * qJD(1) + t462;
t323 = t363 + (t362 * t427 + (-t361 - t382) * t424) * qJD(3);
t322 = qJD(1) * t353 - t386 * t397 + t429;
t321 = t386 * t396 + (-t354 + t432) * qJD(1) + t456;
t320 = -t353 * t396 + t354 * t397 + t431;
t319 = qJD(1) * t364 + t334 * t395 - t344 * t367 - t387 * t397 + t429;
t318 = -t335 * t395 + t344 * t366 + t387 * t396 + (t365 + t432) * qJD(1) + t456;
t317 = -t334 * t366 + t335 * t367 - t364 * t396 - t365 * t397 + t431;
t1 = t367 * ((-t328 * t471 + t377 * t330 + t378 * t332) * t367 + (-t329 * t471 + t331 * t377 + t333 * t378) * t366 + (-t341 * t471 + t342 * t377 + t343 * t378) * t395) / 0.2e1 + t366 * ((t328 * t470 + t330 * t379 + t332 * t380) * t367 + (t329 * t470 + t379 * t331 + t380 * t333) * t366 + (t341 * t470 + t342 * t379 + t343 * t380) * t395) / 0.2e1 + t397 * (t486 * t424 + t430 * t427) / 0.2e1 + t396 * (t430 * t424 - t486 * t427) / 0.2e1 + t395 * ((t328 * t367 + t329 * t366 + t341 * t395) * t406 + ((-t330 * t422 + t332 * t425) * t367 + (-t331 * t422 + t333 * t425) * t366 + (-t342 * t422 + t343 * t425) * t395) * t407) / 0.2e1 + m(7) * (t317 ^ 2 + t318 ^ 2 + t319 ^ 2) / 0.2e1 + m(6) * (t320 ^ 2 + t321 ^ 2 + t322 ^ 2) / 0.2e1 + m(5) * (t323 ^ 2 + t324 ^ 2 + t325 ^ 2) / 0.2e1 + m(4) * (t326 ^ 2 + t327 ^ 2 + t336 ^ 2) / 0.2e1 + m(3) * (t345 ^ 2 + t346 ^ 2) / 0.2e1 + ((t492 * t424 ^ 2 + (t487 * t427 + (-t488 + t491) * t424) * t427) * qJD(3) + (t424 * t490 - t427 * t489) * qJD(1)) * t416 / 0.2e1 + ((t491 * t427 ^ 2 + (t488 * t424 + (-t487 + t492) * t427) * t424) * qJD(3) + (t424 * t489 + t427 * t490) * qJD(1)) * t417 / 0.2e1 + (Icges(2,3) + Icges(3,1) + m(2) * (t402 ^ 2 + t404 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((-t349 * t406 + t351 * t407) * t397 + (-t350 * t406 + t352 * t407) * t396 + ((-t357 * t409 + t359 * t410 - t423 * t371 + t426 * t373) * t427 + (-t358 * t409 + t360 * t410 - t423 * t372 + t374 * t426) * t424) * qJD(3) + (-t406 * t384 + t407 * t385 - t409 * t389 + t410 * t390 - t423 * t399 + t426 * t400) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
