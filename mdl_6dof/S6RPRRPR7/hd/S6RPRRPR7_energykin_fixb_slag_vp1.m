% Calculate kinetic energy for
% S6RPRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 05:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR7_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR7_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR7_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR7_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:19:50
% EndTime: 2019-03-09 05:19:52
% DurationCPUTime: 1.90s
% Computational Cost: add. (1204->246), mult. (1330->378), div. (0->0), fcn. (1208->10), ass. (0->137)
t492 = Icges(5,3) + Icges(6,3);
t420 = qJ(3) + qJ(4);
t409 = pkin(10) + t420;
t406 = sin(t409);
t407 = cos(t409);
t417 = sin(t420);
t418 = cos(t420);
t491 = Icges(5,5) * t417 + Icges(6,5) * t406 + Icges(5,6) * t418 + Icges(6,6) * t407;
t423 = sin(qJ(1));
t426 = cos(qJ(1));
t474 = Icges(6,4) * t406;
t447 = Icges(6,2) * t407 + t474;
t349 = Icges(6,6) * t426 + t423 * t447;
t350 = Icges(6,6) * t423 - t426 * t447;
t473 = Icges(6,4) * t407;
t450 = Icges(6,1) * t406 + t473;
t351 = Icges(6,5) * t426 + t423 * t450;
t352 = Icges(6,5) * t423 - t426 * t450;
t476 = Icges(5,4) * t417;
t448 = Icges(5,2) * t418 + t476;
t357 = Icges(5,6) * t426 + t423 * t448;
t358 = Icges(5,6) * t423 - t426 * t448;
t475 = Icges(5,4) * t418;
t451 = Icges(5,1) * t417 + t475;
t359 = Icges(5,5) * t426 + t423 * t451;
t360 = Icges(5,5) * t423 - t426 * t451;
t384 = -Icges(6,2) * t406 + t473;
t385 = Icges(6,1) * t407 - t474;
t390 = -Icges(5,2) * t417 + t475;
t391 = Icges(5,1) * t418 - t476;
t414 = qJD(3) * t423;
t396 = qJD(4) * t423 + t414;
t415 = qJD(3) * t426;
t397 = qJD(4) * t426 + t415;
t490 = qJD(1) * (t384 * t407 + t385 * t406 + t390 * t418 + t391 * t417) + t396 * (t350 * t407 + t352 * t406 + t358 * t418 + t360 * t417) + t397 * (t349 * t407 + t351 * t406 + t357 * t418 + t359 * t417);
t489 = (t491 * t423 + t426 * t492) * t397 + (t423 * t492 - t491 * t426) * t396 + (Icges(5,5) * t418 + Icges(6,5) * t407 - Icges(5,6) * t417 - Icges(6,6) * t406) * qJD(1);
t422 = sin(qJ(3));
t482 = pkin(3) * t422;
t481 = pkin(4) * t418;
t478 = Icges(4,4) * t422;
t425 = cos(qJ(3));
t477 = Icges(4,4) * t425;
t472 = t407 * t423;
t471 = t407 * t426;
t421 = sin(qJ(6));
t470 = t421 * t423;
t469 = t421 * t426;
t424 = cos(qJ(6));
t468 = t423 * t424;
t467 = t424 * t426;
t393 = qJD(1) * (pkin(1) * t426 + qJ(2) * t423);
t466 = qJD(1) * t426 * pkin(7) + t393;
t416 = qJD(2) * t423;
t462 = pkin(3) * qJD(3) * t425;
t465 = t423 * t462 + t416;
t463 = qJD(6) * t407;
t461 = pkin(4) * t417;
t401 = pkin(1) * t423 - qJ(2) * t426;
t460 = -pkin(7) * t423 - t401;
t459 = qJD(5) * t426 + t396 * t481 + t465;
t381 = pkin(8) * t423 - t426 * t482;
t458 = -t381 + t460;
t457 = pkin(5) * t406 - pkin(9) * t407;
t382 = pkin(8) * t426 + t423 * t482;
t456 = t381 * t415 - t382 * t414;
t455 = rSges(4,1) * t422 + rSges(4,2) * t425;
t454 = rSges(5,1) * t417 + rSges(5,2) * t418;
t453 = rSges(6,1) * t406 + rSges(6,2) * t407;
t452 = Icges(4,1) * t422 + t477;
t449 = Icges(4,2) * t425 + t478;
t446 = Icges(4,5) * t422 + Icges(4,6) * t425;
t370 = Icges(4,6) * t426 + t423 * t449;
t372 = Icges(4,5) * t426 + t423 * t452;
t439 = -t370 * t425 - t372 * t422;
t371 = Icges(4,6) * t423 - t426 * t449;
t373 = Icges(4,5) * t423 - t426 * t452;
t438 = t371 * t425 + t373 * t422;
t399 = -Icges(4,2) * t422 + t477;
t400 = Icges(4,1) * t425 - t478;
t435 = t399 * t425 + t400 * t422;
t339 = qJ(5) * t423 - t426 * t461;
t434 = -t339 + t458;
t433 = t397 * t339 + t456;
t430 = qJD(1) * t382 + (-qJD(2) - t462) * t426 + t466;
t340 = qJ(5) * t426 + t423 * t461;
t429 = qJD(1) * t340 + qJD(5) * t423 + t430;
t404 = rSges(2,1) * t426 - rSges(2,2) * t423;
t403 = rSges(4,1) * t425 - rSges(4,2) * t422;
t402 = rSges(2,1) * t423 + rSges(2,2) * t426;
t398 = Icges(4,5) * t425 - Icges(4,6) * t422;
t395 = qJD(6) * t406 + qJD(1);
t392 = rSges(5,1) * t418 - rSges(5,2) * t417;
t388 = pkin(5) * t407 + pkin(9) * t406;
t386 = rSges(6,1) * t407 - rSges(6,2) * t406;
t380 = -t406 * t467 + t470;
t379 = t406 * t469 + t468;
t378 = t406 * t468 + t469;
t377 = -t406 * t470 + t467;
t376 = rSges(4,3) * t423 - t426 * t455;
t375 = rSges(4,3) * t426 + t423 * t455;
t369 = Icges(4,3) * t423 - t426 * t446;
t368 = Icges(4,3) * t426 + t423 * t446;
t367 = -t423 * t463 + t397;
t366 = t426 * t463 + t396;
t365 = t457 * t426;
t364 = t457 * t423;
t362 = rSges(5,3) * t423 - t426 * t454;
t361 = rSges(5,3) * t426 + t423 * t454;
t354 = rSges(6,3) * t423 - t426 * t453;
t353 = rSges(6,3) * t426 + t423 * t453;
t346 = t393 - qJD(2) * t426 + qJD(1) * (-rSges(3,2) * t426 + rSges(3,3) * t423);
t345 = t416 + (rSges(3,2) * t423 + rSges(3,3) * t426 - t401) * qJD(1);
t344 = rSges(7,3) * t406 + (rSges(7,1) * t424 - rSges(7,2) * t421) * t407;
t343 = Icges(7,5) * t406 + (Icges(7,1) * t424 - Icges(7,4) * t421) * t407;
t342 = Icges(7,6) * t406 + (Icges(7,4) * t424 - Icges(7,2) * t421) * t407;
t341 = Icges(7,3) * t406 + (Icges(7,5) * t424 - Icges(7,6) * t421) * t407;
t336 = (-t375 * t423 + t376 * t426) * qJD(3);
t335 = rSges(7,1) * t380 + rSges(7,2) * t379 + rSges(7,3) * t471;
t334 = rSges(7,1) * t378 + rSges(7,2) * t377 - rSges(7,3) * t472;
t333 = Icges(7,1) * t380 + Icges(7,4) * t379 + Icges(7,5) * t471;
t332 = Icges(7,1) * t378 + Icges(7,4) * t377 - Icges(7,5) * t472;
t331 = Icges(7,4) * t380 + Icges(7,2) * t379 + Icges(7,6) * t471;
t330 = Icges(7,4) * t378 + Icges(7,2) * t377 - Icges(7,6) * t472;
t329 = Icges(7,5) * t380 + Icges(7,6) * t379 + Icges(7,3) * t471;
t328 = Icges(7,5) * t378 + Icges(7,6) * t377 - Icges(7,3) * t472;
t327 = qJD(1) * t375 + (-qJD(3) * t403 - qJD(2)) * t426 + t466;
t326 = t403 * t414 + t416 + (-t376 + t460) * qJD(1);
t325 = qJD(1) * t361 - t392 * t397 + t430;
t324 = t392 * t396 + (-t362 + t458) * qJD(1) + t465;
t323 = -t361 * t396 + t362 * t397 + t456;
t322 = qJD(1) * t353 + (-t386 - t481) * t397 + t429;
t321 = t386 * t396 + (-t354 + t434) * qJD(1) + t459;
t320 = t354 * t397 + (-t340 - t353) * t396 + t433;
t319 = qJD(1) * t364 + t334 * t395 - t344 * t367 + (-t388 - t481) * t397 + t429;
t318 = -t335 * t395 + t344 * t366 + t388 * t396 + (t365 + t434) * qJD(1) + t459;
t317 = -t334 * t366 + t335 * t367 - t365 * t397 + (-t340 - t364) * t396 + t433;
t1 = ((t426 * t398 + t423 * t435) * qJD(1) + (t426 ^ 2 * t368 + (t438 * t423 + (t369 - t439) * t426) * t423) * qJD(3)) * t415 / 0.2e1 + t366 * ((t328 * t471 + t330 * t379 + t332 * t380) * t367 + (t329 * t471 + t379 * t331 + t380 * t333) * t366 + (t341 * t471 + t342 * t379 + t343 * t380) * t395) / 0.2e1 + t395 * ((t328 * t367 + t329 * t366 + t341 * t395) * t406 + ((-t330 * t421 + t332 * t424) * t367 + (-t331 * t421 + t333 * t424) * t366 + (-t342 * t421 + t343 * t424) * t395) * t407) / 0.2e1 + t367 * ((-t328 * t472 + t377 * t330 + t378 * t332) * t367 + (-t329 * t472 + t331 * t377 + t333 * t378) * t366 + (-t341 * t472 + t342 * t377 + t343 * t378) * t395) / 0.2e1 + ((t423 * t398 - t426 * t435) * qJD(1) + (t423 ^ 2 * t369 + (t439 * t426 + (t368 - t438) * t423) * t426) * qJD(3)) * t414 / 0.2e1 + m(6) * (t320 ^ 2 + t321 ^ 2 + t322 ^ 2) / 0.2e1 + m(7) * (t317 ^ 2 + t318 ^ 2 + t319 ^ 2) / 0.2e1 + m(5) * (t323 ^ 2 + t324 ^ 2 + t325 ^ 2) / 0.2e1 + m(4) * (t326 ^ 2 + t327 ^ 2 + t336 ^ 2) / 0.2e1 + m(3) * (t345 ^ 2 + t346 ^ 2) / 0.2e1 + (t489 * t423 - t490 * t426) * t396 / 0.2e1 + (t490 * t423 + t489 * t426) * t397 / 0.2e1 + (m(2) * (t402 ^ 2 + t404 ^ 2) + Icges(2,3) + Icges(3,1)) * qJD(1) ^ 2 / 0.2e1 + (((-t370 * t422 + t372 * t425) * t426 + (-t371 * t422 + t373 * t425) * t423) * qJD(3) + (-t349 * t406 + t351 * t407 - t357 * t417 + t359 * t418) * t397 + (-t350 * t406 + t352 * t407 - t358 * t417 + t360 * t418) * t396 + (-t406 * t384 + t407 * t385 - t417 * t390 + t418 * t391 - t422 * t399 + t425 * t400) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
