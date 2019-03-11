% Calculate kinetic energy for
% S6RPPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:20:17
% EndTime: 2019-03-09 02:20:18
% DurationCPUTime: 1.29s
% Computational Cost: add. (1726->238), mult. (1337->382), div. (0->0), fcn. (1298->12), ass. (0->126)
t414 = sin(qJ(1));
t456 = t414 * pkin(1);
t411 = cos(pkin(11));
t455 = pkin(3) * t411;
t415 = cos(qJ(5));
t454 = pkin(5) * t415;
t407 = pkin(11) + qJ(4);
t401 = sin(t407);
t452 = Icges(5,4) * t401;
t403 = cos(t407);
t451 = Icges(5,4) * t403;
t408 = qJ(1) + pkin(10);
t402 = sin(t408);
t450 = t401 * t402;
t404 = cos(t408);
t449 = t401 * t404;
t409 = qJ(5) + qJ(6);
t405 = sin(t409);
t448 = t402 * t405;
t406 = cos(t409);
t447 = t402 * t406;
t413 = sin(qJ(5));
t446 = t402 * t413;
t445 = t402 * t415;
t444 = t404 * t405;
t443 = t404 * t406;
t442 = t404 * t413;
t441 = t404 * t415;
t416 = cos(qJ(1));
t400 = qJD(1) * t416 * pkin(1);
t439 = qJD(1) * (pkin(2) * t404 + qJ(3) * t402) + t400;
t396 = qJD(4) * t402;
t437 = qJD(5) * t401;
t381 = t404 * t437 + t396;
t438 = qJD(4) * t404;
t436 = qJD(6) * t401;
t432 = pkin(4) * t403 + pkin(8) * t401;
t373 = t432 * t402;
t374 = t432 * t404;
t435 = t373 * t396 + t374 * t438 + qJD(2);
t434 = -pkin(2) * t402 + qJ(3) * t404 - t456;
t382 = t402 * t437 - t438;
t433 = pkin(7) * t404 - t455 * t402 + t434;
t410 = sin(pkin(11));
t431 = rSges(4,1) * t411 - rSges(4,2) * t410;
t430 = rSges(5,1) * t403 - rSges(5,2) * t401;
t429 = Icges(5,1) * t403 - t452;
t428 = -Icges(5,2) * t401 + t451;
t427 = Icges(5,5) * t403 - Icges(5,6) * t401;
t352 = -Icges(5,6) * t404 + t428 * t402;
t354 = -Icges(5,5) * t404 + t429 * t402;
t426 = t352 * t401 - t354 * t403;
t353 = Icges(5,6) * t402 + t428 * t404;
t355 = Icges(5,5) * t402 + t429 * t404;
t425 = -t353 * t401 + t355 * t403;
t386 = Icges(5,2) * t403 + t452;
t387 = Icges(5,1) * t401 + t451;
t424 = -t386 * t401 + t387 * t403;
t423 = -qJD(3) * t404 + qJD(1) * (pkin(7) * t402 + t455 * t404) + t439;
t422 = pkin(9) * t401 + t454 * t403;
t390 = pkin(4) * t401 - pkin(8) * t403;
t421 = qJD(1) * t374 - t390 * t396 + t423;
t397 = qJD(3) * t402;
t420 = t397 + (-t373 + t433) * qJD(1) - t390 * t438;
t418 = qJD(2) ^ 2;
t395 = -qJD(5) * t403 + qJD(1);
t394 = rSges(2,1) * t416 - rSges(2,2) * t414;
t393 = rSges(2,1) * t414 + rSges(2,2) * t416;
t388 = rSges(5,1) * t401 + rSges(5,2) * t403;
t385 = Icges(5,5) * t401 + Icges(5,6) * t403;
t383 = qJD(1) + (-qJD(5) - qJD(6)) * t403;
t380 = t403 * t441 + t446;
t379 = -t403 * t442 + t445;
t378 = t403 * t445 - t442;
t377 = -t403 * t446 - t441;
t376 = t400 + qJD(1) * (rSges(3,1) * t404 - rSges(3,2) * t402);
t375 = (-rSges(3,1) * t402 - rSges(3,2) * t404 - t456) * qJD(1);
t371 = t403 * t443 + t448;
t370 = -t403 * t444 + t447;
t369 = t403 * t447 - t444;
t368 = -t403 * t448 - t443;
t367 = -rSges(6,3) * t403 + (rSges(6,1) * t415 - rSges(6,2) * t413) * t401;
t366 = -Icges(6,5) * t403 + (Icges(6,1) * t415 - Icges(6,4) * t413) * t401;
t365 = -Icges(6,6) * t403 + (Icges(6,4) * t415 - Icges(6,2) * t413) * t401;
t364 = -Icges(6,3) * t403 + (Icges(6,5) * t415 - Icges(6,6) * t413) * t401;
t363 = -rSges(7,3) * t403 + (rSges(7,1) * t406 - rSges(7,2) * t405) * t401;
t362 = -Icges(7,5) * t403 + (Icges(7,1) * t406 - Icges(7,4) * t405) * t401;
t361 = -Icges(7,6) * t403 + (Icges(7,4) * t406 - Icges(7,2) * t405) * t401;
t360 = -Icges(7,3) * t403 + (Icges(7,5) * t406 - Icges(7,6) * t405) * t401;
t357 = rSges(5,3) * t402 + t430 * t404;
t356 = -rSges(5,3) * t404 + t430 * t402;
t351 = Icges(5,3) * t402 + t427 * t404;
t350 = -Icges(5,3) * t404 + t427 * t402;
t349 = t402 * t436 + t382;
t348 = t404 * t436 + t381;
t347 = -pkin(9) * t403 + t454 * t401;
t344 = qJD(1) * t402 * rSges(4,3) + (qJD(1) * t431 - qJD(3)) * t404 + t439;
t343 = t397 + (t404 * rSges(4,3) - t431 * t402 + t434) * qJD(1);
t342 = rSges(6,1) * t380 + rSges(6,2) * t379 + rSges(6,3) * t449;
t341 = rSges(6,1) * t378 + rSges(6,2) * t377 + rSges(6,3) * t450;
t340 = Icges(6,1) * t380 + Icges(6,4) * t379 + Icges(6,5) * t449;
t339 = Icges(6,1) * t378 + Icges(6,4) * t377 + Icges(6,5) * t450;
t338 = Icges(6,4) * t380 + Icges(6,2) * t379 + Icges(6,6) * t449;
t337 = Icges(6,4) * t378 + Icges(6,2) * t377 + Icges(6,6) * t450;
t336 = Icges(6,5) * t380 + Icges(6,6) * t379 + Icges(6,3) * t449;
t335 = Icges(6,5) * t378 + Icges(6,6) * t377 + Icges(6,3) * t450;
t334 = pkin(5) * t446 + t422 * t404;
t333 = -pkin(5) * t442 + t422 * t402;
t332 = rSges(7,1) * t371 + rSges(7,2) * t370 + rSges(7,3) * t449;
t331 = rSges(7,1) * t369 + rSges(7,2) * t368 + rSges(7,3) * t450;
t330 = Icges(7,1) * t371 + Icges(7,4) * t370 + Icges(7,5) * t449;
t329 = Icges(7,1) * t369 + Icges(7,4) * t368 + Icges(7,5) * t450;
t328 = Icges(7,4) * t371 + Icges(7,2) * t370 + Icges(7,6) * t449;
t327 = Icges(7,4) * t369 + Icges(7,2) * t368 + Icges(7,6) * t450;
t326 = Icges(7,5) * t371 + Icges(7,6) * t370 + Icges(7,3) * t449;
t325 = Icges(7,5) * t369 + Icges(7,6) * t368 + Icges(7,3) * t450;
t324 = qJD(2) + (t356 * t402 + t357 * t404) * qJD(4);
t323 = qJD(1) * t357 - t388 * t396 + t423;
t322 = -t388 * t438 + t397 + (-t356 + t433) * qJD(1);
t321 = t341 * t381 - t342 * t382 + t435;
t320 = t342 * t395 - t367 * t381 + t421;
t319 = -t341 * t395 + t367 * t382 + t420;
t318 = t332 * t383 + t334 * t395 - t347 * t381 - t348 * t363 + t421;
t317 = -t331 * t383 - t333 * t395 + t347 * t382 + t349 * t363 + t420;
t316 = t331 * t348 - t332 * t349 + t333 * t381 - t334 * t382 + t435;
t1 = ((t402 * t385 + t424 * t404) * qJD(1) + (t402 ^ 2 * t351 + (t426 * t404 + (-t350 + t425) * t402) * t404) * qJD(4)) * t396 / 0.2e1 - ((-t404 * t385 + t424 * t402) * qJD(1) + (t404 ^ 2 * t350 + (t425 * t402 + (-t351 + t426) * t404) * t402) * qJD(4)) * t438 / 0.2e1 + qJD(1) * ((t403 * t386 + t401 * t387) * qJD(1) + ((t353 * t403 + t355 * t401) * t402 - (t352 * t403 + t354 * t401) * t404) * qJD(4)) / 0.2e1 + m(6) * (t319 ^ 2 + t320 ^ 2 + t321 ^ 2) / 0.2e1 + m(7) * (t316 ^ 2 + t317 ^ 2 + t318 ^ 2) / 0.2e1 + m(5) * (t322 ^ 2 + t323 ^ 2 + t324 ^ 2) / 0.2e1 + m(3) * (t375 ^ 2 + t376 ^ 2 + t418) / 0.2e1 + m(4) * (t343 ^ 2 + t344 ^ 2 + t418) / 0.2e1 + t348 * ((t326 * t449 + t370 * t328 + t371 * t330) * t348 + (t325 * t449 + t327 * t370 + t329 * t371) * t349 + (t360 * t449 + t361 * t370 + t362 * t371) * t383) / 0.2e1 + t349 * ((t326 * t450 + t328 * t368 + t330 * t369) * t348 + (t325 * t450 + t368 * t327 + t369 * t329) * t349 + (t360 * t450 + t361 * t368 + t362 * t369) * t383) / 0.2e1 + t383 * ((-t325 * t349 - t326 * t348 - t360 * t383) * t403 + ((-t328 * t405 + t330 * t406) * t348 + (-t327 * t405 + t329 * t406) * t349 + (-t361 * t405 + t362 * t406) * t383) * t401) / 0.2e1 + t381 * ((t336 * t449 + t379 * t338 + t380 * t340) * t381 + (t335 * t449 + t337 * t379 + t339 * t380) * t382 + (t364 * t449 + t365 * t379 + t366 * t380) * t395) / 0.2e1 + t382 * ((t336 * t450 + t338 * t377 + t340 * t378) * t381 + (t335 * t450 + t377 * t337 + t378 * t339) * t382 + (t364 * t450 + t365 * t377 + t366 * t378) * t395) / 0.2e1 + t395 * ((-t335 * t382 - t336 * t381 - t364 * t395) * t403 + ((-t338 * t413 + t340 * t415) * t381 + (-t337 * t413 + t339 * t415) * t382 + (-t365 * t413 + t366 * t415) * t395) * t401) / 0.2e1 + (Icges(4,2) * t411 ^ 2 + (Icges(4,1) * t410 + 0.2e1 * Icges(4,4) * t411) * t410 + Icges(2,3) + Icges(3,3) + m(2) * (t393 ^ 2 + t394 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
