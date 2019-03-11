% Calculate kinetic energy for
% S6RPPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRR8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR8_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR8_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR8_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR8_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:35:23
% EndTime: 2019-03-09 02:35:24
% DurationCPUTime: 1.31s
% Computational Cost: add. (1181->239), mult. (1341->383), div. (0->0), fcn. (1304->10), ass. (0->122)
t406 = sin(qJ(1));
t408 = cos(qJ(1));
t446 = qJ(3) * qJD(1) * t408 + qJD(3) * t406;
t400 = pkin(10) + qJ(4);
t391 = sin(t400);
t392 = cos(t400);
t407 = cos(qJ(5));
t442 = pkin(5) * t407;
t445 = -pkin(9) * t392 + t391 * t442;
t402 = sin(pkin(10));
t443 = pkin(3) * t402;
t440 = Icges(5,4) * t391;
t439 = Icges(5,4) * t392;
t438 = t392 * t406;
t437 = t392 * t408;
t401 = qJ(5) + qJ(6);
t398 = sin(t401);
t436 = t398 * t406;
t435 = t398 * t408;
t399 = cos(t401);
t434 = t399 * t406;
t433 = t399 * t408;
t405 = sin(qJ(5));
t432 = t405 * t406;
t431 = t405 * t408;
t430 = t406 * t407;
t429 = t407 * t408;
t393 = qJD(4) * t406;
t426 = qJD(5) * t392;
t375 = t408 * t426 + t393;
t397 = qJD(2) * t406;
t427 = qJD(3) * t408 + t397;
t394 = qJD(4) * t408;
t387 = qJD(5) * t391 + qJD(1);
t382 = qJD(1) * (pkin(1) * t408 + qJ(2) * t406);
t425 = -qJD(2) * t408 + t382;
t424 = qJD(1) * (pkin(7) * t408 + t406 * t443) + t382 + t446;
t423 = pkin(4) * t391 - pkin(8) * t392;
t384 = pkin(1) * t406 - qJ(2) * t408;
t422 = t408 * t443 - t384 + (-pkin(7) - qJ(3)) * t406;
t366 = t423 * t406;
t367 = t423 * t408;
t421 = -t366 * t393 - t367 * t394;
t403 = cos(pkin(10));
t420 = rSges(4,1) * t402 + rSges(4,2) * t403;
t419 = rSges(5,1) * t391 + rSges(5,2) * t392;
t418 = Icges(5,1) * t391 + t439;
t417 = Icges(5,2) * t392 + t440;
t416 = Icges(5,5) * t391 + Icges(5,6) * t392;
t353 = Icges(5,6) * t408 + t406 * t417;
t355 = Icges(5,5) * t408 + t406 * t418;
t415 = -t353 * t392 - t355 * t391;
t354 = Icges(5,6) * t406 - t408 * t417;
t356 = Icges(5,5) * t406 - t408 * t418;
t414 = t354 * t392 + t356 * t391;
t378 = -Icges(5,2) * t391 + t439;
t379 = Icges(5,1) * t392 - t440;
t413 = t378 * t392 + t379 * t391;
t381 = pkin(4) * t392 + pkin(8) * t391;
t412 = t381 * t393 + (t367 + t422) * qJD(1) + t427;
t411 = qJD(1) * t366 + (-qJD(4) * t381 - qJD(2)) * t408 + t424;
t386 = rSges(2,1) * t408 - rSges(2,2) * t406;
t385 = rSges(2,1) * t406 + rSges(2,2) * t408;
t380 = rSges(5,1) * t392 - rSges(5,2) * t391;
t377 = Icges(5,5) * t392 - Icges(5,6) * t391;
t376 = -t406 * t426 + t394;
t374 = qJD(6) * t391 + t387;
t373 = -t391 * t429 + t432;
t372 = t391 * t431 + t430;
t371 = t391 * t430 + t431;
t370 = -t391 * t432 + t429;
t363 = -t391 * t433 + t436;
t362 = t391 * t435 + t434;
t361 = t391 * t434 + t435;
t360 = -t391 * t436 + t433;
t358 = rSges(5,3) * t406 - t408 * t419;
t357 = rSges(5,3) * t408 + t406 * t419;
t352 = Icges(5,3) * t406 - t408 * t416;
t351 = Icges(5,3) * t408 + t406 * t416;
t350 = t394 + (-qJD(5) - qJD(6)) * t438;
t349 = qJD(6) * t437 + t375;
t348 = rSges(6,3) * t391 + (rSges(6,1) * t407 - rSges(6,2) * t405) * t392;
t347 = Icges(6,5) * t391 + (Icges(6,1) * t407 - Icges(6,4) * t405) * t392;
t346 = Icges(6,6) * t391 + (Icges(6,4) * t407 - Icges(6,2) * t405) * t392;
t345 = Icges(6,3) * t391 + (Icges(6,5) * t407 - Icges(6,6) * t405) * t392;
t344 = qJD(1) * (-rSges(3,2) * t408 + rSges(3,3) * t406) + t425;
t343 = t397 + (rSges(3,2) * t406 + rSges(3,3) * t408 - t384) * qJD(1);
t342 = rSges(7,3) * t391 + (rSges(7,1) * t399 - rSges(7,2) * t398) * t392;
t341 = Icges(7,5) * t391 + (Icges(7,1) * t399 - Icges(7,4) * t398) * t392;
t340 = Icges(7,6) * t391 + (Icges(7,4) * t399 - Icges(7,2) * t398) * t392;
t339 = Icges(7,3) * t391 + (Icges(7,5) * t399 - Icges(7,6) * t398) * t392;
t338 = pkin(9) * t391 + t392 * t442;
t337 = rSges(6,1) * t373 + rSges(6,2) * t372 + rSges(6,3) * t437;
t336 = rSges(6,1) * t371 + rSges(6,2) * t370 - rSges(6,3) * t438;
t335 = Icges(6,1) * t373 + Icges(6,4) * t372 + Icges(6,5) * t437;
t334 = Icges(6,1) * t371 + Icges(6,4) * t370 - Icges(6,5) * t438;
t333 = Icges(6,4) * t373 + Icges(6,2) * t372 + Icges(6,6) * t437;
t332 = Icges(6,4) * t371 + Icges(6,2) * t370 - Icges(6,6) * t438;
t331 = Icges(6,5) * t373 + Icges(6,6) * t372 + Icges(6,3) * t437;
t330 = Icges(6,5) * t371 + Icges(6,6) * t370 - Icges(6,3) * t438;
t329 = qJD(1) * (rSges(4,3) * t408 + t406 * t420) + t425 + t446;
t328 = (-t384 + t420 * t408 + (-rSges(4,3) - qJ(3)) * t406) * qJD(1) + t427;
t327 = pkin(5) * t432 - t408 * t445;
t326 = pkin(5) * t431 + t406 * t445;
t325 = rSges(7,1) * t363 + rSges(7,2) * t362 + rSges(7,3) * t437;
t324 = rSges(7,1) * t361 + rSges(7,2) * t360 - rSges(7,3) * t438;
t323 = Icges(7,1) * t363 + Icges(7,4) * t362 + Icges(7,5) * t437;
t322 = Icges(7,1) * t361 + Icges(7,4) * t360 - Icges(7,5) * t438;
t321 = Icges(7,4) * t363 + Icges(7,2) * t362 + Icges(7,6) * t437;
t320 = Icges(7,4) * t361 + Icges(7,2) * t360 - Icges(7,6) * t438;
t319 = Icges(7,5) * t363 + Icges(7,6) * t362 + Icges(7,3) * t437;
t318 = Icges(7,5) * t361 + Icges(7,6) * t360 - Icges(7,3) * t438;
t317 = (-t357 * t406 + t358 * t408) * qJD(4);
t316 = qJD(1) * t357 + (-qJD(4) * t380 - qJD(2)) * t408 + t424;
t315 = t380 * t393 + (-t358 + t422) * qJD(1) + t427;
t314 = -t336 * t375 + t337 * t376 + t421;
t313 = t336 * t387 - t348 * t376 + t411;
t312 = -t337 * t387 + t348 * t375 + t412;
t311 = t324 * t374 + t326 * t387 - t338 * t376 - t342 * t350 + t411;
t310 = -t325 * t374 - t327 * t387 + t338 * t375 + t342 * t349 + t412;
t309 = -t324 * t349 + t325 * t350 - t326 * t375 + t327 * t376 + t421;
t1 = ((t408 * t377 + t406 * t413) * qJD(1) + (t408 ^ 2 * t351 + (t414 * t406 + (t352 - t415) * t408) * t406) * qJD(4)) * t394 / 0.2e1 + m(7) * (t309 ^ 2 + t310 ^ 2 + t311 ^ 2) / 0.2e1 + m(6) * (t312 ^ 2 + t313 ^ 2 + t314 ^ 2) / 0.2e1 + m(5) * (t315 ^ 2 + t316 ^ 2 + t317 ^ 2) / 0.2e1 + m(3) * (t343 ^ 2 + t344 ^ 2) / 0.2e1 + m(4) * (t328 ^ 2 + t329 ^ 2) / 0.2e1 + t376 * ((-t330 * t438 + t370 * t332 + t371 * t334) * t376 + (-t331 * t438 + t333 * t370 + t335 * t371) * t375 + (-t345 * t438 + t346 * t370 + t347 * t371) * t387) / 0.2e1 + t374 * ((t318 * t350 + t319 * t349 + t339 * t374) * t391 + ((-t320 * t398 + t322 * t399) * t350 + (-t321 * t398 + t323 * t399) * t349 + (-t340 * t398 + t341 * t399) * t374) * t392) / 0.2e1 + t350 * ((-t318 * t438 + t360 * t320 + t361 * t322) * t350 + (-t319 * t438 + t321 * t360 + t323 * t361) * t349 + (-t339 * t438 + t340 * t360 + t341 * t361) * t374) / 0.2e1 + t349 * ((t318 * t437 + t320 * t362 + t322 * t363) * t350 + (t319 * t437 + t362 * t321 + t363 * t323) * t349 + (t339 * t437 + t340 * t362 + t341 * t363) * t374) / 0.2e1 + t375 * ((t330 * t437 + t332 * t372 + t334 * t373) * t376 + (t331 * t437 + t372 * t333 + t373 * t335) * t375 + (t345 * t437 + t346 * t372 + t347 * t373) * t387) / 0.2e1 + t387 * ((t330 * t376 + t331 * t375 + t345 * t387) * t391 + ((-t332 * t405 + t334 * t407) * t376 + (-t333 * t405 + t335 * t407) * t375 + (-t346 * t405 + t347 * t407) * t387) * t392) / 0.2e1 + ((t406 * t377 - t408 * t413) * qJD(1) + (t406 ^ 2 * t352 + (t415 * t408 + (t351 - t414) * t406) * t408) * qJD(4)) * t393 / 0.2e1 + qJD(1) * ((-t391 * t378 + t392 * t379) * qJD(1) + ((-t353 * t391 + t355 * t392) * t408 + (-t354 * t391 + t356 * t392) * t406) * qJD(4)) / 0.2e1 + (Icges(4,1) * t403 ^ 2 + (-0.2e1 * Icges(4,4) * t403 + Icges(4,2) * t402) * t402 + m(2) * (t385 ^ 2 + t386 ^ 2) + Icges(2,3) + Icges(3,1)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
