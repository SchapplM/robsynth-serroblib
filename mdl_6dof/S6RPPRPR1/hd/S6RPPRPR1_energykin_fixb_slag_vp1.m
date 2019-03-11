% Calculate kinetic energy for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
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
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRPR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:39:00
% EndTime: 2019-03-09 01:39:01
% DurationCPUTime: 1.39s
% Computational Cost: add. (1630->238), mult. (1277->371), div. (0->0), fcn. (1238->12), ass. (0->121)
t402 = sin(qJ(1));
t446 = t402 * pkin(1);
t399 = cos(pkin(10));
t444 = pkin(3) * t399;
t398 = cos(pkin(11));
t443 = pkin(5) * t398;
t394 = pkin(10) + qJ(4);
t388 = sin(t394);
t442 = Icges(5,4) * t388;
t391 = cos(t394);
t441 = Icges(5,4) * t391;
t395 = qJ(1) + pkin(9);
t389 = sin(t395);
t440 = t388 * t389;
t392 = cos(t395);
t439 = t388 * t392;
t396 = sin(pkin(11));
t438 = t389 * t396;
t437 = t391 * t389;
t436 = t391 * t392;
t435 = t392 * t396;
t434 = t392 * t398;
t403 = cos(qJ(1));
t386 = qJD(1) * t403 * pkin(1);
t431 = qJD(1) * (pkin(2) * t392 + qJ(3) * t389) + t386;
t383 = qJD(3) * t389;
t427 = qJD(5) * t388;
t430 = t392 * t427 + t383;
t429 = qJD(4) * t389;
t428 = qJD(4) * t392;
t426 = qJD(6) * t388;
t423 = -pkin(2) * t389 + qJ(3) * t392 - t446;
t375 = pkin(4) * t388 - qJ(5) * t391;
t422 = qJD(4) * (pkin(8) * t391 - t388 * t443 - t375);
t421 = qJD(4) * (rSges(6,3) * t391 - (rSges(6,1) * t398 - rSges(6,2) * t396) * t388 - t375);
t420 = pkin(7) * t392 - t389 * t444 + t423;
t397 = sin(pkin(10));
t419 = rSges(4,1) * t399 - rSges(4,2) * t397;
t418 = rSges(5,1) * t391 - rSges(5,2) * t388;
t417 = pkin(4) * t391 + qJ(5) * t388;
t416 = Icges(5,1) * t391 - t442;
t415 = -Icges(5,2) * t388 + t441;
t414 = Icges(5,5) * t391 - Icges(5,6) * t388;
t342 = -Icges(5,6) * t392 + t389 * t415;
t345 = -Icges(5,5) * t392 + t389 * t416;
t413 = t342 * t388 - t345 * t391;
t343 = Icges(5,6) * t389 + t392 * t415;
t346 = Icges(5,5) * t389 + t392 * t416;
t412 = -t343 * t388 + t346 * t391;
t373 = Icges(5,2) * t391 + t442;
t374 = Icges(5,1) * t388 + t441;
t411 = -t373 * t388 + t374 * t391;
t361 = t417 * t389;
t410 = -t361 + t420;
t409 = -qJD(3) * t392 + qJD(1) * (pkin(7) * t389 + t392 * t444) + t431;
t362 = t417 * t392;
t408 = -qJD(5) * t391 + t361 * t429 + t362 * t428 + qJD(2);
t407 = qJD(1) * t362 + t389 * t427 + t409;
t406 = pkin(8) * t388 + t391 * t443;
t404 = qJD(2) ^ 2;
t393 = pkin(11) + qJ(6);
t390 = cos(t393);
t387 = sin(t393);
t382 = -qJD(6) * t391 + qJD(1);
t381 = rSges(2,1) * t403 - rSges(2,2) * t402;
t380 = rSges(2,1) * t402 + rSges(2,2) * t403;
t376 = rSges(5,1) * t388 + rSges(5,2) * t391;
t372 = Icges(5,5) * t388 + Icges(5,6) * t391;
t370 = t389 * t426 - t428;
t369 = t392 * t426 + t429;
t368 = t386 + qJD(1) * (rSges(3,1) * t392 - rSges(3,2) * t389);
t367 = (-rSges(3,1) * t389 - rSges(3,2) * t392 - t446) * qJD(1);
t366 = t391 * t434 + t438;
t365 = t389 * t398 - t391 * t435;
t364 = t398 * t437 - t435;
t363 = -t396 * t437 - t434;
t359 = t387 * t389 + t390 * t436;
t358 = -t387 * t436 + t389 * t390;
t357 = -t387 * t392 + t390 * t437;
t356 = -t387 * t437 - t390 * t392;
t354 = -Icges(6,5) * t391 + (Icges(6,1) * t398 - Icges(6,4) * t396) * t388;
t353 = -Icges(6,6) * t391 + (Icges(6,4) * t398 - Icges(6,2) * t396) * t388;
t352 = -Icges(6,3) * t391 + (Icges(6,5) * t398 - Icges(6,6) * t396) * t388;
t351 = rSges(5,3) * t389 + t392 * t418;
t350 = -rSges(5,3) * t392 + t389 * t418;
t349 = -rSges(7,3) * t391 + (rSges(7,1) * t390 - rSges(7,2) * t387) * t388;
t344 = -Icges(7,5) * t391 + (Icges(7,1) * t390 - Icges(7,4) * t387) * t388;
t341 = -Icges(7,6) * t391 + (Icges(7,4) * t390 - Icges(7,2) * t387) * t388;
t340 = Icges(5,3) * t389 + t392 * t414;
t339 = -Icges(5,3) * t392 + t389 * t414;
t338 = -Icges(7,3) * t391 + (Icges(7,5) * t390 - Icges(7,6) * t387) * t388;
t334 = qJD(1) * t389 * rSges(4,3) + (qJD(1) * t419 - qJD(3)) * t392 + t431;
t333 = t383 + (t392 * rSges(4,3) - t389 * t419 + t423) * qJD(1);
t332 = rSges(6,1) * t366 + rSges(6,2) * t365 + rSges(6,3) * t439;
t331 = rSges(6,1) * t364 + rSges(6,2) * t363 + rSges(6,3) * t440;
t330 = Icges(6,1) * t366 + Icges(6,4) * t365 + Icges(6,5) * t439;
t329 = Icges(6,1) * t364 + Icges(6,4) * t363 + Icges(6,5) * t440;
t328 = Icges(6,4) * t366 + Icges(6,2) * t365 + Icges(6,6) * t439;
t327 = Icges(6,4) * t364 + Icges(6,2) * t363 + Icges(6,6) * t440;
t326 = Icges(6,5) * t366 + Icges(6,6) * t365 + Icges(6,3) * t439;
t325 = Icges(6,5) * t364 + Icges(6,6) * t363 + Icges(6,3) * t440;
t324 = pkin(5) * t438 + t392 * t406;
t323 = -pkin(5) * t435 + t389 * t406;
t322 = rSges(7,1) * t359 + rSges(7,2) * t358 + rSges(7,3) * t439;
t321 = rSges(7,1) * t357 + rSges(7,2) * t356 + rSges(7,3) * t440;
t320 = Icges(7,1) * t359 + Icges(7,4) * t358 + Icges(7,5) * t439;
t319 = Icges(7,1) * t357 + Icges(7,4) * t356 + Icges(7,5) * t440;
t318 = Icges(7,4) * t359 + Icges(7,2) * t358 + Icges(7,6) * t439;
t317 = Icges(7,4) * t357 + Icges(7,2) * t356 + Icges(7,6) * t440;
t316 = Icges(7,5) * t359 + Icges(7,6) * t358 + Icges(7,3) * t439;
t315 = Icges(7,5) * t357 + Icges(7,6) * t356 + Icges(7,3) * t440;
t314 = qJD(2) + (t350 * t389 + t351 * t392) * qJD(4);
t313 = qJD(1) * t351 - t376 * t429 + t409;
t312 = -t376 * t428 + t383 + (-t350 + t420) * qJD(1);
t311 = (t331 * t389 + t332 * t392) * qJD(4) + t408;
t310 = qJD(1) * t332 + t389 * t421 + t407;
t309 = t392 * t421 + (-t331 + t410) * qJD(1) + t430;
t308 = qJD(1) * t324 + t322 * t382 - t349 * t369 + t389 * t422 + t407;
t307 = -t321 * t382 + t349 * t370 + t392 * t422 + (-t323 + t410) * qJD(1) + t430;
t306 = t321 * t369 - t322 * t370 + (t323 * t389 + t324 * t392) * qJD(4) + t408;
t1 = t369 * ((t316 * t439 + t358 * t318 + t359 * t320) * t369 + (t315 * t439 + t317 * t358 + t319 * t359) * t370 + (t338 * t439 + t341 * t358 + t344 * t359) * t382) / 0.2e1 + t370 * ((t316 * t440 + t318 * t356 + t320 * t357) * t369 + (t315 * t440 + t356 * t317 + t357 * t319) * t370 + (t338 * t440 + t341 * t356 + t344 * t357) * t382) / 0.2e1 + t382 * ((-t315 * t370 - t316 * t369 - t338 * t382) * t391 + ((-t318 * t387 + t320 * t390) * t369 + (-t317 * t387 + t319 * t390) * t370 + (-t341 * t387 + t344 * t390) * t382) * t388) / 0.2e1 + m(3) * (t367 ^ 2 + t368 ^ 2 + t404) / 0.2e1 + m(4) * (t333 ^ 2 + t334 ^ 2 + t404) / 0.2e1 + m(5) * (t312 ^ 2 + t313 ^ 2 + t314 ^ 2) / 0.2e1 + m(6) * (t309 ^ 2 + t310 ^ 2 + t311 ^ 2) / 0.2e1 + m(7) * (t306 ^ 2 + t307 ^ 2 + t308 ^ 2) / 0.2e1 + (((t343 * t391 + t346 * t388) * t389 - (t342 * t391 + t345 * t388) * t392 + (t325 * t392 - t326 * t389) * t391 + ((-t328 * t396 + t330 * t398) * t389 - (-t327 * t396 + t329 * t398) * t392) * t388) * qJD(4) + ((t373 - t352) * t391 + (-t353 * t396 + t354 * t398 + t374) * t388) * qJD(1)) * qJD(1) / 0.2e1 + (((-t325 * t439 - t327 * t365 - t329 * t366 + t413 * t392) * t392 + ((-t339 + t412) * t392 + t326 * t439 + t365 * t328 + t366 * t330 + t389 * t340) * t389) * qJD(4) + (t352 * t439 + t353 * t365 + t354 * t366 + t389 * t372 + t392 * t411) * qJD(1)) * t429 / 0.2e1 - (((-t325 * t440 - t363 * t327 - t364 * t329 + t392 * t339) * t392 + (t326 * t440 + t328 * t363 + t330 * t364 + (-t340 + t413) * t392 + t412 * t389) * t389) * qJD(4) + (t352 * t440 + t353 * t363 + t354 * t364 - t392 * t372 + t411 * t389) * qJD(1)) * t428 / 0.2e1 + (m(2) * (t380 ^ 2 + t381 ^ 2) + Icges(2,3) + Icges(3,3) + Icges(4,2) * t399 ^ 2 + (Icges(4,1) * t397 + 0.2e1 * Icges(4,4) * t399) * t397) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
