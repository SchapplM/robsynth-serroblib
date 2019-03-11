% Calculate kinetic energy for
% S6RPPRRR4
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
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:25:22
% EndTime: 2019-03-09 02:25:23
% DurationCPUTime: 1.37s
% Computational Cost: add. (1220->233), mult. (2306->381), div. (0->0), fcn. (2766->10), ass. (0->117)
t441 = cos(qJ(1));
t440 = sin(qJ(1));
t403 = cos(qJ(5));
t439 = pkin(5) * t403;
t437 = cos(pkin(10));
t436 = sin(pkin(10));
t402 = sin(qJ(4));
t435 = Icges(5,4) * t402;
t404 = cos(qJ(4));
t434 = Icges(5,4) * t404;
t378 = -t440 * t436 - t441 * t437;
t401 = sin(qJ(5));
t433 = t378 * t401;
t432 = t378 * t402;
t379 = t441 * t436 - t440 * t437;
t431 = t379 * t401;
t430 = t379 * t402;
t429 = t404 * t378;
t428 = t404 * t379;
t376 = qJD(4) * t379;
t425 = qJD(5) * t402;
t358 = -t378 * t425 + t376;
t377 = qJD(4) * t378;
t359 = -t379 * t425 - t377;
t427 = qJD(4) * (-rSges(5,1) * t402 - rSges(5,2) * t404);
t426 = qJD(4) * (-pkin(4) * t402 + pkin(8) * t404);
t424 = qJD(6) * t402;
t390 = qJD(5) * t404 + qJD(1);
t386 = t440 * pkin(1) - t441 * qJ(2);
t423 = -t440 * pkin(2) - t386;
t422 = -pkin(4) * t404 - pkin(8) * t402;
t421 = -qJD(2) * t441 + qJD(1) * (t441 * pkin(1) + t440 * qJ(2));
t420 = -rSges(5,1) * t404 + rSges(5,2) * t402;
t419 = -Icges(5,1) * t404 + t435;
t418 = Icges(5,2) * t402 - t434;
t417 = -Icges(5,5) * t404 + Icges(5,6) * t402;
t340 = -Icges(5,6) * t378 + t418 * t379;
t342 = -Icges(5,5) * t378 + t419 * t379;
t416 = -t340 * t402 + t342 * t404;
t341 = Icges(5,6) * t379 + t418 * t378;
t343 = Icges(5,5) * t379 + t419 * t378;
t415 = t341 * t402 - t343 * t404;
t383 = -Icges(5,2) * t404 - t435;
t384 = -Icges(5,1) * t402 - t434;
t414 = t383 * t402 - t384 * t404;
t413 = pkin(3) * t379 + pkin(7) * t378 + t423;
t412 = qJD(1) * t441 * pkin(2) + t421;
t356 = t422 * t379;
t357 = t422 * t378;
t411 = t356 * t376 + t357 * t377 - qJD(3);
t410 = qJD(1) * (-pkin(3) * t378 + pkin(7) * t379) + t412;
t409 = -pkin(9) * t402 - t439 * t404;
t397 = qJD(2) * t440;
t408 = t397 + (-t356 + t413) * qJD(1) - t378 * t426;
t407 = qJD(1) * t357 - t379 * t426 + t410;
t400 = qJ(5) + qJ(6);
t399 = cos(t400);
t398 = sin(t400);
t388 = t441 * rSges(2,1) - t440 * rSges(2,2);
t387 = t440 * rSges(2,1) + t441 * rSges(2,2);
t382 = -Icges(5,5) * t402 - Icges(5,6) * t404;
t380 = qJD(6) * t404 + t390;
t373 = rSges(6,3) * t404 + (-rSges(6,1) * t403 + rSges(6,2) * t401) * t402;
t372 = Icges(6,5) * t404 + (-Icges(6,1) * t403 + Icges(6,4) * t401) * t402;
t371 = Icges(6,6) * t404 + (-Icges(6,4) * t403 + Icges(6,2) * t401) * t402;
t370 = Icges(6,3) * t404 + (-Icges(6,5) * t403 + Icges(6,6) * t401) * t402;
t369 = rSges(7,3) * t404 + (-rSges(7,1) * t399 + rSges(7,2) * t398) * t402;
t368 = Icges(7,5) * t404 + (-Icges(7,1) * t399 + Icges(7,4) * t398) * t402;
t367 = Icges(7,6) * t404 + (-Icges(7,4) * t399 + Icges(7,2) * t398) * t402;
t366 = Icges(7,3) * t404 + (-Icges(7,5) * t399 + Icges(7,6) * t398) * t402;
t365 = pkin(9) * t404 - t439 * t402;
t363 = qJD(1) * (t441 * rSges(3,1) + t440 * rSges(3,3)) + t421;
t362 = t397 + (-t440 * rSges(3,1) + t441 * rSges(3,3) - t386) * qJD(1);
t355 = -t403 * t429 + t431;
t354 = t379 * t403 + t401 * t429;
t353 = -t403 * t428 - t433;
t352 = -t378 * t403 + t401 * t428;
t350 = t379 * t398 - t399 * t429;
t349 = t379 * t399 + t398 * t429;
t348 = -t378 * t398 - t399 * t428;
t347 = -t378 * t399 + t398 * t428;
t345 = rSges(5,3) * t379 + t420 * t378;
t344 = -rSges(5,3) * t378 + t420 * t379;
t339 = Icges(5,3) * t379 + t417 * t378;
t338 = -Icges(5,3) * t378 + t417 * t379;
t337 = -t379 * t424 + t359;
t336 = -t378 * t424 + t358;
t335 = qJD(1) * (-rSges(4,1) * t378 - rSges(4,2) * t379) + t412;
t334 = t397 + (t379 * rSges(4,1) - t378 * rSges(4,2) + t423) * qJD(1);
t333 = pkin(5) * t431 + t409 * t378;
t332 = -pkin(5) * t433 + t409 * t379;
t331 = rSges(6,1) * t355 + rSges(6,2) * t354 - rSges(6,3) * t432;
t330 = rSges(6,1) * t353 + rSges(6,2) * t352 - rSges(6,3) * t430;
t329 = Icges(6,1) * t355 + Icges(6,4) * t354 - Icges(6,5) * t432;
t328 = Icges(6,1) * t353 + Icges(6,4) * t352 - Icges(6,5) * t430;
t327 = Icges(6,4) * t355 + Icges(6,2) * t354 - Icges(6,6) * t432;
t326 = Icges(6,4) * t353 + Icges(6,2) * t352 - Icges(6,6) * t430;
t325 = Icges(6,5) * t355 + Icges(6,6) * t354 - Icges(6,3) * t432;
t324 = Icges(6,5) * t353 + Icges(6,6) * t352 - Icges(6,3) * t430;
t323 = rSges(7,1) * t350 + rSges(7,2) * t349 - rSges(7,3) * t432;
t322 = rSges(7,1) * t348 + rSges(7,2) * t347 - rSges(7,3) * t430;
t321 = Icges(7,1) * t350 + Icges(7,4) * t349 - Icges(7,5) * t432;
t320 = Icges(7,1) * t348 + Icges(7,4) * t347 - Icges(7,5) * t430;
t319 = Icges(7,4) * t350 + Icges(7,2) * t349 - Icges(7,6) * t432;
t318 = Icges(7,4) * t348 + Icges(7,2) * t347 - Icges(7,6) * t430;
t317 = Icges(7,5) * t350 + Icges(7,6) * t349 - Icges(7,3) * t432;
t316 = Icges(7,5) * t348 + Icges(7,6) * t347 - Icges(7,3) * t430;
t315 = qJD(1) * t345 - t379 * t427 + t410;
t314 = -t378 * t427 + t397 + (-t344 + t413) * qJD(1);
t313 = -qJD(3) + (t344 * t379 + t345 * t378) * qJD(4);
t312 = t390 * t331 - t358 * t373 + t407;
t311 = -t390 * t330 + t359 * t373 + t408;
t310 = t330 * t358 - t331 * t359 + t411;
t309 = t380 * t323 + t390 * t333 - t336 * t369 - t358 * t365 + t407;
t308 = -t380 * t322 - t390 * t332 + t337 * t369 + t359 * t365 + t408;
t307 = t322 * t336 - t323 * t337 + t332 * t358 - t333 * t359 + t411;
t1 = ((t414 * t378 + t379 * t382) * qJD(1) + (t379 ^ 2 * t339 + (t416 * t378 + (-t338 + t415) * t379) * t378) * qJD(4)) * t376 / 0.2e1 + m(5) * (t313 ^ 2 + t314 ^ 2 + t315 ^ 2) / 0.2e1 + m(3) * (t362 ^ 2 + t363 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 + t334 ^ 2 + t335 ^ 2) / 0.2e1 + m(6) * (t310 ^ 2 + t311 ^ 2 + t312 ^ 2) / 0.2e1 + m(7) * (t307 ^ 2 + t308 ^ 2 + t309 ^ 2) / 0.2e1 - ((-t378 * t382 + t414 * t379) * qJD(1) + (t378 ^ 2 * t338 + (t415 * t379 + (-t339 + t416) * t378) * t379) * qJD(4)) * t377 / 0.2e1 + qJD(1) * ((-t404 * t383 - t402 * t384) * qJD(1) + ((-t341 * t404 - t343 * t402) * t379 - (-t340 * t404 - t342 * t402) * t378) * qJD(4)) / 0.2e1 + t358 * ((-t325 * t432 + t354 * t327 + t355 * t329) * t358 + (-t324 * t432 + t326 * t354 + t328 * t355) * t359 + (t354 * t371 + t355 * t372 - t370 * t432) * t390) / 0.2e1 + t359 * ((-t325 * t430 + t327 * t352 + t329 * t353) * t358 + (-t324 * t430 + t352 * t326 + t353 * t328) * t359 + (t352 * t371 + t353 * t372 - t370 * t430) * t390) / 0.2e1 + t390 * ((t324 * t359 + t325 * t358 + t370 * t390) * t404 + ((t327 * t401 - t329 * t403) * t358 + (t326 * t401 - t328 * t403) * t359 + (t371 * t401 - t372 * t403) * t390) * t402) / 0.2e1 + t337 * ((-t317 * t430 + t319 * t347 + t321 * t348) * t336 + (-t316 * t430 + t347 * t318 + t348 * t320) * t337 + (t347 * t367 + t348 * t368 - t366 * t430) * t380) / 0.2e1 + t380 * ((t316 * t337 + t317 * t336 + t366 * t380) * t404 + ((t319 * t398 - t321 * t399) * t336 + (t318 * t398 - t320 * t399) * t337 + (t367 * t398 - t368 * t399) * t380) * t402) / 0.2e1 + t336 * ((-t317 * t432 + t349 * t319 + t350 * t321) * t336 + (-t316 * t432 + t318 * t349 + t320 * t350) * t337 + (t349 * t367 + t350 * t368 - t366 * t432) * t380) / 0.2e1 + (m(2) * (t387 ^ 2 + t388 ^ 2) + Icges(2,3) + Icges(3,2) + Icges(4,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
