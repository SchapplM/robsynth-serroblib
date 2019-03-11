% Calculate kinetic energy for
% S6RPPRRR1
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
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:17:55
% EndTime: 2019-03-09 02:17:56
% DurationCPUTime: 1.13s
% Computational Cost: add. (1591->213), mult. (1087->342), div. (0->0), fcn. (992->12), ass. (0->121)
t395 = sin(qJ(1));
t441 = t395 * pkin(1);
t392 = cos(pkin(11));
t439 = t392 * pkin(3);
t389 = pkin(11) + qJ(4);
t382 = sin(t389);
t438 = Icges(5,4) * t382;
t384 = cos(t389);
t437 = Icges(5,4) * t384;
t386 = qJ(5) + t389;
t378 = sin(t386);
t436 = Icges(6,4) * t378;
t379 = cos(t386);
t435 = Icges(6,4) * t379;
t390 = qJ(1) + pkin(10);
t383 = sin(t390);
t434 = t378 * t383;
t385 = cos(t390);
t433 = t378 * t385;
t394 = sin(qJ(6));
t432 = t383 * t394;
t396 = cos(qJ(6));
t431 = t383 * t396;
t430 = t385 * t394;
t429 = t385 * t396;
t397 = cos(qJ(1));
t381 = qJD(1) * t397 * pkin(1);
t427 = qJD(1) * (pkin(2) * t385 + qJ(3) * t383) + t381;
t426 = pkin(4) * t384;
t376 = qJD(4) * t383;
t364 = qJD(5) * t383 + t376;
t424 = qJD(4) * t385;
t423 = qJD(6) * t378;
t422 = pkin(4) * qJD(4) * t382;
t324 = -pkin(8) * t385 + t383 * t426;
t325 = pkin(8) * t383 + t385 * t426;
t421 = t324 * t376 + t325 * t424 + qJD(2);
t420 = -pkin(2) * t383 + qJ(3) * t385 - t441;
t365 = (-qJD(4) - qJD(5)) * t385;
t419 = pkin(7) * t385 - t383 * t439 + t420;
t418 = pkin(5) * t379 + pkin(9) * t378;
t391 = sin(pkin(11));
t417 = rSges(4,1) * t392 - rSges(4,2) * t391;
t416 = rSges(5,1) * t384 - rSges(5,2) * t382;
t415 = rSges(6,1) * t379 - rSges(6,2) * t378;
t414 = Icges(5,1) * t384 - t438;
t413 = Icges(6,1) * t379 - t436;
t412 = -Icges(5,2) * t382 + t437;
t411 = -Icges(6,2) * t378 + t435;
t410 = Icges(5,5) * t384 - Icges(5,6) * t382;
t409 = Icges(6,5) * t379 - Icges(6,6) * t378;
t338 = -Icges(5,6) * t385 + t383 * t412;
t340 = -Icges(5,5) * t385 + t383 * t414;
t408 = t338 * t382 - t340 * t384;
t339 = Icges(5,6) * t383 + t385 * t412;
t341 = Icges(5,5) * t383 + t385 * t414;
t407 = -t339 * t382 + t341 * t384;
t367 = Icges(5,2) * t384 + t438;
t368 = Icges(5,1) * t382 + t437;
t406 = -t367 * t382 + t368 * t384;
t405 = -t324 + t419;
t377 = qJD(3) * t383;
t404 = -t385 * t422 + t377;
t403 = -qJD(3) * t385 + qJD(1) * (pkin(7) * t383 + t385 * t439) + t427;
t402 = (Icges(6,5) * t378 + Icges(6,6) * t379) * qJD(1) + (-Icges(6,3) * t385 + t383 * t409) * t365 + (Icges(6,3) * t383 + t385 * t409) * t364;
t401 = qJD(1) * t325 - t383 * t422 + t403;
t330 = -Icges(6,6) * t385 + t383 * t411;
t331 = Icges(6,6) * t383 + t385 * t411;
t332 = -Icges(6,5) * t385 + t383 * t413;
t333 = Icges(6,5) * t383 + t385 * t413;
t359 = Icges(6,2) * t379 + t436;
t360 = Icges(6,1) * t378 + t435;
t400 = (-t331 * t378 + t333 * t379) * t364 + (-t330 * t378 + t332 * t379) * t365 + (-t359 * t378 + t360 * t379) * qJD(1);
t398 = qJD(2) ^ 2;
t374 = rSges(2,1) * t397 - rSges(2,2) * t395;
t373 = rSges(2,1) * t395 + rSges(2,2) * t397;
t372 = -qJD(6) * t379 + qJD(1);
t369 = rSges(5,1) * t382 + rSges(5,2) * t384;
t366 = Icges(5,5) * t382 + Icges(5,6) * t384;
t362 = pkin(5) * t378 - pkin(9) * t379;
t361 = rSges(6,1) * t378 + rSges(6,2) * t379;
t357 = t381 + qJD(1) * (rSges(3,1) * t385 - rSges(3,2) * t383);
t356 = (-rSges(3,1) * t383 - rSges(3,2) * t385 - t441) * qJD(1);
t355 = t379 * t429 + t432;
t354 = -t379 * t430 + t431;
t353 = t379 * t431 - t430;
t352 = -t379 * t432 - t429;
t351 = t418 * t385;
t350 = t418 * t383;
t349 = t383 * t423 + t365;
t348 = t385 * t423 + t364;
t347 = -rSges(7,3) * t379 + (rSges(7,1) * t396 - rSges(7,2) * t394) * t378;
t346 = -Icges(7,5) * t379 + (Icges(7,1) * t396 - Icges(7,4) * t394) * t378;
t345 = -Icges(7,6) * t379 + (Icges(7,4) * t396 - Icges(7,2) * t394) * t378;
t344 = -Icges(7,3) * t379 + (Icges(7,5) * t396 - Icges(7,6) * t394) * t378;
t343 = rSges(5,3) * t383 + t385 * t416;
t342 = -rSges(5,3) * t385 + t383 * t416;
t337 = Icges(5,3) * t383 + t385 * t410;
t336 = -Icges(5,3) * t385 + t383 * t410;
t335 = rSges(6,3) * t383 + t385 * t415;
t334 = -rSges(6,3) * t385 + t383 * t415;
t320 = qJD(1) * t383 * rSges(4,3) + (qJD(1) * t417 - qJD(3)) * t385 + t427;
t319 = t377 + (t385 * rSges(4,3) - t383 * t417 + t420) * qJD(1);
t318 = rSges(7,1) * t355 + rSges(7,2) * t354 + rSges(7,3) * t433;
t317 = rSges(7,1) * t353 + rSges(7,2) * t352 + rSges(7,3) * t434;
t316 = Icges(7,1) * t355 + Icges(7,4) * t354 + Icges(7,5) * t433;
t315 = Icges(7,1) * t353 + Icges(7,4) * t352 + Icges(7,5) * t434;
t314 = Icges(7,4) * t355 + Icges(7,2) * t354 + Icges(7,6) * t433;
t313 = Icges(7,4) * t353 + Icges(7,2) * t352 + Icges(7,6) * t434;
t312 = Icges(7,5) * t355 + Icges(7,6) * t354 + Icges(7,3) * t433;
t311 = Icges(7,5) * t353 + Icges(7,6) * t352 + Icges(7,3) * t434;
t310 = qJD(2) + (t342 * t383 + t343 * t385) * qJD(4);
t309 = qJD(1) * t343 - t369 * t376 + t403;
t308 = -t369 * t424 + t377 + (-t342 + t419) * qJD(1);
t307 = qJD(1) * t335 - t361 * t364 + t401;
t306 = t361 * t365 + (-t334 + t405) * qJD(1) + t404;
t305 = t334 * t364 - t335 * t365 + t421;
t304 = qJD(1) * t351 + t318 * t372 - t347 * t348 - t362 * t364 + t401;
t303 = -t317 * t372 + t347 * t349 + t362 * t365 + (-t350 + t405) * qJD(1) + t404;
t302 = t317 * t348 - t318 * t349 + t350 * t364 - t351 * t365 + t421;
t1 = ((t383 * t366 + t385 * t406) * qJD(1) + (t383 ^ 2 * t337 + (t408 * t385 + (-t336 + t407) * t383) * t385) * qJD(4)) * t376 / 0.2e1 + t349 * ((t312 * t434 + t314 * t352 + t316 * t353) * t348 + (t311 * t434 + t352 * t313 + t353 * t315) * t349 + (t344 * t434 + t345 * t352 + t346 * t353) * t372) / 0.2e1 + t372 * ((-t311 * t349 - t312 * t348 - t344 * t372) * t379 + ((-t314 * t394 + t316 * t396) * t348 + (-t313 * t394 + t315 * t396) * t349 + (-t345 * t394 + t346 * t396) * t372) * t378) / 0.2e1 + t348 * ((t312 * t433 + t354 * t314 + t355 * t316) * t348 + (t311 * t433 + t313 * t354 + t315 * t355) * t349 + (t344 * t433 + t345 * t354 + t346 * t355) * t372) / 0.2e1 + t364 * (t402 * t383 + t400 * t385) / 0.2e1 + t365 * (t400 * t383 - t402 * t385) / 0.2e1 - ((-t385 * t366 + t383 * t406) * qJD(1) + (t385 ^ 2 * t336 + (t407 * t383 + (-t337 + t408) * t385) * t383) * qJD(4)) * t424 / 0.2e1 + m(7) * (t302 ^ 2 + t303 ^ 2 + t304 ^ 2) / 0.2e1 + m(6) * (t305 ^ 2 + t306 ^ 2 + t307 ^ 2) / 0.2e1 + m(4) * (t319 ^ 2 + t320 ^ 2 + t398) / 0.2e1 + m(5) * (t308 ^ 2 + t309 ^ 2 + t310 ^ 2) / 0.2e1 + m(3) * (t356 ^ 2 + t357 ^ 2 + t398) / 0.2e1 + ((t331 * t379 + t333 * t378) * t364 + (t330 * t379 + t332 * t378) * t365 + ((t339 * t384 + t341 * t382) * t383 - (t338 * t384 + t340 * t382) * t385) * qJD(4) + (t379 * t359 + t378 * t360 + t384 * t367 + t382 * t368) * qJD(1)) * qJD(1) / 0.2e1 + (Icges(4,2) * t392 ^ 2 + (Icges(4,1) * t391 + 0.2e1 * Icges(4,4) * t392) * t391 + Icges(3,3) + Icges(2,3) + m(2) * (t373 ^ 2 + t374 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
