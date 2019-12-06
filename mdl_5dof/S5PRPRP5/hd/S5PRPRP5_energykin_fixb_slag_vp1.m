% Calculate kinetic energy for
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRP5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP5_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP5_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP5_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:37:33
% EndTime: 2019-12-05 15:37:34
% DurationCPUTime: 1.34s
% Computational Cost: add. (887->161), mult. (1406->258), div. (0->0), fcn. (1445->8), ass. (0->90)
t413 = Icges(5,1) + Icges(6,1);
t412 = Icges(5,4) - Icges(6,5);
t411 = Icges(6,4) + Icges(5,5);
t410 = Icges(5,2) + Icges(6,3);
t409 = Icges(6,6) - Icges(5,6);
t408 = -Icges(5,3) - Icges(6,2);
t352 = cos(pkin(7));
t393 = t352 ^ 2;
t350 = sin(pkin(7));
t394 = t350 ^ 2;
t395 = t393 + t394;
t407 = t395 * qJD(2);
t406 = rSges(6,1) + pkin(4);
t405 = rSges(6,3) + qJ(5);
t348 = pkin(8) + qJ(4);
t346 = sin(t348);
t347 = cos(t348);
t355 = cos(qJ(2));
t385 = t350 * t355;
t325 = t346 * t385 + t347 * t352;
t326 = -t346 * t352 + t347 * t385;
t354 = sin(qJ(2));
t386 = t350 * t354;
t404 = t410 * t325 - t412 * t326 + t409 * t386;
t383 = t352 * t355;
t327 = t346 * t383 - t350 * t347;
t328 = t346 * t350 + t347 * t383;
t384 = t352 * t354;
t403 = t410 * t327 - t412 * t328 + t409 * t384;
t402 = t409 * t325 + t411 * t326 - t408 * t386;
t401 = t409 * t327 + t411 * t328 - t408 * t384;
t400 = -t412 * t325 + t413 * t326 + t411 * t386;
t399 = -t412 * t327 + t413 * t328 + t411 * t384;
t398 = t409 * t355 + (-t410 * t346 + t412 * t347) * t354;
t397 = t408 * t355 + (t409 * t346 + t411 * t347) * t354;
t396 = t411 * t355 + (t412 * t346 - t413 * t347) * t354;
t392 = qJD(2) ^ 2;
t351 = cos(pkin(8));
t389 = pkin(3) * t351;
t349 = sin(pkin(8));
t388 = t349 * t350;
t387 = t349 * t352;
t381 = rSges(6,2) * t386 + t405 * t325 + t406 * t326;
t380 = -rSges(6,2) * t384 - t405 * t327 - t406 * t328;
t379 = -rSges(6,2) * t355 + (t405 * t346 + t406 * t347) * t354;
t378 = qJD(2) * t350;
t377 = qJD(2) * t352;
t376 = qJD(3) * t354;
t375 = qJD(4) * t354;
t374 = qJD(4) * t355;
t340 = pkin(2) * t354 - qJ(3) * t355;
t370 = qJD(2) * (pkin(6) * t355 - t354 * t389 - t340);
t369 = qJD(2) * (rSges(4,3) * t355 - (rSges(4,1) * t351 - rSges(4,2) * t349) * t354 - t340);
t366 = Icges(3,5) * t355 - Icges(3,6) * t354;
t363 = -qJD(3) * t355 + qJD(1) + (pkin(2) * t355 + qJ(3) * t354) * t407;
t343 = t350 * t376;
t362 = t350 * t370 + t343;
t344 = t352 * t376;
t361 = t352 * t370 + t344;
t357 = pkin(6) * t354 + t355 * t389;
t358 = (-pkin(3) * t387 + t350 * t357) * t378 + (pkin(3) * t388 + t352 * t357) * t377 + t363;
t341 = rSges(3,1) * t354 + rSges(3,2) * t355;
t338 = t350 * t375 - t377;
t337 = t352 * t375 + t378;
t336 = t351 * t383 + t388;
t335 = -t349 * t383 + t350 * t351;
t334 = t351 * t385 - t387;
t333 = -t349 * t385 - t351 * t352;
t320 = Icges(3,3) * t350 + t352 * t366;
t319 = -Icges(3,3) * t352 + t350 * t366;
t318 = -rSges(5,3) * t355 + (rSges(5,1) * t347 - rSges(5,2) * t346) * t354;
t307 = Icges(4,1) * t336 + Icges(4,4) * t335 + Icges(4,5) * t384;
t306 = Icges(4,1) * t334 + Icges(4,4) * t333 + Icges(4,5) * t386;
t305 = Icges(4,4) * t336 + Icges(4,2) * t335 + Icges(4,6) * t384;
t304 = Icges(4,4) * t334 + Icges(4,2) * t333 + Icges(4,6) * t386;
t303 = Icges(4,5) * t336 + Icges(4,6) * t335 + Icges(4,3) * t384;
t302 = Icges(4,5) * t334 + Icges(4,6) * t333 + Icges(4,3) * t386;
t301 = t352 * t369 + t344;
t300 = t350 * t369 + t343;
t297 = qJD(1) + (rSges(3,1) * t355 - rSges(3,2) * t354) * t407;
t296 = rSges(5,1) * t328 - rSges(5,2) * t327 + rSges(5,3) * t384;
t294 = rSges(5,1) * t326 - rSges(5,2) * t325 + rSges(5,3) * t386;
t280 = (t350 * (rSges(4,1) * t334 + rSges(4,2) * t333 + rSges(4,3) * t386) + t352 * (rSges(4,1) * t336 + rSges(4,2) * t335 + rSges(4,3) * t384)) * qJD(2) + t363;
t279 = t294 * t374 + t318 * t338 + t361;
t278 = -t296 * t374 - t318 * t337 + t362;
t277 = qJD(5) * t327 + t338 * t379 + t374 * t381 + t361;
t276 = qJD(5) * t325 - t337 * t379 + t374 * t380 + t362;
t275 = t294 * t337 - t296 * t338 + t358;
t274 = qJD(5) * t346 * t354 + t337 * t381 + t338 * t380 + t358;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t395 * t392 * t341 ^ 2 + t297 ^ 2) / 0.2e1 + m(4) * (t280 ^ 2 + t300 ^ 2 + t301 ^ 2) / 0.2e1 + m(5) * (t275 ^ 2 + t278 ^ 2 + t279 ^ 2) / 0.2e1 + m(6) * (t274 ^ 2 + t276 ^ 2 + t277 ^ 2) / 0.2e1 + ((t398 * t327 + t396 * t328 - t397 * t384) * t374 + (t404 * t327 + t400 * t328 + t402 * t384) * t338 + (t403 * t327 + t399 * t328 + t401 * t384) * t337) * t337 / 0.2e1 + ((t398 * t325 + t396 * t326 - t397 * t386) * t374 + (t404 * t325 + t400 * t326 + t402 * t386) * t338 + (t403 * t325 + t399 * t326 + t401 * t386) * t337) * t338 / 0.2e1 - ((-t401 * t337 - t402 * t338 + t397 * t374) * t355 + ((t398 * t346 + t396 * t347) * t374 + (t404 * t346 + t400 * t347) * t338 + (t403 * t346 + t399 * t347) * t337) * t354) * t374 / 0.2e1 + (t394 * t320 + (t303 * t384 + t335 * t305 + t336 * t307) * t350 + (-t302 * t384 - t304 * t335 - t306 * t336 - t350 * t319) * t352) * t350 * t392 / 0.2e1 - (t393 * t319 - (t302 * t386 + t333 * t304 + t334 * t306) * t352 + (t303 * t386 + t305 * t333 + t307 * t334 - t352 * t320) * t350) * t352 * t392 / 0.2e1;
T = t1;
