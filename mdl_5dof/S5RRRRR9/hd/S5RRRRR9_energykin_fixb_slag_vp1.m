% Calculate kinetic energy for
% S5RRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR9_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR9_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR9_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR9_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:27:42
% EndTime: 2019-12-31 22:27:43
% DurationCPUTime: 1.99s
% Computational Cost: add. (1315->283), mult. (1923->467), div. (0->0), fcn. (1932->10), ass. (0->140)
t378 = cos(qJ(3));
t416 = t378 * pkin(3);
t376 = sin(qJ(2));
t414 = Icges(3,4) * t376;
t379 = cos(qJ(2));
t413 = Icges(3,4) * t379;
t375 = sin(qJ(3));
t377 = sin(qJ(1));
t412 = t375 * t377;
t380 = cos(qJ(1));
t411 = t375 * t380;
t410 = t376 * t377;
t409 = t376 * t380;
t408 = t379 * t377;
t407 = t379 * t380;
t374 = qJ(3) + qJ(4);
t397 = pkin(2) * t379 + pkin(7) * t376;
t344 = t397 * t377;
t345 = t397 * t380;
t368 = qJD(2) * t377;
t403 = qJD(2) * t380;
t406 = t344 * t368 + t345 * t403;
t370 = cos(t374);
t405 = pkin(4) * t370;
t402 = qJD(3) * t376;
t346 = t380 * t402 + t368;
t401 = qJD(4) * t376;
t400 = qJD(5) * t376;
t399 = -qJD(3) - qJD(4);
t318 = t380 * t401 + t346;
t369 = sin(t374);
t398 = pkin(4) * t369;
t347 = t377 * t402 - t403;
t319 = t377 * t401 + t347;
t396 = rSges(3,1) * t379 - rSges(3,2) * t376;
t395 = Icges(3,1) * t379 - t414;
t394 = -Icges(3,2) * t376 + t413;
t393 = Icges(3,5) * t379 - Icges(3,6) * t376;
t324 = -Icges(3,6) * t380 + t377 * t394;
t327 = -Icges(3,5) * t380 + t377 * t395;
t392 = t324 * t376 - t327 * t379;
t325 = Icges(3,6) * t377 + t380 * t394;
t328 = Icges(3,5) * t377 + t380 * t395;
t391 = -t325 * t376 + t328 * t379;
t353 = Icges(3,2) * t379 + t414;
t354 = Icges(3,1) * t376 + t413;
t390 = -t353 * t376 + t354 * t379;
t387 = pkin(8) * t376 + t416 * t379;
t299 = -pkin(3) * t411 + t377 * t387;
t300 = pkin(3) * t412 + t380 * t387;
t389 = t346 * t299 - t300 * t347 + t406;
t350 = qJD(1) * (pkin(1) * t380 + pkin(6) * t377);
t358 = pkin(2) * t376 - pkin(7) * t379;
t388 = qJD(1) * t345 - t358 * t368 + t350;
t359 = pkin(1) * t377 - pkin(6) * t380;
t386 = (-t344 - t359) * qJD(1) - t358 * t403;
t385 = pkin(9) * t376 + t405 * t379;
t309 = -pkin(8) * t379 + t416 * t376;
t364 = -qJD(3) * t379 + qJD(1);
t384 = t364 * t300 - t309 * t346 + t388;
t383 = -t299 * t364 + t347 * t309 + t386;
t371 = qJ(5) + t374;
t366 = cos(t371);
t365 = sin(t371);
t357 = rSges(2,1) * t380 - rSges(2,2) * t377;
t356 = rSges(2,1) * t377 + rSges(2,2) * t380;
t355 = rSges(3,1) * t376 + rSges(3,2) * t379;
t352 = Icges(3,5) * t376 + Icges(3,6) * t379;
t349 = t399 * t379 + qJD(1);
t343 = t378 * t407 + t412;
t342 = -t375 * t407 + t377 * t378;
t341 = t378 * t408 - t411;
t340 = -t375 * t408 - t378 * t380;
t339 = qJD(1) + (-qJD(5) + t399) * t379;
t335 = t369 * t377 + t370 * t407;
t334 = -t369 * t407 + t370 * t377;
t333 = -t369 * t380 + t370 * t408;
t332 = -t369 * t408 - t370 * t380;
t331 = rSges(3,3) * t377 + t380 * t396;
t330 = -rSges(3,3) * t380 + t377 * t396;
t329 = -rSges(4,3) * t379 + (rSges(4,1) * t378 - rSges(4,2) * t375) * t376;
t326 = -Icges(4,5) * t379 + (Icges(4,1) * t378 - Icges(4,4) * t375) * t376;
t323 = -Icges(4,6) * t379 + (Icges(4,4) * t378 - Icges(4,2) * t375) * t376;
t322 = Icges(3,3) * t377 + t380 * t393;
t321 = -Icges(3,3) * t380 + t377 * t393;
t320 = -Icges(4,3) * t379 + (Icges(4,5) * t378 - Icges(4,6) * t375) * t376;
t317 = t365 * t377 + t366 * t407;
t316 = -t365 * t407 + t366 * t377;
t315 = -t365 * t380 + t366 * t408;
t314 = -t365 * t408 - t366 * t380;
t313 = -rSges(5,3) * t379 + (rSges(5,1) * t370 - rSges(5,2) * t369) * t376;
t312 = -Icges(5,5) * t379 + (Icges(5,1) * t370 - Icges(5,4) * t369) * t376;
t311 = -Icges(5,6) * t379 + (Icges(5,4) * t370 - Icges(5,2) * t369) * t376;
t310 = -Icges(5,3) * t379 + (Icges(5,5) * t370 - Icges(5,6) * t369) * t376;
t308 = -rSges(6,3) * t379 + (rSges(6,1) * t366 - rSges(6,2) * t365) * t376;
t307 = -Icges(6,5) * t379 + (Icges(6,1) * t366 - Icges(6,4) * t365) * t376;
t306 = -Icges(6,6) * t379 + (Icges(6,4) * t366 - Icges(6,2) * t365) * t376;
t305 = -Icges(6,3) * t379 + (Icges(6,5) * t366 - Icges(6,6) * t365) * t376;
t304 = t377 * t400 + t319;
t303 = t380 * t400 + t318;
t301 = -pkin(9) * t379 + t405 * t376;
t298 = rSges(4,1) * t343 + rSges(4,2) * t342 + rSges(4,3) * t409;
t297 = rSges(4,1) * t341 + rSges(4,2) * t340 + rSges(4,3) * t410;
t296 = Icges(4,1) * t343 + Icges(4,4) * t342 + Icges(4,5) * t409;
t295 = Icges(4,1) * t341 + Icges(4,4) * t340 + Icges(4,5) * t410;
t294 = Icges(4,4) * t343 + Icges(4,2) * t342 + Icges(4,6) * t409;
t293 = Icges(4,4) * t341 + Icges(4,2) * t340 + Icges(4,6) * t410;
t292 = Icges(4,5) * t343 + Icges(4,6) * t342 + Icges(4,3) * t409;
t291 = Icges(4,5) * t341 + Icges(4,6) * t340 + Icges(4,3) * t410;
t290 = qJD(1) * t331 - t355 * t368 + t350;
t289 = -t355 * t403 + (-t330 - t359) * qJD(1);
t288 = (t330 * t377 + t331 * t380) * qJD(2);
t286 = rSges(5,1) * t335 + rSges(5,2) * t334 + rSges(5,3) * t409;
t285 = rSges(5,1) * t333 + rSges(5,2) * t332 + rSges(5,3) * t410;
t284 = Icges(5,1) * t335 + Icges(5,4) * t334 + Icges(5,5) * t409;
t283 = Icges(5,1) * t333 + Icges(5,4) * t332 + Icges(5,5) * t410;
t282 = Icges(5,4) * t335 + Icges(5,2) * t334 + Icges(5,6) * t409;
t281 = Icges(5,4) * t333 + Icges(5,2) * t332 + Icges(5,6) * t410;
t280 = Icges(5,5) * t335 + Icges(5,6) * t334 + Icges(5,3) * t409;
t279 = Icges(5,5) * t333 + Icges(5,6) * t332 + Icges(5,3) * t410;
t278 = rSges(6,1) * t317 + rSges(6,2) * t316 + rSges(6,3) * t409;
t277 = rSges(6,1) * t315 + rSges(6,2) * t314 + rSges(6,3) * t410;
t275 = Icges(6,1) * t317 + Icges(6,4) * t316 + Icges(6,5) * t409;
t274 = Icges(6,1) * t315 + Icges(6,4) * t314 + Icges(6,5) * t410;
t273 = Icges(6,4) * t317 + Icges(6,2) * t316 + Icges(6,6) * t409;
t272 = Icges(6,4) * t315 + Icges(6,2) * t314 + Icges(6,6) * t410;
t271 = Icges(6,5) * t317 + Icges(6,6) * t316 + Icges(6,3) * t409;
t270 = Icges(6,5) * t315 + Icges(6,6) * t314 + Icges(6,3) * t410;
t269 = t377 * t398 + t380 * t385;
t268 = t377 * t385 - t380 * t398;
t267 = t298 * t364 - t329 * t346 + t388;
t266 = -t297 * t364 + t329 * t347 + t386;
t265 = t297 * t346 - t298 * t347 + t406;
t264 = t286 * t349 - t313 * t318 + t384;
t263 = -t285 * t349 + t313 * t319 + t383;
t262 = t285 * t318 - t286 * t319 + t389;
t261 = t269 * t349 + t278 * t339 - t301 * t318 - t303 * t308 + t384;
t260 = -t268 * t349 - t277 * t339 + t301 * t319 + t304 * t308 + t383;
t259 = t268 * t318 - t269 * t319 + t277 * t303 - t278 * t304 + t389;
t1 = m(3) * (t288 ^ 2 + t289 ^ 2 + t290 ^ 2) / 0.2e1 + ((t377 * t352 + t380 * t390) * qJD(1) + (t377 ^ 2 * t322 + (t392 * t380 + (-t321 + t391) * t377) * t380) * qJD(2)) * t368 / 0.2e1 - ((-t380 * t352 + t377 * t390) * qJD(1) + (t380 ^ 2 * t321 + (t391 * t377 + (-t322 + t392) * t380) * t377) * qJD(2)) * t403 / 0.2e1 + qJD(1) * ((t379 * t353 + t376 * t354) * qJD(1) + ((t325 * t379 + t328 * t376) * t377 - (t324 * t379 + t327 * t376) * t380) * qJD(2)) / 0.2e1 + m(4) * (t265 ^ 2 + t266 ^ 2 + t267 ^ 2) / 0.2e1 + t346 * ((t292 * t409 + t294 * t342 + t296 * t343) * t346 + (t291 * t409 + t293 * t342 + t295 * t343) * t347 + (t320 * t409 + t323 * t342 + t326 * t343) * t364) / 0.2e1 + t347 * ((t292 * t410 + t294 * t340 + t296 * t341) * t346 + (t291 * t410 + t293 * t340 + t295 * t341) * t347 + (t320 * t410 + t323 * t340 + t326 * t341) * t364) / 0.2e1 + t364 * ((-t291 * t347 - t292 * t346 - t320 * t364) * t379 + ((-t294 * t375 + t296 * t378) * t346 + (-t293 * t375 + t295 * t378) * t347 + (-t323 * t375 + t326 * t378) * t364) * t376) / 0.2e1 + m(5) * (t262 ^ 2 + t263 ^ 2 + t264 ^ 2) / 0.2e1 + t318 * ((t280 * t409 + t282 * t334 + t284 * t335) * t318 + (t279 * t409 + t281 * t334 + t283 * t335) * t319 + (t310 * t409 + t311 * t334 + t312 * t335) * t349) / 0.2e1 + t319 * ((t280 * t410 + t282 * t332 + t284 * t333) * t318 + (t279 * t410 + t281 * t332 + t283 * t333) * t319 + (t310 * t410 + t311 * t332 + t312 * t333) * t349) / 0.2e1 + t349 * ((-t279 * t319 - t280 * t318 - t310 * t349) * t379 + ((-t282 * t369 + t284 * t370) * t318 + (-t281 * t369 + t283 * t370) * t319 + (-t311 * t369 + t312 * t370) * t349) * t376) / 0.2e1 + m(6) * (t259 ^ 2 + t260 ^ 2 + t261 ^ 2) / 0.2e1 + t303 * ((t271 * t409 + t273 * t316 + t275 * t317) * t303 + (t270 * t409 + t272 * t316 + t274 * t317) * t304 + (t305 * t409 + t306 * t316 + t307 * t317) * t339) / 0.2e1 + t304 * ((t271 * t410 + t273 * t314 + t275 * t315) * t303 + (t270 * t410 + t272 * t314 + t274 * t315) * t304 + (t305 * t410 + t306 * t314 + t307 * t315) * t339) / 0.2e1 + t339 * ((-t270 * t304 - t271 * t303 - t305 * t339) * t379 + ((-t273 * t365 + t275 * t366) * t303 + (-t272 * t365 + t274 * t366) * t304 + (-t306 * t365 + t307 * t366) * t339) * t376) / 0.2e1 + (m(2) * (t356 ^ 2 + t357 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
