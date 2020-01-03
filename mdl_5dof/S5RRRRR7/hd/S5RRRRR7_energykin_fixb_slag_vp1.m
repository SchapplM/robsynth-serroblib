% Calculate kinetic energy for
% S5RRRRR7
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
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR7_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR7_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR7_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR7_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:21:20
% EndTime: 2019-12-31 22:21:22
% DurationCPUTime: 1.44s
% Computational Cost: add. (1241->233), mult. (1321->384), div. (0->0), fcn. (1218->10), ass. (0->137)
t359 = qJ(2) + qJ(3);
t354 = sin(t359);
t417 = pkin(3) * t354;
t364 = cos(qJ(2));
t415 = t364 * pkin(2);
t361 = sin(qJ(2));
t413 = Icges(3,4) * t361;
t412 = Icges(3,4) * t364;
t411 = Icges(4,4) * t354;
t355 = cos(t359);
t410 = Icges(4,4) * t355;
t356 = qJ(4) + t359;
t348 = sin(t356);
t409 = Icges(5,4) * t348;
t349 = cos(t356);
t408 = Icges(5,4) * t349;
t362 = sin(qJ(1));
t407 = t348 * t362;
t365 = cos(qJ(1));
t406 = t348 * t365;
t360 = sin(qJ(5));
t405 = t360 * t362;
t404 = t360 * t365;
t363 = cos(qJ(5));
t403 = t362 * t363;
t402 = t363 * t365;
t300 = -pkin(7) * t365 + t362 * t415;
t301 = pkin(7) * t362 + t365 * t415;
t353 = qJD(2) * t362;
t397 = qJD(2) * t365;
t401 = t300 * t353 + t301 * t397;
t347 = pkin(1) * t362 - pkin(6) * t365;
t400 = -t300 - t347;
t399 = pkin(3) * t355;
t338 = qJD(3) * t362 + t353;
t396 = qJD(5) * t348;
t395 = -qJD(2) - qJD(3);
t394 = pkin(2) * qJD(2) * t361;
t281 = -pkin(8) * t365 + t362 * t399;
t393 = -t281 + t400;
t330 = qJD(4) * t362 + t338;
t392 = t365 * t394;
t391 = pkin(4) * t349 + pkin(9) * t348;
t390 = rSges(3,1) * t364 - rSges(3,2) * t361;
t389 = rSges(4,1) * t355 - rSges(4,2) * t354;
t388 = rSges(5,1) * t349 - rSges(5,2) * t348;
t331 = (-qJD(4) + t395) * t365;
t387 = Icges(3,1) * t364 - t413;
t386 = Icges(4,1) * t355 - t411;
t385 = Icges(5,1) * t349 - t409;
t384 = -Icges(3,2) * t361 + t412;
t383 = -Icges(4,2) * t354 + t410;
t382 = -Icges(5,2) * t348 + t408;
t381 = Icges(3,5) * t364 - Icges(3,6) * t361;
t380 = Icges(4,5) * t355 - Icges(4,6) * t354;
t379 = Icges(5,5) * t349 - Icges(5,6) * t348;
t314 = -Icges(3,6) * t365 + t362 * t384;
t316 = -Icges(3,5) * t365 + t362 * t387;
t378 = t314 * t361 - t316 * t364;
t315 = Icges(3,6) * t362 + t365 * t384;
t317 = Icges(3,5) * t362 + t365 * t387;
t377 = -t315 * t361 + t317 * t364;
t342 = Icges(3,2) * t364 + t413;
t343 = Icges(3,1) * t361 + t412;
t376 = -t342 * t361 + t343 * t364;
t339 = t395 * t365;
t375 = t339 * t417 - t392;
t282 = pkin(8) * t362 + t365 * t399;
t374 = t338 * t281 - t282 * t339 + t401;
t337 = qJD(1) * (pkin(1) * t365 + pkin(6) * t362);
t373 = qJD(1) * t301 - t362 * t394 + t337;
t372 = (Icges(5,5) * t348 + Icges(5,6) * t349) * qJD(1) + (-Icges(5,3) * t365 + t362 * t379) * t331 + (Icges(5,3) * t362 + t365 * t379) * t330;
t371 = (Icges(4,5) * t354 + Icges(4,6) * t355) * qJD(1) + (-Icges(4,3) * t365 + t362 * t380) * t339 + (Icges(4,3) * t362 + t365 * t380) * t338;
t370 = qJD(1) * t282 - t338 * t417 + t373;
t292 = -Icges(5,6) * t365 + t362 * t382;
t293 = Icges(5,6) * t362 + t365 * t382;
t294 = -Icges(5,5) * t365 + t362 * t385;
t295 = Icges(5,5) * t362 + t365 * t385;
t326 = Icges(5,2) * t349 + t409;
t327 = Icges(5,1) * t348 + t408;
t369 = (-t293 * t348 + t295 * t349) * t330 + (-t292 * t348 + t294 * t349) * t331 + (-t326 * t348 + t327 * t349) * qJD(1);
t304 = -Icges(4,6) * t365 + t362 * t383;
t305 = Icges(4,6) * t362 + t365 * t383;
t306 = -Icges(4,5) * t365 + t362 * t386;
t307 = Icges(4,5) * t362 + t365 * t386;
t333 = Icges(4,2) * t355 + t411;
t334 = Icges(4,1) * t354 + t410;
t368 = (-t305 * t354 + t307 * t355) * t338 + (-t304 * t354 + t306 * t355) * t339 + (-t333 * t354 + t334 * t355) * qJD(1);
t346 = rSges(2,1) * t365 - rSges(2,2) * t362;
t345 = rSges(2,1) * t362 + rSges(2,2) * t365;
t344 = rSges(3,1) * t361 + rSges(3,2) * t364;
t341 = Icges(3,5) * t361 + Icges(3,6) * t364;
t340 = -qJD(5) * t349 + qJD(1);
t335 = rSges(4,1) * t354 + rSges(4,2) * t355;
t329 = pkin(4) * t348 - pkin(9) * t349;
t328 = rSges(5,1) * t348 + rSges(5,2) * t349;
t323 = t349 * t402 + t405;
t322 = -t349 * t404 + t403;
t321 = t349 * t403 - t404;
t320 = -t349 * t405 - t402;
t319 = rSges(3,3) * t362 + t365 * t390;
t318 = -rSges(3,3) * t365 + t362 * t390;
t313 = Icges(3,3) * t362 + t365 * t381;
t312 = -Icges(3,3) * t365 + t362 * t381;
t311 = t391 * t365;
t310 = t391 * t362;
t309 = rSges(4,3) * t362 + t365 * t389;
t308 = -rSges(4,3) * t365 + t362 * t389;
t299 = rSges(5,3) * t362 + t365 * t388;
t298 = -rSges(5,3) * t365 + t362 * t388;
t297 = t362 * t396 + t331;
t296 = t365 * t396 + t330;
t288 = -rSges(6,3) * t349 + (rSges(6,1) * t363 - rSges(6,2) * t360) * t348;
t285 = -Icges(6,5) * t349 + (Icges(6,1) * t363 - Icges(6,4) * t360) * t348;
t284 = -Icges(6,6) * t349 + (Icges(6,4) * t363 - Icges(6,2) * t360) * t348;
t283 = -Icges(6,3) * t349 + (Icges(6,5) * t363 - Icges(6,6) * t360) * t348;
t279 = qJD(1) * t319 - t344 * t353 + t337;
t278 = -t344 * t397 + (-t318 - t347) * qJD(1);
t276 = (t318 * t362 + t319 * t365) * qJD(2);
t275 = rSges(6,1) * t323 + rSges(6,2) * t322 + rSges(6,3) * t406;
t274 = rSges(6,1) * t321 + rSges(6,2) * t320 + rSges(6,3) * t407;
t273 = Icges(6,1) * t323 + Icges(6,4) * t322 + Icges(6,5) * t406;
t272 = Icges(6,1) * t321 + Icges(6,4) * t320 + Icges(6,5) * t407;
t271 = Icges(6,4) * t323 + Icges(6,2) * t322 + Icges(6,6) * t406;
t270 = Icges(6,4) * t321 + Icges(6,2) * t320 + Icges(6,6) * t407;
t269 = Icges(6,5) * t323 + Icges(6,6) * t322 + Icges(6,3) * t406;
t268 = Icges(6,5) * t321 + Icges(6,6) * t320 + Icges(6,3) * t407;
t267 = qJD(1) * t309 - t335 * t338 + t373;
t266 = -t392 + t335 * t339 + (-t308 + t400) * qJD(1);
t265 = t308 * t338 - t309 * t339 + t401;
t264 = qJD(1) * t299 - t328 * t330 + t370;
t263 = t328 * t331 + (-t298 + t393) * qJD(1) + t375;
t262 = t298 * t330 - t299 * t331 + t374;
t261 = qJD(1) * t311 + t275 * t340 - t288 * t296 - t329 * t330 + t370;
t260 = -t274 * t340 + t288 * t297 + t329 * t331 + (-t310 + t393) * qJD(1) + t375;
t259 = t274 * t296 - t275 * t297 + t310 * t330 - t311 * t331 + t374;
t1 = m(3) * (t276 ^ 2 + t278 ^ 2 + t279 ^ 2) / 0.2e1 + ((t362 * t341 + t365 * t376) * qJD(1) + (t362 ^ 2 * t313 + (t378 * t365 + (-t312 + t377) * t362) * t365) * qJD(2)) * t353 / 0.2e1 - ((-t365 * t341 + t362 * t376) * qJD(1) + (t365 ^ 2 * t312 + (t377 * t362 + (-t313 + t378) * t365) * t362) * qJD(2)) * t397 / 0.2e1 + m(4) * (t265 ^ 2 + t266 ^ 2 + t267 ^ 2) / 0.2e1 + t338 * (t371 * t362 + t368 * t365) / 0.2e1 + t339 * (t368 * t362 - t371 * t365) / 0.2e1 + m(5) * (t262 ^ 2 + t263 ^ 2 + t264 ^ 2) / 0.2e1 + t330 * (t372 * t362 + t369 * t365) / 0.2e1 + t331 * (t369 * t362 - t372 * t365) / 0.2e1 + m(6) * (t259 ^ 2 + t260 ^ 2 + t261 ^ 2) / 0.2e1 + t296 * ((t269 * t406 + t322 * t271 + t323 * t273) * t296 + (t268 * t406 + t270 * t322 + t272 * t323) * t297 + (t283 * t406 + t284 * t322 + t285 * t323) * t340) / 0.2e1 + t297 * ((t269 * t407 + t271 * t320 + t273 * t321) * t296 + (t268 * t407 + t320 * t270 + t321 * t272) * t297 + (t283 * t407 + t284 * t320 + t285 * t321) * t340) / 0.2e1 + t340 * ((-t268 * t297 - t269 * t296 - t283 * t340) * t349 + ((-t271 * t360 + t273 * t363) * t296 + (-t270 * t360 + t272 * t363) * t297 + (-t284 * t360 + t285 * t363) * t340) * t348) / 0.2e1 + (m(2) * (t345 ^ 2 + t346 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t315 * t364 + t317 * t361) * t362 - (t314 * t364 + t316 * t361) * t365) * qJD(2) + (t305 * t355 + t307 * t354) * t338 + (t304 * t355 + t306 * t354) * t339 + (t293 * t349 + t295 * t348) * t330 + (t292 * t349 + t294 * t348) * t331 + (t349 * t326 + t348 * t327 + t355 * t333 + t354 * t334 + t364 * t342 + t361 * t343) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
