% Calculate kinetic energy for
% S6RPPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:28:14
% EndTime: 2019-03-09 02:28:15
% DurationCPUTime: 1.03s
% Computational Cost: add. (709->207), mult. (1037->332), div. (0->0), fcn. (944->8), ass. (0->112)
t359 = sin(qJ(4));
t404 = pkin(4) * t359;
t360 = sin(qJ(1));
t403 = pkin(7) * t360;
t400 = Icges(5,4) * t359;
t362 = cos(qJ(4));
t399 = Icges(5,4) * t362;
t357 = qJ(4) + qJ(5);
t355 = sin(t357);
t398 = Icges(6,4) * t355;
t356 = cos(t357);
t397 = Icges(6,4) * t356;
t396 = t356 * t360;
t363 = cos(qJ(1));
t395 = t356 * t363;
t358 = sin(qJ(6));
t394 = t358 * t360;
t393 = t358 * t363;
t361 = cos(qJ(6));
t392 = t360 * t361;
t391 = t361 * t363;
t351 = qJD(4) * t363;
t338 = qJD(5) * t363 + t351;
t354 = qJD(2) * t360;
t390 = qJD(3) * t363 + t354;
t389 = qJD(4) * t360;
t388 = qJD(6) * t356;
t387 = pkin(4) * qJD(4) * t362;
t386 = t363 * t387 + t390;
t385 = -qJD(2) * t363 + qJD(1) * (pkin(1) * t363 + qJ(2) * t360);
t337 = (-qJD(4) - qJD(5)) * t360;
t384 = pkin(5) * t355 - pkin(9) * t356;
t383 = rSges(5,1) * t359 + rSges(5,2) * t362;
t382 = rSges(6,1) * t355 + rSges(6,2) * t356;
t381 = Icges(5,1) * t359 + t399;
t380 = Icges(6,1) * t355 + t397;
t379 = Icges(5,2) * t362 + t400;
t378 = Icges(6,2) * t356 + t398;
t377 = Icges(5,5) * t359 + Icges(5,6) * t362;
t376 = Icges(6,5) * t355 + Icges(6,6) * t356;
t314 = Icges(5,6) * t363 + t379 * t360;
t316 = Icges(5,5) * t363 + t381 * t360;
t375 = t314 * t362 + t316 * t359;
t315 = -Icges(5,6) * t360 + t379 * t363;
t317 = -Icges(5,5) * t360 + t381 * t363;
t374 = -t315 * t362 - t317 * t359;
t340 = -Icges(5,2) * t359 + t399;
t341 = Icges(5,1) * t362 - t400;
t373 = t340 * t362 + t341 * t359;
t372 = qJD(1) * t363 * qJ(3) + qJD(3) * t360 + t385;
t342 = pkin(1) * t360 - qJ(2) * t363;
t371 = -pkin(7) * t363 - qJ(3) * t360 - t342;
t325 = -pkin(8) * t360 + t363 * t404;
t326 = pkin(8) * t363 + t360 * t404;
t370 = (-t325 * t363 - t326 * t360) * qJD(4);
t369 = -t326 + t371;
t368 = qJD(1) * t325 + t360 * t387 + t372;
t367 = qJD(1) * (Icges(6,5) * t356 - Icges(6,6) * t355) + (Icges(6,3) * t363 + t376 * t360) * t338 + (-Icges(6,3) * t360 + t376 * t363) * t337;
t306 = Icges(6,6) * t363 + t378 * t360;
t307 = -Icges(6,6) * t360 + t378 * t363;
t308 = Icges(6,5) * t363 + t380 * t360;
t309 = -Icges(6,5) * t360 + t380 * t363;
t332 = -Icges(6,2) * t355 + t397;
t333 = Icges(6,1) * t356 - t398;
t366 = (t307 * t356 + t309 * t355) * t337 + (t306 * t356 + t308 * t355) * t338 + (t332 * t356 + t333 * t355) * qJD(1);
t346 = qJD(6) * t355 + qJD(1);
t345 = rSges(2,1) * t363 - rSges(2,2) * t360;
t344 = rSges(5,1) * t362 - rSges(5,2) * t359;
t343 = rSges(2,1) * t360 + rSges(2,2) * t363;
t339 = Icges(5,5) * t362 - Icges(5,6) * t359;
t335 = pkin(5) * t356 + pkin(9) * t355;
t334 = rSges(6,1) * t356 - rSges(6,2) * t355;
t330 = t355 * t391 - t394;
t329 = -t355 * t393 - t392;
t328 = t355 * t392 + t393;
t327 = -t355 * t394 + t391;
t324 = t384 * t363;
t323 = t384 * t360;
t322 = -rSges(5,3) * t360 + t383 * t363;
t321 = rSges(5,3) * t363 + t383 * t360;
t320 = -t360 * t388 + t338;
t319 = -t363 * t388 + t337;
t313 = -Icges(5,3) * t360 + t377 * t363;
t312 = Icges(5,3) * t363 + t377 * t360;
t311 = -rSges(6,3) * t360 + t382 * t363;
t310 = rSges(6,3) * t363 + t382 * t360;
t303 = rSges(7,3) * t355 + (rSges(7,1) * t361 - rSges(7,2) * t358) * t356;
t302 = Icges(7,5) * t355 + (Icges(7,1) * t361 - Icges(7,4) * t358) * t356;
t301 = Icges(7,6) * t355 + (Icges(7,4) * t361 - Icges(7,2) * t358) * t356;
t300 = Icges(7,3) * t355 + (Icges(7,5) * t361 - Icges(7,6) * t358) * t356;
t299 = qJD(1) * (-rSges(3,2) * t363 + rSges(3,3) * t360) + t385;
t298 = t354 + (rSges(3,2) * t360 + rSges(3,3) * t363 - t342) * qJD(1);
t297 = qJD(1) * (rSges(4,2) * t360 + rSges(4,3) * t363) + t372;
t296 = (t363 * rSges(4,2) - t342 + (-rSges(4,3) - qJ(3)) * t360) * qJD(1) + t390;
t295 = rSges(7,1) * t330 + rSges(7,2) * t329 - rSges(7,3) * t395;
t294 = rSges(7,1) * t328 + rSges(7,2) * t327 - rSges(7,3) * t396;
t293 = Icges(7,1) * t330 + Icges(7,4) * t329 - Icges(7,5) * t395;
t292 = Icges(7,1) * t328 + Icges(7,4) * t327 - Icges(7,5) * t396;
t291 = Icges(7,4) * t330 + Icges(7,2) * t329 - Icges(7,6) * t395;
t290 = Icges(7,4) * t328 + Icges(7,2) * t327 - Icges(7,6) * t396;
t289 = Icges(7,5) * t330 + Icges(7,6) * t329 - Icges(7,3) * t395;
t288 = Icges(7,5) * t328 + Icges(7,6) * t327 - Icges(7,3) * t396;
t287 = (-t321 * t360 - t322 * t363) * qJD(4);
t286 = t344 * t389 + (t322 - t403) * qJD(1) + t372;
t285 = t344 * t351 + (-t321 + t371) * qJD(1) + t390;
t284 = -t334 * t337 + (t311 - t403) * qJD(1) + t368;
t283 = t334 * t338 + (-t310 + t369) * qJD(1) + t386;
t282 = t310 * t337 - t311 * t338 + t370;
t281 = t295 * t346 - t303 * t319 - t335 * t337 + (t324 - t403) * qJD(1) + t368;
t280 = -t294 * t346 + t303 * t320 + t335 * t338 + (-t323 + t369) * qJD(1) + t386;
t279 = t294 * t319 - t295 * t320 + t323 * t337 - t324 * t338 + t370;
t1 = m(5) * (t285 ^ 2 + t286 ^ 2 + t287 ^ 2) / 0.2e1 + m(6) * (t282 ^ 2 + t283 ^ 2 + t284 ^ 2) / 0.2e1 + m(7) * (t279 ^ 2 + t280 ^ 2 + t281 ^ 2) / 0.2e1 + ((t363 * t339 + t373 * t360) * qJD(1) + (t363 ^ 2 * t312 + (t374 * t360 + (-t313 + t375) * t363) * t360) * qJD(4)) * t351 / 0.2e1 + t337 * (-t367 * t360 + t366 * t363) / 0.2e1 + t338 * (t366 * t360 + t367 * t363) / 0.2e1 + t319 * ((-t289 * t395 + t329 * t291 + t330 * t293) * t319 + (-t288 * t395 + t290 * t329 + t292 * t330) * t320 + (-t300 * t395 + t301 * t329 + t302 * t330) * t346) / 0.2e1 + t320 * ((-t289 * t396 + t291 * t327 + t293 * t328) * t319 + (-t288 * t396 + t327 * t290 + t328 * t292) * t320 + (-t300 * t396 + t301 * t327 + t302 * t328) * t346) / 0.2e1 + t346 * ((t288 * t320 + t289 * t319 + t300 * t346) * t355 + ((-t291 * t358 + t293 * t361) * t319 + (-t290 * t358 + t292 * t361) * t320 + (-t301 * t358 + t302 * t361) * t346) * t356) / 0.2e1 + m(3) * (t298 ^ 2 + t299 ^ 2) / 0.2e1 + m(4) * (t296 ^ 2 + t297 ^ 2) / 0.2e1 - ((-t360 * t339 + t373 * t363) * qJD(1) + (t360 ^ 2 * t313 + (t375 * t363 + (-t312 + t374) * t360) * t363) * qJD(4)) * t389 / 0.2e1 + ((-(-t315 * t359 + t317 * t362) * t360 + (-t314 * t359 + t316 * t362) * t363) * qJD(4) + (-t307 * t355 + t309 * t356) * t337 + (-t306 * t355 + t308 * t356) * t338 + (-t355 * t332 + t356 * t333 - t359 * t340 + t362 * t341) * qJD(1)) * qJD(1) / 0.2e1 + (Icges(2,3) + Icges(3,1) + Icges(4,1) + m(2) * (t343 ^ 2 + t345 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
