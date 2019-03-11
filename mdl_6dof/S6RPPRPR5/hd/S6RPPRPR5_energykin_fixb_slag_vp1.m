% Calculate kinetic energy for
% S6RPPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
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
% Datum: 2019-03-09 01:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRPR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:48:30
% EndTime: 2019-03-09 01:48:31
% DurationCPUTime: 1.37s
% Computational Cost: add. (684->234), mult. (1251->359), div. (0->0), fcn. (1214->8), ass. (0->111)
t359 = cos(pkin(9));
t400 = pkin(5) * t359;
t399 = pkin(7) * qJD(1);
t361 = sin(qJ(4));
t398 = Icges(5,4) * t361;
t363 = cos(qJ(4));
t397 = Icges(5,4) * t363;
t358 = sin(pkin(9));
t362 = sin(qJ(1));
t396 = t358 * t362;
t364 = cos(qJ(1));
t395 = t358 * t364;
t394 = t361 * t362;
t393 = t361 * t364;
t392 = t362 * t363;
t391 = t363 * t364;
t356 = qJD(2) * t362;
t389 = qJD(3) * t364 + t356;
t388 = qJD(4) * t362;
t387 = qJD(4) * t364;
t386 = qJD(5) * t363;
t385 = qJD(6) * t363;
t345 = pkin(4) * t363 + qJ(5) * t361;
t384 = t345 * t387 + t389;
t381 = -qJD(2) * t364 + qJD(1) * (pkin(1) * t364 + qJ(2) * t362);
t380 = rSges(5,1) * t361 + rSges(5,2) * t363;
t379 = pkin(4) * t361 - qJ(5) * t363;
t378 = Icges(5,1) * t361 + t397;
t377 = Icges(5,2) * t363 + t398;
t376 = Icges(5,5) * t361 + Icges(5,6) * t363;
t322 = Icges(5,6) * t364 + t362 * t377;
t324 = Icges(5,5) * t364 + t362 * t378;
t375 = t322 * t363 + t324 * t361;
t323 = -Icges(5,6) * t362 + t364 * t377;
t325 = -Icges(5,5) * t362 + t364 * t378;
t374 = -t323 * t363 - t325 * t361;
t341 = -Icges(5,2) * t361 + t397;
t342 = Icges(5,1) * t363 - t398;
t373 = t341 * t363 + t342 * t361;
t372 = qJD(4) * (pkin(8) * t361 + t400 * t363) - t386;
t371 = qJD(4) * (rSges(6,3) * t361 + (rSges(6,1) * t359 - rSges(6,2) * t358) * t363) - t386;
t370 = qJD(1) * t364 * qJ(3) + qJD(3) * t362 + t381;
t343 = pkin(1) * t362 - qJ(2) * t364;
t369 = -pkin(7) * t364 - qJ(3) * t362 - t343;
t333 = t379 * t362;
t368 = -t333 + t369;
t334 = t379 * t364;
t367 = qJD(1) * t334 + t345 * t388 + t370;
t366 = -pkin(8) * t363 + t400 * t361;
t357 = pkin(9) + qJ(6);
t353 = qJD(5) * t361;
t352 = cos(t357);
t351 = sin(t357);
t348 = qJD(6) * t361 + qJD(1);
t347 = rSges(2,1) * t364 - rSges(2,2) * t362;
t346 = rSges(5,1) * t363 - rSges(5,2) * t361;
t344 = rSges(2,1) * t362 + rSges(2,2) * t364;
t340 = Icges(5,5) * t363 - Icges(5,6) * t361;
t338 = -t362 * t385 + t387;
t337 = -t364 * t385 - t388;
t332 = t359 * t393 - t396;
t331 = -t358 * t393 - t359 * t362;
t330 = t359 * t394 + t395;
t329 = -t358 * t394 + t359 * t364;
t327 = -rSges(5,3) * t362 + t364 * t380;
t326 = rSges(5,3) * t364 + t362 * t380;
t321 = -Icges(5,3) * t362 + t364 * t376;
t320 = Icges(5,3) * t364 + t362 * t376;
t319 = -t351 * t362 + t352 * t393;
t318 = -t351 * t393 - t352 * t362;
t317 = t351 * t364 + t352 * t394;
t316 = -t351 * t394 + t352 * t364;
t314 = Icges(6,5) * t361 + (Icges(6,1) * t359 - Icges(6,4) * t358) * t363;
t313 = Icges(6,6) * t361 + (Icges(6,4) * t359 - Icges(6,2) * t358) * t363;
t312 = Icges(6,3) * t361 + (Icges(6,5) * t359 - Icges(6,6) * t358) * t363;
t311 = rSges(7,3) * t361 + (rSges(7,1) * t352 - rSges(7,2) * t351) * t363;
t310 = Icges(7,5) * t361 + (Icges(7,1) * t352 - Icges(7,4) * t351) * t363;
t309 = Icges(7,6) * t361 + (Icges(7,4) * t352 - Icges(7,2) * t351) * t363;
t308 = Icges(7,3) * t361 + (Icges(7,5) * t352 - Icges(7,6) * t351) * t363;
t306 = qJD(1) * (-rSges(3,2) * t364 + rSges(3,3) * t362) + t381;
t305 = t356 + (rSges(3,2) * t362 + rSges(3,3) * t364 - t343) * qJD(1);
t304 = qJD(1) * (rSges(4,2) * t362 + rSges(4,3) * t364) + t370;
t303 = (t364 * rSges(4,2) - t343 + (-rSges(4,3) - qJ(3)) * t362) * qJD(1) + t389;
t302 = -pkin(5) * t396 + t364 * t366;
t301 = pkin(5) * t395 + t362 * t366;
t300 = rSges(6,1) * t332 + rSges(6,2) * t331 - rSges(6,3) * t391;
t299 = rSges(6,1) * t330 + rSges(6,2) * t329 - rSges(6,3) * t392;
t298 = Icges(6,1) * t332 + Icges(6,4) * t331 - Icges(6,5) * t391;
t297 = Icges(6,1) * t330 + Icges(6,4) * t329 - Icges(6,5) * t392;
t296 = Icges(6,4) * t332 + Icges(6,2) * t331 - Icges(6,6) * t391;
t295 = Icges(6,4) * t330 + Icges(6,2) * t329 - Icges(6,6) * t392;
t294 = Icges(6,5) * t332 + Icges(6,6) * t331 - Icges(6,3) * t391;
t293 = Icges(6,5) * t330 + Icges(6,6) * t329 - Icges(6,3) * t392;
t292 = (-t326 * t362 - t327 * t364) * qJD(4);
t291 = rSges(7,1) * t319 + rSges(7,2) * t318 - rSges(7,3) * t391;
t290 = rSges(7,1) * t317 + rSges(7,2) * t316 - rSges(7,3) * t392;
t289 = Icges(7,1) * t319 + Icges(7,4) * t318 - Icges(7,5) * t391;
t288 = Icges(7,1) * t317 + Icges(7,4) * t316 - Icges(7,5) * t392;
t287 = Icges(7,4) * t319 + Icges(7,2) * t318 - Icges(7,6) * t391;
t286 = Icges(7,4) * t317 + Icges(7,2) * t316 - Icges(7,6) * t392;
t285 = Icges(7,5) * t319 + Icges(7,6) * t318 - Icges(7,3) * t391;
t284 = Icges(7,5) * t317 + Icges(7,6) * t316 - Icges(7,3) * t392;
t283 = t346 * t388 + (-pkin(7) * t362 + t327) * qJD(1) + t370;
t282 = t346 * t387 + (-t326 + t369) * qJD(1) + t389;
t281 = t353 + ((-t300 - t334) * t364 + (-t299 - t333) * t362) * qJD(4);
t280 = qJD(1) * t300 + (t371 - t399) * t362 + t367;
t279 = t371 * t364 + (-t299 + t368) * qJD(1) + t384;
t278 = qJD(1) * t302 + t291 * t348 - t311 * t337 + (t372 - t399) * t362 + t367;
t277 = -t290 * t348 + t311 * t338 + t372 * t364 + (-t301 + t368) * qJD(1) + t384;
t276 = t290 * t337 - t291 * t338 + t353 + ((-t302 - t334) * t364 + (-t301 - t333) * t362) * qJD(4);
t1 = t337 * ((-t285 * t391 + t318 * t287 + t319 * t289) * t337 + (-t284 * t391 + t286 * t318 + t288 * t319) * t338 + (-t308 * t391 + t309 * t318 + t310 * t319) * t348) / 0.2e1 + t338 * ((-t285 * t392 + t287 * t316 + t289 * t317) * t337 + (-t284 * t392 + t316 * t286 + t317 * t288) * t338 + (-t308 * t392 + t309 * t316 + t310 * t317) * t348) / 0.2e1 + t348 * ((t284 * t338 + t285 * t337 + t308 * t348) * t361 + ((-t287 * t351 + t289 * t352) * t337 + (-t286 * t351 + t288 * t352) * t338 + (-t309 * t351 + t310 * t352) * t348) * t363) / 0.2e1 + m(7) * (t276 ^ 2 + t277 ^ 2 + t278 ^ 2) / 0.2e1 + m(6) * (t279 ^ 2 + t280 ^ 2 + t281 ^ 2) / 0.2e1 + m(3) * (t305 ^ 2 + t306 ^ 2) / 0.2e1 + m(4) * (t303 ^ 2 + t304 ^ 2) / 0.2e1 + m(5) * (t282 ^ 2 + t283 ^ 2 + t292 ^ 2) / 0.2e1 + (((t293 * t364 - t294 * t362) * t361 + (-(-t296 * t358 + t298 * t359) * t362 + (-t295 * t358 + t297 * t359) * t364) * t363 - (-t323 * t361 + t325 * t363) * t362 + (-t322 * t361 + t363 * t324) * t364) * qJD(4) + ((-t313 * t358 + t314 * t359 + t342) * t363 + (t312 - t341) * t361) * qJD(1)) * qJD(1) / 0.2e1 - (((-t293 * t391 + t295 * t331 + t297 * t332 + t375 * t364) * t364 + ((-t320 + t374) * t364 + t294 * t391 - t296 * t331 - t298 * t332 + t321 * t362) * t362) * qJD(4) + (-t312 * t391 + t313 * t331 + t314 * t332 - t362 * t340 + t364 * t373) * qJD(1)) * t388 / 0.2e1 + (((-t293 * t392 + t295 * t329 + t297 * t330 + t320 * t364) * t364 + (t294 * t392 - t296 * t329 - t298 * t330 + (-t321 + t375) * t364 + t374 * t362) * t362) * qJD(4) + (-t312 * t392 + t313 * t329 + t314 * t330 + t364 * t340 + t362 * t373) * qJD(1)) * t387 / 0.2e1 + (Icges(4,1) + Icges(2,3) + m(2) * (t344 ^ 2 + t347 ^ 2) + Icges(3,1)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
