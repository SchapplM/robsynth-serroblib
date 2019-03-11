% Calculate kinetic energy for
% S6RPPRRR6
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
% Datum: 2019-03-09 02:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:30:30
% EndTime: 2019-03-09 02:30:32
% DurationCPUTime: 1.27s
% Computational Cost: add. (720->231), mult. (1311->374), div. (0->0), fcn. (1274->8), ass. (0->115)
t368 = sin(qJ(1));
t405 = pkin(7) * t368;
t369 = cos(qJ(5));
t404 = pkin(5) * t369;
t367 = sin(qJ(4));
t402 = Icges(5,4) * t367;
t370 = cos(qJ(4));
t401 = Icges(5,4) * t370;
t366 = sin(qJ(5));
t400 = t366 * t368;
t371 = cos(qJ(1));
t399 = t366 * t371;
t398 = t367 * t368;
t397 = t367 * t371;
t396 = t368 * t369;
t395 = t368 * t370;
t394 = t369 * t371;
t393 = t370 * t371;
t362 = qJD(2) * t368;
t392 = qJD(3) * t371 + t362;
t391 = qJD(4) * t368;
t359 = qJD(4) * t371;
t390 = qJD(5) * t370;
t355 = qJD(5) * t367 + qJD(1);
t389 = -qJD(2) * t371 + qJD(1) * (pkin(1) * t371 + qJ(2) * t368);
t388 = (-qJD(5) - qJD(6)) * t370;
t387 = pkin(4) * t367 - pkin(8) * t370;
t386 = rSges(5,1) * t367 + rSges(5,2) * t370;
t385 = Icges(5,1) * t367 + t401;
t384 = Icges(5,2) * t370 + t402;
t383 = Icges(5,5) * t367 + Icges(5,6) * t370;
t322 = Icges(5,6) * t371 + t368 * t384;
t325 = Icges(5,5) * t371 + t368 * t385;
t382 = t322 * t370 + t325 * t367;
t323 = -Icges(5,6) * t368 + t371 * t384;
t326 = -Icges(5,5) * t368 + t371 * t385;
t381 = -t323 * t370 - t326 * t367;
t348 = -Icges(5,2) * t367 + t401;
t349 = Icges(5,1) * t370 - t402;
t380 = t348 * t370 + t349 * t367;
t379 = qJD(1) * t371 * qJ(3) + qJD(3) * t368 + t389;
t350 = pkin(1) * t368 - qJ(2) * t371;
t378 = -pkin(7) * t371 - qJ(3) * t368 - t350;
t339 = t387 * t368;
t340 = t387 * t371;
t377 = (-t339 * t368 - t340 * t371) * qJD(4);
t376 = -pkin(9) * t370 + t367 * t404;
t354 = pkin(4) * t370 + pkin(8) * t367;
t375 = t354 * t391 + t379 + (t340 - t405) * qJD(1);
t374 = t354 * t359 + (-t339 + t378) * qJD(1) + t392;
t365 = qJ(5) + qJ(6);
t364 = cos(t365);
t363 = sin(t365);
t353 = rSges(2,1) * t371 - rSges(2,2) * t368;
t352 = rSges(5,1) * t370 - rSges(5,2) * t367;
t351 = rSges(2,1) * t368 + rSges(2,2) * t371;
t347 = Icges(5,5) * t370 - Icges(5,6) * t367;
t345 = qJD(6) * t367 + t355;
t344 = -t368 * t390 + t359;
t343 = -t371 * t390 - t391;
t338 = t367 * t394 - t400;
t337 = -t366 * t397 - t396;
t336 = t367 * t396 + t399;
t335 = -t366 * t398 + t394;
t333 = -t363 * t368 + t364 * t397;
t332 = -t363 * t397 - t364 * t368;
t331 = t363 * t371 + t364 * t398;
t330 = -t363 * t398 + t364 * t371;
t329 = -rSges(5,3) * t368 + t371 * t386;
t328 = rSges(6,3) * t367 + (rSges(6,1) * t369 - rSges(6,2) * t366) * t370;
t327 = rSges(5,3) * t371 + t368 * t386;
t324 = Icges(6,5) * t367 + (Icges(6,1) * t369 - Icges(6,4) * t366) * t370;
t321 = Icges(6,6) * t367 + (Icges(6,4) * t369 - Icges(6,2) * t366) * t370;
t320 = -Icges(5,3) * t368 + t371 * t383;
t319 = Icges(5,3) * t371 + t368 * t383;
t318 = Icges(6,3) * t367 + (Icges(6,5) * t369 - Icges(6,6) * t366) * t370;
t317 = t368 * t388 + t359;
t316 = t371 * t388 - t391;
t315 = rSges(7,3) * t367 + (rSges(7,1) * t364 - rSges(7,2) * t363) * t370;
t314 = Icges(7,5) * t367 + (Icges(7,1) * t364 - Icges(7,4) * t363) * t370;
t313 = Icges(7,6) * t367 + (Icges(7,4) * t364 - Icges(7,2) * t363) * t370;
t312 = Icges(7,3) * t367 + (Icges(7,5) * t364 - Icges(7,6) * t363) * t370;
t311 = pkin(9) * t367 + t370 * t404;
t310 = qJD(1) * (-rSges(3,2) * t371 + rSges(3,3) * t368) + t389;
t309 = t362 + (rSges(3,2) * t368 + rSges(3,3) * t371 - t350) * qJD(1);
t308 = qJD(1) * (rSges(4,2) * t368 + rSges(4,3) * t371) + t379;
t307 = (t371 * rSges(4,2) - t350 + (-rSges(4,3) - qJ(3)) * t368) * qJD(1) + t392;
t306 = -pkin(5) * t400 + t371 * t376;
t305 = pkin(5) * t399 + t368 * t376;
t304 = rSges(6,1) * t338 + rSges(6,2) * t337 - rSges(6,3) * t393;
t303 = rSges(6,1) * t336 + rSges(6,2) * t335 - rSges(6,3) * t395;
t302 = Icges(6,1) * t338 + Icges(6,4) * t337 - Icges(6,5) * t393;
t301 = Icges(6,1) * t336 + Icges(6,4) * t335 - Icges(6,5) * t395;
t300 = Icges(6,4) * t338 + Icges(6,2) * t337 - Icges(6,6) * t393;
t299 = Icges(6,4) * t336 + Icges(6,2) * t335 - Icges(6,6) * t395;
t298 = Icges(6,5) * t338 + Icges(6,6) * t337 - Icges(6,3) * t393;
t297 = Icges(6,5) * t336 + Icges(6,6) * t335 - Icges(6,3) * t395;
t296 = (-t327 * t368 - t329 * t371) * qJD(4);
t295 = rSges(7,1) * t333 + rSges(7,2) * t332 - rSges(7,3) * t393;
t294 = rSges(7,1) * t331 + rSges(7,2) * t330 - rSges(7,3) * t395;
t293 = Icges(7,1) * t333 + Icges(7,4) * t332 - Icges(7,5) * t393;
t292 = Icges(7,1) * t331 + Icges(7,4) * t330 - Icges(7,5) * t395;
t291 = Icges(7,4) * t333 + Icges(7,2) * t332 - Icges(7,6) * t393;
t290 = Icges(7,4) * t331 + Icges(7,2) * t330 - Icges(7,6) * t395;
t289 = Icges(7,5) * t333 + Icges(7,6) * t332 - Icges(7,3) * t393;
t288 = Icges(7,5) * t331 + Icges(7,6) * t330 - Icges(7,3) * t395;
t287 = t352 * t391 + (t329 - t405) * qJD(1) + t379;
t286 = t352 * t359 + (-t327 + t378) * qJD(1) + t392;
t285 = t304 * t355 - t328 * t343 + t375;
t284 = -t303 * t355 + t328 * t344 + t374;
t283 = t303 * t343 - t304 * t344 + t377;
t282 = t295 * t345 + t306 * t355 - t311 * t343 - t315 * t316 + t375;
t281 = -t294 * t345 - t305 * t355 + t311 * t344 + t315 * t317 + t374;
t280 = t294 * t316 - t295 * t317 + t305 * t343 - t306 * t344 + t377;
t1 = m(7) * (t280 ^ 2 + t281 ^ 2 + t282 ^ 2) / 0.2e1 + m(6) * (t283 ^ 2 + t284 ^ 2 + t285 ^ 2) / 0.2e1 + m(5) * (t286 ^ 2 + t287 ^ 2 + t296 ^ 2) / 0.2e1 + m(3) * (t309 ^ 2 + t310 ^ 2) / 0.2e1 + m(4) * (t307 ^ 2 + t308 ^ 2) / 0.2e1 + t316 * ((-t289 * t393 + t332 * t291 + t333 * t293) * t316 + (-t288 * t393 + t290 * t332 + t292 * t333) * t317 + (-t312 * t393 + t313 * t332 + t314 * t333) * t345) / 0.2e1 + t317 * ((-t289 * t395 + t291 * t330 + t293 * t331) * t316 + (-t288 * t395 + t330 * t290 + t331 * t292) * t317 + (-t312 * t395 + t313 * t330 + t314 * t331) * t345) / 0.2e1 + t345 * ((t288 * t317 + t289 * t316 + t312 * t345) * t367 + ((-t291 * t363 + t293 * t364) * t316 + (-t290 * t363 + t292 * t364) * t317 + (-t313 * t363 + t314 * t364) * t345) * t370) / 0.2e1 + t343 * ((-t298 * t393 + t337 * t300 + t338 * t302) * t343 + (-t297 * t393 + t299 * t337 + t301 * t338) * t344 + (-t318 * t393 + t321 * t337 + t324 * t338) * t355) / 0.2e1 + t344 * ((-t298 * t395 + t300 * t335 + t302 * t336) * t343 + (-t297 * t395 + t335 * t299 + t336 * t301) * t344 + (-t318 * t395 + t321 * t335 + t324 * t336) * t355) / 0.2e1 + t355 * ((t297 * t344 + t298 * t343 + t318 * t355) * t367 + ((-t300 * t366 + t302 * t369) * t343 + (-t299 * t366 + t301 * t369) * t344 + (-t321 * t366 + t324 * t369) * t355) * t370) / 0.2e1 + ((t371 * t347 + t368 * t380) * qJD(1) + (t371 ^ 2 * t319 + (t381 * t368 + (-t320 + t382) * t371) * t368) * qJD(4)) * t359 / 0.2e1 + qJD(1) * ((-t367 * t348 + t370 * t349) * qJD(1) + (-(-t323 * t367 + t326 * t370) * t368 + (-t322 * t367 + t325 * t370) * t371) * qJD(4)) / 0.2e1 - ((-t368 * t347 + t371 * t380) * qJD(1) + (t368 ^ 2 * t320 + (t382 * t371 + (-t319 + t381) * t368) * t371) * qJD(4)) * t391 / 0.2e1 + (m(2) * (t351 ^ 2 + t353 ^ 2) + Icges(2,3) + Icges(3,1) + Icges(4,1)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
