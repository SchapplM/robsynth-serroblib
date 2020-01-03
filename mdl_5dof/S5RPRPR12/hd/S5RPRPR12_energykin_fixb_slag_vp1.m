% Calculate kinetic energy for
% S5RPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR12_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR12_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR12_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR12_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR12_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:29:24
% EndTime: 2019-12-31 18:29:26
% DurationCPUTime: 1.32s
% Computational Cost: add. (1094->224), mult. (1242->358), div. (0->0), fcn. (1224->10), ass. (0->114)
t347 = cos(pkin(8));
t390 = pkin(2) * t347;
t346 = cos(pkin(9));
t389 = pkin(4) * t346;
t343 = pkin(8) + qJ(3);
t338 = sin(t343);
t388 = Icges(4,4) * t338;
t340 = cos(t343);
t387 = Icges(4,4) * t340;
t350 = sin(qJ(1));
t386 = t338 * t350;
t351 = cos(qJ(1));
t385 = t338 * t351;
t384 = t340 * t350;
t383 = t340 * t351;
t344 = sin(pkin(9));
t382 = t344 * t350;
t381 = t344 * t351;
t380 = t346 * t350;
t379 = t346 * t351;
t331 = pkin(1) * t350 - qJ(2) * t351;
t376 = pkin(6) * t351 - t350 * t390 - t331;
t341 = qJD(2) * t350;
t372 = qJD(4) * t338;
t375 = t351 * t372 + t341;
t374 = qJD(3) * t350;
t373 = qJD(3) * t351;
t371 = qJD(5) * t338;
t361 = pkin(3) * t340 + qJ(4) * t338;
t315 = t361 * t350;
t370 = -t315 + t376;
t326 = pkin(3) * t338 - qJ(4) * t340;
t367 = qJD(3) * (pkin(7) * t340 - t338 * t389 - t326);
t366 = qJD(3) * (t340 * rSges(5,3) - (rSges(5,1) * t346 - rSges(5,2) * t344) * t338 - t326);
t328 = qJD(1) * (pkin(1) * t351 + qJ(2) * t350);
t365 = -qJD(2) * t351 + qJD(1) * (pkin(6) * t350 + t351 * t390) + t328;
t316 = t361 * t351;
t364 = -qJD(4) * t340 + t315 * t374 + t316 * t373;
t345 = sin(pkin(8));
t363 = rSges(3,1) * t347 - rSges(3,2) * t345;
t362 = rSges(4,1) * t340 - rSges(4,2) * t338;
t360 = Icges(4,1) * t340 - t388;
t359 = -Icges(4,2) * t338 + t387;
t358 = Icges(4,5) * t340 - Icges(4,6) * t338;
t302 = -Icges(4,6) * t351 + t350 * t359;
t304 = -Icges(4,5) * t351 + t350 * t360;
t357 = t302 * t338 - t304 * t340;
t303 = Icges(4,6) * t350 + t351 * t359;
t305 = Icges(4,5) * t350 + t351 * t360;
t356 = -t303 * t338 + t305 * t340;
t324 = Icges(4,2) * t340 + t388;
t325 = Icges(4,1) * t338 + t387;
t355 = -t324 * t338 + t325 * t340;
t354 = qJD(1) * t316 + t350 * t372 + t365;
t353 = pkin(7) * t338 + t340 * t389;
t342 = pkin(9) + qJ(5);
t339 = cos(t342);
t337 = sin(t342);
t334 = -qJD(5) * t340 + qJD(1);
t333 = rSges(2,1) * t351 - rSges(2,2) * t350;
t332 = rSges(2,1) * t350 + rSges(2,2) * t351;
t327 = rSges(4,1) * t338 + rSges(4,2) * t340;
t323 = Icges(4,5) * t338 + Icges(4,6) * t340;
t322 = t350 * t371 - t373;
t321 = t351 * t371 + t374;
t320 = t340 * t379 + t382;
t319 = -t340 * t381 + t380;
t318 = t340 * t380 - t381;
t317 = -t340 * t382 - t379;
t313 = t337 * t350 + t339 * t383;
t312 = -t337 * t383 + t339 * t350;
t311 = -t337 * t351 + t339 * t384;
t310 = -t337 * t384 - t339 * t351;
t309 = t350 * rSges(4,3) + t351 * t362;
t308 = -t351 * rSges(4,3) + t350 * t362;
t301 = Icges(4,3) * t350 + t351 * t358;
t300 = -Icges(4,3) * t351 + t350 * t358;
t296 = -Icges(5,5) * t340 + (Icges(5,1) * t346 - Icges(5,4) * t344) * t338;
t295 = -Icges(5,6) * t340 + (Icges(5,4) * t346 - Icges(5,2) * t344) * t338;
t294 = -Icges(5,3) * t340 + (Icges(5,5) * t346 - Icges(5,6) * t344) * t338;
t293 = -t340 * rSges(6,3) + (rSges(6,1) * t339 - rSges(6,2) * t337) * t338;
t292 = -Icges(6,5) * t340 + (Icges(6,1) * t339 - Icges(6,4) * t337) * t338;
t291 = -Icges(6,6) * t340 + (Icges(6,4) * t339 - Icges(6,2) * t337) * t338;
t290 = -Icges(6,3) * t340 + (Icges(6,5) * t339 - Icges(6,6) * t337) * t338;
t288 = qJD(1) * t350 * rSges(3,3) + t328 + (qJD(1) * t363 - qJD(2)) * t351;
t287 = t341 + (t351 * rSges(3,3) - t350 * t363 - t331) * qJD(1);
t286 = rSges(5,1) * t320 + rSges(5,2) * t319 + rSges(5,3) * t385;
t285 = rSges(5,1) * t318 + rSges(5,2) * t317 + rSges(5,3) * t386;
t284 = Icges(5,1) * t320 + Icges(5,4) * t319 + Icges(5,5) * t385;
t283 = Icges(5,1) * t318 + Icges(5,4) * t317 + Icges(5,5) * t386;
t282 = Icges(5,4) * t320 + Icges(5,2) * t319 + Icges(5,6) * t385;
t281 = Icges(5,4) * t318 + Icges(5,2) * t317 + Icges(5,6) * t386;
t280 = Icges(5,5) * t320 + Icges(5,6) * t319 + Icges(5,3) * t385;
t279 = Icges(5,5) * t318 + Icges(5,6) * t317 + Icges(5,3) * t386;
t278 = pkin(4) * t382 + t351 * t353;
t277 = -pkin(4) * t381 + t350 * t353;
t276 = (t308 * t350 + t309 * t351) * qJD(3);
t275 = rSges(6,1) * t313 + rSges(6,2) * t312 + rSges(6,3) * t385;
t274 = rSges(6,1) * t311 + rSges(6,2) * t310 + rSges(6,3) * t386;
t273 = Icges(6,1) * t313 + Icges(6,4) * t312 + Icges(6,5) * t385;
t272 = Icges(6,1) * t311 + Icges(6,4) * t310 + Icges(6,5) * t386;
t271 = Icges(6,4) * t313 + Icges(6,2) * t312 + Icges(6,6) * t385;
t270 = Icges(6,4) * t311 + Icges(6,2) * t310 + Icges(6,6) * t386;
t269 = Icges(6,5) * t313 + Icges(6,6) * t312 + Icges(6,3) * t385;
t268 = Icges(6,5) * t311 + Icges(6,6) * t310 + Icges(6,3) * t386;
t267 = qJD(1) * t309 - t327 * t374 + t365;
t266 = -t327 * t373 + t341 + (-t308 + t376) * qJD(1);
t265 = (t285 * t350 + t286 * t351) * qJD(3) + t364;
t264 = qJD(1) * t286 + t350 * t366 + t354;
t263 = t351 * t366 + (-t285 + t370) * qJD(1) + t375;
t262 = qJD(1) * t278 + t275 * t334 - t293 * t321 + t350 * t367 + t354;
t261 = -t274 * t334 + t293 * t322 + t351 * t367 + (-t277 + t370) * qJD(1) + t375;
t260 = t274 * t321 - t275 * t322 + (t277 * t350 + t278 * t351) * qJD(3) + t364;
t1 = m(3) * (t287 ^ 2 + t288 ^ 2) / 0.2e1 + m(4) * (t266 ^ 2 + t267 ^ 2 + t276 ^ 2) / 0.2e1 + m(5) * (t263 ^ 2 + t264 ^ 2 + t265 ^ 2) / 0.2e1 + m(6) * (t260 ^ 2 + t261 ^ 2 + t262 ^ 2) / 0.2e1 + t321 * ((t269 * t385 + t312 * t271 + t313 * t273) * t321 + (t268 * t385 + t270 * t312 + t272 * t313) * t322 + (t290 * t385 + t291 * t312 + t292 * t313) * t334) / 0.2e1 + t322 * ((t269 * t386 + t271 * t310 + t273 * t311) * t321 + (t268 * t386 + t310 * t270 + t311 * t272) * t322 + (t290 * t386 + t291 * t310 + t292 * t311) * t334) / 0.2e1 + t334 * ((-t268 * t322 - t269 * t321 - t290 * t334) * t340 + ((-t271 * t337 + t273 * t339) * t321 + (-t270 * t337 + t272 * t339) * t322 + (-t291 * t337 + t292 * t339) * t334) * t338) / 0.2e1 + (((t340 * t303 + t305 * t338) * t350 - (t302 * t340 + t304 * t338) * t351 + (t279 * t351 - t280 * t350) * t340 + ((-t282 * t344 + t284 * t346) * t350 - (-t281 * t344 + t283 * t346) * t351) * t338) * qJD(3) + ((t324 - t294) * t340 + (-t295 * t344 + t296 * t346 + t325) * t338) * qJD(1)) * qJD(1) / 0.2e1 + (((-t279 * t385 - t281 * t319 - t283 * t320 + t357 * t351) * t351 + ((-t300 + t356) * t351 + t280 * t385 + t319 * t282 + t320 * t284 + t350 * t301) * t350) * qJD(3) + (t294 * t385 + t319 * t295 + t320 * t296 + t350 * t323 + t351 * t355) * qJD(1)) * t374 / 0.2e1 - (((-t279 * t386 - t317 * t281 - t318 * t283 + t351 * t300) * t351 + ((-t301 + t357) * t351 + t280 * t386 + t282 * t317 + t284 * t318 + t356 * t350) * t350) * qJD(3) + (t294 * t386 + t295 * t317 + t296 * t318 - t351 * t323 + t350 * t355) * qJD(1)) * t373 / 0.2e1 + (m(2) * (t332 ^ 2 + t333 ^ 2) + Icges(2,3) + Icges(3,2) * t347 ^ 2 + (Icges(3,1) * t345 + 0.2e1 * Icges(3,4) * t347) * t345) * qJD(1) ^ 2 / 0.2e1;
T = t1;
