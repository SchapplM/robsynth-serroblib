% Calculate kinetic energy for
% S5PPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR3_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR3_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR3_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:21
% EndTime: 2019-12-05 15:16:22
% DurationCPUTime: 1.26s
% Computational Cost: add. (1067->227), mult. (2254->368), div. (0->0), fcn. (2612->10), ass. (0->106)
t364 = cos(qJ(4));
t387 = t364 * pkin(4);
t358 = sin(pkin(9));
t363 = sin(qJ(3));
t385 = t358 * t363;
t365 = cos(qJ(3));
t384 = t358 * t365;
t359 = sin(pkin(8));
t383 = t359 * t358;
t360 = cos(pkin(9));
t382 = t360 * t359;
t362 = sin(qJ(4));
t381 = t360 * t362;
t361 = cos(pkin(8));
t380 = t361 * t358;
t379 = t361 * t363;
t378 = t361 * t365;
t341 = t363 * t382 + t378;
t376 = qJD(3) * t358;
t349 = t359 * t376;
t329 = qJD(4) * t341 + t349;
t343 = -t359 * t365 + t360 * t379;
t350 = t361 * t376;
t330 = qJD(4) * t343 + t350;
t377 = qJD(2) * t361;
t375 = qJD(3) * t360;
t374 = t362 * t383;
t373 = t362 * t380;
t342 = t365 * t382 - t379;
t319 = t342 * pkin(3) + t341 * pkin(6);
t347 = (pkin(3) * t365 + pkin(6) * t363) * t358;
t354 = qJD(2) * t359;
t372 = t319 * t375 + t347 * t349 + t354;
t348 = qJD(4) * t385 - t375;
t344 = t359 * t363 + t360 * t378;
t320 = t344 * pkin(3) + t343 * pkin(6);
t370 = t319 * t350 - t320 * t349 + qJD(1);
t369 = (-t320 * t360 - t347 * t380) * qJD(3) - t377;
t368 = qJD(1) ^ 2;
t357 = qJ(4) + qJ(5);
t356 = cos(t357);
t355 = sin(t357);
t346 = t364 * t384 - t381;
t345 = -t360 * t364 - t362 * t384;
t338 = -t360 * t355 + t356 * t384;
t337 = -t355 * t384 - t360 * t356;
t336 = -t360 * rSges(4,3) + (rSges(4,1) * t365 - rSges(4,2) * t363) * t358;
t335 = -Icges(4,5) * t360 + (Icges(4,1) * t365 - Icges(4,4) * t363) * t358;
t334 = -Icges(4,6) * t360 + (Icges(4,4) * t365 - Icges(4,2) * t363) * t358;
t333 = -Icges(4,3) * t360 + (Icges(4,5) * t365 - Icges(4,6) * t363) * t358;
t332 = qJD(5) * t385 + t348;
t328 = t344 * t364 + t373;
t327 = -t344 * t362 + t364 * t380;
t326 = t342 * t364 + t374;
t325 = -t342 * t362 + t364 * t383;
t324 = t344 * t356 + t355 * t380;
t323 = -t344 * t355 + t356 * t380;
t322 = t342 * t356 + t355 * t383;
t321 = -t342 * t355 + t356 * t383;
t317 = -pkin(4) * t381 + (pkin(7) * t363 + t365 * t387) * t358;
t316 = t346 * rSges(5,1) + t345 * rSges(5,2) + rSges(5,3) * t385;
t315 = Icges(5,1) * t346 + Icges(5,4) * t345 + Icges(5,5) * t385;
t314 = Icges(5,4) * t346 + Icges(5,2) * t345 + Icges(5,6) * t385;
t313 = Icges(5,5) * t346 + Icges(5,6) * t345 + Icges(5,3) * t385;
t311 = t344 * rSges(4,1) - t343 * rSges(4,2) + rSges(4,3) * t380;
t310 = t342 * rSges(4,1) - t341 * rSges(4,2) + rSges(4,3) * t383;
t309 = Icges(4,1) * t344 - Icges(4,4) * t343 + Icges(4,5) * t380;
t308 = Icges(4,1) * t342 - Icges(4,4) * t341 + Icges(4,5) * t383;
t307 = Icges(4,4) * t344 - Icges(4,2) * t343 + Icges(4,6) * t380;
t306 = Icges(4,4) * t342 - Icges(4,2) * t341 + Icges(4,6) * t383;
t305 = Icges(4,5) * t344 - Icges(4,6) * t343 + Icges(4,3) * t380;
t304 = Icges(4,5) * t342 - Icges(4,6) * t341 + Icges(4,3) * t383;
t303 = qJD(5) * t343 + t330;
t302 = qJD(5) * t341 + t329;
t301 = t338 * rSges(6,1) + t337 * rSges(6,2) + rSges(6,3) * t385;
t300 = Icges(6,1) * t338 + Icges(6,4) * t337 + Icges(6,5) * t385;
t299 = Icges(6,4) * t338 + Icges(6,2) * t337 + Icges(6,6) * t385;
t298 = Icges(6,5) * t338 + Icges(6,6) * t337 + Icges(6,3) * t385;
t297 = t328 * rSges(5,1) + t327 * rSges(5,2) + t343 * rSges(5,3);
t296 = t326 * rSges(5,1) + t325 * rSges(5,2) + t341 * rSges(5,3);
t295 = Icges(5,1) * t328 + Icges(5,4) * t327 + Icges(5,5) * t343;
t294 = Icges(5,1) * t326 + Icges(5,4) * t325 + Icges(5,5) * t341;
t293 = Icges(5,4) * t328 + Icges(5,2) * t327 + Icges(5,6) * t343;
t292 = Icges(5,4) * t326 + Icges(5,2) * t325 + Icges(5,6) * t341;
t291 = Icges(5,5) * t328 + Icges(5,6) * t327 + Icges(5,3) * t343;
t290 = Icges(5,5) * t326 + Icges(5,6) * t325 + Icges(5,3) * t341;
t289 = -t377 + (-t311 * t360 - t336 * t380) * qJD(3);
t288 = t354 + (t310 * t360 + t336 * t383) * qJD(3);
t287 = pkin(4) * t373 + pkin(7) * t343 + t344 * t387;
t286 = pkin(4) * t374 + pkin(7) * t341 + t342 * t387;
t285 = t324 * rSges(6,1) + t323 * rSges(6,2) + t343 * rSges(6,3);
t284 = t322 * rSges(6,1) + t321 * rSges(6,2) + t341 * rSges(6,3);
t283 = Icges(6,1) * t324 + Icges(6,4) * t323 + Icges(6,5) * t343;
t282 = Icges(6,1) * t322 + Icges(6,4) * t321 + Icges(6,5) * t341;
t281 = Icges(6,4) * t324 + Icges(6,2) * t323 + Icges(6,6) * t343;
t280 = Icges(6,4) * t322 + Icges(6,2) * t321 + Icges(6,6) * t341;
t279 = Icges(6,5) * t324 + Icges(6,6) * t323 + Icges(6,3) * t343;
t278 = Icges(6,5) * t322 + Icges(6,6) * t321 + Icges(6,3) * t341;
t277 = qJD(1) + (t310 * t361 - t311 * t359) * t376;
t276 = t348 * t297 - t330 * t316 + t369;
t275 = -t348 * t296 + t329 * t316 + t372;
t274 = t330 * t296 - t329 * t297 + t370;
t273 = t332 * t285 + t348 * t287 - t303 * t301 - t330 * t317 + t369;
t272 = -t332 * t284 - t348 * t286 + t302 * t301 + t329 * t317 + t372;
t271 = t303 * t284 - t302 * t285 + t330 * t286 - t329 * t287 + t370;
t1 = m(2) * t368 / 0.2e1 + m(3) * (t368 + (t359 ^ 2 + t361 ^ 2) * qJD(2) ^ 2) / 0.2e1 + m(4) * (t277 ^ 2 + t288 ^ 2 + t289 ^ 2) / 0.2e1 + m(5) * (t274 ^ 2 + t275 ^ 2 + t276 ^ 2) / 0.2e1 + t330 * ((t343 * t291 + t327 * t293 + t328 * t295) * t330 + (t343 * t290 + t327 * t292 + t328 * t294) * t329 + (t343 * t313 + t327 * t314 + t328 * t315) * t348) / 0.2e1 + t329 * ((t341 * t291 + t325 * t293 + t326 * t295) * t330 + (t341 * t290 + t325 * t292 + t326 * t294) * t329 + (t341 * t313 + t325 * t314 + t326 * t315) * t348) / 0.2e1 + t348 * ((t291 * t385 + t345 * t293 + t346 * t295) * t330 + (t290 * t385 + t345 * t292 + t346 * t294) * t329 + (t313 * t385 + t345 * t314 + t346 * t315) * t348) / 0.2e1 + m(6) * (t271 ^ 2 + t272 ^ 2 + t273 ^ 2) / 0.2e1 + t303 * ((t343 * t279 + t323 * t281 + t324 * t283) * t303 + (t343 * t278 + t323 * t280 + t324 * t282) * t302 + (t343 * t298 + t323 * t299 + t324 * t300) * t332) / 0.2e1 + t302 * ((t341 * t279 + t321 * t281 + t322 * t283) * t303 + (t341 * t278 + t321 * t280 + t322 * t282) * t302 + (t341 * t298 + t321 * t299 + t322 * t300) * t332) / 0.2e1 + t332 * ((t279 * t385 + t337 * t281 + t338 * t283) * t303 + (t278 * t385 + t337 * t280 + t338 * t282) * t302 + (t298 * t385 + t337 * t299 + t338 * t300) * t332) / 0.2e1 + (-t360 * (t360 ^ 2 * t333 + (((-t307 * t363 + t309 * t365) * t361 + (-t306 * t363 + t308 * t365) * t359) * t358 + (-t304 * t359 - t305 * t361 + t334 * t363 - t335 * t365) * t360) * t358) / 0.2e1 + (t361 * ((t305 * t380 - t343 * t307 + t344 * t309) * t380 + (t304 * t380 - t343 * t306 + t344 * t308) * t383 - (t333 * t380 - t343 * t334 + t344 * t335) * t360) + t359 * ((t305 * t383 - t341 * t307 + t342 * t309) * t380 + (t304 * t383 - t341 * t306 + t342 * t308) * t383 - (t333 * t383 - t341 * t334 + t342 * t335) * t360)) * t358 / 0.2e1) * qJD(3) ^ 2;
T = t1;
