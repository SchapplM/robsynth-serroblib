% Calculate kinetic energy for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:37:46
% EndTime: 2019-12-05 18:37:47
% DurationCPUTime: 1.47s
% Computational Cost: add. (1124->209), mult. (1103->322), div. (0->0), fcn. (944->10), ass. (0->127)
t419 = Icges(4,3) + Icges(5,3);
t353 = qJ(2) + qJ(3);
t344 = pkin(9) + t353;
t339 = sin(t344);
t340 = cos(t344);
t348 = sin(t353);
t349 = cos(t353);
t418 = Icges(4,5) * t349 + Icges(5,5) * t340 - Icges(4,6) * t348 - Icges(5,6) * t339;
t355 = sin(qJ(1));
t357 = cos(qJ(1));
t402 = Icges(5,4) * t340;
t377 = -Icges(5,2) * t339 + t402;
t285 = -Icges(5,6) * t357 + t377 * t355;
t286 = Icges(5,6) * t355 + t377 * t357;
t403 = Icges(5,4) * t339;
t381 = Icges(5,1) * t340 - t403;
t287 = -Icges(5,5) * t357 + t381 * t355;
t288 = Icges(5,5) * t355 + t381 * t357;
t404 = Icges(4,4) * t349;
t378 = -Icges(4,2) * t348 + t404;
t296 = -Icges(4,6) * t357 + t378 * t355;
t297 = Icges(4,6) * t355 + t378 * t357;
t405 = Icges(4,4) * t348;
t382 = Icges(4,1) * t349 - t405;
t298 = -Icges(4,5) * t357 + t382 * t355;
t299 = Icges(4,5) * t355 + t382 * t357;
t316 = Icges(5,2) * t340 + t403;
t317 = Icges(5,1) * t339 + t402;
t323 = Icges(4,2) * t349 + t405;
t324 = Icges(4,1) * t348 + t404;
t347 = qJD(2) * t355;
t328 = qJD(3) * t355 + t347;
t392 = -qJD(2) - qJD(3);
t329 = t392 * t357;
t417 = (-t285 * t339 + t287 * t340 - t296 * t348 + t298 * t349) * t329 + (-t286 * t339 + t288 * t340 - t297 * t348 + t299 * t349) * t328 + (-t316 * t339 + t317 * t340 - t323 * t348 + t324 * t349) * qJD(1);
t416 = (t418 * t355 - t419 * t357) * t329 + (t419 * t355 + t418 * t357) * t328 + (Icges(4,5) * t348 + Icges(5,5) * t339 + Icges(4,6) * t349 + Icges(5,6) * t340) * qJD(1);
t412 = pkin(3) * t348;
t411 = pkin(4) * t339;
t356 = cos(qJ(2));
t409 = t356 * pkin(2);
t354 = sin(qJ(2));
t407 = Icges(3,4) * t354;
t406 = Icges(3,4) * t356;
t342 = qJ(5) + t344;
t337 = sin(t342);
t401 = Icges(6,4) * t337;
t338 = cos(t342);
t400 = Icges(6,4) * t338;
t292 = -pkin(7) * t357 + t409 * t355;
t293 = pkin(7) * t355 + t409 * t357;
t393 = qJD(2) * t357;
t399 = t292 * t347 + t293 * t393;
t336 = pkin(1) * t355 - pkin(6) * t357;
t398 = -t292 - t336;
t397 = pkin(4) * t340;
t396 = pkin(3) * t349;
t391 = pkin(2) * qJD(2) * t354;
t271 = -qJ(4) * t357 + t396 * t355;
t390 = t328 * t271 + t399;
t389 = -t271 + t398;
t388 = t357 * t391;
t387 = rSges(3,1) * t356 - rSges(3,2) * t354;
t386 = rSges(4,1) * t349 - rSges(4,2) * t348;
t385 = rSges(5,1) * t340 - rSges(5,2) * t339;
t384 = rSges(6,1) * t338 - rSges(6,2) * t337;
t383 = Icges(3,1) * t356 - t407;
t380 = Icges(6,1) * t338 - t401;
t379 = -Icges(3,2) * t354 + t406;
t376 = -Icges(6,2) * t337 + t400;
t375 = Icges(3,5) * t356 - Icges(3,6) * t354;
t372 = Icges(6,5) * t338 - Icges(6,6) * t337;
t305 = -Icges(3,6) * t357 + t379 * t355;
t307 = -Icges(3,5) * t357 + t383 * t355;
t371 = t305 * t354 - t307 * t356;
t306 = Icges(3,6) * t355 + t379 * t357;
t308 = Icges(3,5) * t355 + t383 * t357;
t370 = -t306 * t354 + t308 * t356;
t331 = Icges(3,2) * t356 + t407;
t332 = Icges(3,1) * t354 + t406;
t369 = -t331 * t354 + t332 * t356;
t327 = qJD(1) * (pkin(1) * t357 + pkin(6) * t355);
t368 = qJD(1) * t293 - t355 * t391 + t327;
t367 = qJD(4) * t355 + t329 * t412 - t388;
t320 = qJD(5) * t355 + t328;
t321 = (-qJD(5) + t392) * t357;
t366 = qJD(1) * (Icges(6,5) * t337 + Icges(6,6) * t338) + (-Icges(6,3) * t357 + t372 * t355) * t321 + (Icges(6,3) * t355 + t372 * t357) * t320;
t272 = qJ(4) * t355 + t396 * t357;
t363 = qJD(1) * t272 - qJD(4) * t357 + t368;
t275 = -Icges(6,6) * t357 + t376 * t355;
t276 = Icges(6,6) * t355 + t376 * t357;
t277 = -Icges(6,5) * t357 + t380 * t355;
t278 = Icges(6,5) * t355 + t380 * t357;
t310 = Icges(6,2) * t338 + t401;
t311 = Icges(6,1) * t337 + t400;
t362 = (-t276 * t337 + t278 * t338) * t320 + (-t275 * t337 + t277 * t338) * t321 + (-t310 * t337 + t311 * t338) * qJD(1);
t335 = rSges(2,1) * t357 - rSges(2,2) * t355;
t334 = rSges(2,1) * t355 + rSges(2,2) * t357;
t333 = rSges(3,1) * t354 + rSges(3,2) * t356;
t330 = Icges(3,5) * t354 + Icges(3,6) * t356;
t325 = rSges(4,1) * t348 + rSges(4,2) * t349;
t318 = rSges(5,1) * t339 + rSges(5,2) * t340;
t314 = rSges(6,1) * t337 + rSges(6,2) * t338;
t313 = rSges(3,3) * t355 + t387 * t357;
t312 = -rSges(3,3) * t357 + t387 * t355;
t304 = Icges(3,3) * t355 + t375 * t357;
t303 = -Icges(3,3) * t357 + t375 * t355;
t301 = rSges(4,3) * t355 + t386 * t357;
t300 = -rSges(4,3) * t357 + t386 * t355;
t290 = rSges(5,3) * t355 + t385 * t357;
t289 = -rSges(5,3) * t357 + t385 * t355;
t282 = rSges(6,3) * t355 + t384 * t357;
t281 = -rSges(6,3) * t357 + t384 * t355;
t269 = qJD(1) * t313 - t333 * t347 + t327;
t268 = -t333 * t393 + (-t312 - t336) * qJD(1);
t267 = (t312 * t355 + t313 * t357) * qJD(2);
t265 = pkin(8) * t355 + t397 * t357;
t264 = -pkin(8) * t357 + t397 * t355;
t263 = qJD(1) * t301 - t325 * t328 + t368;
t262 = -t388 + t325 * t329 + (-t300 + t398) * qJD(1);
t261 = t300 * t328 - t301 * t329 + t399;
t260 = qJD(1) * t290 + (-t318 - t412) * t328 + t363;
t259 = t318 * t329 + (-t289 + t389) * qJD(1) + t367;
t258 = t289 * t328 + (-t272 - t290) * t329 + t390;
t257 = -t314 * t320 + (-t411 - t412) * t328 + (t265 + t282) * qJD(1) + t363;
t256 = t329 * t411 + t314 * t321 + (-t264 - t281 + t389) * qJD(1) + t367;
t255 = t264 * t328 + t281 * t320 - t282 * t321 + (-t265 - t272) * t329 + t390;
t1 = m(3) * (t267 ^ 2 + t268 ^ 2 + t269 ^ 2) / 0.2e1 + ((t355 * t330 + t369 * t357) * qJD(1) + (t355 ^ 2 * t304 + (t371 * t357 + (-t303 + t370) * t355) * t357) * qJD(2)) * t347 / 0.2e1 - ((-t357 * t330 + t369 * t355) * qJD(1) + (t357 ^ 2 * t303 + (t370 * t355 + (-t304 + t371) * t357) * t355) * qJD(2)) * t393 / 0.2e1 + m(4) * (t261 ^ 2 + t262 ^ 2 + t263 ^ 2) / 0.2e1 + m(5) * (t258 ^ 2 + t259 ^ 2 + t260 ^ 2) / 0.2e1 + m(6) * (t255 ^ 2 + t256 ^ 2 + t257 ^ 2) / 0.2e1 + t320 * (t366 * t355 + t362 * t357) / 0.2e1 + t321 * (t362 * t355 - t366 * t357) / 0.2e1 + (t416 * t355 + t417 * t357) * t328 / 0.2e1 + (t417 * t355 - t416 * t357) * t329 / 0.2e1 + (m(2) * (t334 ^ 2 + t335 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t306 * t356 + t308 * t354) * t355 - (t305 * t356 + t307 * t354) * t357) * qJD(2) + (t276 * t338 + t278 * t337) * t320 + (t275 * t338 + t277 * t337) * t321 + (t285 * t340 + t287 * t339 + t296 * t349 + t298 * t348) * t329 + (t286 * t340 + t288 * t339 + t297 * t349 + t299 * t348) * t328 + (t338 * t310 + t337 * t311 + t340 * t316 + t339 * t317 + t349 * t323 + t348 * t324 + t356 * t331 + t354 * t332) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
