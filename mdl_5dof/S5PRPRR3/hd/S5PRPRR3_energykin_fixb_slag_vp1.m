% Calculate kinetic energy for
% S5PRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR3_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR3_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR3_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:46:35
% EndTime: 2019-12-05 15:46:37
% DurationCPUTime: 1.45s
% Computational Cost: add. (1088->180), mult. (1295->301), div. (0->0), fcn. (1288->10), ass. (0->101)
t421 = Icges(3,3) + Icges(4,3);
t357 = qJ(2) + pkin(9);
t353 = sin(t357);
t354 = cos(t357);
t363 = sin(qJ(2));
t365 = cos(qJ(2));
t420 = Icges(3,5) * t365 + Icges(4,5) * t354 - Icges(3,6) * t363 - Icges(4,6) * t353;
t360 = cos(pkin(8));
t411 = t360 ^ 2;
t359 = sin(pkin(8));
t412 = t359 ^ 2;
t414 = t411 + t412;
t413 = qJD(2) * t414;
t419 = t360 * t359;
t418 = t420 * t359 - t421 * t360;
t417 = t421 * t359 + t420 * t360;
t410 = qJD(2) ^ 2;
t409 = pkin(2) * t363;
t408 = pkin(2) * t365;
t364 = cos(qJ(4));
t407 = pkin(4) * t364;
t404 = t353 * t359;
t403 = t353 * t360;
t358 = qJ(4) + qJ(5);
t355 = sin(t358);
t402 = t355 * t359;
t401 = t355 * t360;
t356 = cos(t358);
t400 = t356 * t359;
t399 = t356 * t360;
t362 = sin(qJ(4));
t398 = t359 * t362;
t397 = t359 * t364;
t396 = t360 * t362;
t395 = t360 * t364;
t352 = qJD(2) * t359;
t392 = qJD(4) * t353;
t340 = t360 * t392 + t352;
t394 = qJD(2) * t360;
t393 = qJD(3) * t360;
t391 = qJD(4) * t354;
t390 = qJD(5) * t353;
t389 = qJD(1) + (-qJ(3) * t360 + t408 * t359) * t352 + (qJ(3) * t359 + t408 * t360) * t394;
t341 = t359 * t392 - t394;
t386 = qJD(2) * (-rSges(4,1) * t353 - rSges(4,2) * t354 - t409);
t385 = qJD(2) * (-pkin(3) * t353 + pkin(6) * t354 - t409);
t384 = t389 + (pkin(3) * t354 + pkin(6) * t353) * t413;
t351 = qJD(3) * t359;
t370 = t360 * t385 + t351;
t369 = pkin(7) * t353 + t407 * t354;
t368 = t359 * t385 - t393;
t348 = rSges(3,1) * t363 + rSges(3,2) * t365;
t342 = (-qJD(4) - qJD(5)) * t354;
t339 = t354 * t395 + t398;
t338 = -t354 * t396 + t397;
t337 = t354 * t397 - t396;
t336 = -t354 * t398 - t395;
t329 = t354 * t399 + t402;
t328 = -t354 * t401 + t400;
t327 = t354 * t400 - t401;
t326 = -t354 * t402 - t399;
t317 = t359 * t390 + t341;
t316 = t360 * t390 + t340;
t315 = -t354 * rSges(5,3) + (rSges(5,1) * t364 - rSges(5,2) * t362) * t353;
t314 = -Icges(5,5) * t354 + (Icges(5,1) * t364 - Icges(5,4) * t362) * t353;
t313 = -Icges(5,6) * t354 + (Icges(5,4) * t364 - Icges(5,2) * t362) * t353;
t312 = -Icges(5,3) * t354 + (Icges(5,5) * t364 - Icges(5,6) * t362) * t353;
t311 = -t354 * rSges(6,3) + (rSges(6,1) * t356 - rSges(6,2) * t355) * t353;
t310 = -Icges(6,5) * t354 + (Icges(6,1) * t356 - Icges(6,4) * t355) * t353;
t309 = -Icges(6,6) * t354 + (Icges(6,4) * t356 - Icges(6,2) * t355) * t353;
t308 = -Icges(6,3) * t354 + (Icges(6,5) * t356 - Icges(6,6) * t355) * t353;
t305 = -pkin(7) * t354 + t407 * t353;
t304 = t360 * t386 + t351;
t303 = t359 * t386 - t393;
t302 = rSges(5,1) * t339 + rSges(5,2) * t338 + rSges(5,3) * t403;
t301 = rSges(5,1) * t337 + rSges(5,2) * t336 + rSges(5,3) * t404;
t300 = Icges(5,1) * t339 + Icges(5,4) * t338 + Icges(5,5) * t403;
t299 = Icges(5,1) * t337 + Icges(5,4) * t336 + Icges(5,5) * t404;
t298 = Icges(5,4) * t339 + Icges(5,2) * t338 + Icges(5,6) * t403;
t297 = Icges(5,4) * t337 + Icges(5,2) * t336 + Icges(5,6) * t404;
t296 = Icges(5,5) * t339 + Icges(5,6) * t338 + Icges(5,3) * t403;
t295 = Icges(5,5) * t337 + Icges(5,6) * t336 + Icges(5,3) * t404;
t294 = pkin(4) * t398 + t369 * t360;
t293 = -pkin(4) * t396 + t369 * t359;
t292 = qJD(1) + (rSges(3,1) * t365 - rSges(3,2) * t363) * t413;
t291 = rSges(6,1) * t329 + rSges(6,2) * t328 + rSges(6,3) * t403;
t290 = rSges(6,1) * t327 + rSges(6,2) * t326 + rSges(6,3) * t404;
t289 = Icges(6,1) * t329 + Icges(6,4) * t328 + Icges(6,5) * t403;
t288 = Icges(6,1) * t327 + Icges(6,4) * t326 + Icges(6,5) * t404;
t287 = Icges(6,4) * t329 + Icges(6,2) * t328 + Icges(6,6) * t403;
t286 = Icges(6,4) * t327 + Icges(6,2) * t326 + Icges(6,6) * t404;
t285 = Icges(6,5) * t329 + Icges(6,6) * t328 + Icges(6,3) * t403;
t284 = Icges(6,5) * t327 + Icges(6,6) * t326 + Icges(6,3) * t404;
t283 = t389 + (rSges(4,1) * t354 - rSges(4,2) * t353) * t413;
t282 = t301 * t391 + t315 * t341 + t370;
t281 = -t302 * t391 - t315 * t340 + t368;
t280 = t301 * t340 - t302 * t341 + t384;
t279 = -t290 * t342 + t293 * t391 + t305 * t341 + t311 * t317 + t370;
t278 = t291 * t342 - t294 * t391 - t305 * t340 - t311 * t316 + t368;
t277 = t290 * t316 - t291 * t317 + t293 * t340 - t294 * t341 + t384;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t348 ^ 2 * t410 * t414 + t292 ^ 2) / 0.2e1 + m(4) * (t283 ^ 2 + t303 ^ 2 + t304 ^ 2) / 0.2e1 + m(5) * (t280 ^ 2 + t281 ^ 2 + t282 ^ 2) / 0.2e1 + t340 * ((t296 * t403 + t338 * t298 + t339 * t300) * t340 + (t295 * t403 + t297 * t338 + t299 * t339) * t341 - (t312 * t403 + t313 * t338 + t314 * t339) * t391) / 0.2e1 + t341 * ((t296 * t404 + t298 * t336 + t300 * t337) * t340 + (t295 * t404 + t336 * t297 + t337 * t299) * t341 - (t312 * t404 + t313 * t336 + t314 * t337) * t391) / 0.2e1 - ((-t295 * t341 - t296 * t340 + t312 * t391) * t354 + ((-t298 * t362 + t300 * t364) * t340 + (-t297 * t362 + t299 * t364) * t341 - (-t313 * t362 + t314 * t364) * t391) * t353) * t391 / 0.2e1 + m(6) * (t277 ^ 2 + t278 ^ 2 + t279 ^ 2) / 0.2e1 + t316 * ((t285 * t403 + t328 * t287 + t329 * t289) * t316 + (t284 * t403 + t286 * t328 + t288 * t329) * t317 + (t308 * t403 + t309 * t328 + t310 * t329) * t342) / 0.2e1 + t317 * ((t285 * t404 + t287 * t326 + t289 * t327) * t316 + (t284 * t404 + t326 * t286 + t327 * t288) * t317 + (t308 * t404 + t309 * t326 + t310 * t327) * t342) / 0.2e1 + t342 * ((-t284 * t317 - t285 * t316 - t308 * t342) * t354 + ((-t287 * t355 + t289 * t356) * t316 + (-t286 * t355 + t288 * t356) * t317 + (-t309 * t355 + t310 * t356) * t342) * t353) / 0.2e1 + (t417 * t412 - t418 * t419) * t359 * t410 / 0.2e1 - (t418 * t411 - t417 * t419) * t360 * t410 / 0.2e1;
T = t1;
