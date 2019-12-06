% Calculate kinetic energy for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR7_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR7_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR7_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:59:42
% EndTime: 2019-12-05 15:59:44
% DurationCPUTime: 1.52s
% Computational Cost: add. (673->180), mult. (1290->305), div. (0->0), fcn. (1283->8), ass. (0->94)
t409 = Icges(4,1) + Icges(3,3);
t358 = sin(qJ(2));
t360 = cos(qJ(2));
t408 = (-Icges(4,4) + Icges(3,5)) * t360 + (Icges(4,5) - Icges(3,6)) * t358;
t356 = cos(pkin(8));
t399 = t356 ^ 2;
t355 = sin(pkin(8));
t400 = t355 ^ 2;
t402 = t399 + t400;
t401 = qJD(2) * t402;
t407 = t356 * t355;
t406 = t408 * t355 - t409 * t356;
t405 = t409 * t355 + t408 * t356;
t398 = qJD(2) ^ 2;
t359 = cos(qJ(4));
t397 = pkin(4) * t359;
t395 = t355 * t358;
t394 = t355 * t360;
t393 = t356 * t358;
t392 = t356 * t360;
t357 = sin(qJ(4));
t391 = t357 * t358;
t390 = t358 * t359;
t351 = qJD(2) * t355;
t386 = qJD(4) * t360;
t339 = t356 * t386 + t351;
t389 = qJD(2) * t356;
t388 = qJD(3) * t358;
t387 = qJD(4) * t358;
t385 = qJD(5) * t360;
t343 = pkin(2) * t358 - qJ(3) * t360;
t382 = qJD(2) * (rSges(4,2) * t358 + rSges(4,3) * t360 - t343);
t340 = t355 * t386 - t389;
t381 = qJD(2) * (-pkin(6) * t358 - t343);
t370 = -qJD(3) * t360 + qJD(1) + (pkin(2) * t360 + qJ(3) * t358) * t401;
t366 = pkin(4) * t391 + pkin(7) * t360;
t348 = t355 * t388;
t365 = t355 * t381 + t348;
t349 = t356 * t388;
t364 = t356 * t381 + t349;
t363 = (-pkin(3) * t356 + pkin(6) * t394) * t351 + (pkin(3) * t355 + pkin(6) * t392) * t389 + t370;
t354 = qJ(4) + qJ(5);
t353 = cos(t354);
t352 = sin(t354);
t345 = rSges(3,1) * t358 + rSges(3,2) * t360;
t341 = (qJD(4) + qJD(5)) * t358;
t338 = t355 * t391 - t356 * t359;
t337 = t355 * t390 + t356 * t357;
t336 = t355 * t359 + t356 * t391;
t335 = -t355 * t357 + t356 * t390;
t332 = -pkin(4) * t357 * t360 + pkin(7) * t358;
t331 = rSges(5,3) * t358 + (-rSges(5,1) * t357 - rSges(5,2) * t359) * t360;
t330 = Icges(5,5) * t358 + (-Icges(5,1) * t357 - Icges(5,4) * t359) * t360;
t329 = Icges(5,6) * t358 + (-Icges(5,4) * t357 - Icges(5,2) * t359) * t360;
t328 = Icges(5,3) * t358 + (-Icges(5,5) * t357 - Icges(5,6) * t359) * t360;
t327 = t352 * t395 - t353 * t356;
t326 = t352 * t356 + t353 * t395;
t325 = t352 * t393 + t353 * t355;
t324 = -t352 * t355 + t353 * t393;
t309 = t355 * t385 + t340;
t308 = t356 * t385 + t339;
t307 = rSges(6,3) * t358 + (-rSges(6,1) * t352 - rSges(6,2) * t353) * t360;
t306 = Icges(6,5) * t358 + (-Icges(6,1) * t352 - Icges(6,4) * t353) * t360;
t305 = Icges(6,6) * t358 + (-Icges(6,4) * t352 - Icges(6,2) * t353) * t360;
t304 = Icges(6,3) * t358 + (-Icges(6,5) * t352 - Icges(6,6) * t353) * t360;
t303 = t356 * t382 + t349;
t302 = t355 * t382 + t348;
t301 = t355 * t366 - t356 * t397;
t300 = t355 * t397 + t356 * t366;
t299 = rSges(5,1) * t338 + rSges(5,2) * t337 + rSges(5,3) * t394;
t298 = rSges(5,1) * t336 + rSges(5,2) * t335 + rSges(5,3) * t392;
t297 = Icges(5,1) * t338 + Icges(5,4) * t337 + Icges(5,5) * t394;
t296 = Icges(5,1) * t336 + Icges(5,4) * t335 + Icges(5,5) * t392;
t295 = Icges(5,4) * t338 + Icges(5,2) * t337 + Icges(5,6) * t394;
t294 = Icges(5,4) * t336 + Icges(5,2) * t335 + Icges(5,6) * t392;
t293 = Icges(5,5) * t338 + Icges(5,6) * t337 + Icges(5,3) * t394;
t292 = Icges(5,5) * t336 + Icges(5,6) * t335 + Icges(5,3) * t392;
t291 = rSges(6,1) * t327 + rSges(6,2) * t326 + rSges(6,3) * t394;
t290 = rSges(6,1) * t325 + rSges(6,2) * t324 + rSges(6,3) * t392;
t289 = Icges(6,1) * t327 + Icges(6,4) * t326 + Icges(6,5) * t394;
t288 = Icges(6,1) * t325 + Icges(6,4) * t324 + Icges(6,5) * t392;
t287 = Icges(6,4) * t327 + Icges(6,2) * t326 + Icges(6,6) * t394;
t286 = Icges(6,4) * t325 + Icges(6,2) * t324 + Icges(6,6) * t392;
t285 = Icges(6,5) * t327 + Icges(6,6) * t326 + Icges(6,3) * t394;
t284 = Icges(6,5) * t325 + Icges(6,6) * t324 + Icges(6,3) * t392;
t283 = qJD(1) + (rSges(3,1) * t360 - rSges(3,2) * t358) * t401;
t282 = t370 + (-rSges(4,2) * t360 + rSges(4,3) * t358) * t401;
t281 = -t299 * t387 + t331 * t340 + t364;
t280 = t298 * t387 - t331 * t339 + t365;
t279 = -t298 * t340 + t299 * t339 + t363;
t278 = -t291 * t341 - t301 * t387 + t307 * t309 + t332 * t340 + t364;
t277 = t290 * t341 + t300 * t387 - t307 * t308 - t332 * t339 + t365;
t276 = -t290 * t309 + t291 * t308 - t300 * t340 + t301 * t339 + t363;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t402 * t398 * t345 ^ 2 + t283 ^ 2) / 0.2e1 + m(4) * (t282 ^ 2 + t302 ^ 2 + t303 ^ 2) / 0.2e1 + m(5) * (t279 ^ 2 + t280 ^ 2 + t281 ^ 2) / 0.2e1 + t339 * ((t292 * t392 + t335 * t294 + t336 * t296) * t339 + (t293 * t392 + t295 * t335 + t297 * t336) * t340 + (t328 * t392 + t329 * t335 + t330 * t336) * t387) / 0.2e1 + t340 * ((t292 * t394 + t294 * t337 + t296 * t338) * t339 + (t293 * t394 + t337 * t295 + t338 * t297) * t340 + (t328 * t394 + t329 * t337 + t330 * t338) * t387) / 0.2e1 + ((t292 * t339 + t293 * t340 + t328 * t387) * t358 + ((-t294 * t359 - t296 * t357) * t339 + (-t295 * t359 - t297 * t357) * t340 + (-t329 * t359 - t330 * t357) * t387) * t360) * t387 / 0.2e1 + m(6) * (t276 ^ 2 + t277 ^ 2 + t278 ^ 2) / 0.2e1 + t308 * ((t284 * t392 + t324 * t286 + t325 * t288) * t308 + (t285 * t392 + t287 * t324 + t289 * t325) * t309 + (t304 * t392 + t305 * t324 + t306 * t325) * t341) / 0.2e1 + t309 * ((t284 * t394 + t286 * t326 + t288 * t327) * t308 + (t285 * t394 + t326 * t287 + t327 * t289) * t309 + (t304 * t394 + t305 * t326 + t306 * t327) * t341) / 0.2e1 + t341 * ((t284 * t308 + t285 * t309 + t304 * t341) * t358 + ((-t286 * t353 - t288 * t352) * t308 + (-t287 * t353 - t289 * t352) * t309 + (-t305 * t353 - t306 * t352) * t341) * t360) / 0.2e1 + (t405 * t400 - t406 * t407) * t355 * t398 / 0.2e1 - (t406 * t399 - t405 * t407) * t356 * t398 / 0.2e1;
T = t1;
