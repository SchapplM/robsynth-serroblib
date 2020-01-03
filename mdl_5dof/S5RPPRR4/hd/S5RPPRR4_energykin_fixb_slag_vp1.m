% Calculate kinetic energy for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:30:24
% EndTime: 2020-01-03 11:30:25
% DurationCPUTime: 1.25s
% Computational Cost: add. (950->195), mult. (1152->322), div. (0->0), fcn. (1161->10), ass. (0->98)
t344 = sin(pkin(8));
t346 = cos(pkin(8));
t342 = pkin(9) + qJ(4);
t338 = cos(t342);
t349 = cos(qJ(1));
t371 = t349 * t338;
t337 = sin(t342);
t348 = sin(qJ(1));
t379 = t348 * t337;
t317 = -t346 * t379 - t371;
t372 = t349 * t337;
t378 = t348 * t338;
t318 = t346 * t378 - t372;
t319 = t346 * t372 - t378;
t320 = -t346 * t371 - t379;
t369 = t349 * t344;
t376 = t348 * t344;
t353 = (Icges(5,5) * t318 + Icges(5,6) * t317 + Icges(5,3) * t376) * t348 - (Icges(5,5) * t320 + Icges(5,6) * t319 - Icges(5,3) * t369) * t349;
t388 = t353 * t344;
t365 = pkin(4) * t338;
t387 = pkin(7) * t344 + t365 * t346;
t345 = cos(pkin(9));
t386 = pkin(3) * t345 * t346 + pkin(6) * t344;
t385 = t346 ^ 2;
t339 = qJ(5) + t342;
t334 = sin(t339);
t381 = t348 * t334;
t335 = cos(t339);
t380 = t348 * t335;
t343 = sin(pkin(9));
t377 = t348 * t343;
t375 = t348 * t345;
t374 = t349 * t334;
t373 = t349 * t335;
t370 = t349 * t343;
t368 = t349 * t345;
t330 = -t349 * pkin(1) - t348 * qJ(2);
t354 = pkin(2) * t346 + qJ(3) * t344;
t366 = t354 * t349 - t330;
t363 = qJD(2) * t349;
t362 = qJD(3) * t346;
t361 = qJD(3) * t349;
t360 = qJD(4) * t344;
t359 = qJD(4) * t348;
t358 = qJD(4) + qJD(5);
t357 = pkin(4) * t337;
t328 = qJD(1) * (t348 * pkin(1) - t349 * qJ(2));
t356 = t328 + (qJD(1) * t354 - qJD(2)) * t348;
t355 = rSges(3,1) * t346 - rSges(3,2) * t344;
t352 = qJD(1) * (-pkin(3) * t370 + t386 * t348) + t356;
t332 = qJD(3) * t376;
t351 = t332 + (pkin(3) * t377 + t386 * t349 + t366) * qJD(1);
t333 = -qJD(4) * t346 + qJD(1);
t331 = -t349 * rSges(2,1) + t348 * rSges(2,2);
t329 = t348 * rSges(2,1) + t349 * rSges(2,2);
t326 = -t358 * t346 + qJD(1);
t323 = t358 * t369;
t322 = t358 * t376;
t316 = -t346 * t373 - t381;
t315 = t346 * t374 - t380;
t314 = t346 * t380 - t374;
t313 = -t346 * t381 - t373;
t312 = -t346 * rSges(5,3) + (rSges(5,1) * t338 - rSges(5,2) * t337) * t344;
t311 = -Icges(5,5) * t346 + (Icges(5,1) * t338 - Icges(5,4) * t337) * t344;
t310 = -Icges(5,6) * t346 + (Icges(5,4) * t338 - Icges(5,2) * t337) * t344;
t309 = -Icges(5,3) * t346 + (Icges(5,5) * t338 - Icges(5,6) * t337) * t344;
t308 = -t346 * rSges(6,3) + (rSges(6,1) * t335 - rSges(6,2) * t334) * t344;
t307 = -Icges(6,5) * t346 + (Icges(6,1) * t335 - Icges(6,4) * t334) * t344;
t306 = -Icges(6,6) * t346 + (Icges(6,4) * t335 - Icges(6,2) * t334) * t344;
t305 = -Icges(6,3) * t346 + (Icges(6,5) * t335 - Icges(6,6) * t334) * t344;
t304 = -t363 + (t348 * rSges(3,3) + t355 * t349 - t330) * qJD(1);
t303 = -qJD(1) * t349 * rSges(3,3) + t328 + (qJD(1) * t355 - qJD(2)) * t348;
t302 = -pkin(7) * t346 + t365 * t344;
t299 = t320 * rSges(5,1) + t319 * rSges(5,2) - rSges(5,3) * t369;
t298 = t318 * rSges(5,1) + t317 * rSges(5,2) + rSges(5,3) * t376;
t297 = Icges(5,1) * t320 + Icges(5,4) * t319 - Icges(5,5) * t369;
t296 = Icges(5,1) * t318 + Icges(5,4) * t317 + Icges(5,5) * t376;
t295 = Icges(5,4) * t320 + Icges(5,2) * t319 - Icges(5,6) * t369;
t294 = Icges(5,4) * t318 + Icges(5,2) * t317 + Icges(5,6) * t376;
t291 = t316 * rSges(6,1) + t315 * rSges(6,2) - rSges(6,3) * t369;
t290 = t314 * rSges(6,1) + t313 * rSges(6,2) + rSges(6,3) * t376;
t289 = Icges(6,1) * t316 + Icges(6,4) * t315 - Icges(6,5) * t369;
t288 = Icges(6,1) * t314 + Icges(6,4) * t313 + Icges(6,5) * t376;
t287 = Icges(6,4) * t316 + Icges(6,2) * t315 - Icges(6,6) * t369;
t286 = Icges(6,4) * t314 + Icges(6,2) * t313 + Icges(6,6) * t376;
t285 = Icges(6,5) * t316 + Icges(6,6) * t315 - Icges(6,3) * t369;
t284 = Icges(6,5) * t314 + Icges(6,6) * t313 + Icges(6,3) * t376;
t283 = -t357 * t348 - t387 * t349;
t282 = t387 * t348 - t357 * t349;
t281 = -t363 + t332 + (-(-t346 * t368 - t377) * rSges(4,1) - (t346 * t370 - t375) * rSges(4,2) + rSges(4,3) * t369 + t366) * qJD(1);
t280 = -t344 * t361 + qJD(1) * ((t346 * t375 - t370) * rSges(4,1) + (-t346 * t377 - t368) * rSges(4,2) + rSges(4,3) * t376) + t356;
t279 = -t362 + (t298 * t349 + t299 * t348) * t360;
t278 = -t333 * t299 + (-t312 * t360 - qJD(2)) * t349 + t351;
t277 = t333 * t298 + (-t312 * t359 - t361) * t344 + t352;
t276 = -t362 + t323 * t290 + t322 * t291 + (t282 * t349 + t283 * t348) * t360;
t275 = -t333 * t283 - t326 * t291 - t323 * t308 + (-t302 * t360 - qJD(2)) * t349 + t351;
t274 = t333 * t282 + t326 * t290 - t322 * t308 + (-t302 * t359 - t361) * t344 + t352;
t1 = m(3) * (t303 ^ 2 + t304 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 * t385 + t280 ^ 2 + t281 ^ 2) / 0.2e1 + m(5) * (t277 ^ 2 + t278 ^ 2 + t279 ^ 2) / 0.2e1 + t333 * ((-t346 * t309 + (-t310 * t337 + t311 * t338) * t344) * t333 + (((-t294 * t337 + t296 * t338) * t348 - (-t295 * t337 + t297 * t338) * t349) * t344 - t353 * t346) * t360) / 0.2e1 + t344 * ((t309 * t376 + t317 * t310 + t318 * t311) * t333 + (-(t317 * t295 + t318 * t297) * t349 + (t317 * t294 + t318 * t296 + t388) * t348) * t360) * t359 / 0.2e1 - t349 * ((-t309 * t369 + t319 * t310 + t320 * t311) * t333 + ((t319 * t294 + t320 * t296) * t348 + (-t319 * t295 - t320 * t297 - t388) * t349) * t360) * t360 / 0.2e1 + m(6) * (t274 ^ 2 + t275 ^ 2 + t276 ^ 2) / 0.2e1 + t326 * ((-t284 * t322 + t285 * t323 - t305 * t326) * t346 + ((-t306 * t334 + t307 * t335) * t326 + (-t286 * t334 + t288 * t335) * t322 - (-t287 * t334 + t289 * t335) * t323) * t344) / 0.2e1 + t322 * ((t305 * t376 + t313 * t306 + t314 * t307) * t326 + (t284 * t376 + t313 * t286 + t314 * t288) * t322 - (t285 * t376 + t313 * t287 + t314 * t289) * t323) / 0.2e1 - t323 * ((-t305 * t369 + t315 * t306 + t316 * t307) * t326 + (-t284 * t369 + t315 * t286 + t316 * t288) * t322 - (-t285 * t369 + t315 * t287 + t316 * t289) * t323) / 0.2e1 + (m(2) * (t329 ^ 2 + t331 ^ 2) + Icges(2,3) + (Icges(3,2) + Icges(4,3)) * t385 + ((Icges(3,1) + Icges(4,1) * t345 ^ 2 + (-0.2e1 * Icges(4,4) * t345 + Icges(4,2) * t343) * t343) * t344 + 0.2e1 * (-Icges(4,5) * t345 + Icges(4,6) * t343 + Icges(3,4)) * t346) * t344) * qJD(1) ^ 2 / 0.2e1;
T = t1;
