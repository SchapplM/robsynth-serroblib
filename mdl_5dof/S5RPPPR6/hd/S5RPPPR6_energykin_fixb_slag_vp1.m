% Calculate kinetic energy for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPPR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR6_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR6_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR6_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:36
% EndTime: 2019-12-31 17:47:36
% DurationCPUTime: 0.55s
% Computational Cost: add. (395->136), mult. (968->216), div. (0->0), fcn. (1069->8), ass. (0->70)
t330 = cos(pkin(7));
t360 = t330 ^ 2;
t329 = cos(pkin(8));
t357 = t329 * t330;
t332 = sin(qJ(1));
t356 = t330 * t332;
t333 = cos(qJ(5));
t355 = t330 * t333;
t334 = cos(qJ(1));
t354 = t330 * t334;
t327 = sin(pkin(8));
t353 = t332 * t327;
t352 = t332 * t329;
t351 = t334 * t327;
t350 = t334 * t329;
t317 = t332 * pkin(1) - t334 * qJ(2);
t328 = sin(pkin(7));
t337 = pkin(2) * t330 + qJ(3) * t328;
t349 = -t337 * t332 - t317;
t326 = qJD(2) * t332;
t346 = qJD(3) * t328;
t348 = t334 * t346 + t326;
t347 = qJD(1) * t332;
t345 = qJD(4) * t330;
t307 = -t328 * t350 + t353;
t344 = qJD(5) * t307;
t309 = t328 * t352 + t351;
t343 = qJD(5) * t309;
t314 = qJD(1) * (t334 * pkin(1) + t332 * qJ(2));
t342 = qJD(1) * t337 * t334 + t332 * t346 + t314;
t341 = t334 * pkin(3) - qJ(4) * t356 + t349;
t340 = t334 * t345 + t348;
t315 = -qJD(3) * t330 + qJD(4) * t328;
t339 = rSges(3,1) * t330 - rSges(3,2) * t328;
t338 = rSges(4,2) * t330 - rSges(4,3) * t328;
t336 = -qJD(2) * t334 + qJD(1) * (t332 * pkin(3) + qJ(4) * t354) + t332 * t345 + t342;
t331 = sin(qJ(5));
t319 = t334 * rSges(2,1) - t332 * rSges(2,2);
t318 = t332 * rSges(2,1) + t334 * rSges(2,2);
t316 = qJD(5) * t357 + qJD(1);
t310 = t328 * t353 - t350;
t308 = t328 * t351 + t352;
t306 = -t327 * t355 + t328 * t331;
t305 = t330 * t327 * t331 + t328 * t333;
t303 = t310 * t333 + t331 * t356;
t302 = -t310 * t331 + t332 * t355;
t301 = t308 * t333 + t331 * t354;
t300 = -t308 * t331 + t333 * t354;
t299 = rSges(3,3) * t347 + t314 + (qJD(1) * t339 - qJD(2)) * t334;
t298 = t326 + (t334 * rSges(3,3) - t339 * t332 - t317) * qJD(1);
t297 = t306 * rSges(6,1) + t305 * rSges(6,2) + rSges(6,3) * t357;
t296 = Icges(6,1) * t306 + Icges(6,4) * t305 + Icges(6,5) * t357;
t295 = Icges(6,4) * t306 + Icges(6,2) * t305 + Icges(6,6) * t357;
t294 = Icges(6,5) * t306 + Icges(6,6) * t305 + Icges(6,3) * t357;
t293 = rSges(4,1) * t347 + (-qJD(1) * t338 - qJD(2)) * t334 + t342;
t292 = (t334 * rSges(4,1) + t338 * t332 + t349) * qJD(1) + t348;
t291 = t303 * rSges(6,1) + t302 * rSges(6,2) - t309 * rSges(6,3);
t290 = t301 * rSges(6,1) + t300 * rSges(6,2) + t307 * rSges(6,3);
t289 = Icges(6,1) * t303 + Icges(6,4) * t302 - Icges(6,5) * t309;
t288 = Icges(6,1) * t301 + Icges(6,4) * t300 + Icges(6,5) * t307;
t287 = Icges(6,4) * t303 + Icges(6,2) * t302 - Icges(6,6) * t309;
t286 = Icges(6,4) * t301 + Icges(6,2) * t300 + Icges(6,6) * t307;
t285 = Icges(6,5) * t303 + Icges(6,6) * t302 - Icges(6,3) * t309;
t284 = Icges(6,5) * t301 + Icges(6,6) * t300 + Icges(6,3) * t307;
t283 = qJD(1) * (t308 * rSges(5,1) - t307 * rSges(5,2) + rSges(5,3) * t354) + t336;
t282 = (-t310 * rSges(5,1) - t309 * rSges(5,2) - rSges(5,3) * t356 + t341) * qJD(1) + t340;
t281 = (t290 * t309 + t291 * t307) * qJD(5) + t315;
t280 = qJD(1) * (t308 * pkin(4) + t307 * pkin(6)) - t297 * t344 + t316 * t290 + t336;
t279 = -t297 * t343 - t316 * t291 + (-t310 * pkin(4) + t309 * pkin(6) + t341) * qJD(1) + t340;
t1 = m(3) * (t298 ^ 2 + t299 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 * t360 + t292 ^ 2 + t293 ^ 2) / 0.2e1 + m(5) * (t282 ^ 2 + t283 ^ 2 + t315 ^ 2) / 0.2e1 + m(6) * (t279 ^ 2 + t280 ^ 2 + t281 ^ 2) / 0.2e1 + ((t307 * t294 + t300 * t295 + t301 * t296) * t316 + ((t307 * t284 + t300 * t286 + t301 * t288) * t307 - (t307 * t285 + t300 * t287 + t301 * t289) * t309) * qJD(5)) * t344 / 0.2e1 - ((-t309 * t294 + t302 * t295 + t303 * t296) * t316 + ((-t309 * t284 + t302 * t286 + t303 * t288) * t307 - (-t309 * t285 + t302 * t287 + t303 * t289) * t309) * qJD(5)) * t343 / 0.2e1 + t316 * ((t294 * t357 + t305 * t295 + t306 * t296) * t316 + ((t284 * t357 + t305 * t286 + t306 * t288) * t307 - (t285 * t357 + t305 * t287 + t306 * t289) * t309) * qJD(5)) / 0.2e1 + (m(2) * (t318 ^ 2 + t319 ^ 2) + Icges(2,3) + ((Icges(3,1) + Icges(4,2) + Icges(5,3)) * t328 + 0.2e1 * (-Icges(5,5) * t327 - Icges(5,6) * t329 + Icges(3,4) + Icges(4,6)) * t330) * t328 + (t329 ^ 2 * Icges(5,2) + (Icges(5,1) * t327 + 0.2e1 * Icges(5,4) * t329) * t327 + Icges(3,2) + Icges(4,3)) * t360) * qJD(1) ^ 2 / 0.2e1;
T = t1;
