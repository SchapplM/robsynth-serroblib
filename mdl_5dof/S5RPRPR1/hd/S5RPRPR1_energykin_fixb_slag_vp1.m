% Calculate kinetic energy for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:01
% EndTime: 2019-12-05 17:47:02
% DurationCPUTime: 1.14s
% Computational Cost: add. (605->170), mult. (766->260), div. (0->0), fcn. (636->8), ass. (0->99)
t381 = Icges(4,3) + Icges(5,3);
t318 = qJ(3) + pkin(8);
t311 = sin(t318);
t312 = cos(t318);
t320 = sin(qJ(3));
t322 = cos(qJ(3));
t380 = Icges(4,5) * t320 + Icges(5,5) * t311 + Icges(4,6) * t322 + Icges(5,6) * t312;
t321 = sin(qJ(1));
t323 = cos(qJ(1));
t379 = t380 * t321 + t381 * t323;
t378 = t381 * t321 - t380 * t323;
t377 = Icges(4,5) * t322 + Icges(5,5) * t312 - Icges(4,6) * t320 - Icges(5,6) * t311;
t361 = Icges(5,4) * t312;
t293 = -Icges(5,2) * t311 + t361;
t362 = Icges(5,4) * t311;
t294 = Icges(5,1) * t312 - t362;
t363 = Icges(4,4) * t322;
t301 = -Icges(4,2) * t320 + t363;
t364 = Icges(4,4) * t320;
t302 = Icges(4,1) * t322 - t364;
t376 = t293 * t312 + t294 * t311 + t301 * t322 + t302 * t320;
t339 = Icges(5,2) * t312 + t362;
t271 = Icges(5,6) * t321 - t339 * t323;
t342 = Icges(5,1) * t311 + t361;
t273 = Icges(5,5) * t321 - t342 * t323;
t340 = Icges(4,2) * t322 + t364;
t281 = Icges(4,6) * t321 - t340 * t323;
t343 = Icges(4,1) * t320 + t363;
t283 = Icges(4,5) * t321 - t343 * t323;
t375 = t271 * t312 + t273 * t311 + t281 * t322 + t283 * t320;
t270 = Icges(5,6) * t323 + t339 * t321;
t272 = Icges(5,5) * t323 + t342 * t321;
t280 = Icges(4,6) * t323 + t340 * t321;
t282 = Icges(4,5) * t323 + t343 * t321;
t374 = -t270 * t312 - t272 * t311 - t280 * t322 - t282 * t320;
t313 = qJ(5) + t318;
t309 = cos(t313);
t308 = sin(t313);
t360 = Icges(6,4) * t308;
t338 = Icges(6,2) * t309 + t360;
t262 = Icges(6,6) * t323 + t338 * t321;
t263 = Icges(6,6) * t321 - t338 * t323;
t359 = Icges(6,4) * t309;
t341 = Icges(6,1) * t308 + t359;
t264 = Icges(6,5) * t323 + t341 * t321;
t265 = Icges(6,5) * t321 - t341 * t323;
t289 = -Icges(6,2) * t308 + t359;
t290 = Icges(6,1) * t309 - t360;
t354 = qJD(3) + qJD(5);
t298 = t354 * t321;
t299 = t354 * t323;
t373 = (t262 * t309 + t264 * t308) * t299 + (t263 * t309 + t265 * t308) * t298 + (t289 * t309 + t290 * t308) * qJD(1);
t369 = pkin(3) * t320;
t368 = pkin(3) * t322;
t367 = pkin(4) * t312;
t296 = qJD(1) * (pkin(1) * t323 + qJ(2) * t321);
t358 = qJD(1) * t323 * pkin(6) + t296;
t356 = qJD(3) * t321;
t355 = qJD(3) * t323;
t316 = qJD(2) * t321;
t353 = qJD(4) * t323 + t356 * t368 + t316;
t350 = pkin(4) * t311;
t303 = pkin(1) * t321 - qJ(2) * t323;
t349 = -pkin(6) * t321 - t303;
t287 = qJ(4) * t323 + t321 * t369;
t348 = qJD(1) * t287 + qJD(4) * t321 + t358;
t286 = qJ(4) * t321 - t323 * t369;
t347 = -t286 + t349;
t346 = rSges(4,1) * t320 + rSges(4,2) * t322;
t345 = rSges(5,1) * t311 + rSges(5,2) * t312;
t344 = rSges(6,1) * t308 + rSges(6,2) * t309;
t335 = Icges(6,5) * t308 + Icges(6,6) * t309;
t325 = (Icges(6,5) * t309 - Icges(6,6) * t308) * qJD(1) + (Icges(6,3) * t323 + t335 * t321) * t299 + (Icges(6,3) * t321 - t335 * t323) * t298;
t306 = rSges(2,1) * t323 - rSges(2,2) * t321;
t305 = rSges(4,1) * t322 - rSges(4,2) * t320;
t304 = rSges(2,1) * t321 + rSges(2,2) * t323;
t295 = rSges(5,1) * t312 - rSges(5,2) * t311;
t291 = rSges(6,1) * t309 - rSges(6,2) * t308;
t285 = rSges(4,3) * t321 - t346 * t323;
t284 = rSges(4,3) * t323 + t346 * t321;
t276 = t286 * t355;
t275 = rSges(5,3) * t321 - t345 * t323;
t274 = rSges(5,3) * t323 + t345 * t321;
t267 = rSges(6,3) * t321 - t344 * t323;
t266 = rSges(6,3) * t323 + t344 * t321;
t259 = t296 - qJD(2) * t323 + qJD(1) * (-rSges(3,2) * t323 + rSges(3,3) * t321);
t258 = t316 + (rSges(3,2) * t321 + rSges(3,3) * t323 - t303) * qJD(1);
t257 = pkin(7) * t323 + t350 * t321;
t256 = pkin(7) * t321 - t350 * t323;
t255 = (-t284 * t321 + t285 * t323) * qJD(3);
t254 = qJD(1) * t284 + (-qJD(3) * t305 - qJD(2)) * t323 + t358;
t253 = t305 * t356 + t316 + (-t285 + t349) * qJD(1);
t252 = qJD(1) * t274 + (-qJD(2) + (-t295 - t368) * qJD(3)) * t323 + t348;
t251 = t295 * t356 + (-t275 + t347) * qJD(1) + t353;
t250 = t276 + (t275 * t323 + (-t274 - t287) * t321) * qJD(3);
t249 = -t291 * t299 + (t257 + t266) * qJD(1) + (-qJD(2) + (-t367 - t368) * qJD(3)) * t323 + t348;
t248 = t356 * t367 + t291 * t298 + (-t256 - t267 + t347) * qJD(1) + t353;
t247 = -t266 * t298 + t267 * t299 + t276 + (t256 * t323 + (-t257 - t287) * t321) * qJD(3);
t1 = m(3) * (t258 ^ 2 + t259 ^ 2) / 0.2e1 + m(4) * (t253 ^ 2 + t254 ^ 2 + t255 ^ 2) / 0.2e1 + m(5) * (t250 ^ 2 + t251 ^ 2 + t252 ^ 2) / 0.2e1 + m(6) * (t247 ^ 2 + t248 ^ 2 + t249 ^ 2) / 0.2e1 + t299 * (t373 * t321 + t325 * t323) / 0.2e1 + t298 * (t325 * t321 - t373 * t323) / 0.2e1 + ((t378 * t321 ^ 2 + (t374 * t323 + (-t375 + t379) * t321) * t323) * qJD(3) + (t377 * t321 - t376 * t323) * qJD(1)) * t356 / 0.2e1 + ((t379 * t323 ^ 2 + (t375 * t321 + (-t374 + t378) * t323) * t321) * qJD(3) + (t376 * t321 + t377 * t323) * qJD(1)) * t355 / 0.2e1 + (m(2) * (t304 ^ 2 + t306 ^ 2) + Icges(2,3) + Icges(3,1)) * qJD(1) ^ 2 / 0.2e1 + ((-t262 * t308 + t264 * t309) * t299 + (-t263 * t308 + t265 * t309) * t298 + ((-t270 * t311 + t272 * t312 - t280 * t320 + t282 * t322) * t323 + (-t271 * t311 + t273 * t312 - t281 * t320 + t283 * t322) * t321) * qJD(3) + (-t308 * t289 + t309 * t290 - t311 * t293 + t312 * t294 - t320 * t301 + t322 * t302) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
