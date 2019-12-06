% Calculate kinetic energy for
% S5RRRPR3
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
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR3_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR3_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:42:17
% EndTime: 2019-12-05 18:42:18
% DurationCPUTime: 1.03s
% Computational Cost: add. (1014->166), mult. (764->258), div. (0->0), fcn. (634->10), ass. (0->103)
t382 = Icges(4,3) + Icges(5,3);
t316 = qJ(3) + pkin(9);
t308 = sin(t316);
t309 = cos(t316);
t319 = sin(qJ(3));
t321 = cos(qJ(3));
t381 = Icges(4,5) * t321 + Icges(5,5) * t309 - Icges(4,6) * t319 - Icges(5,6) * t308;
t317 = qJ(1) + qJ(2);
t311 = sin(t317);
t312 = cos(t317);
t380 = -t381 * t311 + t382 * t312;
t379 = t382 * t311 + t381 * t312;
t378 = Icges(4,5) * t319 + Icges(5,5) * t308 + Icges(4,6) * t321 + Icges(5,6) * t309;
t362 = Icges(5,4) * t308;
t289 = Icges(5,2) * t309 + t362;
t361 = Icges(5,4) * t309;
t290 = Icges(5,1) * t308 + t361;
t364 = Icges(4,4) * t319;
t298 = Icges(4,2) * t321 + t364;
t363 = Icges(4,4) * t321;
t299 = Icges(4,1) * t319 + t363;
t377 = t289 * t308 - t290 * t309 + t298 * t319 - t299 * t321;
t340 = -Icges(5,2) * t308 + t361;
t268 = Icges(5,6) * t311 + t340 * t312;
t343 = Icges(5,1) * t309 - t362;
t270 = Icges(5,5) * t311 + t343 * t312;
t341 = -Icges(4,2) * t319 + t363;
t276 = Icges(4,6) * t311 + t341 * t312;
t344 = Icges(4,1) * t321 - t364;
t278 = Icges(4,5) * t311 + t344 * t312;
t376 = t268 * t308 - t270 * t309 + t276 * t319 - t278 * t321;
t267 = Icges(5,6) * t312 - t340 * t311;
t269 = Icges(5,5) * t312 - t343 * t311;
t275 = Icges(4,6) * t312 - t341 * t311;
t277 = Icges(4,5) * t312 - t344 * t311;
t375 = -t267 * t308 + t269 * t309 - t275 * t319 + t277 * t321;
t310 = qJ(5) + t316;
t305 = sin(t310);
t306 = cos(t310);
t359 = Icges(6,4) * t306;
t339 = -Icges(6,2) * t305 + t359;
t257 = Icges(6,6) * t312 - t339 * t311;
t258 = Icges(6,6) * t311 + t339 * t312;
t360 = Icges(6,4) * t305;
t342 = Icges(6,1) * t306 - t360;
t259 = Icges(6,5) * t312 - t342 * t311;
t260 = Icges(6,5) * t311 + t342 * t312;
t285 = Icges(6,2) * t306 + t360;
t286 = Icges(6,1) * t305 + t359;
t353 = qJD(3) + qJD(5);
t292 = t353 * t311;
t293 = t353 * t312;
t315 = qJD(1) + qJD(2);
t374 = (t285 * t305 - t286 * t306) * t315 + (t257 * t305 - t259 * t306) * t293 + (t258 * t305 - t260 * t306) * t292;
t369 = pkin(3) * t319;
t368 = pkin(4) * t308;
t367 = t321 * pkin(3);
t365 = pkin(1) * qJD(1);
t262 = qJ(4) * t311 + t367 * t312;
t294 = pkin(2) * t312 + pkin(7) * t311;
t358 = -t262 - t294;
t357 = pkin(4) * t309;
t355 = qJD(3) * t311;
t354 = qJD(3) * t312;
t320 = sin(qJ(1));
t352 = t320 * t365;
t322 = cos(qJ(1));
t351 = t322 * t365;
t348 = t315 * (-pkin(2) * t311 + pkin(7) * t312) - t352;
t347 = rSges(4,1) * t321 - rSges(4,2) * t319;
t346 = rSges(5,1) * t309 - rSges(5,2) * t308;
t345 = rSges(6,1) * t306 - rSges(6,2) * t305;
t336 = Icges(6,5) * t306 - Icges(6,6) * t305;
t326 = qJD(4) * t312 + t355 * t369 - t351;
t261 = qJ(4) * t312 - t367 * t311;
t325 = qJD(4) * t311 + t315 * t261 + t348;
t324 = (Icges(6,3) * t312 - t336 * t311) * t293 + (Icges(6,3) * t311 + t336 * t312) * t292 + (Icges(6,5) * t305 + Icges(6,6) * t306) * t315;
t302 = rSges(2,1) * t322 - rSges(2,2) * t320;
t301 = -rSges(2,1) * t320 - rSges(2,2) * t322;
t300 = rSges(4,1) * t319 + rSges(4,2) * t321;
t291 = rSges(5,1) * t308 + rSges(5,2) * t309;
t287 = rSges(6,1) * t305 + rSges(6,2) * t306;
t282 = -t351 - t315 * (rSges(3,1) * t312 - rSges(3,2) * t311);
t281 = -t352 + t315 * (-rSges(3,1) * t311 - rSges(3,2) * t312);
t280 = rSges(4,3) * t311 + t347 * t312;
t279 = rSges(4,3) * t312 - t347 * t311;
t272 = rSges(5,3) * t311 + t346 * t312;
t271 = rSges(5,3) * t312 - t346 * t311;
t264 = rSges(6,3) * t311 + t345 * t312;
t263 = rSges(6,3) * t312 - t345 * t311;
t253 = t262 * t354;
t252 = pkin(8) * t311 + t357 * t312;
t251 = pkin(8) * t312 - t357 * t311;
t250 = (-t279 * t311 + t280 * t312) * qJD(3);
t249 = -t351 + t300 * t355 + (-t280 - t294) * t315;
t248 = t279 * t315 - t300 * t354 + t348;
t247 = t291 * t355 + (-t272 + t358) * t315 + t326;
t246 = t271 * t315 + (-t291 - t369) * t354 + t325;
t245 = t253 + (t272 * t312 + (-t261 - t271) * t311) * qJD(3);
t244 = t355 * t368 + t287 * t292 + (-t252 - t264 + t358) * t315 + t326;
t243 = -t287 * t293 + (t251 + t263) * t315 + (-t368 - t369) * t354 + t325;
t242 = -t263 * t292 + t264 * t293 + t253 + (t252 * t312 + (-t251 - t261) * t311) * qJD(3);
t1 = m(3) * (t281 ^ 2 + t282 ^ 2) / 0.2e1 + t315 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t248 ^ 2 + t249 ^ 2 + t250 ^ 2) / 0.2e1 + m(5) * (t245 ^ 2 + t246 ^ 2 + t247 ^ 2) / 0.2e1 + m(6) * (t242 ^ 2 + t243 ^ 2 + t244 ^ 2) / 0.2e1 + t293 * (t374 * t311 + t324 * t312) / 0.2e1 + t292 * (t324 * t311 - t374 * t312) / 0.2e1 + (m(2) * (t301 ^ 2 + t302 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t378 * t311 - t377 * t312) * t315 + (t379 * t311 ^ 2 + (t375 * t312 + (-t376 + t380) * t311) * t312) * qJD(3)) * t355 / 0.2e1 + ((t377 * t311 + t378 * t312) * t315 + (t380 * t312 ^ 2 + (t376 * t311 + (-t375 + t379) * t312) * t311) * qJD(3)) * t354 / 0.2e1 + ((t257 * t306 + t259 * t305) * t293 + (t258 * t306 + t260 * t305) * t292 + ((t267 * t309 + t269 * t308 + t275 * t321 + t277 * t319) * t312 + (t268 * t309 + t270 * t308 + t276 * t321 + t278 * t319) * t311) * qJD(3) + (t306 * t285 + t305 * t286 + t309 * t289 + t308 * t290 + t321 * t298 + t319 * t299) * t315) * t315 / 0.2e1;
T = t1;
