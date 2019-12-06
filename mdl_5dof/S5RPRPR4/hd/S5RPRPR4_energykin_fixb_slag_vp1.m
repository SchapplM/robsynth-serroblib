% Calculate kinetic energy for
% S5RPRPR4
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
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:53:07
% EndTime: 2019-12-05 17:53:08
% DurationCPUTime: 1.11s
% Computational Cost: add. (985->169), mult. (765->255), div. (0->0), fcn. (634->10), ass. (0->101)
t382 = Icges(4,3) + Icges(5,3);
t316 = qJ(3) + pkin(9);
t309 = sin(t316);
t311 = cos(t316);
t319 = sin(qJ(3));
t321 = cos(qJ(3));
t381 = Icges(4,5) * t321 + Icges(5,5) * t311 - Icges(4,6) * t319 - Icges(5,6) * t309;
t317 = qJ(1) + pkin(8);
t310 = sin(t317);
t312 = cos(t317);
t380 = -t381 * t310 + t382 * t312;
t379 = t382 * t310 + t381 * t312;
t378 = Icges(4,5) * t319 + Icges(5,5) * t309 + Icges(4,6) * t321 + Icges(5,6) * t311;
t361 = Icges(5,4) * t309;
t292 = Icges(5,2) * t311 + t361;
t360 = Icges(5,4) * t311;
t293 = Icges(5,1) * t309 + t360;
t363 = Icges(4,4) * t319;
t299 = Icges(4,2) * t321 + t363;
t362 = Icges(4,4) * t321;
t300 = Icges(4,1) * t319 + t362;
t377 = t292 * t309 - t293 * t311 + t299 * t319 - t300 * t321;
t338 = -Icges(5,2) * t309 + t360;
t269 = Icges(5,6) * t310 + t338 * t312;
t341 = Icges(5,1) * t311 - t361;
t271 = Icges(5,5) * t310 + t341 * t312;
t339 = -Icges(4,2) * t319 + t362;
t277 = Icges(4,6) * t310 + t339 * t312;
t342 = Icges(4,1) * t321 - t363;
t279 = Icges(4,5) * t310 + t342 * t312;
t376 = t269 * t309 - t271 * t311 + t277 * t319 - t279 * t321;
t268 = Icges(5,6) * t312 - t338 * t310;
t270 = Icges(5,5) * t312 - t341 * t310;
t276 = Icges(4,6) * t312 - t339 * t310;
t278 = Icges(4,5) * t312 - t342 * t310;
t375 = -t268 * t309 + t270 * t311 - t276 * t319 + t278 * t321;
t313 = qJ(5) + t316;
t306 = sin(t313);
t307 = cos(t313);
t358 = Icges(6,4) * t307;
t337 = -Icges(6,2) * t306 + t358;
t260 = Icges(6,6) * t312 - t337 * t310;
t261 = Icges(6,6) * t310 + t337 * t312;
t359 = Icges(6,4) * t306;
t340 = Icges(6,1) * t307 - t359;
t262 = Icges(6,5) * t312 - t340 * t310;
t263 = Icges(6,5) * t310 + t340 * t312;
t285 = Icges(6,2) * t307 + t359;
t286 = Icges(6,1) * t306 + t358;
t351 = qJD(3) + qJD(5);
t289 = t351 * t310;
t290 = t351 * t312;
t374 = (t285 * t306 - t286 * t307) * qJD(1) + (t260 * t306 - t262 * t307) * t290 + (t261 * t306 - t263 * t307) * t289;
t320 = sin(qJ(1));
t370 = pkin(1) * t320;
t322 = cos(qJ(1));
t369 = pkin(1) * t322;
t368 = pkin(3) * t319;
t367 = pkin(4) * t309;
t365 = t321 * pkin(3);
t357 = pkin(4) * t311;
t354 = qJD(3) * t310;
t356 = qJD(4) * t312 + t354 * t368;
t353 = qJD(3) * t312;
t257 = qJ(4) * t310 + t365 * t312;
t352 = t257 * t353 + qJD(2);
t256 = qJ(4) * t312 - t365 * t310;
t288 = qJD(1) * (-pkin(2) * t310 + pkin(6) * t312);
t350 = qJD(1) * t256 + qJD(4) * t310 + t288;
t347 = -pkin(2) * t312 - pkin(6) * t310 - t369;
t346 = -t257 + t347;
t345 = rSges(4,1) * t321 - rSges(4,2) * t319;
t344 = rSges(5,1) * t311 - rSges(5,2) * t309;
t343 = rSges(6,1) * t307 - rSges(6,2) * t306;
t334 = Icges(6,5) * t307 - Icges(6,6) * t306;
t324 = (Icges(6,5) * t306 + Icges(6,6) * t307) * qJD(1) + (Icges(6,3) * t312 - t334 * t310) * t290 + (Icges(6,3) * t310 + t334 * t312) * t289;
t303 = rSges(2,1) * t322 - rSges(2,2) * t320;
t302 = -rSges(2,1) * t320 - rSges(2,2) * t322;
t301 = rSges(4,1) * t319 + rSges(4,2) * t321;
t294 = rSges(5,1) * t309 + rSges(5,2) * t311;
t287 = rSges(6,1) * t306 + rSges(6,2) * t307;
t283 = (-rSges(3,1) * t312 + rSges(3,2) * t310 - t369) * qJD(1);
t282 = (-rSges(3,1) * t310 - rSges(3,2) * t312 - t370) * qJD(1);
t281 = rSges(4,3) * t310 + t345 * t312;
t280 = rSges(4,3) * t312 - t345 * t310;
t273 = rSges(5,3) * t310 + t344 * t312;
t272 = rSges(5,3) * t312 - t344 * t310;
t265 = rSges(6,3) * t310 + t343 * t312;
t264 = rSges(6,3) * t312 - t343 * t310;
t253 = pkin(7) * t310 + t357 * t312;
t252 = pkin(7) * t312 - t357 * t310;
t251 = t301 * t354 + (-t281 + t347) * qJD(1);
t250 = -t301 * t353 + t288 + (t280 - t370) * qJD(1);
t249 = qJD(2) + (-t280 * t310 + t281 * t312) * qJD(3);
t248 = t294 * t354 + (-t273 + t346) * qJD(1) + t356;
t247 = (-t294 - t368) * t353 + (t272 - t370) * qJD(1) + t350;
t246 = (t273 * t312 + (-t256 - t272) * t310) * qJD(3) + t352;
t245 = t354 * t367 + t287 * t289 + (-t253 - t265 + t346) * qJD(1) + t356;
t244 = -t287 * t290 + (-t367 - t368) * t353 + (t252 + t264 - t370) * qJD(1) + t350;
t243 = -t264 * t289 + t265 * t290 + (t253 * t312 + (-t252 - t256) * t310) * qJD(3) + t352;
t1 = m(3) * (qJD(2) ^ 2 + t282 ^ 2 + t283 ^ 2) / 0.2e1 + m(4) * (t249 ^ 2 + t250 ^ 2 + t251 ^ 2) / 0.2e1 + m(5) * (t246 ^ 2 + t247 ^ 2 + t248 ^ 2) / 0.2e1 + m(6) * (t243 ^ 2 + t244 ^ 2 + t245 ^ 2) / 0.2e1 + t290 * (t374 * t310 + t324 * t312) / 0.2e1 + t289 * (t324 * t310 - t374 * t312) / 0.2e1 + ((t379 * t310 ^ 2 + (t375 * t312 + (-t376 + t380) * t310) * t312) * qJD(3) + (t378 * t310 - t377 * t312) * qJD(1)) * t354 / 0.2e1 + ((t380 * t312 ^ 2 + (t376 * t310 + (-t375 + t379) * t312) * t310) * qJD(3) + (t377 * t310 + t378 * t312) * qJD(1)) * t353 / 0.2e1 + (m(2) * (t302 ^ 2 + t303 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1 + ((t260 * t307 + t262 * t306) * t290 + (t261 * t307 + t263 * t306) * t289 + ((t268 * t311 + t270 * t309 + t276 * t321 + t278 * t319) * t312 + (t269 * t311 + t271 * t309 + t277 * t321 + t279 * t319) * t310) * qJD(3) + (t307 * t285 + t306 * t286 + t311 * t292 + t309 * t293 + t321 * t299 + t319 * t300) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
