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
% m [6x1]
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
% Datum: 2022-01-23 09:23
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:22:22
% EndTime: 2022-01-23 09:22:23
% DurationCPUTime: 1.05s
% Computational Cost: add. (985->165), mult. (765->257), div. (0->0), fcn. (634->10), ass. (0->101)
t369 = Icges(4,3) + Icges(5,3);
t306 = qJ(3) + pkin(9);
t299 = sin(t306);
t301 = cos(t306);
t309 = sin(qJ(3));
t311 = cos(qJ(3));
t368 = Icges(4,5) * t311 + Icges(5,5) * t301 - Icges(4,6) * t309 - Icges(5,6) * t299;
t307 = qJ(1) + pkin(8);
t300 = sin(t307);
t302 = cos(t307);
t367 = t368 * t300 - t369 * t302;
t366 = t369 * t300 + t368 * t302;
t365 = Icges(4,5) * t309 + Icges(5,5) * t299 + Icges(4,6) * t311 + Icges(5,6) * t301;
t351 = Icges(5,4) * t299;
t283 = Icges(5,2) * t301 + t351;
t350 = Icges(5,4) * t301;
t284 = Icges(5,1) * t299 + t350;
t353 = Icges(4,4) * t309;
t289 = Icges(4,2) * t311 + t353;
t352 = Icges(4,4) * t311;
t290 = Icges(4,1) * t309 + t352;
t364 = -t283 * t299 + t284 * t301 - t289 * t309 + t290 * t311;
t328 = -Icges(5,2) * t299 + t350;
t260 = Icges(5,6) * t300 + t328 * t302;
t331 = Icges(5,1) * t301 - t351;
t262 = Icges(5,5) * t300 + t331 * t302;
t329 = -Icges(4,2) * t309 + t352;
t268 = Icges(4,6) * t300 + t329 * t302;
t332 = Icges(4,1) * t311 - t353;
t270 = Icges(4,5) * t300 + t332 * t302;
t363 = -t260 * t299 + t262 * t301 - t268 * t309 + t270 * t311;
t259 = -Icges(5,6) * t302 + t328 * t300;
t261 = -Icges(5,5) * t302 + t331 * t300;
t267 = -Icges(4,6) * t302 + t329 * t300;
t269 = -Icges(4,5) * t302 + t332 * t300;
t362 = t259 * t299 - t261 * t301 + t267 * t309 - t269 * t311;
t358 = pkin(3) * t309;
t310 = sin(qJ(1));
t357 = t310 * pkin(1);
t355 = t311 * pkin(3);
t303 = qJ(5) + t306;
t295 = sin(t303);
t349 = Icges(6,4) * t295;
t296 = cos(t303);
t348 = Icges(6,4) * t296;
t312 = cos(qJ(1));
t298 = qJD(1) * t312 * pkin(1);
t347 = qJD(1) * (t302 * pkin(2) + t300 * pkin(6)) + t298;
t346 = pkin(4) * t301;
t344 = qJD(3) * t300;
t343 = qJD(3) * t302;
t342 = qJD(3) + qJD(5);
t247 = -qJ(4) * t302 + t355 * t300;
t248 = qJ(4) * t300 + t355 * t302;
t341 = t247 * t344 + t248 * t343 + qJD(2);
t338 = -t300 * pkin(2) + t302 * pkin(6) - t357;
t337 = -t247 + t338;
t336 = rSges(4,1) * t311 - rSges(4,2) * t309;
t335 = rSges(5,1) * t301 - rSges(5,2) * t299;
t334 = rSges(6,1) * t296 - rSges(6,2) * t295;
t333 = qJD(3) * (-t299 * rSges(5,1) - t301 * rSges(5,2) - t358);
t330 = Icges(6,1) * t296 - t349;
t327 = -Icges(6,2) * t295 + t348;
t324 = Icges(6,5) * t296 - Icges(6,6) * t295;
t317 = qJD(1) * t248 - qJD(4) * t302 + t347;
t316 = qJD(3) * (-pkin(4) * t299 - t358);
t280 = t342 * t300;
t281 = t342 * t302;
t315 = qJD(1) * (Icges(6,5) * t295 + Icges(6,6) * t296) - (-Icges(6,3) * t302 + t324 * t300) * t281 + (Icges(6,3) * t300 + t324 * t302) * t280;
t251 = -Icges(6,6) * t302 + t327 * t300;
t252 = Icges(6,6) * t300 + t327 * t302;
t253 = -Icges(6,5) * t302 + t330 * t300;
t254 = Icges(6,5) * t300 + t330 * t302;
t276 = Icges(6,2) * t296 + t349;
t277 = Icges(6,1) * t295 + t348;
t314 = (-t252 * t295 + t254 * t296) * t280 - (-t251 * t295 + t253 * t296) * t281 + (-t276 * t295 + t277 * t296) * qJD(1);
t294 = qJD(4) * t300;
t293 = t312 * rSges(2,1) - t310 * rSges(2,2);
t292 = t310 * rSges(2,1) + t312 * rSges(2,2);
t291 = t309 * rSges(4,1) + t311 * rSges(4,2);
t278 = t295 * rSges(6,1) + t296 * rSges(6,2);
t274 = t298 + qJD(1) * (t302 * rSges(3,1) - t300 * rSges(3,2));
t273 = (-t300 * rSges(3,1) - t302 * rSges(3,2) - t357) * qJD(1);
t272 = t300 * rSges(4,3) + t336 * t302;
t271 = -t302 * rSges(4,3) + t336 * t300;
t264 = t300 * rSges(5,3) + t335 * t302;
t263 = -t302 * rSges(5,3) + t335 * t300;
t256 = t300 * rSges(6,3) + t334 * t302;
t255 = -t302 * rSges(6,3) + t334 * t300;
t243 = pkin(7) * t300 + t346 * t302;
t242 = -pkin(7) * t302 + t346 * t300;
t241 = qJD(1) * t272 - t291 * t344 + t347;
t240 = -t291 * t343 + (-t271 + t338) * qJD(1);
t239 = qJD(2) + (t271 * t300 + t272 * t302) * qJD(3);
t238 = qJD(1) * t264 + t300 * t333 + t317;
t237 = t294 + t302 * t333 + (-t263 + t337) * qJD(1);
t236 = (t263 * t300 + t264 * t302) * qJD(3) + t341;
t235 = -t280 * t278 + t300 * t316 + (t243 + t256) * qJD(1) + t317;
t234 = -t281 * t278 + t294 + t302 * t316 + (-t242 - t255 + t337) * qJD(1);
t233 = t280 * t255 + t281 * t256 + (t242 * t300 + t243 * t302) * qJD(3) + t341;
t1 = m(3) * (qJD(2) ^ 2 + t273 ^ 2 + t274 ^ 2) / 0.2e1 + m(4) * (t239 ^ 2 + t240 ^ 2 + t241 ^ 2) / 0.2e1 + m(5) * (t236 ^ 2 + t237 ^ 2 + t238 ^ 2) / 0.2e1 + m(6) * (t233 ^ 2 + t234 ^ 2 + t235 ^ 2) / 0.2e1 + t280 * (t315 * t300 + t314 * t302) / 0.2e1 - t281 * (t314 * t300 - t315 * t302) / 0.2e1 + ((t366 * t300 ^ 2 + (t362 * t302 + (t363 - t367) * t300) * t302) * qJD(3) + (t365 * t300 + t364 * t302) * qJD(1)) * t344 / 0.2e1 - ((t367 * t302 ^ 2 + (t363 * t300 + (t362 - t366) * t302) * t300) * qJD(3) + (t364 * t300 - t365 * t302) * qJD(1)) * t343 / 0.2e1 + (m(2) * (t292 ^ 2 + t293 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1 + ((t252 * t296 + t254 * t295) * t280 - (t251 * t296 + t253 * t295) * t281 + ((-t259 * t301 - t261 * t299 - t267 * t311 - t269 * t309) * t302 + (t260 * t301 + t262 * t299 + t268 * t311 + t270 * t309) * t300) * qJD(3) + (t296 * t276 + t295 * t277 + t301 * t283 + t299 * t284 + t311 * t289 + t309 * t290) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
