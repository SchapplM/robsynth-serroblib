% Calculate kinetic energy for
% S5RPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR13_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR13_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR13_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR13_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR13_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:14:31
% EndTime: 2019-12-31 19:14:32
% DurationCPUTime: 1.28s
% Computational Cost: add. (693->219), mult. (1266->360), div. (0->0), fcn. (1248->8), ass. (0->110)
t334 = sin(qJ(3));
t337 = cos(qJ(3));
t336 = cos(qJ(4));
t366 = pkin(4) * t336;
t368 = -pkin(8) * t337 + t334 * t366;
t364 = Icges(4,4) * t334;
t363 = Icges(4,4) * t337;
t333 = sin(qJ(4));
t335 = sin(qJ(1));
t362 = t333 * t335;
t338 = cos(qJ(1));
t361 = t333 * t338;
t360 = t334 * t335;
t359 = t334 * t338;
t358 = t335 * t336;
t357 = t335 * t337;
t356 = t336 * t338;
t355 = t337 * t338;
t313 = qJD(1) * (pkin(1) * t338 + qJ(2) * t335);
t354 = qJD(1) * t338 * pkin(6) + t313;
t327 = qJD(3) * t335;
t353 = qJD(4) * t337;
t310 = t338 * t353 + t327;
t328 = qJD(3) * t338;
t323 = qJD(4) * t334 + qJD(1);
t317 = pkin(1) * t335 - qJ(2) * t338;
t352 = -pkin(6) * t335 - t317;
t351 = pkin(3) * t334 - pkin(7) * t337;
t307 = t351 * t335;
t308 = t351 * t338;
t350 = -t307 * t327 - t308 * t328;
t349 = rSges(4,1) * t334 + rSges(4,2) * t337;
t348 = Icges(4,1) * t334 + t363;
t347 = Icges(4,2) * t337 + t364;
t346 = Icges(4,5) * t334 + Icges(4,6) * t337;
t289 = Icges(4,6) * t338 + t335 * t347;
t292 = Icges(4,5) * t338 + t335 * t348;
t345 = -t289 * t337 - t292 * t334;
t290 = Icges(4,6) * t335 - t338 * t347;
t293 = Icges(4,5) * t335 - t338 * t348;
t344 = t290 * t337 + t293 * t334;
t315 = -Icges(4,2) * t334 + t363;
t316 = Icges(4,1) * t337 - t364;
t343 = t315 * t337 + t316 * t334;
t321 = pkin(3) * t337 + pkin(7) * t334;
t329 = qJD(2) * t335;
t342 = t321 * t327 + t329 + (t308 + t352) * qJD(1);
t341 = qJD(1) * t307 + (-qJD(3) * t321 - qJD(2)) * t338 + t354;
t332 = qJ(4) + qJ(5);
t331 = cos(t332);
t330 = sin(t332);
t320 = rSges(2,1) * t338 - rSges(2,2) * t335;
t319 = rSges(4,1) * t337 - rSges(4,2) * t334;
t318 = rSges(2,1) * t335 + rSges(2,2) * t338;
t314 = Icges(4,5) * t337 - Icges(4,6) * t334;
t312 = qJD(5) * t334 + t323;
t311 = -t335 * t353 + t328;
t306 = -t334 * t356 + t362;
t305 = t333 * t359 + t358;
t304 = t334 * t358 + t361;
t303 = -t333 * t360 + t356;
t300 = t330 * t335 - t331 * t359;
t299 = t330 * t359 + t331 * t335;
t298 = t330 * t338 + t331 * t360;
t297 = -t330 * t360 + t331 * t338;
t296 = rSges(4,3) * t335 - t338 * t349;
t295 = rSges(5,3) * t334 + (rSges(5,1) * t336 - rSges(5,2) * t333) * t337;
t294 = rSges(4,3) * t338 + t335 * t349;
t291 = Icges(5,5) * t334 + (Icges(5,1) * t336 - Icges(5,4) * t333) * t337;
t288 = Icges(5,6) * t334 + (Icges(5,4) * t336 - Icges(5,2) * t333) * t337;
t287 = Icges(4,3) * t335 - t338 * t346;
t286 = Icges(4,3) * t338 + t335 * t346;
t285 = Icges(5,3) * t334 + (Icges(5,5) * t336 - Icges(5,6) * t333) * t337;
t284 = t328 + (-qJD(4) - qJD(5)) * t357;
t283 = qJD(5) * t355 + t310;
t282 = rSges(6,3) * t334 + (rSges(6,1) * t331 - rSges(6,2) * t330) * t337;
t281 = Icges(6,5) * t334 + (Icges(6,1) * t331 - Icges(6,4) * t330) * t337;
t280 = Icges(6,6) * t334 + (Icges(6,4) * t331 - Icges(6,2) * t330) * t337;
t279 = Icges(6,3) * t334 + (Icges(6,5) * t331 - Icges(6,6) * t330) * t337;
t278 = pkin(8) * t334 + t337 * t366;
t277 = t313 - qJD(2) * t338 + qJD(1) * (-rSges(3,2) * t338 + rSges(3,3) * t335);
t276 = t329 + (rSges(3,2) * t335 + rSges(3,3) * t338 - t317) * qJD(1);
t275 = pkin(4) * t362 - t368 * t338;
t274 = pkin(4) * t361 + t368 * t335;
t273 = rSges(5,1) * t306 + rSges(5,2) * t305 + rSges(5,3) * t355;
t272 = rSges(5,1) * t304 + rSges(5,2) * t303 - rSges(5,3) * t357;
t271 = Icges(5,1) * t306 + Icges(5,4) * t305 + Icges(5,5) * t355;
t270 = Icges(5,1) * t304 + Icges(5,4) * t303 - Icges(5,5) * t357;
t269 = Icges(5,4) * t306 + Icges(5,2) * t305 + Icges(5,6) * t355;
t268 = Icges(5,4) * t304 + Icges(5,2) * t303 - Icges(5,6) * t357;
t267 = Icges(5,5) * t306 + Icges(5,6) * t305 + Icges(5,3) * t355;
t266 = Icges(5,5) * t304 + Icges(5,6) * t303 - Icges(5,3) * t357;
t265 = (-t294 * t335 + t296 * t338) * qJD(3);
t264 = rSges(6,1) * t300 + rSges(6,2) * t299 + rSges(6,3) * t355;
t263 = rSges(6,1) * t298 + rSges(6,2) * t297 - rSges(6,3) * t357;
t262 = Icges(6,1) * t300 + Icges(6,4) * t299 + Icges(6,5) * t355;
t261 = Icges(6,1) * t298 + Icges(6,4) * t297 - Icges(6,5) * t357;
t260 = Icges(6,4) * t300 + Icges(6,2) * t299 + Icges(6,6) * t355;
t259 = Icges(6,4) * t298 + Icges(6,2) * t297 - Icges(6,6) * t357;
t258 = Icges(6,5) * t300 + Icges(6,6) * t299 + Icges(6,3) * t355;
t257 = Icges(6,5) * t298 + Icges(6,6) * t297 - Icges(6,3) * t357;
t256 = qJD(1) * t294 + (-qJD(3) * t319 - qJD(2)) * t338 + t354;
t255 = t319 * t327 + t329 + (-t296 + t352) * qJD(1);
t254 = t272 * t323 - t295 * t311 + t341;
t253 = -t273 * t323 + t295 * t310 + t342;
t252 = -t272 * t310 + t273 * t311 + t350;
t251 = t263 * t312 + t274 * t323 - t278 * t311 - t282 * t284 + t341;
t250 = -t264 * t312 - t275 * t323 + t278 * t310 + t282 * t283 + t342;
t249 = -t263 * t283 + t264 * t284 - t274 * t310 + t275 * t311 + t350;
t1 = m(3) * (t276 ^ 2 + t277 ^ 2) / 0.2e1 + m(4) * (t255 ^ 2 + t256 ^ 2 + t265 ^ 2) / 0.2e1 + ((t338 * t314 + t335 * t343) * qJD(1) + (t338 ^ 2 * t286 + (t344 * t335 + (t287 - t345) * t338) * t335) * qJD(3)) * t328 / 0.2e1 + ((t335 * t314 - t338 * t343) * qJD(1) + (t335 ^ 2 * t287 + (t345 * t338 + (t286 - t344) * t335) * t338) * qJD(3)) * t327 / 0.2e1 + qJD(1) * ((-t334 * t315 + t337 * t316) * qJD(1) + ((-t289 * t334 + t292 * t337) * t338 + (-t290 * t334 + t293 * t337) * t335) * qJD(3)) / 0.2e1 + m(5) * (t252 ^ 2 + t253 ^ 2 + t254 ^ 2) / 0.2e1 + t311 * ((-t266 * t357 + t303 * t268 + t304 * t270) * t311 + (-t267 * t357 + t269 * t303 + t271 * t304) * t310 + (-t285 * t357 + t288 * t303 + t291 * t304) * t323) / 0.2e1 + t310 * ((t266 * t355 + t268 * t305 + t270 * t306) * t311 + (t267 * t355 + t305 * t269 + t306 * t271) * t310 + (t285 * t355 + t288 * t305 + t291 * t306) * t323) / 0.2e1 + t323 * ((t266 * t311 + t267 * t310 + t285 * t323) * t334 + ((-t268 * t333 + t270 * t336) * t311 + (-t269 * t333 + t271 * t336) * t310 + (-t288 * t333 + t291 * t336) * t323) * t337) / 0.2e1 + m(6) * (t249 ^ 2 + t250 ^ 2 + t251 ^ 2) / 0.2e1 + t284 * ((-t257 * t357 + t297 * t259 + t298 * t261) * t284 + (-t258 * t357 + t260 * t297 + t262 * t298) * t283 + (-t279 * t357 + t280 * t297 + t281 * t298) * t312) / 0.2e1 + t283 * ((t257 * t355 + t259 * t299 + t261 * t300) * t284 + (t258 * t355 + t299 * t260 + t300 * t262) * t283 + (t279 * t355 + t280 * t299 + t281 * t300) * t312) / 0.2e1 + t312 * ((t257 * t284 + t258 * t283 + t279 * t312) * t334 + ((-t259 * t330 + t261 * t331) * t284 + (-t260 * t330 + t262 * t331) * t283 + (-t280 * t330 + t281 * t331) * t312) * t337) / 0.2e1 + (m(2) * (t318 ^ 2 + t320 ^ 2) + Icges(2,3) + Icges(3,1)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
