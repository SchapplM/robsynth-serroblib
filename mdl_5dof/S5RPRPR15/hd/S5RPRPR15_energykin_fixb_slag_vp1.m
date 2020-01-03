% Calculate kinetic energy for
% S5RPRPR15
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
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR15_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR15_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR15_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR15_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR15_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:36:36
% EndTime: 2019-12-31 18:36:38
% DurationCPUTime: 1.36s
% Computational Cost: add. (657->223), mult. (1206->348), div. (0->0), fcn. (1188->8), ass. (0->107)
t329 = sin(qJ(3));
t331 = cos(qJ(3));
t327 = cos(pkin(8));
t363 = pkin(4) * t327;
t366 = -pkin(7) * t331 + t363 * t329;
t362 = Icges(4,4) * t329;
t361 = Icges(4,4) * t331;
t326 = sin(pkin(8));
t330 = sin(qJ(1));
t360 = t326 * t330;
t332 = cos(qJ(1));
t359 = t326 * t332;
t358 = t329 * t330;
t357 = t329 * t332;
t356 = t330 * t331;
t355 = t331 * t332;
t340 = pkin(3) * t329 - qJ(4) * t331;
t304 = t340 * t332;
t349 = qJD(3) * t332;
t353 = qJD(4) * t329 - t304 * t349;
t314 = pkin(3) * t331 + qJ(4) * t329;
t324 = qJD(2) * t330;
t350 = qJD(3) * t330;
t352 = t314 * t350 + t324;
t308 = qJD(1) * (pkin(1) * t332 + qJ(2) * t330);
t351 = qJD(1) * t332 * pkin(6) + t308;
t348 = qJD(4) * t331;
t347 = qJD(5) * t331;
t312 = pkin(1) * t330 - qJ(2) * t332;
t344 = -pkin(6) * t330 - t312;
t303 = t340 * t330;
t343 = qJD(1) * t303 + t332 * t348 + t351;
t342 = t304 + t344;
t341 = rSges(4,1) * t329 + rSges(4,2) * t331;
t339 = Icges(4,1) * t329 + t361;
t338 = Icges(4,2) * t331 + t362;
t337 = Icges(4,5) * t329 + Icges(4,6) * t331;
t291 = Icges(4,6) * t332 + t338 * t330;
t293 = Icges(4,5) * t332 + t339 * t330;
t336 = -t291 * t331 - t293 * t329;
t292 = Icges(4,6) * t330 - t338 * t332;
t294 = Icges(4,5) * t330 - t339 * t332;
t335 = t292 * t331 + t294 * t329;
t310 = -Icges(4,2) * t329 + t361;
t311 = Icges(4,1) * t331 - t362;
t334 = t310 * t331 + t311 * t329;
t325 = pkin(8) + qJ(5);
t322 = cos(t325);
t321 = sin(t325);
t318 = qJD(5) * t329 + qJD(1);
t316 = rSges(2,1) * t332 - rSges(2,2) * t330;
t315 = rSges(4,1) * t331 - rSges(4,2) * t329;
t313 = rSges(2,1) * t330 + rSges(2,2) * t332;
t309 = Icges(4,5) * t331 - Icges(4,6) * t329;
t307 = -t330 * t347 + t349;
t306 = t332 * t347 + t350;
t302 = -t327 * t357 + t360;
t301 = t326 * t357 + t327 * t330;
t300 = t327 * t358 + t359;
t299 = -t326 * t358 + t327 * t332;
t297 = rSges(4,3) * t330 - t341 * t332;
t296 = rSges(4,3) * t332 + t341 * t330;
t290 = Icges(4,3) * t330 - t337 * t332;
t289 = Icges(4,3) * t332 + t337 * t330;
t288 = t321 * t330 - t322 * t357;
t287 = t321 * t357 + t322 * t330;
t286 = t321 * t332 + t322 * t358;
t285 = -t321 * t358 + t322 * t332;
t284 = rSges(5,3) * t329 + (rSges(5,1) * t327 - rSges(5,2) * t326) * t331;
t283 = Icges(5,5) * t329 + (Icges(5,1) * t327 - Icges(5,4) * t326) * t331;
t282 = Icges(5,6) * t329 + (Icges(5,4) * t327 - Icges(5,2) * t326) * t331;
t281 = Icges(5,3) * t329 + (Icges(5,5) * t327 - Icges(5,6) * t326) * t331;
t280 = rSges(6,3) * t329 + (rSges(6,1) * t322 - rSges(6,2) * t321) * t331;
t279 = Icges(6,5) * t329 + (Icges(6,1) * t322 - Icges(6,4) * t321) * t331;
t278 = Icges(6,6) * t329 + (Icges(6,4) * t322 - Icges(6,2) * t321) * t331;
t277 = Icges(6,3) * t329 + (Icges(6,5) * t322 - Icges(6,6) * t321) * t331;
t276 = pkin(7) * t329 + t363 * t331;
t275 = t308 - qJD(2) * t332 + qJD(1) * (-rSges(3,2) * t332 + rSges(3,3) * t330);
t274 = t324 + (rSges(3,2) * t330 + rSges(3,3) * t332 - t312) * qJD(1);
t273 = pkin(4) * t360 - t332 * t366;
t272 = pkin(4) * t359 + t330 * t366;
t271 = rSges(5,1) * t302 + rSges(5,2) * t301 + rSges(5,3) * t355;
t270 = rSges(5,1) * t300 + rSges(5,2) * t299 - rSges(5,3) * t356;
t269 = Icges(5,1) * t302 + Icges(5,4) * t301 + Icges(5,5) * t355;
t268 = Icges(5,1) * t300 + Icges(5,4) * t299 - Icges(5,5) * t356;
t267 = Icges(5,4) * t302 + Icges(5,2) * t301 + Icges(5,6) * t355;
t266 = Icges(5,4) * t300 + Icges(5,2) * t299 - Icges(5,6) * t356;
t265 = Icges(5,5) * t302 + Icges(5,6) * t301 + Icges(5,3) * t355;
t264 = Icges(5,5) * t300 + Icges(5,6) * t299 - Icges(5,3) * t356;
t263 = (-t296 * t330 + t297 * t332) * qJD(3);
t262 = rSges(6,1) * t288 + rSges(6,2) * t287 + rSges(6,3) * t355;
t261 = rSges(6,1) * t286 + rSges(6,2) * t285 - rSges(6,3) * t356;
t260 = Icges(6,1) * t288 + Icges(6,4) * t287 + Icges(6,5) * t355;
t259 = Icges(6,1) * t286 + Icges(6,4) * t285 - Icges(6,5) * t356;
t258 = Icges(6,4) * t288 + Icges(6,2) * t287 + Icges(6,6) * t355;
t257 = Icges(6,4) * t286 + Icges(6,2) * t285 - Icges(6,6) * t356;
t256 = Icges(6,5) * t288 + Icges(6,6) * t287 + Icges(6,3) * t355;
t255 = Icges(6,5) * t286 + Icges(6,6) * t285 - Icges(6,3) * t356;
t254 = qJD(1) * t296 + (-qJD(3) * t315 - qJD(2)) * t332 + t351;
t253 = t315 * t350 + t324 + (-t297 + t344) * qJD(1);
t252 = qJD(1) * t270 + (-qJD(2) + (-t284 - t314) * qJD(3)) * t332 + t343;
t251 = (qJD(3) * t284 - t348) * t330 + (-t271 + t342) * qJD(1) + t352;
t250 = (t271 * t332 + (-t270 - t303) * t330) * qJD(3) + t353;
t249 = qJD(1) * t272 + t261 * t318 - t280 * t307 + (-qJD(2) + (-t276 - t314) * qJD(3)) * t332 + t343;
t248 = -t262 * t318 + t280 * t306 + (qJD(3) * t276 - t348) * t330 + (-t273 + t342) * qJD(1) + t352;
t247 = -t261 * t306 + t262 * t307 + (t273 * t332 + (-t272 - t303) * t330) * qJD(3) + t353;
t1 = m(3) * (t274 ^ 2 + t275 ^ 2) / 0.2e1 + m(4) * (t253 ^ 2 + t254 ^ 2 + t263 ^ 2) / 0.2e1 + m(5) * (t250 ^ 2 + t251 ^ 2 + t252 ^ 2) / 0.2e1 + m(6) * (t247 ^ 2 + t248 ^ 2 + t249 ^ 2) / 0.2e1 + t307 * ((-t255 * t356 + t257 * t285 + t259 * t286) * t307 + (-t256 * t356 + t258 * t285 + t260 * t286) * t306 + (-t277 * t356 + t278 * t285 + t279 * t286) * t318) / 0.2e1 + t306 * ((t255 * t355 + t257 * t287 + t259 * t288) * t307 + (t256 * t355 + t258 * t287 + t260 * t288) * t306 + (t277 * t355 + t278 * t287 + t279 * t288) * t318) / 0.2e1 + t318 * ((t255 * t307 + t256 * t306 + t277 * t318) * t329 + ((-t257 * t321 + t259 * t322) * t307 + (-t258 * t321 + t260 * t322) * t306 + (-t278 * t321 + t279 * t322) * t318) * t331) / 0.2e1 + (((-t291 * t329 + t293 * t331) * t332 + (-t292 * t329 + t294 * t331) * t330 + (t264 * t332 + t265 * t330) * t329 + ((-t266 * t326 + t268 * t327) * t332 + (-t267 * t326 + t269 * t327) * t330) * t331) * qJD(3) + ((-t282 * t326 + t283 * t327 + t311) * t331 + (-t310 + t281) * t329) * qJD(1)) * qJD(1) / 0.2e1 + (((t264 * t355 + t266 * t301 + t268 * t302 + t336 * t332) * t332 + ((t289 - t335) * t332 + t265 * t355 + t267 * t301 + t269 * t302 + t290 * t330) * t330) * qJD(3) + (t281 * t355 + t282 * t301 + t283 * t302 + t330 * t309 - t334 * t332) * qJD(1)) * t350 / 0.2e1 + (((-t264 * t356 + t266 * t299 + t268 * t300 + t289 * t332) * t332 + ((t290 - t336) * t332 - t265 * t356 + t267 * t299 + t269 * t300 + t335 * t330) * t330) * qJD(3) + (-t281 * t356 + t282 * t299 + t283 * t300 + t332 * t309 + t334 * t330) * qJD(1)) * t349 / 0.2e1 + (m(2) * (t313 ^ 2 + t316 ^ 2) + Icges(2,3) + Icges(3,1)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
