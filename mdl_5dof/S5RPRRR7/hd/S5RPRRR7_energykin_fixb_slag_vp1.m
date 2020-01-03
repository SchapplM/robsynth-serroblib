% Calculate kinetic energy for
% S5RPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR7_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR7_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR7_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR7_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:03:10
% EndTime: 2019-12-31 19:03:11
% DurationCPUTime: 1.26s
% Computational Cost: add. (1185->217), mult. (1257->359), div. (0->0), fcn. (1238->10), ass. (0->113)
t334 = sin(qJ(1));
t369 = pkin(1) * t334;
t335 = cos(qJ(4));
t368 = pkin(4) * t335;
t333 = sin(qJ(3));
t366 = Icges(4,4) * t333;
t336 = cos(qJ(3));
t365 = Icges(4,4) * t336;
t330 = qJ(1) + pkin(9);
t326 = sin(t330);
t332 = sin(qJ(4));
t364 = t326 * t332;
t363 = t326 * t333;
t327 = cos(t330);
t362 = t327 * t332;
t361 = t327 * t333;
t331 = qJ(4) + qJ(5);
t328 = sin(t331);
t360 = t328 * t336;
t329 = cos(t331);
t359 = t329 * t336;
t358 = t332 * t336;
t357 = t335 * t336;
t337 = cos(qJ(1));
t325 = qJD(1) * t337 * pkin(1);
t356 = qJD(1) * (pkin(2) * t327 + pkin(6) * t326) + t325;
t322 = qJD(3) * t326;
t354 = qJD(4) * t333;
t308 = t327 * t354 + t322;
t355 = qJD(3) * t327;
t353 = qJD(5) * t333;
t350 = pkin(3) * t336 + pkin(7) * t333;
t306 = t350 * t326;
t307 = t350 * t327;
t352 = t306 * t322 + t307 * t355 + qJD(2);
t351 = -pkin(2) * t326 + pkin(6) * t327 - t369;
t309 = t326 * t354 - t355;
t349 = rSges(4,1) * t336 - rSges(4,2) * t333;
t348 = Icges(4,1) * t336 - t366;
t347 = -Icges(4,2) * t333 + t365;
t346 = Icges(4,5) * t336 - Icges(4,6) * t333;
t278 = -Icges(4,6) * t327 + t326 * t347;
t280 = -Icges(4,5) * t327 + t326 * t348;
t345 = t278 * t333 - t280 * t336;
t279 = Icges(4,6) * t326 + t327 * t347;
t281 = Icges(4,5) * t326 + t327 * t348;
t344 = -t279 * t333 + t281 * t336;
t314 = Icges(4,2) * t336 + t366;
t315 = Icges(4,1) * t333 + t365;
t343 = -t314 * t333 + t315 * t336;
t321 = pkin(3) * t333 - pkin(7) * t336;
t342 = qJD(1) * t307 - t321 * t322 + t356;
t341 = pkin(8) * t333 + t336 * t368;
t340 = (-t306 + t351) * qJD(1) - t321 * t355;
t323 = -qJD(4) * t336 + qJD(1);
t320 = rSges(2,1) * t337 - rSges(2,2) * t334;
t319 = rSges(2,1) * t334 + rSges(2,2) * t337;
t318 = rSges(4,1) * t333 + rSges(4,2) * t336;
t313 = Icges(4,5) * t333 + Icges(4,6) * t336;
t312 = qJD(1) + (-qJD(4) - qJD(5)) * t336;
t305 = -rSges(5,3) * t336 + (rSges(5,1) * t335 - rSges(5,2) * t332) * t333;
t304 = -Icges(5,5) * t336 + (Icges(5,1) * t335 - Icges(5,4) * t332) * t333;
t303 = -Icges(5,6) * t336 + (Icges(5,4) * t335 - Icges(5,2) * t332) * t333;
t302 = -Icges(5,3) * t336 + (Icges(5,5) * t335 - Icges(5,6) * t332) * t333;
t301 = t327 * t357 + t364;
t300 = t326 * t335 - t327 * t358;
t299 = t326 * t357 - t362;
t298 = -t326 * t358 - t327 * t335;
t296 = t325 + qJD(1) * (rSges(3,1) * t327 - rSges(3,2) * t326);
t295 = (-rSges(3,1) * t326 - rSges(3,2) * t327 - t369) * qJD(1);
t294 = -rSges(6,3) * t336 + (rSges(6,1) * t329 - rSges(6,2) * t328) * t333;
t293 = -Icges(6,5) * t336 + (Icges(6,1) * t329 - Icges(6,4) * t328) * t333;
t292 = -Icges(6,6) * t336 + (Icges(6,4) * t329 - Icges(6,2) * t328) * t333;
t291 = -Icges(6,3) * t336 + (Icges(6,5) * t329 - Icges(6,6) * t328) * t333;
t290 = t326 * t328 + t327 * t359;
t289 = t326 * t329 - t327 * t360;
t288 = t326 * t359 - t327 * t328;
t287 = -t326 * t360 - t327 * t329;
t286 = -pkin(8) * t336 + t333 * t368;
t283 = rSges(4,3) * t326 + t327 * t349;
t282 = -rSges(4,3) * t327 + t326 * t349;
t277 = Icges(4,3) * t326 + t327 * t346;
t276 = -Icges(4,3) * t327 + t326 * t346;
t275 = t326 * t353 + t309;
t274 = t327 * t353 + t308;
t273 = pkin(4) * t364 + t327 * t341;
t272 = -pkin(4) * t362 + t326 * t341;
t271 = rSges(5,1) * t301 + rSges(5,2) * t300 + rSges(5,3) * t361;
t270 = rSges(5,1) * t299 + rSges(5,2) * t298 + rSges(5,3) * t363;
t269 = Icges(5,1) * t301 + Icges(5,4) * t300 + Icges(5,5) * t361;
t268 = Icges(5,1) * t299 + Icges(5,4) * t298 + Icges(5,5) * t363;
t267 = Icges(5,4) * t301 + Icges(5,2) * t300 + Icges(5,6) * t361;
t266 = Icges(5,4) * t299 + Icges(5,2) * t298 + Icges(5,6) * t363;
t265 = Icges(5,5) * t301 + Icges(5,6) * t300 + Icges(5,3) * t361;
t264 = Icges(5,5) * t299 + Icges(5,6) * t298 + Icges(5,3) * t363;
t263 = rSges(6,1) * t290 + rSges(6,2) * t289 + rSges(6,3) * t361;
t262 = rSges(6,1) * t288 + rSges(6,2) * t287 + rSges(6,3) * t363;
t261 = Icges(6,1) * t290 + Icges(6,4) * t289 + Icges(6,5) * t361;
t260 = Icges(6,1) * t288 + Icges(6,4) * t287 + Icges(6,5) * t363;
t259 = Icges(6,4) * t290 + Icges(6,2) * t289 + Icges(6,6) * t361;
t258 = Icges(6,4) * t288 + Icges(6,2) * t287 + Icges(6,6) * t363;
t257 = Icges(6,5) * t290 + Icges(6,6) * t289 + Icges(6,3) * t361;
t256 = Icges(6,5) * t288 + Icges(6,6) * t287 + Icges(6,3) * t363;
t255 = qJD(1) * t283 - t318 * t322 + t356;
t254 = -t318 * t355 + (-t282 + t351) * qJD(1);
t253 = qJD(2) + (t282 * t326 + t283 * t327) * qJD(3);
t252 = t271 * t323 - t305 * t308 + t342;
t251 = -t270 * t323 + t305 * t309 + t340;
t250 = t270 * t308 - t271 * t309 + t352;
t249 = t263 * t312 + t273 * t323 - t274 * t294 - t286 * t308 + t342;
t248 = -t262 * t312 - t272 * t323 + t275 * t294 + t286 * t309 + t340;
t247 = t262 * t274 - t263 * t275 + t272 * t308 - t273 * t309 + t352;
t1 = m(3) * (qJD(2) ^ 2 + t295 ^ 2 + t296 ^ 2) / 0.2e1 + m(4) * (t253 ^ 2 + t254 ^ 2 + t255 ^ 2) / 0.2e1 + ((t326 * t313 + t327 * t343) * qJD(1) + (t326 ^ 2 * t277 + (t345 * t327 + (-t276 + t344) * t326) * t327) * qJD(3)) * t322 / 0.2e1 - ((-t327 * t313 + t326 * t343) * qJD(1) + (t327 ^ 2 * t276 + (t344 * t326 + (-t277 + t345) * t327) * t326) * qJD(3)) * t355 / 0.2e1 + qJD(1) * ((t336 * t314 + t333 * t315) * qJD(1) + ((t279 * t336 + t281 * t333) * t326 - (t278 * t336 + t280 * t333) * t327) * qJD(3)) / 0.2e1 + m(5) * (t250 ^ 2 + t251 ^ 2 + t252 ^ 2) / 0.2e1 + t308 * ((t265 * t361 + t300 * t267 + t301 * t269) * t308 + (t264 * t361 + t300 * t266 + t301 * t268) * t309 + (t300 * t303 + t301 * t304 + t302 * t361) * t323) / 0.2e1 + t309 * ((t265 * t363 + t267 * t298 + t269 * t299) * t308 + (t264 * t363 + t298 * t266 + t299 * t268) * t309 + (t298 * t303 + t299 * t304 + t302 * t363) * t323) / 0.2e1 + t323 * ((-t264 * t309 - t265 * t308 - t302 * t323) * t336 + ((-t267 * t332 + t269 * t335) * t308 + (-t266 * t332 + t268 * t335) * t309 + (-t303 * t332 + t304 * t335) * t323) * t333) / 0.2e1 + m(6) * (t247 ^ 2 + t248 ^ 2 + t249 ^ 2) / 0.2e1 + t274 * ((t257 * t361 + t289 * t259 + t290 * t261) * t274 + (t256 * t361 + t258 * t289 + t260 * t290) * t275 + (t289 * t292 + t290 * t293 + t291 * t361) * t312) / 0.2e1 + t275 * ((t257 * t363 + t259 * t287 + t261 * t288) * t274 + (t256 * t363 + t287 * t258 + t288 * t260) * t275 + (t287 * t292 + t288 * t293 + t291 * t363) * t312) / 0.2e1 + t312 * ((-t256 * t275 - t257 * t274 - t291 * t312) * t336 + ((-t259 * t328 + t261 * t329) * t274 + (-t258 * t328 + t260 * t329) * t275 + (-t292 * t328 + t293 * t329) * t312) * t333) / 0.2e1 + (m(2) * (t319 ^ 2 + t320 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
