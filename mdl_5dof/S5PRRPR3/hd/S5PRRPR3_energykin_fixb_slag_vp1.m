% Calculate kinetic energy for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR3_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR3_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR3_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:01
% EndTime: 2019-12-05 16:19:02
% DurationCPUTime: 1.00s
% Computational Cost: add. (973->157), mult. (739->249), div. (0->0), fcn. (622->8), ass. (0->97)
t362 = Icges(4,3) + Icges(5,3);
t305 = qJ(3) + pkin(9);
t298 = sin(t305);
t300 = cos(t305);
t307 = sin(qJ(3));
t308 = cos(qJ(3));
t361 = Icges(4,5) * t308 + Icges(5,5) * t300 - Icges(4,6) * t307 - Icges(5,6) * t298;
t304 = pkin(8) + qJ(2);
t297 = sin(t304);
t299 = cos(t304);
t360 = t361 * t297 - t362 * t299;
t359 = t362 * t297 + t361 * t299;
t358 = Icges(4,5) * t307 + Icges(5,5) * t298 + Icges(4,6) * t308 + Icges(5,6) * t300;
t346 = Icges(5,4) * t298;
t282 = Icges(5,2) * t300 + t346;
t345 = Icges(5,4) * t300;
t283 = Icges(5,1) * t298 + t345;
t348 = Icges(4,4) * t307;
t290 = Icges(4,2) * t308 + t348;
t347 = Icges(4,4) * t308;
t291 = Icges(4,1) * t307 + t347;
t357 = -t282 * t298 + t283 * t300 - t290 * t307 + t291 * t308;
t324 = -Icges(5,2) * t298 + t345;
t261 = Icges(5,6) * t297 + t299 * t324;
t327 = Icges(5,1) * t300 - t346;
t263 = Icges(5,5) * t297 + t299 * t327;
t325 = -Icges(4,2) * t307 + t347;
t269 = Icges(4,6) * t297 + t299 * t325;
t328 = Icges(4,1) * t308 - t348;
t271 = Icges(4,5) * t297 + t299 * t328;
t356 = -t261 * t298 + t263 * t300 - t269 * t307 + t271 * t308;
t260 = -Icges(5,6) * t299 + t297 * t324;
t262 = -Icges(5,5) * t299 + t297 * t327;
t268 = -Icges(4,6) * t299 + t297 * t325;
t270 = -Icges(4,5) * t299 + t297 * t328;
t355 = t260 * t298 - t262 * t300 + t268 * t307 - t270 * t308;
t352 = pkin(3) * t307;
t350 = t308 * pkin(3);
t301 = qJ(5) + t305;
t294 = sin(t301);
t344 = Icges(6,4) * t294;
t295 = cos(t301);
t343 = Icges(6,4) * t295;
t248 = -qJ(4) * t299 + t297 * t350;
t287 = pkin(2) * t297 - pkin(6) * t299;
t342 = -t248 - t287;
t341 = pkin(4) * t300;
t339 = qJD(3) * t297;
t338 = qJD(3) * t299;
t337 = qJD(3) + qJD(5);
t249 = qJ(4) * t297 + t299 * t350;
t336 = t248 * t339 + t249 * t338 + qJD(1);
t278 = qJD(2) * (pkin(2) * t299 + pkin(6) * t297);
t333 = qJD(2) * t249 - qJD(4) * t299 + t278;
t332 = rSges(4,1) * t308 - rSges(4,2) * t307;
t331 = rSges(5,1) * t300 - rSges(5,2) * t298;
t330 = rSges(6,1) * t295 - rSges(6,2) * t294;
t329 = qJD(3) * (-rSges(5,1) * t298 - rSges(5,2) * t300 - t352);
t326 = Icges(6,1) * t295 - t344;
t323 = -Icges(6,2) * t294 + t343;
t320 = Icges(6,5) * t295 - Icges(6,6) * t294;
t313 = qJD(3) * (-pkin(4) * t298 - t352);
t279 = t337 * t297;
t280 = t337 * t299;
t312 = qJD(2) * (Icges(6,5) * t294 + Icges(6,6) * t295) - (-Icges(6,3) * t299 + t297 * t320) * t280 + (Icges(6,3) * t297 + t299 * t320) * t279;
t252 = -Icges(6,6) * t299 + t297 * t323;
t253 = Icges(6,6) * t297 + t299 * t323;
t254 = -Icges(6,5) * t299 + t297 * t326;
t255 = Icges(6,5) * t297 + t299 * t326;
t275 = Icges(6,2) * t295 + t344;
t276 = Icges(6,1) * t294 + t343;
t311 = (-t253 * t294 + t255 * t295) * t279 - (-t252 * t294 + t254 * t295) * t280 + (-t275 * t294 + t276 * t295) * qJD(2);
t310 = qJD(1) ^ 2;
t309 = qJD(2) ^ 2;
t293 = qJD(4) * t297;
t292 = rSges(4,1) * t307 + rSges(4,2) * t308;
t286 = rSges(3,1) * t299 - rSges(3,2) * t297;
t284 = rSges(3,1) * t297 + rSges(3,2) * t299;
t277 = rSges(6,1) * t294 + rSges(6,2) * t295;
t273 = rSges(4,3) * t297 + t299 * t332;
t272 = -rSges(4,3) * t299 + t297 * t332;
t265 = rSges(5,3) * t297 + t299 * t331;
t264 = -rSges(5,3) * t299 + t297 * t331;
t257 = rSges(6,3) * t297 + t299 * t330;
t256 = -rSges(6,3) * t299 + t297 * t330;
t244 = pkin(7) * t297 + t299 * t341;
t243 = -pkin(7) * t299 + t297 * t341;
t242 = qJD(2) * t273 - t292 * t339 + t278;
t241 = -t292 * t338 + (-t272 - t287) * qJD(2);
t240 = qJD(1) + (t272 * t297 + t273 * t299) * qJD(3);
t239 = qJD(2) * t265 + t297 * t329 + t333;
t238 = t293 + t299 * t329 + (-t264 + t342) * qJD(2);
t237 = (t264 * t297 + t265 * t299) * qJD(3) + t336;
t236 = -t277 * t279 + t297 * t313 + (t244 + t257) * qJD(2) + t333;
t235 = -t277 * t280 + t293 + t299 * t313 + (-t243 - t256 + t342) * qJD(2);
t234 = t256 * t279 + t257 * t280 + (t243 * t297 + t244 * t299) * qJD(3) + t336;
t1 = m(2) * t310 / 0.2e1 + m(3) * (t310 + (t284 ^ 2 + t286 ^ 2) * t309) / 0.2e1 + t309 * Icges(3,3) / 0.2e1 + m(4) * (t240 ^ 2 + t241 ^ 2 + t242 ^ 2) / 0.2e1 + m(5) * (t237 ^ 2 + t238 ^ 2 + t239 ^ 2) / 0.2e1 + m(6) * (t234 ^ 2 + t235 ^ 2 + t236 ^ 2) / 0.2e1 + t279 * (t312 * t297 + t311 * t299) / 0.2e1 - t280 * (t311 * t297 - t312 * t299) / 0.2e1 + ((t359 * t297 ^ 2 + (t355 * t299 + (t356 - t360) * t297) * t299) * qJD(3) + (t358 * t297 + t357 * t299) * qJD(2)) * t339 / 0.2e1 - ((t360 * t299 ^ 2 + (t356 * t297 + (t355 - t359) * t299) * t297) * qJD(3) + (t357 * t297 - t358 * t299) * qJD(2)) * t338 / 0.2e1 + ((t253 * t295 + t255 * t294) * t279 - (t252 * t295 + t254 * t294) * t280 + ((-t260 * t300 - t262 * t298 - t268 * t308 - t270 * t307) * t299 + (t261 * t300 + t263 * t298 + t269 * t308 + t271 * t307) * t297) * qJD(3) + (t295 * t275 + t294 * t276 + t300 * t282 + t298 * t283 + t308 * t290 + t307 * t291) * qJD(2)) * qJD(2) / 0.2e1;
T = t1;
