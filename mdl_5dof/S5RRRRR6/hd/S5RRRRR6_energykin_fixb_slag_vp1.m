% Calculate kinetic energy for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2020-01-03 12:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR6_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR6_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:14:45
% EndTime: 2020-01-03 12:14:46
% DurationCPUTime: 0.90s
% Computational Cost: add. (1062->167), mult. (788->278), div. (0->0), fcn. (658->10), ass. (0->108)
t311 = qJ(1) + qJ(2);
t303 = sin(t311);
t305 = cos(t311);
t310 = qJ(3) + qJ(4);
t306 = qJ(5) + t310;
t297 = sin(t306);
t298 = cos(t306);
t353 = Icges(6,4) * t298;
t333 = -Icges(6,2) * t297 + t353;
t251 = -Icges(6,6) * t305 + t303 * t333;
t252 = -Icges(6,6) * t303 - t305 * t333;
t354 = Icges(6,4) * t297;
t336 = Icges(6,1) * t298 - t354;
t253 = -Icges(6,5) * t305 + t303 * t336;
t254 = -Icges(6,5) * t303 - t305 * t336;
t346 = -qJD(3) - qJD(4);
t344 = -qJD(5) + t346;
t275 = t344 * t303;
t276 = t344 * t305;
t279 = Icges(6,2) * t298 + t354;
t280 = Icges(6,1) * t297 + t353;
t308 = qJD(1) + qJD(2);
t366 = (t279 * t297 - t280 * t298) * t308 + (t251 * t297 - t253 * t298) * t276 + (t252 * t297 - t254 * t298) * t275;
t302 = sin(t310);
t304 = cos(t310);
t355 = Icges(5,4) * t304;
t334 = -Icges(5,2) * t302 + t355;
t259 = -Icges(5,6) * t305 + t303 * t334;
t260 = -Icges(5,6) * t303 - t305 * t334;
t356 = Icges(5,4) * t302;
t337 = Icges(5,1) * t304 - t356;
t261 = -Icges(5,5) * t305 + t303 * t337;
t262 = -Icges(5,5) * t303 - t305 * t337;
t282 = t346 * t303;
t283 = t346 * t305;
t285 = Icges(5,2) * t304 + t356;
t286 = Icges(5,1) * t302 + t355;
t365 = (t285 * t302 - t286 * t304) * t308 + (t259 * t302 - t261 * t304) * t283 + (t260 * t302 - t262 * t304) * t282;
t362 = pkin(4) * t302;
t314 = cos(qJ(3));
t361 = t314 * pkin(3);
t359 = pkin(1) * qJD(1);
t312 = sin(qJ(3));
t358 = Icges(4,4) * t312;
t357 = Icges(4,4) * t314;
t248 = -pkin(8) * t303 - t305 * t361;
t288 = -t305 * pkin(2) - t303 * pkin(7);
t352 = -t248 - t288;
t313 = sin(qJ(1));
t300 = t313 * t359;
t351 = t308 * (t303 * pkin(2) - t305 * pkin(7)) + t300;
t350 = pkin(4) * t304;
t348 = qJD(3) * t303;
t347 = qJD(3) * t305;
t345 = pkin(3) * qJD(3) * t312;
t247 = -pkin(8) * t305 + t303 * t361;
t343 = t308 * t247 + t305 * t345 + t351;
t342 = t247 * t348 - t248 * t347;
t341 = rSges(4,1) * t314 - rSges(4,2) * t312;
t340 = rSges(5,1) * t304 - rSges(5,2) * t302;
t339 = rSges(6,1) * t298 - rSges(6,2) * t297;
t338 = Icges(4,1) * t314 - t358;
t335 = -Icges(4,2) * t312 + t357;
t332 = Icges(4,5) * t314 - Icges(4,6) * t312;
t331 = Icges(5,5) * t304 - Icges(5,6) * t302;
t330 = Icges(6,5) * t298 - Icges(6,6) * t297;
t267 = -Icges(4,6) * t305 + t303 * t335;
t269 = -Icges(4,5) * t305 + t303 * t338;
t325 = -t267 * t312 + t269 * t314;
t268 = -Icges(4,6) * t303 - t305 * t335;
t270 = -Icges(4,5) * t303 - t305 * t338;
t324 = t268 * t312 - t270 * t314;
t292 = Icges(4,2) * t314 + t358;
t293 = Icges(4,1) * t312 + t357;
t321 = t292 * t312 - t293 * t314;
t315 = cos(qJ(1));
t301 = t315 * t359;
t320 = -t303 * t345 + t301;
t319 = -(-Icges(6,3) * t305 + t303 * t330) * t276 - (-Icges(6,3) * t303 - t305 * t330) * t275 - (Icges(6,5) * t297 + Icges(6,6) * t298) * t308;
t318 = -(-Icges(5,3) * t305 + t303 * t331) * t283 - (-Icges(5,3) * t303 - t305 * t331) * t282 - (Icges(5,5) * t302 + Icges(5,6) * t304) * t308;
t296 = -t315 * rSges(2,1) + t313 * rSges(2,2);
t295 = t313 * rSges(2,1) + t315 * rSges(2,2);
t294 = t312 * rSges(4,1) + t314 * rSges(4,2);
t291 = Icges(4,5) * t312 + Icges(4,6) * t314;
t287 = t302 * rSges(5,1) + t304 * rSges(5,2);
t281 = t297 * rSges(6,1) + t298 * rSges(6,2);
t274 = t301 - t308 * (-t305 * rSges(3,1) + t303 * rSges(3,2));
t273 = t300 + t308 * (t303 * rSges(3,1) + t305 * rSges(3,2));
t272 = -t303 * rSges(4,3) - t305 * t341;
t271 = -t305 * rSges(4,3) + t303 * t341;
t266 = -Icges(4,3) * t303 - t305 * t332;
t265 = -Icges(4,3) * t305 + t303 * t332;
t264 = -t303 * rSges(5,3) - t305 * t340;
t263 = -t305 * rSges(5,3) + t303 * t340;
t256 = -t303 * rSges(6,3) - t305 * t339;
t255 = -t305 * rSges(6,3) + t303 * t339;
t244 = -pkin(9) * t303 - t305 * t350;
t243 = -pkin(9) * t305 + t303 * t350;
t242 = (t271 * t303 - t272 * t305) * qJD(3);
t241 = -t294 * t348 + t301 + (-t272 - t288) * t308;
t240 = t308 * t271 + t294 * t347 + t351;
t239 = t282 * t287 + (-t264 + t352) * t308 + t320;
t238 = t308 * t263 - t283 * t287 + t343;
t237 = -t282 * t263 + t283 * t264 + t342;
t236 = t282 * t362 + t275 * t281 + (-t244 - t256 + t352) * t308 + t320;
t235 = -t283 * t362 - t276 * t281 + (t243 + t255) * t308 + t343;
t234 = -t282 * t243 + t283 * t244 - t275 * t255 + t276 * t256 + t342;
t1 = m(3) * (t273 ^ 2 + t274 ^ 2) / 0.2e1 + t308 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t240 ^ 2 + t241 ^ 2 + t242 ^ 2) / 0.2e1 - ((-t305 * t291 - t303 * t321) * t308 + (t305 ^ 2 * t265 + (t324 * t303 + (t266 - t325) * t305) * t303) * qJD(3)) * t347 / 0.2e1 - ((-t303 * t291 + t305 * t321) * t308 + (t303 ^ 2 * t266 + (t325 * t305 + (t265 - t324) * t303) * t305) * qJD(3)) * t348 / 0.2e1 + m(5) * (t237 ^ 2 + t238 ^ 2 + t239 ^ 2) / 0.2e1 + t283 * (-t365 * t303 + t318 * t305) / 0.2e1 + t282 * (t318 * t303 + t365 * t305) / 0.2e1 + m(6) * (t234 ^ 2 + t235 ^ 2 + t236 ^ 2) / 0.2e1 + t276 * (-t366 * t303 + t319 * t305) / 0.2e1 + t275 * (t319 * t303 + t366 * t305) / 0.2e1 + (m(2) * (t295 ^ 2 + t296 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((-(t314 * t267 + t312 * t269) * t305 - (t314 * t268 + t312 * t270) * t303) * qJD(3) + (t304 * t259 + t302 * t261) * t283 + (t304 * t260 + t302 * t262) * t282 + (t298 * t251 + t297 * t253) * t276 + (t298 * t252 + t297 * t254) * t275 + (t298 * t279 + t297 * t280 + t304 * t285 + t302 * t286 + t314 * t292 + t312 * t293) * t308) * t308 / 0.2e1;
T = t1;
