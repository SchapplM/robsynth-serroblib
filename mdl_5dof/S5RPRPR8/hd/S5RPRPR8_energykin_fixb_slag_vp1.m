% Calculate kinetic energy for
% S5RPRPR8
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
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR8_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR8_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR8_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR8_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:21:16
% EndTime: 2019-12-31 18:21:18
% DurationCPUTime: 1.41s
% Computational Cost: add. (1125->217), mult. (1197->346), div. (0->0), fcn. (1178->10), ass. (0->111)
t331 = sin(qJ(1));
t370 = pkin(1) * t331;
t328 = cos(pkin(9));
t368 = pkin(4) * t328;
t330 = sin(qJ(3));
t367 = Icges(4,4) * t330;
t332 = cos(qJ(3));
t366 = Icges(4,4) * t332;
t326 = qJ(1) + pkin(8);
t322 = sin(t326);
t327 = sin(pkin(9));
t365 = t322 * t327;
t364 = t322 * t330;
t363 = t322 * t332;
t324 = cos(t326);
t362 = t324 * t327;
t361 = t324 * t330;
t360 = t324 * t332;
t359 = t327 * t332;
t358 = t328 * t332;
t333 = cos(qJ(1));
t320 = qJD(1) * t333 * pkin(1);
t356 = qJD(1) * (pkin(2) * t324 + pkin(6) * t322) + t320;
t355 = qJD(3) * t322;
t354 = qJD(3) * t324;
t353 = qJD(4) * t330;
t352 = qJD(5) * t330;
t349 = -pkin(2) * t322 + pkin(6) * t324 - t370;
t314 = pkin(3) * t330 - qJ(4) * t332;
t348 = qJD(3) * (pkin(7) * t332 - t330 * t368 - t314);
t347 = qJD(3) * (rSges(5,3) * t332 - (rSges(5,1) * t328 - rSges(5,2) * t327) * t330 - t314);
t343 = pkin(3) * t332 + qJ(4) * t330;
t304 = t343 * t324;
t346 = qJD(1) * t304 + t322 * t353 + t356;
t303 = t343 * t322;
t345 = -t303 + t349;
t344 = rSges(4,1) * t332 - rSges(4,2) * t330;
t342 = Icges(4,1) * t332 - t367;
t341 = -Icges(4,2) * t330 + t366;
t340 = Icges(4,5) * t332 - Icges(4,6) * t330;
t275 = -Icges(4,6) * t324 + t322 * t341;
t277 = -Icges(4,5) * t324 + t322 * t342;
t339 = t275 * t330 - t277 * t332;
t276 = Icges(4,6) * t322 + t324 * t341;
t278 = Icges(4,5) * t322 + t324 * t342;
t338 = -t276 * t330 + t278 * t332;
t310 = Icges(4,2) * t332 + t367;
t311 = Icges(4,1) * t330 + t366;
t337 = -t310 * t330 + t311 * t332;
t336 = -qJD(4) * t332 + t303 * t355 + t304 * t354 + qJD(2);
t335 = pkin(7) * t330 + t332 * t368;
t325 = pkin(9) + qJ(5);
t323 = cos(t325);
t321 = sin(t325);
t318 = -qJD(5) * t332 + qJD(1);
t317 = rSges(2,1) * t333 - rSges(2,2) * t331;
t316 = rSges(2,1) * t331 + rSges(2,2) * t333;
t315 = rSges(4,1) * t330 + rSges(4,2) * t332;
t313 = t324 * t353;
t309 = Icges(4,5) * t330 + Icges(4,6) * t332;
t306 = t322 * t352 - t354;
t305 = t324 * t352 + t355;
t301 = -Icges(5,5) * t332 + (Icges(5,1) * t328 - Icges(5,4) * t327) * t330;
t300 = -Icges(5,6) * t332 + (Icges(5,4) * t328 - Icges(5,2) * t327) * t330;
t299 = -Icges(5,3) * t332 + (Icges(5,5) * t328 - Icges(5,6) * t327) * t330;
t298 = t324 * t358 + t365;
t297 = t322 * t328 - t324 * t359;
t296 = t322 * t358 - t362;
t295 = -t322 * t359 - t324 * t328;
t293 = t320 + qJD(1) * (rSges(3,1) * t324 - rSges(3,2) * t322);
t292 = (-rSges(3,1) * t322 - rSges(3,2) * t324 - t370) * qJD(1);
t291 = -rSges(6,3) * t332 + (rSges(6,1) * t323 - rSges(6,2) * t321) * t330;
t290 = -Icges(6,5) * t332 + (Icges(6,1) * t323 - Icges(6,4) * t321) * t330;
t289 = -Icges(6,6) * t332 + (Icges(6,4) * t323 - Icges(6,2) * t321) * t330;
t288 = -Icges(6,3) * t332 + (Icges(6,5) * t323 - Icges(6,6) * t321) * t330;
t287 = t321 * t322 + t323 * t360;
t286 = -t321 * t360 + t322 * t323;
t285 = -t321 * t324 + t323 * t363;
t284 = -t321 * t363 - t323 * t324;
t282 = rSges(4,3) * t322 + t324 * t344;
t281 = -rSges(4,3) * t324 + t322 * t344;
t274 = Icges(4,3) * t322 + t324 * t340;
t273 = -Icges(4,3) * t324 + t322 * t340;
t272 = pkin(4) * t365 + t324 * t335;
t271 = -pkin(4) * t362 + t322 * t335;
t270 = rSges(5,1) * t298 + rSges(5,2) * t297 + rSges(5,3) * t361;
t269 = rSges(5,1) * t296 + rSges(5,2) * t295 + rSges(5,3) * t364;
t268 = Icges(5,1) * t298 + Icges(5,4) * t297 + Icges(5,5) * t361;
t267 = Icges(5,1) * t296 + Icges(5,4) * t295 + Icges(5,5) * t364;
t266 = Icges(5,4) * t298 + Icges(5,2) * t297 + Icges(5,6) * t361;
t265 = Icges(5,4) * t296 + Icges(5,2) * t295 + Icges(5,6) * t364;
t264 = Icges(5,5) * t298 + Icges(5,6) * t297 + Icges(5,3) * t361;
t263 = Icges(5,5) * t296 + Icges(5,6) * t295 + Icges(5,3) * t364;
t262 = rSges(6,1) * t287 + rSges(6,2) * t286 + rSges(6,3) * t361;
t261 = rSges(6,1) * t285 + rSges(6,2) * t284 + rSges(6,3) * t364;
t260 = Icges(6,1) * t287 + Icges(6,4) * t286 + Icges(6,5) * t361;
t259 = Icges(6,1) * t285 + Icges(6,4) * t284 + Icges(6,5) * t364;
t258 = Icges(6,4) * t287 + Icges(6,2) * t286 + Icges(6,6) * t361;
t257 = Icges(6,4) * t285 + Icges(6,2) * t284 + Icges(6,6) * t364;
t256 = Icges(6,5) * t287 + Icges(6,6) * t286 + Icges(6,3) * t361;
t255 = Icges(6,5) * t285 + Icges(6,6) * t284 + Icges(6,3) * t364;
t254 = qJD(1) * t282 - t315 * t355 + t356;
t253 = -t315 * t354 + (-t281 + t349) * qJD(1);
t252 = qJD(2) + (t281 * t322 + t282 * t324) * qJD(3);
t251 = qJD(1) * t270 + t322 * t347 + t346;
t250 = t313 + t324 * t347 + (-t269 + t345) * qJD(1);
t249 = (t269 * t322 + t270 * t324) * qJD(3) + t336;
t248 = qJD(1) * t272 + t262 * t318 - t291 * t305 + t322 * t348 + t346;
t247 = -t261 * t318 + t291 * t306 + t313 + t324 * t348 + (-t271 + t345) * qJD(1);
t246 = t261 * t305 - t262 * t306 + (t271 * t322 + t272 * t324) * qJD(3) + t336;
t1 = m(3) * (qJD(2) ^ 2 + t292 ^ 2 + t293 ^ 2) / 0.2e1 + m(4) * (t252 ^ 2 + t253 ^ 2 + t254 ^ 2) / 0.2e1 + m(5) * (t249 ^ 2 + t250 ^ 2 + t251 ^ 2) / 0.2e1 + m(6) * (t246 ^ 2 + t247 ^ 2 + t248 ^ 2) / 0.2e1 + t305 * ((t256 * t361 + t258 * t286 + t260 * t287) * t305 + (t255 * t361 + t257 * t286 + t259 * t287) * t306 + (t286 * t289 + t287 * t290 + t288 * t361) * t318) / 0.2e1 + t306 * ((t256 * t364 + t258 * t284 + t260 * t285) * t305 + (t255 * t364 + t257 * t284 + t259 * t285) * t306 + (t284 * t289 + t285 * t290 + t288 * t364) * t318) / 0.2e1 + t318 * ((-t255 * t306 - t256 * t305 - t288 * t318) * t332 + ((-t258 * t321 + t260 * t323) * t305 + (-t257 * t321 + t259 * t323) * t306 + (-t289 * t321 + t290 * t323) * t318) * t330) / 0.2e1 + (((t332 * t276 + t278 * t330) * t322 - (t275 * t332 + t277 * t330) * t324 + (t263 * t324 - t264 * t322) * t332 + ((-t266 * t327 + t268 * t328) * t322 - (-t265 * t327 + t267 * t328) * t324) * t330) * qJD(3) + ((t310 - t299) * t332 + (-t300 * t327 + t301 * t328 + t311) * t330) * qJD(1)) * qJD(1) / 0.2e1 + (((-t263 * t361 - t265 * t297 - t267 * t298 + t339 * t324) * t324 + ((-t273 + t338) * t324 + t264 * t361 + t266 * t297 + t268 * t298 + t274 * t322) * t322) * qJD(3) + (t297 * t300 + t298 * t301 + t299 * t361 + t322 * t309 + t324 * t337) * qJD(1)) * t355 / 0.2e1 - (((-t263 * t364 - t265 * t295 - t296 * t267 + t273 * t324) * t324 + ((-t274 + t339) * t324 + t264 * t364 + t266 * t295 + t268 * t296 + t338 * t322) * t322) * qJD(3) + (t295 * t300 + t296 * t301 + t299 * t364 - t324 * t309 + t322 * t337) * qJD(1)) * t354 / 0.2e1 + (m(2) * (t316 ^ 2 + t317 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
