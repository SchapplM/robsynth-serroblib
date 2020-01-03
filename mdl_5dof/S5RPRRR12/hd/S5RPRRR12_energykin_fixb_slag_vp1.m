% Calculate kinetic energy for
% S5RPRRR12
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
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR12_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR12_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR12_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR12_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR12_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:12:27
% EndTime: 2019-12-31 19:12:29
% DurationCPUTime: 1.12s
% Computational Cost: add. (682->194), mult. (992->318), div. (0->0), fcn. (918->8), ass. (0->108)
t324 = sin(qJ(1));
t327 = cos(qJ(1));
t321 = qJ(3) + qJ(4);
t320 = cos(t321);
t319 = sin(t321);
t361 = Icges(5,4) * t319;
t340 = Icges(5,2) * t320 + t361;
t270 = Icges(5,6) * t327 + t340 * t324;
t271 = Icges(5,6) * t324 - t340 * t327;
t360 = Icges(5,4) * t320;
t342 = Icges(5,1) * t319 + t360;
t272 = Icges(5,5) * t327 + t342 * t324;
t273 = Icges(5,5) * t324 - t342 * t327;
t297 = -Icges(5,2) * t319 + t360;
t298 = Icges(5,1) * t320 - t361;
t316 = qJD(3) * t324;
t302 = qJD(4) * t324 + t316;
t317 = qJD(3) * t327;
t303 = qJD(4) * t327 + t317;
t368 = (t270 * t320 + t272 * t319) * t303 + (t271 * t320 + t273 * t319) * t302 + (t297 * t320 + t298 * t319) * qJD(1);
t323 = sin(qJ(3));
t366 = pkin(3) * t323;
t363 = Icges(4,4) * t323;
t326 = cos(qJ(3));
t362 = Icges(4,4) * t326;
t359 = t320 * t324;
t358 = t320 * t327;
t322 = sin(qJ(5));
t357 = t324 * t322;
t325 = cos(qJ(5));
t356 = t324 * t325;
t355 = t327 * t322;
t354 = t327 * t325;
t301 = qJD(1) * (t327 * pkin(1) + t324 * qJ(2));
t353 = qJD(1) * t327 * pkin(6) + t301;
t318 = qJD(2) * t324;
t350 = pkin(3) * qJD(3) * t326;
t352 = t324 * t350 + t318;
t351 = qJD(5) * t320;
t307 = t324 * pkin(1) - t327 * qJ(2);
t349 = -t324 * pkin(6) - t307;
t290 = pkin(7) * t324 - t327 * t366;
t348 = -t290 + t349;
t347 = pkin(4) * t319 - pkin(8) * t320;
t291 = pkin(7) * t327 + t324 * t366;
t346 = t290 * t317 - t291 * t316;
t345 = rSges(4,1) * t323 + rSges(4,2) * t326;
t344 = rSges(5,1) * t319 + rSges(5,2) * t320;
t343 = Icges(4,1) * t323 + t362;
t341 = Icges(4,2) * t326 + t363;
t339 = Icges(4,5) * t323 + Icges(4,6) * t326;
t338 = Icges(5,5) * t319 + Icges(5,6) * t320;
t279 = Icges(4,6) * t327 + t341 * t324;
t281 = Icges(4,5) * t327 + t343 * t324;
t335 = -t279 * t326 - t281 * t323;
t280 = Icges(4,6) * t324 - t341 * t327;
t282 = Icges(4,5) * t324 - t343 * t327;
t334 = t280 * t326 + t282 * t323;
t305 = -Icges(4,2) * t323 + t362;
t306 = Icges(4,1) * t326 - t363;
t332 = t305 * t326 + t306 * t323;
t331 = (Icges(5,5) * t320 - Icges(5,6) * t319) * qJD(1) + (Icges(5,3) * t327 + t338 * t324) * t303 + (Icges(5,3) * t324 - t338 * t327) * t302;
t330 = qJD(1) * t291 + (-qJD(2) - t350) * t327 + t353;
t311 = qJD(5) * t319 + qJD(1);
t310 = t327 * rSges(2,1) - t324 * rSges(2,2);
t309 = t326 * rSges(4,1) - t323 * rSges(4,2);
t308 = t324 * rSges(2,1) + t327 * rSges(2,2);
t304 = Icges(4,5) * t326 - Icges(4,6) * t323;
t300 = t320 * pkin(4) + t319 * pkin(8);
t299 = t320 * rSges(5,1) - t319 * rSges(5,2);
t295 = -t319 * t354 + t357;
t294 = t319 * t355 + t356;
t293 = t319 * t356 + t355;
t292 = -t319 * t357 + t354;
t289 = t347 * t327;
t288 = t347 * t324;
t287 = t324 * rSges(4,3) - t345 * t327;
t286 = t327 * rSges(4,3) + t345 * t324;
t285 = -t324 * t351 + t303;
t284 = t327 * t351 + t302;
t278 = Icges(4,3) * t324 - t339 * t327;
t277 = Icges(4,3) * t327 + t339 * t324;
t275 = t324 * rSges(5,3) - t344 * t327;
t274 = t327 * rSges(5,3) + t344 * t324;
t267 = t319 * rSges(6,3) + (rSges(6,1) * t325 - rSges(6,2) * t322) * t320;
t266 = Icges(6,5) * t319 + (Icges(6,1) * t325 - Icges(6,4) * t322) * t320;
t265 = Icges(6,6) * t319 + (Icges(6,4) * t325 - Icges(6,2) * t322) * t320;
t264 = Icges(6,3) * t319 + (Icges(6,5) * t325 - Icges(6,6) * t322) * t320;
t263 = t301 - qJD(2) * t327 + qJD(1) * (-t327 * rSges(3,2) + t324 * rSges(3,3));
t262 = t318 + (t324 * rSges(3,2) + t327 * rSges(3,3) - t307) * qJD(1);
t261 = t295 * rSges(6,1) + t294 * rSges(6,2) + rSges(6,3) * t358;
t260 = t293 * rSges(6,1) + t292 * rSges(6,2) - rSges(6,3) * t359;
t259 = Icges(6,1) * t295 + Icges(6,4) * t294 + Icges(6,5) * t358;
t258 = Icges(6,1) * t293 + Icges(6,4) * t292 - Icges(6,5) * t359;
t257 = Icges(6,4) * t295 + Icges(6,2) * t294 + Icges(6,6) * t358;
t256 = Icges(6,4) * t293 + Icges(6,2) * t292 - Icges(6,6) * t359;
t255 = Icges(6,5) * t295 + Icges(6,6) * t294 + Icges(6,3) * t358;
t254 = Icges(6,5) * t293 + Icges(6,6) * t292 - Icges(6,3) * t359;
t253 = (-t286 * t324 + t287 * t327) * qJD(3);
t252 = qJD(1) * t286 + (-qJD(3) * t309 - qJD(2)) * t327 + t353;
t251 = t309 * t316 + t318 + (-t287 + t349) * qJD(1);
t250 = qJD(1) * t274 - t303 * t299 + t330;
t249 = t302 * t299 + (-t275 + t348) * qJD(1) + t352;
t248 = -t302 * t274 + t303 * t275 + t346;
t247 = qJD(1) * t288 + t311 * t260 - t285 * t267 - t303 * t300 + t330;
t246 = -t311 * t261 + t284 * t267 + t302 * t300 + (t289 + t348) * qJD(1) + t352;
t245 = -t284 * t260 + t285 * t261 - t302 * t288 - t303 * t289 + t346;
t1 = m(3) * (t262 ^ 2 + t263 ^ 2) / 0.2e1 + m(4) * (t251 ^ 2 + t252 ^ 2 + t253 ^ 2) / 0.2e1 + ((t327 * t304 + t332 * t324) * qJD(1) + (t327 ^ 2 * t277 + (t334 * t324 + (t278 - t335) * t327) * t324) * qJD(3)) * t317 / 0.2e1 + ((t324 * t304 - t332 * t327) * qJD(1) + (t324 ^ 2 * t278 + (t335 * t327 + (t277 - t334) * t324) * t327) * qJD(3)) * t316 / 0.2e1 + m(5) * (t248 ^ 2 + t249 ^ 2 + t250 ^ 2) / 0.2e1 + t303 * (t368 * t324 + t331 * t327) / 0.2e1 + t302 * (t331 * t324 - t368 * t327) / 0.2e1 + m(6) * (t245 ^ 2 + t246 ^ 2 + t247 ^ 2) / 0.2e1 + t285 * ((-t254 * t359 + t292 * t256 + t293 * t258) * t285 + (-t255 * t359 + t292 * t257 + t293 * t259) * t284 + (-t264 * t359 + t292 * t265 + t293 * t266) * t311) / 0.2e1 + t284 * ((t254 * t358 + t294 * t256 + t295 * t258) * t285 + (t255 * t358 + t294 * t257 + t295 * t259) * t284 + (t264 * t358 + t294 * t265 + t295 * t266) * t311) / 0.2e1 + t311 * ((t254 * t285 + t255 * t284 + t264 * t311) * t319 + ((-t256 * t322 + t258 * t325) * t285 + (-t257 * t322 + t259 * t325) * t284 + (-t265 * t322 + t266 * t325) * t311) * t320) / 0.2e1 + (((-t323 * t279 + t326 * t281) * t327 + (-t323 * t280 + t326 * t282) * t324) * qJD(3) + (-t319 * t270 + t320 * t272) * t303 + (-t319 * t271 + t320 * t273) * t302 + (-t319 * t297 + t320 * t298 - t323 * t305 + t326 * t306) * qJD(1)) * qJD(1) / 0.2e1 + (m(2) * (t308 ^ 2 + t310 ^ 2) + Icges(2,3) + Icges(3,1)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
