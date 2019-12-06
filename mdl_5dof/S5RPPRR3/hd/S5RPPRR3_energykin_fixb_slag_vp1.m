% Calculate kinetic energy for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-05 17:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR3_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR3_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:41:36
% EndTime: 2019-12-05 17:41:37
% DurationCPUTime: 0.69s
% Computational Cost: add. (783->137), mult. (559->220), div. (0->0), fcn. (452->10), ass. (0->83)
t296 = qJ(1) + pkin(8);
t289 = sin(t296);
t291 = cos(t296);
t295 = pkin(9) + qJ(4);
t292 = qJ(5) + t295;
t285 = sin(t292);
t286 = cos(t292);
t331 = Icges(6,4) * t286;
t313 = -Icges(6,2) * t285 + t331;
t252 = Icges(6,6) * t291 - t289 * t313;
t253 = Icges(6,6) * t289 + t291 * t313;
t332 = Icges(6,4) * t285;
t315 = Icges(6,1) * t286 - t332;
t254 = Icges(6,5) * t291 - t289 * t315;
t255 = Icges(6,5) * t289 + t291 * t315;
t269 = Icges(6,2) * t286 + t332;
t270 = Icges(6,1) * t285 + t331;
t324 = qJD(4) + qJD(5);
t273 = t324 * t289;
t274 = t324 * t291;
t340 = (t269 * t285 - t270 * t286) * qJD(1) + (t252 * t285 - t254 * t286) * t274 + (t253 * t285 - t255 * t286) * t273;
t300 = sin(qJ(1));
t338 = t300 * pkin(1);
t301 = cos(qJ(1));
t337 = t301 * pkin(1);
t298 = cos(pkin(9));
t335 = t298 * pkin(3);
t288 = sin(t295);
t334 = Icges(5,4) * t288;
t290 = cos(t295);
t333 = Icges(5,4) * t290;
t329 = qJD(1) * (-pkin(2) * t289 + qJ(3) * t291) + qJD(3) * t289;
t328 = pkin(4) * t290;
t326 = qJD(4) * t289;
t325 = qJD(4) * t291;
t323 = pkin(4) * qJD(4) * t288;
t322 = qJD(1) * (pkin(6) * t291 - t289 * t335) + t329;
t321 = -pkin(2) * t291 - qJ(3) * t289 - t337;
t320 = -pkin(6) * t289 - t291 * t335 + t321;
t297 = sin(pkin(9));
t319 = -rSges(4,1) * t298 + rSges(4,2) * t297;
t318 = rSges(5,1) * t290 - rSges(5,2) * t288;
t317 = rSges(6,1) * t286 - rSges(6,2) * t285;
t316 = Icges(5,1) * t290 - t334;
t314 = -Icges(5,2) * t288 + t333;
t312 = Icges(5,5) * t290 - Icges(5,6) * t288;
t311 = Icges(6,5) * t286 - Icges(6,6) * t285;
t260 = Icges(5,6) * t291 - t289 * t314;
t262 = Icges(5,5) * t291 - t289 * t316;
t308 = -t260 * t288 + t262 * t290;
t261 = Icges(5,6) * t289 + t291 * t314;
t263 = Icges(5,5) * t289 + t291 * t316;
t307 = t261 * t288 - t263 * t290;
t276 = Icges(5,2) * t290 + t334;
t277 = Icges(5,1) * t288 + t333;
t305 = t276 * t288 - t277 * t290;
t304 = (Icges(6,5) * t285 + Icges(6,6) * t286) * qJD(1) + (Icges(6,3) * t291 - t289 * t311) * t274 + (Icges(6,3) * t289 + t291 * t311) * t273;
t302 = qJD(2) ^ 2;
t284 = qJD(3) * t291;
t282 = rSges(2,1) * t301 - rSges(2,2) * t300;
t281 = -rSges(2,1) * t300 - rSges(2,2) * t301;
t278 = rSges(5,1) * t288 + rSges(5,2) * t290;
t275 = Icges(5,5) * t288 + Icges(5,6) * t290;
t271 = rSges(6,1) * t285 + rSges(6,2) * t286;
t267 = (-rSges(3,1) * t291 + rSges(3,2) * t289 - t337) * qJD(1);
t266 = (-rSges(3,1) * t289 - rSges(3,2) * t291 - t338) * qJD(1);
t265 = rSges(5,3) * t289 + t291 * t318;
t264 = rSges(5,3) * t291 - t289 * t318;
t259 = Icges(5,3) * t289 + t291 * t312;
t258 = Icges(5,3) * t291 - t289 * t312;
t257 = rSges(6,3) * t289 + t291 * t317;
t256 = rSges(6,3) * t291 - t289 * t317;
t247 = pkin(7) * t289 + t291 * t328;
t246 = pkin(7) * t291 - t289 * t328;
t245 = t284 + (-t289 * rSges(4,3) + t291 * t319 + t321) * qJD(1);
t244 = (t291 * rSges(4,3) + t289 * t319 - t338) * qJD(1) + t329;
t243 = qJD(2) + (-t264 * t289 + t265 * t291) * qJD(4);
t242 = t278 * t326 + t284 + (-t265 + t320) * qJD(1);
t241 = -t278 * t325 + (t264 - t338) * qJD(1) + t322;
t240 = t289 * t323 + t271 * t273 + t284 + (-t247 - t257 + t320) * qJD(1);
t239 = -t291 * t323 - t271 * t274 + (t246 + t256 - t338) * qJD(1) + t322;
t238 = -t256 * t273 + t257 * t274 + qJD(2) + (-t246 * t289 + t247 * t291) * qJD(4);
t1 = m(3) * (t266 ^ 2 + t267 ^ 2 + t302) / 0.2e1 + m(4) * (t244 ^ 2 + t245 ^ 2 + t302) / 0.2e1 + m(5) * (t241 ^ 2 + t242 ^ 2 + t243 ^ 2) / 0.2e1 + ((t275 * t291 + t289 * t305) * qJD(1) + (t291 ^ 2 * t258 + (t307 * t289 + (t259 - t308) * t291) * t289) * qJD(4)) * t325 / 0.2e1 + ((t289 * t275 - t291 * t305) * qJD(1) + (t289 ^ 2 * t259 + (t308 * t291 + (t258 - t307) * t289) * t291) * qJD(4)) * t326 / 0.2e1 + m(6) * (t238 ^ 2 + t239 ^ 2 + t240 ^ 2) / 0.2e1 + t274 * (t340 * t289 + t304 * t291) / 0.2e1 + t273 * (t304 * t289 - t340 * t291) / 0.2e1 + (((t260 * t290 + t262 * t288) * t291 + (t261 * t290 + t263 * t288) * t289) * qJD(4) + (t252 * t286 + t254 * t285) * t274 + (t253 * t286 + t255 * t285) * t273 + (t286 * t269 + t285 * t270 + t290 * t276 + t288 * t277) * qJD(1)) * qJD(1) / 0.2e1 + (m(2) * (t281 ^ 2 + t282 ^ 2) + Icges(2,3) + Icges(3,3) + Icges(4,2) * t298 ^ 2 + (Icges(4,1) * t297 + 0.2e1 * Icges(4,4) * t298) * t297) * qJD(1) ^ 2 / 0.2e1;
T = t1;
