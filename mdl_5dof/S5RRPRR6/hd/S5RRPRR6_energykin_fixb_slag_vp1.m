% Calculate kinetic energy for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR6_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR6_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:05:28
% EndTime: 2020-01-03 12:05:29
% DurationCPUTime: 0.89s
% Computational Cost: add. (981->172), mult. (1025->295), div. (0->0), fcn. (1032->10), ass. (0->86)
t301 = sin(pkin(9));
t300 = qJ(1) + qJ(2);
t294 = sin(t300);
t296 = cos(t300);
t305 = cos(qJ(4));
t302 = cos(pkin(9));
t303 = sin(qJ(4));
t320 = t302 * t303;
t279 = -t294 * t320 - t296 * t305;
t319 = t302 * t305;
t321 = t296 * t303;
t280 = t294 * t319 - t321;
t281 = -t294 * t305 + t296 * t320;
t324 = t294 * t303;
t282 = -t296 * t319 - t324;
t323 = t296 * t301;
t326 = t294 * t301;
t309 = (Icges(5,5) * t280 + Icges(5,6) * t279 + Icges(5,3) * t326) * t294 - (Icges(5,5) * t282 + Icges(5,6) * t281 - Icges(5,3) * t323) * t296;
t334 = t301 * t309;
t298 = qJD(1) + qJD(2);
t333 = -qJD(3) + t298 * (rSges(4,1) * t302 - rSges(4,2) * t301);
t329 = pkin(4) * t305;
t332 = pkin(8) * t301 + t302 * t329;
t327 = pkin(1) * qJD(1);
t325 = t294 * t302;
t322 = t296 * t302;
t304 = sin(qJ(1));
t291 = t304 * t327;
t318 = t298 * (pkin(2) * t294 - qJ(3) * t296) + t291;
t317 = qJD(4) * t301;
t316 = qJD(4) + qJD(5);
t314 = pkin(3) * t302 + pkin(7) * t301;
t315 = t298 * t314 * t294 + t318;
t312 = -(-pkin(8) * t302 + t301 * t329) * t317 - qJD(3);
t311 = -(-rSges(5,3) * t302 + (rSges(5,1) * t305 - rSges(5,2) * t303) * t301) * t317 - qJD(3);
t286 = -pkin(2) * t296 - qJ(3) * t294;
t306 = cos(qJ(1));
t292 = t306 * t327;
t310 = t292 + (t314 * t296 - t286) * t298;
t299 = qJ(4) + qJ(5);
t295 = cos(t299);
t293 = sin(t299);
t289 = -t302 * qJD(4) + t298;
t288 = -rSges(2,1) * t306 + rSges(2,2) * t304;
t287 = rSges(2,1) * t304 + rSges(2,2) * t306;
t285 = -t302 * t316 + t298;
t278 = t316 * t323;
t277 = t316 * t326;
t275 = -Icges(5,5) * t302 + (Icges(5,1) * t305 - Icges(5,4) * t303) * t301;
t274 = -Icges(5,6) * t302 + (Icges(5,4) * t305 - Icges(5,2) * t303) * t301;
t273 = -Icges(5,3) * t302 + (Icges(5,5) * t305 - Icges(5,6) * t303) * t301;
t271 = -t293 * t294 - t295 * t322;
t270 = t293 * t322 - t294 * t295;
t269 = -t293 * t296 + t295 * t325;
t268 = -t293 * t325 - t295 * t296;
t267 = t292 - t298 * (-rSges(3,1) * t296 + rSges(3,2) * t294);
t266 = t291 + t298 * (rSges(3,1) * t294 + rSges(3,2) * t296);
t265 = -rSges(6,3) * t302 + (rSges(6,1) * t295 - rSges(6,2) * t293) * t301;
t264 = -Icges(6,5) * t302 + (Icges(6,1) * t295 - Icges(6,4) * t293) * t301;
t263 = -Icges(6,6) * t302 + (Icges(6,4) * t295 - Icges(6,2) * t293) * t301;
t262 = -Icges(6,3) * t302 + (Icges(6,5) * t295 - Icges(6,6) * t293) * t301;
t260 = rSges(5,1) * t282 + rSges(5,2) * t281 - rSges(5,3) * t323;
t259 = rSges(5,1) * t280 + rSges(5,2) * t279 + rSges(5,3) * t326;
t258 = -pkin(4) * t324 - t332 * t296;
t257 = -pkin(4) * t321 + t332 * t294;
t256 = Icges(5,1) * t282 + Icges(5,4) * t281 - Icges(5,5) * t323;
t255 = Icges(5,1) * t280 + Icges(5,4) * t279 + Icges(5,5) * t326;
t254 = Icges(5,4) * t282 + Icges(5,2) * t281 - Icges(5,6) * t323;
t253 = Icges(5,4) * t280 + Icges(5,2) * t279 + Icges(5,6) * t326;
t250 = rSges(6,1) * t271 + rSges(6,2) * t270 - rSges(6,3) * t323;
t249 = rSges(6,1) * t269 + rSges(6,2) * t268 + rSges(6,3) * t326;
t248 = t292 + (rSges(4,3) * t294 - t286) * t298 + t333 * t296;
t247 = -t298 * rSges(4,3) * t296 + t333 * t294 + t318;
t246 = Icges(6,1) * t271 + Icges(6,4) * t270 - Icges(6,5) * t323;
t245 = Icges(6,1) * t269 + Icges(6,4) * t268 + Icges(6,5) * t326;
t244 = Icges(6,4) * t271 + Icges(6,2) * t270 - Icges(6,6) * t323;
t243 = Icges(6,4) * t269 + Icges(6,2) * t268 + Icges(6,6) * t326;
t242 = Icges(6,5) * t271 + Icges(6,6) * t270 - Icges(6,3) * t323;
t241 = Icges(6,5) * t269 + Icges(6,6) * t268 + Icges(6,3) * t326;
t240 = (t259 * t296 + t260 * t294) * t317;
t239 = -t260 * t289 + t296 * t311 + t310;
t238 = t259 * t289 + t294 * t311 + t315;
t237 = -t250 * t285 - t258 * t289 - t265 * t278 + t296 * t312 + t310;
t236 = t249 * t285 + t257 * t289 - t265 * t277 + t294 * t312 + t315;
t235 = t249 * t278 + t250 * t277 + (t257 * t296 + t258 * t294) * t317;
t1 = m(3) * (t266 ^ 2 + t267 ^ 2) / 0.2e1 + m(4) * (t247 ^ 2 + t248 ^ 2) / 0.2e1 + m(5) * (t238 ^ 2 + t239 ^ 2 + t240 ^ 2) / 0.2e1 + t289 * ((-t302 * t273 + (-t274 * t303 + t275 * t305) * t301) * t289 + (((-t253 * t303 + t255 * t305) * t294 - (-t254 * t303 + t256 * t305) * t296) * t301 - t309 * t302) * t317) / 0.2e1 + t294 * ((t273 * t326 + t274 * t279 + t275 * t280) * t289 + (-(t254 * t279 + t256 * t280) * t296 + (t279 * t253 + t280 * t255 + t334) * t294) * t317) * t317 / 0.2e1 - t296 * ((-t273 * t323 + t274 * t281 + t275 * t282) * t289 + ((t253 * t281 + t255 * t282) * t294 + (-t281 * t254 - t282 * t256 - t334) * t296) * t317) * t317 / 0.2e1 + m(6) * (t235 ^ 2 + t236 ^ 2 + t237 ^ 2) / 0.2e1 + t285 * ((-t241 * t277 + t242 * t278 - t262 * t285) * t302 + ((-t263 * t293 + t264 * t295) * t285 + (-t243 * t293 + t245 * t295) * t277 - (-t244 * t293 + t246 * t295) * t278) * t301) / 0.2e1 + t277 * ((t262 * t326 + t263 * t268 + t264 * t269) * t285 + (t241 * t326 + t243 * t268 + t245 * t269) * t277 - (t242 * t326 + t244 * t268 + t246 * t269) * t278) / 0.2e1 - t278 * ((-t262 * t323 + t263 * t270 + t264 * t271) * t285 + (-t241 * t323 + t243 * t270 + t245 * t271) * t277 - (-t242 * t323 + t244 * t270 + t246 * t271) * t278) / 0.2e1 + (Icges(3,3) + Icges(4,2) * t302 ^ 2 + (Icges(4,1) * t301 + 0.2e1 * Icges(4,4) * t302) * t301) * t298 ^ 2 / 0.2e1 + (m(2) * (t287 ^ 2 + t288 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
