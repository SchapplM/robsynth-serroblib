% Calculate kinetic energy for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR7_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR7_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR7_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:02
% EndTime: 2019-12-31 16:36:03
% DurationCPUTime: 0.92s
% Computational Cost: add. (1064->216), mult. (2733->357), div. (0->0), fcn. (3302->10), ass. (0->99)
t337 = qJD(2) ^ 2;
t336 = cos(qJ(3));
t314 = sin(pkin(8));
t315 = sin(pkin(4));
t335 = t314 * t315;
t319 = sin(qJ(3));
t334 = t315 * t319;
t322 = cos(qJ(2));
t333 = t315 * t322;
t316 = cos(pkin(8));
t332 = t316 * t315;
t317 = cos(pkin(4));
t320 = sin(qJ(2));
t331 = t317 * t320;
t330 = t317 * t322;
t303 = t314 * t330 + t316 * t320;
t329 = qJD(2) * t315;
t311 = t314 * t329;
t293 = qJD(3) * t303 + t311;
t313 = qJD(2) * t317;
t327 = t315 * t336;
t326 = t316 * t329;
t301 = t314 * t320 - t316 * t330;
t302 = t314 * t322 + t316 * t331;
t283 = pkin(2) * t302 + pkin(6) * t301;
t304 = -t314 * t331 + t316 * t322;
t284 = pkin(2) * t304 + pkin(6) * t303;
t325 = t283 * t311 + t284 * t326 + qJD(1);
t294 = qJD(3) * t301 - t326;
t308 = -qJD(3) * t333 + t313;
t307 = (pkin(2) * t320 - pkin(6) * t322) * t315;
t324 = t284 * t313 - t307 * t311;
t323 = (-t283 * t317 - t307 * t332) * qJD(2);
t321 = cos(qJ(4));
t318 = sin(qJ(4));
t306 = t317 * t319 + t320 * t327;
t305 = -t317 * t336 + t320 * t334;
t298 = t317 * rSges(3,3) + (rSges(3,1) * t320 + rSges(3,2) * t322) * t315;
t297 = Icges(3,5) * t317 + (Icges(3,1) * t320 + Icges(3,4) * t322) * t315;
t296 = Icges(3,6) * t317 + (Icges(3,4) * t320 + Icges(3,2) * t322) * t315;
t295 = Icges(3,3) * t317 + (Icges(3,5) * t320 + Icges(3,6) * t322) * t315;
t292 = t306 * t321 - t318 * t333;
t291 = -t306 * t318 - t321 * t333;
t290 = t304 * t336 + t314 * t334;
t289 = t304 * t319 - t314 * t327;
t288 = t302 * t336 - t319 * t332;
t287 = t302 * t319 + t316 * t327;
t286 = qJD(4) * t305 + t308;
t285 = pkin(3) * t306 + pkin(7) * t305;
t281 = t306 * rSges(4,1) - t305 * rSges(4,2) - rSges(4,3) * t333;
t280 = Icges(4,1) * t306 - Icges(4,4) * t305 - Icges(4,5) * t333;
t279 = Icges(4,4) * t306 - Icges(4,2) * t305 - Icges(4,6) * t333;
t278 = Icges(4,5) * t306 - Icges(4,6) * t305 - Icges(4,3) * t333;
t275 = rSges(3,1) * t304 - rSges(3,2) * t303 + rSges(3,3) * t335;
t274 = rSges(3,1) * t302 - rSges(3,2) * t301 - rSges(3,3) * t332;
t273 = Icges(3,1) * t304 - Icges(3,4) * t303 + Icges(3,5) * t335;
t272 = Icges(3,1) * t302 - Icges(3,4) * t301 - Icges(3,5) * t332;
t271 = Icges(3,4) * t304 - Icges(3,2) * t303 + Icges(3,6) * t335;
t270 = Icges(3,4) * t302 - Icges(3,2) * t301 - Icges(3,6) * t332;
t269 = Icges(3,5) * t304 - Icges(3,6) * t303 + Icges(3,3) * t335;
t268 = Icges(3,5) * t302 - Icges(3,6) * t301 - Icges(3,3) * t332;
t267 = t290 * t321 + t303 * t318;
t266 = -t290 * t318 + t303 * t321;
t265 = t288 * t321 + t301 * t318;
t264 = -t288 * t318 + t301 * t321;
t263 = qJD(4) * t287 + t294;
t262 = qJD(4) * t289 + t293;
t261 = pkin(3) * t290 + pkin(7) * t289;
t260 = pkin(3) * t288 + pkin(7) * t287;
t259 = (-t274 * t317 - t298 * t332) * qJD(2);
t258 = (t275 * t317 - t298 * t335) * qJD(2);
t257 = rSges(5,1) * t292 + rSges(5,2) * t291 + rSges(5,3) * t305;
t256 = Icges(5,1) * t292 + Icges(5,4) * t291 + Icges(5,5) * t305;
t255 = Icges(5,4) * t292 + Icges(5,2) * t291 + Icges(5,6) * t305;
t254 = Icges(5,5) * t292 + Icges(5,6) * t291 + Icges(5,3) * t305;
t253 = rSges(4,1) * t290 - rSges(4,2) * t289 + rSges(4,3) * t303;
t252 = rSges(4,1) * t288 - rSges(4,2) * t287 + rSges(4,3) * t301;
t251 = Icges(4,1) * t290 - Icges(4,4) * t289 + Icges(4,5) * t303;
t250 = Icges(4,1) * t288 - Icges(4,4) * t287 + Icges(4,5) * t301;
t249 = Icges(4,4) * t290 - Icges(4,2) * t289 + Icges(4,6) * t303;
t248 = Icges(4,4) * t288 - Icges(4,2) * t287 + Icges(4,6) * t301;
t247 = Icges(4,5) * t290 - Icges(4,6) * t289 + Icges(4,3) * t303;
t246 = Icges(4,5) * t288 - Icges(4,6) * t287 + Icges(4,3) * t301;
t245 = qJD(1) + (t274 * t314 + t275 * t316) * t329;
t244 = rSges(5,1) * t267 + rSges(5,2) * t266 + rSges(5,3) * t289;
t243 = rSges(5,1) * t265 + rSges(5,2) * t264 + rSges(5,3) * t287;
t242 = Icges(5,1) * t267 + Icges(5,4) * t266 + Icges(5,5) * t289;
t241 = Icges(5,1) * t265 + Icges(5,4) * t264 + Icges(5,5) * t287;
t240 = Icges(5,4) * t267 + Icges(5,2) * t266 + Icges(5,6) * t289;
t239 = Icges(5,4) * t265 + Icges(5,2) * t264 + Icges(5,6) * t287;
t238 = Icges(5,5) * t267 + Icges(5,6) * t266 + Icges(5,3) * t289;
t237 = Icges(5,5) * t265 + Icges(5,6) * t264 + Icges(5,3) * t287;
t236 = -t252 * t308 + t281 * t294 + t323;
t235 = t253 * t308 - t281 * t293 + t324;
t234 = t252 * t293 - t253 * t294 + t325;
t233 = -t243 * t286 + t257 * t263 - t260 * t308 + t285 * t294 + t323;
t232 = t244 * t286 - t257 * t262 + t261 * t308 - t285 * t293 + t324;
t231 = t243 * t262 - t244 * t263 + t260 * t293 - t261 * t294 + t325;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t245 ^ 2 + t258 ^ 2 + t259 ^ 2) / 0.2e1 - t337 * ((-t269 * t332 - t271 * t301 + t273 * t302) * t335 - (-t268 * t332 - t270 * t301 + t272 * t302) * t332 + (-t295 * t332 - t296 * t301 + t297 * t302) * t317) * t332 / 0.2e1 + m(4) * (t234 ^ 2 + t235 ^ 2 + t236 ^ 2) / 0.2e1 + t293 * ((t303 * t247 - t289 * t249 + t290 * t251) * t293 + (t246 * t303 - t248 * t289 + t250 * t290) * t294 + (t278 * t303 - t279 * t289 + t280 * t290) * t308) / 0.2e1 + t294 * ((t247 * t301 - t249 * t287 + t251 * t288) * t293 + (t301 * t246 - t287 * t248 + t288 * t250) * t294 + (t278 * t301 - t279 * t287 + t280 * t288) * t308) / 0.2e1 + t308 * ((-t247 * t333 - t305 * t249 + t306 * t251) * t293 + (-t246 * t333 - t305 * t248 + t306 * t250) * t294 + (-t278 * t333 - t305 * t279 + t306 * t280) * t308) / 0.2e1 + m(5) * (t231 ^ 2 + t232 ^ 2 + t233 ^ 2) / 0.2e1 + t262 * ((t289 * t238 + t266 * t240 + t267 * t242) * t262 + (t237 * t289 + t239 * t266 + t241 * t267) * t263 + (t254 * t289 + t255 * t266 + t256 * t267) * t286) / 0.2e1 + t263 * ((t238 * t287 + t240 * t264 + t242 * t265) * t262 + (t287 * t237 + t264 * t239 + t265 * t241) * t263 + (t254 * t287 + t255 * t264 + t256 * t265) * t286) / 0.2e1 + t286 * ((t238 * t305 + t240 * t291 + t242 * t292) * t262 + (t237 * t305 + t239 * t291 + t241 * t292) * t263 + (t305 * t254 + t291 * t255 + t292 * t256) * t286) / 0.2e1 + (((t269 * t335 - t271 * t303 + t273 * t304) * t335 - (t268 * t335 - t270 * t303 + t272 * t304) * t332 + (t295 * t335 - t296 * t303 + t297 * t304) * t317) * t335 + t317 * (t317 ^ 2 * t295 + (((t271 * t322 + t273 * t320) * t314 - (t270 * t322 + t272 * t320) * t316) * t315 + (-t268 * t316 + t269 * t314 + t296 * t322 + t297 * t320) * t317) * t315)) * t337 / 0.2e1;
T = t1;
