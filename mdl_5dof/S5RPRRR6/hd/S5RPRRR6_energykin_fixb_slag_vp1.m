% Calculate kinetic energy for
% S5RPRRR6
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
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR6_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR6_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR6_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR6_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:00:55
% EndTime: 2019-12-31 19:00:56
% DurationCPUTime: 1.07s
% Computational Cost: add. (1132->192), mult. (991->317), div. (0->0), fcn. (916->10), ass. (0->108)
t312 = sin(qJ(1));
t352 = pkin(1) * t312;
t314 = cos(qJ(3));
t350 = pkin(3) * t314;
t311 = sin(qJ(3));
t348 = Icges(4,4) * t311;
t347 = Icges(4,4) * t314;
t309 = qJ(3) + qJ(4);
t306 = sin(t309);
t346 = Icges(5,4) * t306;
t307 = cos(t309);
t345 = Icges(5,4) * t307;
t308 = qJ(1) + pkin(9);
t304 = sin(t308);
t344 = t304 * t306;
t305 = cos(t308);
t343 = t305 * t306;
t310 = sin(qJ(5));
t342 = t307 * t310;
t313 = cos(qJ(5));
t341 = t307 * t313;
t315 = cos(qJ(1));
t303 = qJD(1) * t315 * pkin(1);
t340 = qJD(1) * (pkin(2) * t305 + pkin(6) * t304) + t303;
t301 = qJD(3) * t304;
t285 = qJD(4) * t304 + t301;
t339 = qJD(3) * t305;
t338 = qJD(5) * t306;
t337 = pkin(3) * qJD(3) * t311;
t252 = -pkin(7) * t305 + t350 * t304;
t253 = pkin(7) * t304 + t350 * t305;
t336 = t252 * t301 + t253 * t339 + qJD(2);
t335 = -pkin(2) * t304 + pkin(6) * t305 - t352;
t286 = (-qJD(3) - qJD(4)) * t305;
t334 = t305 * t337;
t333 = -t252 + t335;
t332 = pkin(4) * t307 + pkin(8) * t306;
t331 = rSges(4,1) * t314 - rSges(4,2) * t311;
t330 = rSges(5,1) * t307 - rSges(5,2) * t306;
t329 = Icges(4,1) * t314 - t348;
t328 = Icges(5,1) * t307 - t346;
t327 = -Icges(4,2) * t311 + t347;
t326 = -Icges(5,2) * t306 + t345;
t325 = Icges(4,5) * t314 - Icges(4,6) * t311;
t324 = Icges(5,5) * t307 - Icges(5,6) * t306;
t264 = -Icges(4,6) * t305 + t327 * t304;
t266 = -Icges(4,5) * t305 + t329 * t304;
t323 = t264 * t311 - t266 * t314;
t265 = Icges(4,6) * t304 + t327 * t305;
t267 = Icges(4,5) * t304 + t329 * t305;
t322 = -t265 * t311 + t267 * t314;
t294 = Icges(4,2) * t314 + t348;
t295 = Icges(4,1) * t311 + t347;
t321 = -t294 * t311 + t295 * t314;
t320 = qJD(1) * t253 - t304 * t337 + t340;
t319 = (Icges(5,5) * t306 + Icges(5,6) * t307) * qJD(1) + (-Icges(5,3) * t305 + t324 * t304) * t286 + (Icges(5,3) * t304 + t324 * t305) * t285;
t256 = -Icges(5,6) * t305 + t326 * t304;
t257 = Icges(5,6) * t304 + t326 * t305;
t258 = -Icges(5,5) * t305 + t328 * t304;
t259 = Icges(5,5) * t304 + t328 * t305;
t289 = Icges(5,2) * t307 + t346;
t290 = Icges(5,1) * t306 + t345;
t318 = (-t257 * t306 + t259 * t307) * t285 + (-t256 * t306 + t258 * t307) * t286 + (-t289 * t306 + t290 * t307) * qJD(1);
t299 = -qJD(5) * t307 + qJD(1);
t298 = rSges(2,1) * t315 - rSges(2,2) * t312;
t297 = rSges(2,1) * t312 + rSges(2,2) * t315;
t296 = rSges(4,1) * t311 + rSges(4,2) * t314;
t293 = Icges(4,5) * t311 + Icges(4,6) * t314;
t292 = pkin(4) * t306 - pkin(8) * t307;
t291 = rSges(5,1) * t306 + rSges(5,2) * t307;
t283 = t304 * t310 + t305 * t341;
t282 = t304 * t313 - t305 * t342;
t281 = t304 * t341 - t305 * t310;
t280 = -t304 * t342 - t305 * t313;
t279 = t303 + qJD(1) * (rSges(3,1) * t305 - rSges(3,2) * t304);
t278 = (-rSges(3,1) * t304 - rSges(3,2) * t305 - t352) * qJD(1);
t277 = t332 * t305;
t276 = t332 * t304;
t275 = -rSges(6,3) * t307 + (rSges(6,1) * t313 - rSges(6,2) * t310) * t306;
t274 = -Icges(6,5) * t307 + (Icges(6,1) * t313 - Icges(6,4) * t310) * t306;
t273 = -Icges(6,6) * t307 + (Icges(6,4) * t313 - Icges(6,2) * t310) * t306;
t272 = -Icges(6,3) * t307 + (Icges(6,5) * t313 - Icges(6,6) * t310) * t306;
t271 = rSges(4,3) * t304 + t331 * t305;
t270 = -rSges(4,3) * t305 + t331 * t304;
t269 = t304 * t338 + t286;
t268 = t305 * t338 + t285;
t263 = Icges(4,3) * t304 + t325 * t305;
t262 = -Icges(4,3) * t305 + t325 * t304;
t261 = rSges(5,3) * t304 + t330 * t305;
t260 = -rSges(5,3) * t305 + t330 * t304;
t248 = rSges(6,1) * t283 + rSges(6,2) * t282 + rSges(6,3) * t343;
t247 = rSges(6,1) * t281 + rSges(6,2) * t280 + rSges(6,3) * t344;
t246 = Icges(6,1) * t283 + Icges(6,4) * t282 + Icges(6,5) * t343;
t245 = Icges(6,1) * t281 + Icges(6,4) * t280 + Icges(6,5) * t344;
t244 = Icges(6,4) * t283 + Icges(6,2) * t282 + Icges(6,6) * t343;
t243 = Icges(6,4) * t281 + Icges(6,2) * t280 + Icges(6,6) * t344;
t242 = Icges(6,5) * t283 + Icges(6,6) * t282 + Icges(6,3) * t343;
t241 = Icges(6,5) * t281 + Icges(6,6) * t280 + Icges(6,3) * t344;
t240 = qJD(1) * t271 - t296 * t301 + t340;
t239 = -t296 * t339 + (-t270 + t335) * qJD(1);
t238 = qJD(2) + (t270 * t304 + t271 * t305) * qJD(3);
t237 = qJD(1) * t261 - t285 * t291 + t320;
t236 = -t334 + t286 * t291 + (-t260 + t333) * qJD(1);
t235 = t260 * t285 - t261 * t286 + t336;
t234 = qJD(1) * t277 + t248 * t299 - t268 * t275 - t285 * t292 + t320;
t233 = -t334 - t247 * t299 + t269 * t275 + t286 * t292 + (-t276 + t333) * qJD(1);
t232 = t247 * t268 - t248 * t269 + t276 * t285 - t277 * t286 + t336;
t1 = m(3) * (qJD(2) ^ 2 + t278 ^ 2 + t279 ^ 2) / 0.2e1 + m(4) * (t238 ^ 2 + t239 ^ 2 + t240 ^ 2) / 0.2e1 + ((t304 * t293 + t321 * t305) * qJD(1) + (t304 ^ 2 * t263 + (t323 * t305 + (-t262 + t322) * t304) * t305) * qJD(3)) * t301 / 0.2e1 - ((-t305 * t293 + t321 * t304) * qJD(1) + (t305 ^ 2 * t262 + (t322 * t304 + (-t263 + t323) * t305) * t304) * qJD(3)) * t339 / 0.2e1 + m(5) * (t235 ^ 2 + t236 ^ 2 + t237 ^ 2) / 0.2e1 + t285 * (t319 * t304 + t318 * t305) / 0.2e1 + t286 * (t318 * t304 - t319 * t305) / 0.2e1 + m(6) * (t232 ^ 2 + t233 ^ 2 + t234 ^ 2) / 0.2e1 + t268 * ((t242 * t343 + t282 * t244 + t283 * t246) * t268 + (t241 * t343 + t243 * t282 + t245 * t283) * t269 + (t272 * t343 + t273 * t282 + t274 * t283) * t299) / 0.2e1 + t269 * ((t242 * t344 + t244 * t280 + t246 * t281) * t268 + (t241 * t344 + t280 * t243 + t281 * t245) * t269 + (t272 * t344 + t273 * t280 + t274 * t281) * t299) / 0.2e1 + t299 * ((-t241 * t269 - t242 * t268 - t272 * t299) * t307 + ((-t244 * t310 + t246 * t313) * t268 + (-t243 * t310 + t245 * t313) * t269 + (-t273 * t310 + t274 * t313) * t299) * t306) / 0.2e1 + (((t265 * t314 + t267 * t311) * t304 - (t264 * t314 + t266 * t311) * t305) * qJD(3) + (t257 * t307 + t259 * t306) * t285 + (t256 * t307 + t258 * t306) * t286 + (t307 * t289 + t306 * t290 + t314 * t294 + t311 * t295) * qJD(1)) * qJD(1) / 0.2e1 + (m(2) * (t297 ^ 2 + t298 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
