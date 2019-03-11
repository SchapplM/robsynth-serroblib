% Calculate kinetic energy for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPPRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:29:50
% EndTime: 2019-03-09 01:29:50
% DurationCPUTime: 0.62s
% Computational Cost: add. (657->167), mult. (777->265), div. (0->0), fcn. (706->8), ass. (0->85)
t304 = qJ(1) + pkin(9);
t302 = sin(t304);
t336 = pkin(7) * t302;
t307 = sin(qJ(1));
t335 = t307 * pkin(1);
t306 = sin(qJ(5));
t334 = Icges(6,4) * t306;
t309 = cos(qJ(5));
t333 = Icges(6,4) * t309;
t332 = t302 * t309;
t303 = cos(t304);
t331 = t303 * t309;
t305 = sin(qJ(6));
t330 = t305 * t306;
t308 = cos(qJ(6));
t329 = t306 * t308;
t299 = qJD(3) * t302;
t328 = qJD(4) * t303 + t299;
t327 = qJD(5) * t302;
t326 = qJD(5) * t303;
t325 = qJD(6) * t309;
t324 = -pkin(2) * t302 + qJ(3) * t303 - t335;
t323 = pkin(5) * t306 - pkin(8) * t309;
t310 = cos(qJ(1));
t301 = qJD(1) * t310 * pkin(1);
t322 = -qJD(3) * t303 + qJD(1) * (pkin(2) * t303 + qJ(3) * t302) + t301;
t321 = rSges(6,1) * t306 + rSges(6,2) * t309;
t320 = Icges(6,1) * t306 + t333;
t319 = Icges(6,2) * t309 + t334;
t318 = Icges(6,5) * t306 + Icges(6,6) * t309;
t267 = Icges(6,6) * t303 + t319 * t302;
t269 = Icges(6,5) * t303 + t320 * t302;
t317 = t267 * t309 + t269 * t306;
t268 = -Icges(6,6) * t302 + t319 * t303;
t270 = -Icges(6,5) * t302 + t320 * t303;
t316 = -t268 * t309 - t270 * t306;
t290 = -Icges(6,2) * t306 + t333;
t291 = Icges(6,1) * t309 - t334;
t315 = t290 * t309 + t291 * t306;
t314 = qJD(1) * t303 * qJ(4) + qJD(4) * t302 + t322;
t313 = -pkin(7) * t303 - qJ(4) * t302 + t324;
t311 = qJD(2) ^ 2;
t300 = qJD(6) * t306 + qJD(1);
t295 = pkin(5) * t309 + pkin(8) * t306;
t294 = rSges(2,1) * t310 - rSges(2,2) * t307;
t293 = rSges(6,1) * t309 - rSges(6,2) * t306;
t292 = rSges(2,1) * t307 + rSges(2,2) * t310;
t289 = Icges(6,5) * t309 - Icges(6,6) * t306;
t286 = -t302 * t325 + t326;
t285 = -t303 * t325 - t327;
t284 = t323 * t303;
t283 = t323 * t302;
t282 = rSges(7,3) * t306 + (rSges(7,1) * t308 - rSges(7,2) * t305) * t309;
t281 = Icges(7,5) * t306 + (Icges(7,1) * t308 - Icges(7,4) * t305) * t309;
t280 = Icges(7,6) * t306 + (Icges(7,4) * t308 - Icges(7,2) * t305) * t309;
t279 = Icges(7,3) * t306 + (Icges(7,5) * t308 - Icges(7,6) * t305) * t309;
t278 = -t302 * t305 + t303 * t329;
t277 = -t302 * t308 - t303 * t330;
t276 = t302 * t329 + t303 * t305;
t275 = -t302 * t330 + t303 * t308;
t274 = t301 + qJD(1) * (rSges(3,1) * t303 - rSges(3,2) * t302);
t273 = (-rSges(3,1) * t302 - rSges(3,2) * t303 - t335) * qJD(1);
t272 = -rSges(6,3) * t302 + t321 * t303;
t271 = rSges(6,3) * t303 + t321 * t302;
t266 = -Icges(6,3) * t302 + t318 * t303;
t265 = Icges(6,3) * t303 + t318 * t302;
t264 = qJD(1) * (-rSges(4,2) * t303 + rSges(4,3) * t302) + t322;
t263 = t299 + (rSges(4,2) * t302 + rSges(4,3) * t303 + t324) * qJD(1);
t262 = rSges(7,1) * t278 + rSges(7,2) * t277 - rSges(7,3) * t331;
t261 = rSges(7,1) * t276 + rSges(7,2) * t275 - rSges(7,3) * t332;
t260 = Icges(7,1) * t278 + Icges(7,4) * t277 - Icges(7,5) * t331;
t259 = Icges(7,1) * t276 + Icges(7,4) * t275 - Icges(7,5) * t332;
t258 = Icges(7,4) * t278 + Icges(7,2) * t277 - Icges(7,6) * t331;
t257 = Icges(7,4) * t276 + Icges(7,2) * t275 - Icges(7,6) * t332;
t256 = Icges(7,5) * t278 + Icges(7,6) * t277 - Icges(7,3) * t331;
t255 = Icges(7,5) * t276 + Icges(7,6) * t275 - Icges(7,3) * t332;
t254 = qJD(1) * (rSges(5,2) * t302 + rSges(5,3) * t303) + t314;
t253 = (t303 * rSges(5,2) + (-rSges(5,3) - qJ(4)) * t302 + t324) * qJD(1) + t328;
t252 = qJD(2) + (-t271 * t302 - t272 * t303) * qJD(5);
t251 = t293 * t327 + (t272 - t336) * qJD(1) + t314;
t250 = t293 * t326 + (-t271 + t313) * qJD(1) + t328;
t249 = t295 * t327 + t262 * t300 - t282 * t285 + (t284 - t336) * qJD(1) + t314;
t248 = t295 * t326 - t261 * t300 + t282 * t286 + (-t283 + t313) * qJD(1) + t328;
t247 = t261 * t285 - t262 * t286 + qJD(2) + (-t283 * t302 - t284 * t303) * qJD(5);
t1 = m(6) * (t250 ^ 2 + t251 ^ 2 + t252 ^ 2) / 0.2e1 + m(7) * (t247 ^ 2 + t248 ^ 2 + t249 ^ 2) / 0.2e1 + m(3) * (t273 ^ 2 + t274 ^ 2 + t311) / 0.2e1 + m(4) * (t263 ^ 2 + t264 ^ 2 + t311) / 0.2e1 + m(5) * (t253 ^ 2 + t254 ^ 2 + t311) / 0.2e1 + t285 * ((-t256 * t331 + t277 * t258 + t278 * t260) * t285 + (-t255 * t331 + t257 * t277 + t259 * t278) * t286 + (t277 * t280 + t278 * t281 - t279 * t331) * t300) / 0.2e1 + t286 * ((-t256 * t332 + t258 * t275 + t260 * t276) * t285 + (-t255 * t332 + t275 * t257 + t276 * t259) * t286 + (t275 * t280 + t276 * t281 - t279 * t332) * t300) / 0.2e1 + t300 * ((t255 * t286 + t256 * t285 + t279 * t300) * t306 + ((-t258 * t305 + t260 * t308) * t285 + (-t257 * t305 + t259 * t308) * t286 + (-t280 * t305 + t281 * t308) * t300) * t309) / 0.2e1 - ((-t302 * t289 + t315 * t303) * qJD(1) + (t302 ^ 2 * t266 + (t317 * t303 + (-t265 + t316) * t302) * t303) * qJD(5)) * t327 / 0.2e1 + ((t303 * t289 + t315 * t302) * qJD(1) + (t303 ^ 2 * t265 + (t316 * t302 + (-t266 + t317) * t303) * t302) * qJD(5)) * t326 / 0.2e1 + qJD(1) * ((-t306 * t290 + t309 * t291) * qJD(1) + (-(-t268 * t306 + t270 * t309) * t302 + (-t267 * t306 + t269 * t309) * t303) * qJD(5)) / 0.2e1 + (m(2) * (t292 ^ 2 + t294 ^ 2) + Icges(4,1) + Icges(5,1) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
