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
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:35:40
% EndTime: 2019-12-05 18:35:41
% DurationCPUTime: 0.87s
% Computational Cost: add. (981->171), mult. (1025->294), div. (0->0), fcn. (1032->10), ass. (0->89)
t312 = sin(pkin(9));
t311 = qJ(1) + qJ(2);
t305 = sin(t311);
t307 = cos(t311);
t316 = cos(qJ(4));
t313 = cos(pkin(9));
t314 = sin(qJ(4));
t334 = t313 * t314;
t290 = t305 * t334 + t307 * t316;
t333 = t313 * t316;
t335 = t307 * t314;
t291 = -t305 * t333 + t335;
t292 = t305 * t316 - t307 * t334;
t338 = t305 * t314;
t293 = t307 * t333 + t338;
t337 = t307 * t312;
t340 = t305 * t312;
t323 = (Icges(5,5) * t291 + Icges(5,6) * t290 - Icges(5,3) * t340) * t305 - (Icges(5,5) * t293 + Icges(5,6) * t292 + Icges(5,3) * t337) * t307;
t347 = t312 * t323;
t343 = pkin(4) * t316;
t346 = pkin(8) * t312 + t313 * t343;
t341 = pkin(1) * qJD(1);
t339 = t305 * t313;
t336 = t307 * t313;
t332 = qJD(4) * t312;
t309 = qJD(1) + qJD(2);
t331 = qJD(4) + qJD(5);
t315 = sin(qJ(1));
t330 = t315 * t341;
t317 = cos(qJ(1));
t329 = t317 * t341;
t328 = t305 * t332;
t327 = t307 * t332;
t326 = qJD(3) * t307 - t329;
t325 = pkin(3) * t313 + pkin(7) * t312;
t324 = -rSges(4,1) * t313 + rSges(4,2) * t312;
t322 = t309 * (-pkin(2) * t305 + qJ(3) * t307) + qJD(3) * t305 - t330;
t321 = -t309 * t325 * t305 + t322;
t297 = pkin(2) * t307 + qJ(3) * t305;
t320 = (-t325 * t307 - t297) * t309 + t326;
t310 = qJ(4) + qJ(5);
t306 = cos(t310);
t304 = sin(t310);
t300 = -t313 * qJD(4) + t309;
t299 = rSges(2,1) * t317 - rSges(2,2) * t315;
t298 = -rSges(2,1) * t315 - rSges(2,2) * t317;
t296 = -t313 * t331 + t309;
t289 = t331 * t337;
t288 = t331 * t340;
t287 = -rSges(5,3) * t313 + (rSges(5,1) * t316 - rSges(5,2) * t314) * t312;
t286 = -Icges(5,5) * t313 + (Icges(5,1) * t316 - Icges(5,4) * t314) * t312;
t285 = -Icges(5,6) * t313 + (Icges(5,4) * t316 - Icges(5,2) * t314) * t312;
t284 = -Icges(5,3) * t313 + (Icges(5,5) * t316 - Icges(5,6) * t314) * t312;
t282 = t304 * t305 + t306 * t336;
t281 = -t304 * t336 + t305 * t306;
t280 = t304 * t307 - t306 * t339;
t279 = t304 * t339 + t306 * t307;
t278 = -t329 - t309 * (rSges(3,1) * t307 - rSges(3,2) * t305);
t277 = -t330 + t309 * (-rSges(3,1) * t305 - rSges(3,2) * t307);
t276 = -rSges(6,3) * t313 + (rSges(6,1) * t306 - rSges(6,2) * t304) * t312;
t275 = -Icges(6,5) * t313 + (Icges(6,1) * t306 - Icges(6,4) * t304) * t312;
t274 = -Icges(6,6) * t313 + (Icges(6,4) * t306 - Icges(6,2) * t304) * t312;
t273 = -Icges(6,3) * t313 + (Icges(6,5) * t306 - Icges(6,6) * t304) * t312;
t272 = -pkin(8) * t313 + t312 * t343;
t271 = rSges(5,1) * t293 + rSges(5,2) * t292 + rSges(5,3) * t337;
t270 = rSges(5,1) * t291 + rSges(5,2) * t290 - rSges(5,3) * t340;
t269 = pkin(4) * t338 + t346 * t307;
t268 = pkin(4) * t335 - t346 * t305;
t267 = Icges(5,1) * t293 + Icges(5,4) * t292 + Icges(5,5) * t337;
t266 = Icges(5,1) * t291 + Icges(5,4) * t290 - Icges(5,5) * t340;
t265 = Icges(5,4) * t293 + Icges(5,2) * t292 + Icges(5,6) * t337;
t264 = Icges(5,4) * t291 + Icges(5,2) * t290 - Icges(5,6) * t340;
t261 = rSges(6,1) * t282 + rSges(6,2) * t281 + rSges(6,3) * t337;
t260 = rSges(6,1) * t280 + rSges(6,2) * t279 - rSges(6,3) * t340;
t259 = (-rSges(4,3) * t305 + t307 * t324 - t297) * t309 + t326;
t258 = t309 * (rSges(4,3) * t307 + t305 * t324) + t322;
t257 = Icges(6,1) * t282 + Icges(6,4) * t281 + Icges(6,5) * t337;
t256 = Icges(6,1) * t280 + Icges(6,4) * t279 - Icges(6,5) * t340;
t255 = Icges(6,4) * t282 + Icges(6,2) * t281 + Icges(6,6) * t337;
t254 = Icges(6,4) * t280 + Icges(6,2) * t279 - Icges(6,6) * t340;
t253 = Icges(6,5) * t282 + Icges(6,6) * t281 + Icges(6,3) * t337;
t252 = Icges(6,5) * t280 + Icges(6,6) * t279 - Icges(6,3) * t340;
t251 = (-t270 * t307 - t271 * t305) * t332;
t250 = -t271 * t300 + t287 * t327 + t320;
t249 = t270 * t300 + t287 * t328 + t321;
t248 = -t261 * t296 - t269 * t300 + t272 * t327 + t276 * t289 + t320;
t247 = t260 * t296 + t268 * t300 + t272 * t328 + t276 * t288 + t321;
t246 = -t260 * t289 - t261 * t288 + (-t268 * t307 - t269 * t305) * t332;
t1 = m(3) * (t277 ^ 2 + t278 ^ 2) / 0.2e1 + m(4) * (t258 ^ 2 + t259 ^ 2) / 0.2e1 + m(5) * (t249 ^ 2 + t250 ^ 2 + t251 ^ 2) / 0.2e1 + t300 * ((-t284 * t313 + (-t285 * t314 + t286 * t316) * t312) * t300 + ((-(-t264 * t314 + t266 * t316) * t305 + (-t265 * t314 + t267 * t316) * t307) * t312 + t323 * t313) * t332) / 0.2e1 - ((-t284 * t340 + t285 * t290 + t286 * t291) * t300 + ((t265 * t290 + t267 * t291) * t307 + (-t290 * t264 - t291 * t266 + t347) * t305) * t332) * t328 / 0.2e1 + ((t284 * t337 + t285 * t292 + t286 * t293) * t300 + (-(t264 * t292 + t293 * t266) * t305 + (t265 * t292 + t267 * t293 - t347) * t307) * t332) * t327 / 0.2e1 + m(6) * (t246 ^ 2 + t247 ^ 2 + t248 ^ 2) / 0.2e1 + t296 * ((t252 * t288 - t253 * t289 - t273 * t296) * t313 + ((-t274 * t304 + t275 * t306) * t296 - (-t254 * t304 + t256 * t306) * t288 + (-t255 * t304 + t257 * t306) * t289) * t312) / 0.2e1 - t288 * ((-t273 * t340 + t274 * t279 + t275 * t280) * t296 - (-t252 * t340 + t254 * t279 + t256 * t280) * t288 + (-t253 * t340 + t255 * t279 + t257 * t280) * t289) / 0.2e1 + t289 * ((t273 * t337 + t274 * t281 + t275 * t282) * t296 - (t252 * t337 + t254 * t281 + t256 * t282) * t288 + (t253 * t337 + t255 * t281 + t257 * t282) * t289) / 0.2e1 + (Icges(3,3) + Icges(4,2) * t313 ^ 2 + (Icges(4,1) * t312 + 0.2e1 * Icges(4,4) * t313) * t312) * t309 ^ 2 / 0.2e1 + (m(2) * (t298 ^ 2 + t299 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
