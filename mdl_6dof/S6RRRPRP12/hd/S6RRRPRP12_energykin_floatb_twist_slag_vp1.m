% Calculate kinetic energy for
% S6RRRPRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP12_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP12_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRP12_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP12_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP12_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP12_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:51:32
% EndTime: 2019-03-09 17:51:36
% DurationCPUTime: 4.02s
% Computational Cost: add. (2544->340), mult. (5831->474), div. (0->0), fcn. (6966->10), ass. (0->157)
t376 = Icges(4,1) + Icges(5,2);
t375 = Icges(5,1) + Icges(4,3);
t374 = Icges(6,1) + Icges(7,1);
t373 = -Icges(4,4) - Icges(5,6);
t372 = Icges(5,4) - Icges(4,5);
t371 = Icges(6,4) - Icges(7,5);
t370 = Icges(7,4) + Icges(6,5);
t369 = Icges(5,5) - Icges(4,6);
t368 = Icges(4,2) + Icges(5,3);
t367 = Icges(6,2) + Icges(7,3);
t366 = Icges(7,2) + Icges(6,3);
t365 = Icges(6,6) - Icges(7,6);
t364 = rSges(7,1) + pkin(5);
t363 = rSges(7,3) + qJ(6);
t300 = cos(pkin(6));
t303 = sin(qJ(1));
t304 = cos(qJ(2));
t328 = t303 * t304;
t302 = sin(qJ(2));
t305 = cos(qJ(1));
t329 = t302 * t305;
t268 = t300 * t329 + t328;
t299 = sin(pkin(6));
t338 = cos(qJ(3));
t319 = t299 * t338;
t336 = sin(qJ(3));
t247 = t268 * t336 + t305 * t319;
t327 = t304 * t305;
t330 = t302 * t303;
t267 = -t300 * t327 + t330;
t301 = sin(qJ(5));
t337 = cos(qJ(5));
t211 = -t247 * t337 + t267 * t301;
t212 = t247 * t301 + t267 * t337;
t318 = t299 * t336;
t248 = t268 * t338 - t305 * t318;
t362 = t367 * t211 - t371 * t212 - t365 * t248;
t270 = -t300 * t330 + t327;
t249 = t270 * t336 - t303 * t319;
t269 = t300 * t328 + t329;
t213 = -t249 * t337 + t269 * t301;
t214 = t249 * t301 + t269 * t337;
t250 = t270 * t338 + t303 * t318;
t361 = t367 * t213 - t371 * t214 - t365 * t250;
t360 = -t365 * t211 + t370 * t212 + t366 * t248;
t359 = -t365 * t213 + t370 * t214 + t366 * t250;
t358 = -t371 * t211 + t374 * t212 + t370 * t248;
t357 = -t371 * t213 + t374 * t214 + t370 * t250;
t265 = -t300 * t338 + t302 * t318;
t332 = t299 * t304;
t243 = t265 * t337 + t301 * t332;
t244 = t265 * t301 - t332 * t337;
t266 = t300 * t336 + t302 * t319;
t356 = -t367 * t243 - t371 * t244 - t365 * t266;
t355 = t365 * t243 + t370 * t244 + t366 * t266;
t354 = t371 * t243 + t374 * t244 + t370 * t266;
t353 = t368 * t247 + t373 * t248 + t369 * t267;
t352 = t368 * t249 + t373 * t250 + t369 * t269;
t351 = t369 * t247 - t372 * t248 + t375 * t267;
t350 = t369 * t249 - t372 * t250 + t375 * t269;
t349 = t373 * t247 + t376 * t248 - t372 * t267;
t348 = t373 * t249 + t376 * t250 - t372 * t269;
t347 = t368 * t265 + t373 * t266 - t369 * t332;
t346 = t373 * t265 + t376 * t266 + t372 * t332;
t345 = t369 * t265 - t372 * t266 - t375 * t332;
t335 = pkin(8) * t300;
t334 = Icges(2,4) * t303;
t333 = t299 * t303;
t331 = t299 * t305;
t326 = rSges(7,2) * t248 + t363 * t211 + t212 * t364;
t325 = rSges(7,2) * t250 + t363 * t213 + t214 * t364;
t324 = rSges(7,2) * t266 - t363 * t243 + t244 * t364;
t323 = qJD(2) * t299;
t322 = V_base(5) * pkin(7) + V_base(1);
t279 = t303 * t323 + V_base(4);
t296 = V_base(6) + qJD(1);
t246 = qJD(3) * t269 + t279;
t280 = qJD(2) * t300 + t296;
t278 = -t305 * t323 + V_base(5);
t273 = t303 * pkin(1) - pkin(8) * t331;
t317 = -t273 * t296 + V_base(5) * t335 + t322;
t274 = pkin(1) * t305 + pkin(8) * t333;
t316 = V_base(4) * t273 - t274 * V_base(5) + V_base(3);
t245 = qJD(3) * t267 + t278;
t263 = -qJD(3) * t332 + t280;
t315 = t296 * t274 + V_base(2) + (-pkin(7) - t335) * V_base(4);
t237 = pkin(2) * t268 + pkin(9) * t267;
t272 = (pkin(2) * t302 - pkin(9) * t304) * t299;
t314 = -t237 * t280 + t278 * t272 + t317;
t238 = pkin(2) * t270 + pkin(9) * t269;
t313 = t279 * t237 - t238 * t278 + t316;
t235 = pkin(3) * t266 + qJ(4) * t265;
t312 = qJD(4) * t249 + t245 * t235 + t314;
t206 = pkin(3) * t248 + qJ(4) * t247;
t311 = qJD(4) * t265 + t246 * t206 + t313;
t310 = t280 * t238 - t272 * t279 + t315;
t207 = pkin(3) * t250 + qJ(4) * t249;
t309 = qJD(4) * t247 + t263 * t207 + t310;
t215 = pkin(4) * t267 + pkin(10) * t248;
t255 = -pkin(4) * t332 + pkin(10) * t266;
t308 = t245 * t255 + (-t206 - t215) * t263 + t312;
t216 = pkin(4) * t269 + pkin(10) * t250;
t307 = t246 * t215 + (-t207 - t216) * t245 + t311;
t306 = t263 * t216 + (-t235 - t255) * t246 + t309;
t297 = Icges(2,4) * t305;
t288 = rSges(2,1) * t305 - t303 * rSges(2,2);
t287 = t303 * rSges(2,1) + rSges(2,2) * t305;
t286 = Icges(2,1) * t305 - t334;
t285 = Icges(2,1) * t303 + t297;
t284 = -Icges(2,2) * t303 + t297;
t283 = Icges(2,2) * t305 + t334;
t277 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t276 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t275 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t259 = rSges(3,3) * t300 + (rSges(3,1) * t302 + rSges(3,2) * t304) * t299;
t258 = Icges(3,5) * t300 + (Icges(3,1) * t302 + Icges(3,4) * t304) * t299;
t257 = Icges(3,6) * t300 + (Icges(3,4) * t302 + Icges(3,2) * t304) * t299;
t256 = Icges(3,3) * t300 + (Icges(3,5) * t302 + Icges(3,6) * t304) * t299;
t254 = V_base(5) * rSges(2,3) - t287 * t296 + t322;
t253 = t288 * t296 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t251 = t287 * V_base(4) - t288 * V_base(5) + V_base(3);
t236 = qJD(5) * t266 + t263;
t234 = rSges(3,1) * t270 - rSges(3,2) * t269 + rSges(3,3) * t333;
t233 = t268 * rSges(3,1) - t267 * rSges(3,2) - rSges(3,3) * t331;
t232 = Icges(3,1) * t270 - Icges(3,4) * t269 + Icges(3,5) * t333;
t231 = Icges(3,1) * t268 - Icges(3,4) * t267 - Icges(3,5) * t331;
t230 = Icges(3,4) * t270 - Icges(3,2) * t269 + Icges(3,6) * t333;
t229 = Icges(3,4) * t268 - Icges(3,2) * t267 - Icges(3,6) * t331;
t228 = Icges(3,5) * t270 - Icges(3,6) * t269 + Icges(3,3) * t333;
t227 = Icges(3,5) * t268 - Icges(3,6) * t267 - Icges(3,3) * t331;
t226 = rSges(4,1) * t266 - rSges(4,2) * t265 - rSges(4,3) * t332;
t225 = -rSges(5,1) * t332 - rSges(5,2) * t266 + rSges(5,3) * t265;
t209 = qJD(5) * t250 + t246;
t208 = qJD(5) * t248 + t245;
t201 = rSges(4,1) * t250 - rSges(4,2) * t249 + rSges(4,3) * t269;
t200 = rSges(4,1) * t248 - rSges(4,2) * t247 + rSges(4,3) * t267;
t199 = rSges(5,1) * t269 - rSges(5,2) * t250 + rSges(5,3) * t249;
t198 = rSges(5,1) * t267 - rSges(5,2) * t248 + rSges(5,3) * t247;
t184 = rSges(6,1) * t244 + rSges(6,2) * t243 + rSges(6,3) * t266;
t173 = -t233 * t280 + t259 * t278 + t317;
t172 = t234 * t280 - t259 * t279 + t315;
t171 = rSges(6,1) * t214 - rSges(6,2) * t213 + rSges(6,3) * t250;
t169 = rSges(6,1) * t212 - rSges(6,2) * t211 + rSges(6,3) * t248;
t155 = t233 * t279 - t234 * t278 + t316;
t154 = -t200 * t263 + t226 * t245 + t314;
t153 = t201 * t263 - t226 * t246 + t310;
t152 = t200 * t246 - t201 * t245 + t313;
t151 = t225 * t245 + (-t198 - t206) * t263 + t312;
t150 = t199 * t263 + (-t225 - t235) * t246 + t309;
t149 = t198 * t246 + (-t199 - t207) * t245 + t311;
t148 = -t169 * t236 + t184 * t208 + t308;
t147 = t171 * t236 - t184 * t209 + t306;
t146 = t169 * t209 - t171 * t208 + t307;
t145 = qJD(6) * t213 + t208 * t324 - t236 * t326 + t308;
t144 = qJD(6) * t211 - t209 * t324 + t236 * t325 + t306;
t143 = -qJD(6) * t243 - t208 * t325 + t209 * t326 + t307;
t1 = m(5) * (t149 ^ 2 + t150 ^ 2 + t151 ^ 2) / 0.2e1 + m(4) * (t152 ^ 2 + t153 ^ 2 + t154 ^ 2) / 0.2e1 + m(6) * (t146 ^ 2 + t147 ^ 2 + t148 ^ 2) / 0.2e1 + m(7) * (t143 ^ 2 + t144 ^ 2 + t145 ^ 2) / 0.2e1 + t280 * ((t227 * t278 + t228 * t279 + t256 * t280) * t300 + ((t230 * t304 + t232 * t302) * t279 + (t229 * t304 + t231 * t302) * t278 + (t257 * t304 + t258 * t302) * t280) * t299) / 0.2e1 + t279 * ((t228 * t333 - t269 * t230 + t270 * t232) * t279 + (t227 * t333 - t229 * t269 + t231 * t270) * t278 + (t256 * t333 - t257 * t269 + t258 * t270) * t280) / 0.2e1 + t278 * ((-t228 * t331 - t267 * t230 + t268 * t232) * t279 + (-t227 * t331 - t267 * t229 + t268 * t231) * t278 + (-t256 * t331 - t267 * t257 + t268 * t258) * t280) / 0.2e1 + m(1) * (t275 ^ 2 + t276 ^ 2 + t277 ^ 2) / 0.2e1 + m(3) * (t155 ^ 2 + t172 ^ 2 + t173 ^ 2) / 0.2e1 + m(2) * (t251 ^ 2 + t253 ^ 2 + t254 ^ 2) / 0.2e1 + ((t211 * t356 + t212 * t354 + t248 * t355) * t236 + (t211 * t361 + t212 * t357 + t248 * t359) * t209 + (t362 * t211 + t358 * t212 + t360 * t248) * t208) * t208 / 0.2e1 + ((t213 * t356 + t214 * t354 + t250 * t355) * t236 + (t361 * t213 + t357 * t214 + t359 * t250) * t209 + (t213 * t362 + t358 * t214 + t360 * t250) * t208) * t209 / 0.2e1 + ((-t356 * t243 + t354 * t244 + t355 * t266) * t236 + (-t243 * t361 + t244 * t357 + t266 * t359) * t209 + (-t243 * t362 + t358 * t244 + t360 * t266) * t208) * t236 / 0.2e1 + ((t247 * t347 + t248 * t346 + t267 * t345) * t263 + (t247 * t352 + t248 * t348 + t267 * t350) * t246 + (t353 * t247 + t349 * t248 + t351 * t267) * t245) * t245 / 0.2e1 + ((t249 * t347 + t250 * t346 + t269 * t345) * t263 + (t352 * t249 + t348 * t250 + t350 * t269) * t246 + (t249 * t353 + t250 * t349 + t269 * t351) * t245) * t246 / 0.2e1 + ((t347 * t265 + t346 * t266 - t345 * t332) * t263 + (t265 * t352 + t266 * t348 - t332 * t350) * t246 + (t265 * t353 + t266 * t349 - t332 * t351) * t245) * t263 / 0.2e1 + ((-t303 * t283 + t285 * t305 + Icges(1,4)) * V_base(5) + (-t303 * t284 + t305 * t286 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t305 * t283 + t303 * t285 + Icges(1,2)) * V_base(5) + (t284 * t305 + t303 * t286 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t303 + Icges(2,6) * t305) * V_base(5) + (Icges(2,5) * t305 - Icges(2,6) * t303) * V_base(4) + Icges(2,3) * t296 / 0.2e1) * t296;
T  = t1;
