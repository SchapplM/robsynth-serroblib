% Calculate kinetic energy for
% S6PRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRP5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRPRP5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRP5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRP5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:45:22
% EndTime: 2019-03-08 21:45:26
% DurationCPUTime: 4.09s
% Computational Cost: add. (2484->340), mult. (5831->472), div. (0->0), fcn. (6966->10), ass. (0->155)
t377 = Icges(4,1) + Icges(5,2);
t376 = Icges(5,1) + Icges(4,3);
t375 = Icges(6,1) + Icges(7,1);
t374 = -Icges(4,4) - Icges(5,6);
t373 = Icges(5,4) - Icges(4,5);
t372 = Icges(6,4) - Icges(7,5);
t371 = Icges(7,4) + Icges(6,5);
t370 = Icges(5,5) - Icges(4,6);
t369 = Icges(4,2) + Icges(5,3);
t368 = Icges(6,2) + Icges(7,3);
t367 = Icges(7,2) + Icges(6,3);
t366 = Icges(6,6) - Icges(7,6);
t365 = rSges(7,1) + pkin(5);
t364 = rSges(7,3) + qJ(6);
t298 = sin(pkin(10));
t300 = cos(pkin(10));
t304 = cos(qJ(2));
t301 = cos(pkin(6));
t303 = sin(qJ(2));
t329 = t301 * t303;
t264 = t298 * t304 + t300 * t329;
t299 = sin(pkin(6));
t337 = cos(qJ(3));
t318 = t299 * t337;
t335 = sin(qJ(3));
t246 = t264 * t335 + t300 * t318;
t328 = t301 * t304;
t263 = t298 * t303 - t300 * t328;
t302 = sin(qJ(5));
t336 = cos(qJ(5));
t211 = -t246 * t336 + t263 * t302;
t212 = t246 * t302 + t263 * t336;
t317 = t299 * t335;
t247 = t264 * t337 - t300 * t317;
t363 = t368 * t211 - t372 * t212 - t366 * t247;
t266 = -t298 * t329 + t300 * t304;
t248 = t266 * t335 - t298 * t318;
t265 = t298 * t328 + t300 * t303;
t213 = -t248 * t336 + t265 * t302;
t214 = t248 * t302 + t265 * t336;
t249 = t266 * t337 + t298 * t317;
t362 = t368 * t213 - t372 * t214 - t366 * t249;
t361 = -t366 * t211 + t371 * t212 + t367 * t247;
t360 = -t366 * t213 + t371 * t214 + t367 * t249;
t359 = -t372 * t211 + t375 * t212 + t371 * t247;
t358 = -t372 * t213 + t375 * t214 + t371 * t249;
t357 = t369 * t246 + t374 * t247 + t370 * t263;
t356 = t369 * t248 + t374 * t249 + t370 * t265;
t355 = t370 * t246 - t373 * t247 + t376 * t263;
t354 = t370 * t248 - t373 * t249 + t376 * t265;
t353 = t374 * t246 + t377 * t247 - t373 * t263;
t352 = t374 * t248 + t377 * t249 - t373 * t265;
t270 = -t301 * t337 + t303 * t317;
t330 = t299 * t304;
t250 = t270 * t336 + t302 * t330;
t251 = t270 * t302 - t330 * t336;
t271 = t301 * t335 + t303 * t318;
t351 = -t368 * t250 - t372 * t251 - t366 * t271;
t350 = t366 * t250 + t371 * t251 + t367 * t271;
t349 = t372 * t250 + t375 * t251 + t371 * t271;
t348 = t369 * t270 + t374 * t271 - t370 * t330;
t347 = t374 * t270 + t377 * t271 + t373 * t330;
t346 = t370 * t270 - t373 * t271 - t376 * t330;
t334 = pkin(7) * t301;
t333 = Icges(2,4) * t298;
t332 = t298 * t299;
t331 = t299 * t300;
t327 = rSges(7,2) * t247 + t364 * t211 + t212 * t365;
t326 = rSges(7,2) * t249 + t364 * t213 + t214 * t365;
t325 = rSges(7,2) * t271 - t364 * t250 + t251 * t365;
t324 = qJD(2) * t299;
t323 = V_base(5) * qJ(1) + V_base(1);
t319 = qJD(1) + V_base(3);
t279 = t298 * t324 + V_base(4);
t291 = qJD(2) * t301 + V_base(6);
t245 = qJD(3) * t265 + t279;
t278 = -t300 * t324 + V_base(5);
t244 = qJD(3) * t263 + t278;
t267 = -qJD(3) * t330 + t291;
t273 = pkin(1) * t298 - pkin(7) * t331;
t316 = -t273 * V_base(6) + V_base(5) * t334 + t323;
t274 = pkin(1) * t300 + pkin(7) * t332;
t315 = V_base(4) * t273 - t274 * V_base(5) + t319;
t314 = V_base(6) * t274 + V_base(2) + (-qJ(1) - t334) * V_base(4);
t235 = pkin(2) * t264 + pkin(8) * t263;
t272 = (pkin(2) * t303 - pkin(8) * t304) * t299;
t313 = -t235 * t291 + t278 * t272 + t316;
t236 = pkin(2) * t266 + pkin(8) * t265;
t312 = t279 * t235 - t236 * t278 + t315;
t311 = t291 * t236 - t272 * t279 + t314;
t237 = pkin(3) * t271 + qJ(4) * t270;
t310 = qJD(4) * t248 + t244 * t237 + t313;
t205 = pkin(3) * t247 + qJ(4) * t246;
t309 = qJD(4) * t270 + t245 * t205 + t312;
t206 = pkin(3) * t249 + qJ(4) * t248;
t308 = qJD(4) * t246 + t267 * t206 + t311;
t215 = pkin(4) * t263 + pkin(9) * t247;
t253 = -pkin(4) * t330 + t271 * pkin(9);
t307 = t244 * t253 + (-t205 - t215) * t267 + t310;
t216 = pkin(4) * t265 + pkin(9) * t249;
t306 = t245 * t215 + (-t206 - t216) * t244 + t309;
t305 = t267 * t216 + (-t237 - t253) * t245 + t308;
t296 = Icges(2,4) * t300;
t287 = rSges(2,1) * t300 - rSges(2,2) * t298;
t286 = rSges(2,1) * t298 + rSges(2,2) * t300;
t285 = Icges(2,1) * t300 - t333;
t284 = Icges(2,1) * t298 + t296;
t283 = -Icges(2,2) * t298 + t296;
t282 = Icges(2,2) * t300 + t333;
t277 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t276 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t275 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t259 = t301 * rSges(3,3) + (rSges(3,1) * t303 + rSges(3,2) * t304) * t299;
t258 = Icges(3,5) * t301 + (Icges(3,1) * t303 + Icges(3,4) * t304) * t299;
t257 = Icges(3,6) * t301 + (Icges(3,4) * t303 + Icges(3,2) * t304) * t299;
t256 = Icges(3,3) * t301 + (Icges(3,5) * t303 + Icges(3,6) * t304) * t299;
t255 = V_base(5) * rSges(2,3) - t286 * V_base(6) + t323;
t254 = t287 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t243 = t286 * V_base(4) - t287 * V_base(5) + t319;
t238 = qJD(5) * t271 + t267;
t234 = t271 * rSges(4,1) - t270 * rSges(4,2) - rSges(4,3) * t330;
t233 = -rSges(5,1) * t330 - t271 * rSges(5,2) + t270 * rSges(5,3);
t226 = rSges(3,1) * t266 - rSges(3,2) * t265 + rSges(3,3) * t332;
t225 = rSges(3,1) * t264 - rSges(3,2) * t263 - rSges(3,3) * t331;
t224 = Icges(3,1) * t266 - Icges(3,4) * t265 + Icges(3,5) * t332;
t223 = Icges(3,1) * t264 - Icges(3,4) * t263 - Icges(3,5) * t331;
t222 = Icges(3,4) * t266 - Icges(3,2) * t265 + Icges(3,6) * t332;
t221 = Icges(3,4) * t264 - Icges(3,2) * t263 - Icges(3,6) * t331;
t220 = Icges(3,5) * t266 - Icges(3,6) * t265 + Icges(3,3) * t332;
t219 = Icges(3,5) * t264 - Icges(3,6) * t263 - Icges(3,3) * t331;
t209 = qJD(5) * t249 + t245;
t208 = qJD(5) * t247 + t244;
t202 = rSges(6,1) * t251 + rSges(6,2) * t250 + rSges(6,3) * t271;
t192 = rSges(4,1) * t249 - rSges(4,2) * t248 + rSges(4,3) * t265;
t191 = rSges(4,1) * t247 - rSges(4,2) * t246 + rSges(4,3) * t263;
t190 = rSges(5,1) * t265 - rSges(5,2) * t249 + rSges(5,3) * t248;
t189 = rSges(5,1) * t263 - rSges(5,2) * t247 + rSges(5,3) * t246;
t173 = -t225 * t291 + t259 * t278 + t316;
t172 = t226 * t291 - t259 * t279 + t314;
t171 = rSges(6,1) * t214 - rSges(6,2) * t213 + rSges(6,3) * t249;
t169 = rSges(6,1) * t212 - rSges(6,2) * t211 + rSges(6,3) * t247;
t155 = t225 * t279 - t226 * t278 + t315;
t154 = -t191 * t267 + t234 * t244 + t313;
t153 = t192 * t267 - t234 * t245 + t311;
t152 = t191 * t245 - t192 * t244 + t312;
t151 = t233 * t244 + (-t189 - t205) * t267 + t310;
t150 = t190 * t267 + (-t233 - t237) * t245 + t308;
t149 = t189 * t245 + (-t190 - t206) * t244 + t309;
t148 = -t169 * t238 + t202 * t208 + t307;
t147 = t171 * t238 - t202 * t209 + t305;
t146 = t169 * t209 - t171 * t208 + t306;
t145 = qJD(6) * t213 + t208 * t325 - t238 * t327 + t307;
t144 = qJD(6) * t211 - t209 * t325 + t238 * t326 + t305;
t143 = -qJD(6) * t250 - t208 * t326 + t209 * t327 + t306;
t1 = t291 * ((t219 * t278 + t220 * t279 + t256 * t291) * t301 + ((t222 * t304 + t224 * t303) * t279 + (t221 * t304 + t223 * t303) * t278 + (t257 * t304 + t258 * t303) * t291) * t299) / 0.2e1 + t279 * ((t220 * t332 - t265 * t222 + t266 * t224) * t279 + (t219 * t332 - t221 * t265 + t223 * t266) * t278 + (t256 * t332 - t257 * t265 + t258 * t266) * t291) / 0.2e1 + t278 * ((-t220 * t331 - t222 * t263 + t224 * t264) * t279 + (-t219 * t331 - t263 * t221 + t264 * t223) * t278 + (-t256 * t331 - t257 * t263 + t258 * t264) * t291) / 0.2e1 + m(6) * (t146 ^ 2 + t147 ^ 2 + t148 ^ 2) / 0.2e1 + m(7) * (t143 ^ 2 + t144 ^ 2 + t145 ^ 2) / 0.2e1 + m(1) * (t275 ^ 2 + t276 ^ 2 + t277 ^ 2) / 0.2e1 + m(5) * (t149 ^ 2 + t150 ^ 2 + t151 ^ 2) / 0.2e1 + m(4) * (t152 ^ 2 + t153 ^ 2 + t154 ^ 2) / 0.2e1 + m(3) * (t155 ^ 2 + t172 ^ 2 + t173 ^ 2) / 0.2e1 + m(2) * (t243 ^ 2 + t254 ^ 2 + t255 ^ 2) / 0.2e1 + ((t211 * t351 + t212 * t349 + t247 * t350) * t238 + (t211 * t362 + t212 * t358 + t247 * t360) * t209 + (t363 * t211 + t359 * t212 + t361 * t247) * t208) * t208 / 0.2e1 + ((t213 * t351 + t214 * t349 + t249 * t350) * t238 + (t362 * t213 + t358 * t214 + t360 * t249) * t209 + (t213 * t363 + t359 * t214 + t361 * t249) * t208) * t209 / 0.2e1 + ((-t351 * t250 + t349 * t251 + t350 * t271) * t238 + (-t362 * t250 + t251 * t358 + t360 * t271) * t209 + (-t250 * t363 + t359 * t251 + t361 * t271) * t208) * t238 / 0.2e1 + ((t246 * t348 + t247 * t347 + t263 * t346) * t267 + (t246 * t356 + t247 * t352 + t263 * t354) * t245 + (t357 * t246 + t353 * t247 + t355 * t263) * t244) * t244 / 0.2e1 + ((t248 * t348 + t249 * t347 + t265 * t346) * t267 + (t356 * t248 + t352 * t249 + t354 * t265) * t245 + (t248 * t357 + t249 * t353 + t265 * t355) * t244) * t245 / 0.2e1 + ((t348 * t270 + t347 * t271 - t346 * t330) * t267 + (t270 * t356 + t271 * t352 - t354 * t330) * t245 + (t270 * t357 + t271 * t353 - t355 * t330) * t244) * t267 / 0.2e1 + ((-t282 * t298 + t284 * t300 + Icges(1,4)) * V_base(5) + (-t298 * t283 + t300 * t285 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t300 * t282 + t298 * t284 + Icges(1,2)) * V_base(5) + (t283 * t300 + t285 * t298 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t300 - Icges(2,6) * t298 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t298 + Icges(2,6) * t300 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
