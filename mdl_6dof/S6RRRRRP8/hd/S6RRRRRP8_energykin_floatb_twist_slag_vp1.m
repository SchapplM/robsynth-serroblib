% Calculate kinetic energy for
% S6RRRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP8_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRP8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP8_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP8_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP8_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:47:10
% EndTime: 2019-03-10 01:47:14
% DurationCPUTime: 4.51s
% Computational Cost: add. (3564->387), mult. (6139->564), div. (0->0), fcn. (7307->12), ass. (0->177)
t384 = Icges(6,1) + Icges(7,1);
t383 = -Icges(6,4) + Icges(7,5);
t382 = Icges(7,4) + Icges(6,5);
t381 = Icges(6,2) + Icges(7,3);
t380 = Icges(7,2) + Icges(6,3);
t379 = -Icges(6,6) + Icges(7,6);
t378 = rSges(7,1) + pkin(5);
t377 = rSges(7,3) + qJ(6);
t317 = cos(pkin(6));
t321 = sin(qJ(1));
t323 = cos(qJ(2));
t351 = t321 * t323;
t320 = sin(qJ(2));
t324 = cos(qJ(1));
t352 = t320 * t324;
t284 = t317 * t352 + t351;
t349 = qJ(3) + qJ(4);
t313 = sin(t349);
t339 = cos(t349);
t316 = sin(pkin(6));
t355 = t316 * t324;
t253 = t284 * t339 - t313 * t355;
t350 = t323 * t324;
t353 = t320 * t321;
t283 = -t317 * t350 + t353;
t318 = sin(qJ(5));
t364 = cos(qJ(5));
t221 = t253 * t318 - t283 * t364;
t222 = t253 * t364 + t283 * t318;
t338 = t316 * t339;
t252 = t284 * t313 + t324 * t338;
t376 = t221 * t381 + t222 * t383 + t252 * t379;
t286 = -t317 * t353 + t350;
t358 = t316 * t321;
t255 = t286 * t339 + t313 * t358;
t285 = t317 * t351 + t352;
t223 = t255 * t318 - t285 * t364;
t224 = t255 * t364 + t285 * t318;
t254 = t286 * t313 - t321 * t338;
t375 = t223 * t381 + t224 * t383 + t254 * t379;
t374 = t221 * t379 + t222 * t382 + t252 * t380;
t373 = t223 * t379 + t224 * t382 + t254 * t380;
t372 = t383 * t221 + t222 * t384 + t382 * t252;
t371 = t383 * t223 + t224 * t384 + t382 * t254;
t272 = t317 * t313 + t320 * t338;
t356 = t316 * t323;
t250 = t272 * t318 + t356 * t364;
t251 = t272 * t364 - t318 * t356;
t359 = t316 * t320;
t271 = t313 * t359 - t317 * t339;
t370 = t250 * t381 + t251 * t383 + t271 * t379;
t369 = t250 * t379 + t251 * t382 + t271 * t380;
t368 = t383 * t250 + t251 * t384 + t382 * t271;
t363 = pkin(8) * t317;
t322 = cos(qJ(3));
t362 = pkin(3) * t322;
t360 = Icges(2,4) * t321;
t357 = t316 * t322;
t319 = sin(qJ(3));
t354 = t317 * t319;
t348 = rSges(7,2) * t252 + t377 * t221 + t222 * t378;
t347 = rSges(7,2) * t254 + t377 * t223 + t224 * t378;
t346 = rSges(7,2) * t271 + t377 * t250 + t251 * t378;
t345 = qJD(2) * t316;
t344 = V_base(5) * pkin(7) + V_base(1);
t341 = t319 * t358;
t340 = t319 * t355;
t295 = t321 * t345 + V_base(4);
t312 = V_base(6) + qJD(1);
t257 = qJD(3) * t285 + t295;
t297 = qJD(2) * t317 + t312;
t232 = qJD(4) * t285 + t257;
t294 = -t324 * t345 + V_base(5);
t289 = pkin(1) * t321 - pkin(8) * t355;
t337 = -t289 * t312 + V_base(5) * t363 + t344;
t290 = pkin(1) * t324 + pkin(8) * t358;
t336 = V_base(4) * t289 - t290 * V_base(5) + V_base(3);
t256 = qJD(3) * t283 + t294;
t231 = qJD(4) * t283 + t256;
t335 = t312 * t290 + V_base(2) + (-pkin(7) - t363) * V_base(4);
t266 = (-qJD(3) - qJD(4)) * t356 + t297;
t248 = t284 * pkin(2) + t283 * pkin(9);
t288 = (pkin(2) * t320 - pkin(9) * t323) * t316;
t334 = -t248 * t297 + t294 * t288 + t337;
t249 = t286 * pkin(2) + t285 * pkin(9);
t333 = t295 * t248 - t249 * t294 + t336;
t332 = t297 * t249 - t288 * t295 + t335;
t205 = -pkin(3) * t340 + pkin(10) * t283 + t284 * t362;
t245 = pkin(3) * t354 + (-pkin(10) * t323 + t320 * t362) * t316;
t279 = -qJD(3) * t356 + t297;
t331 = -t205 * t279 + t256 * t245 + t334;
t206 = pkin(3) * t341 + pkin(10) * t285 + t286 * t362;
t330 = t257 * t205 - t206 * t256 + t333;
t329 = t279 * t206 - t245 * t257 + t332;
t219 = pkin(4) * t253 + pkin(11) * t252;
t238 = pkin(4) * t272 + pkin(11) * t271;
t328 = -t219 * t266 + t231 * t238 + t331;
t220 = pkin(4) * t255 + pkin(11) * t254;
t327 = t232 * t219 - t220 * t231 + t330;
t326 = t266 * t220 - t232 * t238 + t329;
t314 = Icges(2,4) * t324;
t305 = rSges(2,1) * t324 - rSges(2,2) * t321;
t304 = rSges(2,1) * t321 + rSges(2,2) * t324;
t303 = Icges(2,1) * t324 - t360;
t302 = Icges(2,1) * t321 + t314;
t301 = -Icges(2,2) * t321 + t314;
t300 = Icges(2,2) * t324 + t360;
t293 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t292 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t291 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t282 = t320 * t357 + t354;
t281 = t317 * t322 - t319 * t359;
t270 = rSges(3,3) * t317 + (rSges(3,1) * t320 + rSges(3,2) * t323) * t316;
t269 = Icges(3,5) * t317 + (Icges(3,1) * t320 + Icges(3,4) * t323) * t316;
t268 = Icges(3,6) * t317 + (Icges(3,4) * t320 + Icges(3,2) * t323) * t316;
t267 = Icges(3,3) * t317 + (Icges(3,5) * t320 + Icges(3,6) * t323) * t316;
t265 = V_base(5) * rSges(2,3) - t304 * t312 + t344;
t264 = t305 * t312 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t262 = t304 * V_base(4) - t305 * V_base(5) + V_base(3);
t261 = t286 * t322 + t341;
t260 = -t286 * t319 + t321 * t357;
t259 = t284 * t322 - t340;
t258 = -t284 * t319 - t322 * t355;
t247 = rSges(3,1) * t286 - rSges(3,2) * t285 + rSges(3,3) * t358;
t246 = rSges(3,1) * t284 - rSges(3,2) * t283 - rSges(3,3) * t355;
t244 = Icges(3,1) * t286 - Icges(3,4) * t285 + Icges(3,5) * t358;
t243 = Icges(3,1) * t284 - Icges(3,4) * t283 - Icges(3,5) * t355;
t242 = Icges(3,4) * t286 - Icges(3,2) * t285 + Icges(3,6) * t358;
t241 = Icges(3,4) * t284 - Icges(3,2) * t283 - Icges(3,6) * t355;
t240 = Icges(3,5) * t286 - Icges(3,6) * t285 + Icges(3,3) * t358;
t239 = Icges(3,5) * t284 - Icges(3,6) * t283 - Icges(3,3) * t355;
t237 = rSges(4,1) * t282 + rSges(4,2) * t281 - rSges(4,3) * t356;
t236 = Icges(4,1) * t282 + Icges(4,4) * t281 - Icges(4,5) * t356;
t235 = Icges(4,4) * t282 + Icges(4,2) * t281 - Icges(4,6) * t356;
t234 = Icges(4,5) * t282 + Icges(4,6) * t281 - Icges(4,3) * t356;
t229 = qJD(5) * t271 + t266;
t228 = rSges(5,1) * t272 - rSges(5,2) * t271 - rSges(5,3) * t356;
t227 = Icges(5,1) * t272 - Icges(5,4) * t271 - Icges(5,5) * t356;
t226 = Icges(5,4) * t272 - Icges(5,2) * t271 - Icges(5,6) * t356;
t225 = Icges(5,5) * t272 - Icges(5,6) * t271 - Icges(5,3) * t356;
t216 = rSges(4,1) * t261 + rSges(4,2) * t260 + rSges(4,3) * t285;
t215 = rSges(4,1) * t259 + rSges(4,2) * t258 + rSges(4,3) * t283;
t214 = Icges(4,1) * t261 + Icges(4,4) * t260 + Icges(4,5) * t285;
t213 = Icges(4,1) * t259 + Icges(4,4) * t258 + Icges(4,5) * t283;
t212 = Icges(4,4) * t261 + Icges(4,2) * t260 + Icges(4,6) * t285;
t211 = Icges(4,4) * t259 + Icges(4,2) * t258 + Icges(4,6) * t283;
t210 = Icges(4,5) * t261 + Icges(4,6) * t260 + Icges(4,3) * t285;
t209 = Icges(4,5) * t259 + Icges(4,6) * t258 + Icges(4,3) * t283;
t208 = qJD(5) * t254 + t232;
t207 = qJD(5) * t252 + t231;
t204 = rSges(5,1) * t255 - rSges(5,2) * t254 + rSges(5,3) * t285;
t203 = rSges(5,1) * t253 - rSges(5,2) * t252 + rSges(5,3) * t283;
t202 = Icges(5,1) * t255 - Icges(5,4) * t254 + Icges(5,5) * t285;
t201 = Icges(5,1) * t253 - Icges(5,4) * t252 + Icges(5,5) * t283;
t200 = Icges(5,4) * t255 - Icges(5,2) * t254 + Icges(5,6) * t285;
t199 = Icges(5,4) * t253 - Icges(5,2) * t252 + Icges(5,6) * t283;
t198 = Icges(5,5) * t255 - Icges(5,6) * t254 + Icges(5,3) * t285;
t197 = Icges(5,5) * t253 - Icges(5,6) * t252 + Icges(5,3) * t283;
t195 = rSges(6,1) * t251 - rSges(6,2) * t250 + rSges(6,3) * t271;
t181 = -t246 * t297 + t270 * t294 + t337;
t180 = t247 * t297 - t270 * t295 + t335;
t179 = t246 * t295 - t247 * t294 + t336;
t178 = rSges(6,1) * t224 - rSges(6,2) * t223 + rSges(6,3) * t254;
t176 = rSges(6,1) * t222 - rSges(6,2) * t221 + rSges(6,3) * t252;
t162 = -t215 * t279 + t237 * t256 + t334;
t161 = t216 * t279 - t237 * t257 + t332;
t160 = t215 * t257 - t216 * t256 + t333;
t159 = -t203 * t266 + t228 * t231 + t331;
t158 = t204 * t266 - t228 * t232 + t329;
t157 = t203 * t232 - t204 * t231 + t330;
t156 = -t176 * t229 + t195 * t207 + t328;
t155 = t178 * t229 - t195 * t208 + t326;
t154 = t176 * t208 - t178 * t207 + t327;
t153 = qJD(6) * t223 + t207 * t346 - t229 * t348 + t328;
t152 = qJD(6) * t221 - t208 * t346 + t229 * t347 + t326;
t151 = qJD(6) * t250 - t207 * t347 + t208 * t348 + t327;
t1 = m(1) * (t291 ^ 2 + t292 ^ 2 + t293 ^ 2) / 0.2e1 + t257 * ((t210 * t285 + t212 * t260 + t214 * t261) * t257 + (t209 * t285 + t211 * t260 + t213 * t261) * t256 + (t234 * t285 + t235 * t260 + t236 * t261) * t279) / 0.2e1 + t232 * ((t198 * t285 - t200 * t254 + t202 * t255) * t232 + (t197 * t285 - t199 * t254 + t201 * t255) * t231 + (t225 * t285 - t226 * t254 + t227 * t255) * t266) / 0.2e1 + t256 * ((t210 * t283 + t212 * t258 + t214 * t259) * t257 + (t209 * t283 + t211 * t258 + t213 * t259) * t256 + (t234 * t283 + t235 * t258 + t236 * t259) * t279) / 0.2e1 + t231 * ((t198 * t283 - t200 * t252 + t202 * t253) * t232 + (t197 * t283 - t199 * t252 + t201 * t253) * t231 + (t225 * t283 - t226 * t252 + t227 * t253) * t266) / 0.2e1 + m(2) * (t262 ^ 2 + t264 ^ 2 + t265 ^ 2) / 0.2e1 + t295 * ((t240 * t358 - t242 * t285 + t244 * t286) * t295 + (t239 * t358 - t241 * t285 + t243 * t286) * t294 + (t267 * t358 - t268 * t285 + t269 * t286) * t297) / 0.2e1 + t279 * ((-t210 * t356 + t212 * t281 + t214 * t282) * t257 + (-t209 * t356 + t211 * t281 + t213 * t282) * t256 + (-t234 * t356 + t235 * t281 + t236 * t282) * t279) / 0.2e1 + t266 * ((-t198 * t356 - t200 * t271 + t202 * t272) * t232 + (-t197 * t356 - t199 * t271 + t201 * t272) * t231 + (-t225 * t356 - t226 * t271 + t227 * t272) * t266) / 0.2e1 + t294 * ((-t240 * t355 - t242 * t283 + t244 * t284) * t295 + (-t239 * t355 - t241 * t283 + t243 * t284) * t294 + (-t267 * t355 - t268 * t283 + t269 * t284) * t297) / 0.2e1 + m(3) * (t179 ^ 2 + t180 ^ 2 + t181 ^ 2) / 0.2e1 + m(5) * (t157 ^ 2 + t158 ^ 2 + t159 ^ 2) / 0.2e1 + m(4) * (t160 ^ 2 + t161 ^ 2 + t162 ^ 2) / 0.2e1 + m(6) * (t154 ^ 2 + t155 ^ 2 + t156 ^ 2) / 0.2e1 + m(7) * (t151 ^ 2 + t152 ^ 2 + t153 ^ 2) / 0.2e1 + t297 * ((t239 * t294 + t240 * t295 + t267 * t297) * t317 + ((t242 * t323 + t244 * t320) * t295 + (t241 * t323 + t243 * t320) * t294 + (t268 * t323 + t269 * t320) * t297) * t316) / 0.2e1 + ((t221 * t370 + t222 * t368 + t252 * t369) * t229 + (t221 * t375 + t222 * t371 + t252 * t373) * t208 + (t376 * t221 + t372 * t222 + t374 * t252) * t207) * t207 / 0.2e1 + ((t223 * t370 + t224 * t368 + t254 * t369) * t229 + (t375 * t223 + t371 * t224 + t373 * t254) * t208 + (t223 * t376 + t372 * t224 + t374 * t254) * t207) * t208 / 0.2e1 + ((t370 * t250 + t368 * t251 + t369 * t271) * t229 + (t250 * t375 + t251 * t371 + t271 * t373) * t208 + (t250 * t376 + t372 * t251 + t374 * t271) * t207) * t229 / 0.2e1 + ((-t300 * t321 + t302 * t324 + Icges(1,4)) * V_base(5) + (-t301 * t321 + t303 * t324 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t300 * t324 + t302 * t321 + Icges(1,2)) * V_base(5) + (t301 * t324 + t303 * t321 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t321 + Icges(2,6) * t324) * V_base(5) + (Icges(2,5) * t324 - Icges(2,6) * t321) * V_base(4) + Icges(2,3) * t312 / 0.2e1) * t312;
T  = t1;
