% Calculate kinetic energy for
% S6PRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRPR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPR2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:05:16
% EndTime: 2019-03-08 23:05:20
% DurationCPUTime: 4.78s
% Computational Cost: add. (3627->438), mult. (6052->630), div. (0->0), fcn. (7184->14), ass. (0->191)
t371 = Icges(5,2) + Icges(6,3);
t311 = sin(pkin(11));
t314 = cos(pkin(11));
t320 = cos(qJ(2));
t315 = cos(pkin(6));
t318 = sin(qJ(2));
t348 = t315 * t318;
t271 = t311 * t320 + t314 * t348;
t346 = qJ(3) + qJ(4);
t307 = sin(t346);
t335 = cos(t346);
t312 = sin(pkin(6));
t354 = t312 * t314;
t245 = t271 * t335 - t307 * t354;
t347 = t315 * t320;
t270 = t311 * t318 - t314 * t347;
t310 = sin(pkin(12));
t313 = cos(pkin(12));
t209 = -t245 * t310 + t270 * t313;
t357 = t270 * t310;
t210 = t245 * t313 + t357;
t334 = t312 * t335;
t244 = t271 * t307 + t314 * t334;
t369 = -Icges(5,4) * t245 + Icges(6,5) * t210 - Icges(5,6) * t270 + Icges(6,6) * t209 + t371 * t244;
t273 = -t311 * t348 + t314 * t320;
t355 = t311 * t312;
t247 = t273 * t335 + t307 * t355;
t272 = t311 * t347 + t314 * t318;
t211 = -t247 * t310 + t272 * t313;
t356 = t272 * t310;
t212 = t247 * t313 + t356;
t246 = t273 * t307 - t311 * t334;
t368 = -Icges(5,4) * t247 + Icges(6,5) * t212 - Icges(5,6) * t272 + Icges(6,6) * t211 + t371 * t246;
t265 = t315 * t307 + t318 * t334;
t350 = t312 * t320;
t242 = -t265 * t310 - t313 * t350;
t338 = t310 * t350;
t243 = t265 * t313 - t338;
t352 = t312 * t318;
t264 = t307 * t352 - t315 * t335;
t367 = -Icges(5,4) * t265 + Icges(6,5) * t243 + Icges(5,6) * t350 + Icges(6,6) * t242 + t371 * t264;
t362 = pkin(7) * t315;
t319 = cos(qJ(3));
t361 = pkin(3) * t319;
t360 = pkin(5) * t313;
t358 = Icges(2,4) * t311;
t317 = sin(qJ(3));
t353 = t312 * t317;
t351 = t312 * t319;
t349 = t315 * t317;
t344 = qJD(2) * t312;
t343 = V_base(5) * qJ(1) + V_base(1);
t339 = qJD(1) + V_base(3);
t337 = t311 * t353;
t336 = t314 * t353;
t286 = t311 * t344 + V_base(4);
t297 = qJD(2) * t315 + V_base(6);
t250 = qJD(3) * t272 + t286;
t219 = qJD(4) * t272 + t250;
t285 = -t314 * t344 + V_base(5);
t249 = qJD(3) * t270 + t285;
t280 = pkin(1) * t311 - pkin(7) * t354;
t333 = -t280 * V_base(6) + V_base(5) * t362 + t343;
t281 = pkin(1) * t314 + pkin(7) * t355;
t332 = V_base(4) * t280 - t281 * V_base(5) + t339;
t218 = qJD(4) * t270 + t249;
t258 = (-qJD(3) - qJD(4)) * t350 + t297;
t331 = V_base(6) * t281 + V_base(2) + (-qJ(1) - t362) * V_base(4);
t238 = t271 * pkin(2) + t270 * pkin(8);
t279 = (pkin(2) * t318 - pkin(8) * t320) * t312;
t330 = -t238 * t297 + t285 * t279 + t333;
t239 = t273 * pkin(2) + t272 * pkin(8);
t329 = t286 * t238 - t239 * t285 + t332;
t328 = t297 * t239 - t279 * t286 + t331;
t190 = -pkin(3) * t336 + pkin(9) * t270 + t271 * t361;
t235 = pkin(3) * t349 + (-pkin(9) * t320 + t318 * t361) * t312;
t274 = -qJD(3) * t350 + t297;
t327 = -t190 * t274 + t249 * t235 + t330;
t191 = pkin(3) * t337 + pkin(9) * t272 + t273 * t361;
t326 = t250 * t190 - t191 * t249 + t329;
t325 = t274 * t191 - t235 * t250 + t328;
t233 = pkin(4) * t265 + qJ(5) * t264;
t324 = qJD(5) * t246 + t218 * t233 + t327;
t203 = pkin(4) * t245 + qJ(5) * t244;
t323 = qJD(5) * t264 + t219 * t203 + t326;
t204 = pkin(4) * t247 + qJ(5) * t246;
t322 = qJD(5) * t244 + t258 * t204 + t325;
t309 = pkin(12) + qJ(6);
t306 = Icges(2,4) * t314;
t305 = cos(t309);
t304 = sin(t309);
t295 = rSges(2,1) * t314 - rSges(2,2) * t311;
t294 = rSges(2,1) * t311 + rSges(2,2) * t314;
t293 = Icges(2,1) * t314 - t358;
t292 = Icges(2,1) * t311 + t306;
t291 = -Icges(2,2) * t311 + t306;
t290 = Icges(2,2) * t314 + t358;
t284 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t283 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t282 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t278 = t318 * t351 + t349;
t277 = t315 * t319 - t317 * t352;
t263 = rSges(3,3) * t315 + (rSges(3,1) * t318 + rSges(3,2) * t320) * t312;
t262 = Icges(3,5) * t315 + (Icges(3,1) * t318 + Icges(3,4) * t320) * t312;
t261 = Icges(3,6) * t315 + (Icges(3,4) * t318 + Icges(3,2) * t320) * t312;
t260 = Icges(3,3) * t315 + (Icges(3,5) * t318 + Icges(3,6) * t320) * t312;
t257 = V_base(5) * rSges(2,3) - t294 * V_base(6) + t343;
t256 = t295 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t254 = t273 * t319 + t337;
t253 = -t273 * t317 + t311 * t351;
t252 = t271 * t319 - t336;
t251 = -t271 * t317 - t314 * t351;
t248 = t294 * V_base(4) - t295 * V_base(5) + t339;
t237 = t265 * t305 - t304 * t350;
t236 = -t265 * t304 - t305 * t350;
t234 = rSges(4,1) * t278 + rSges(4,2) * t277 - rSges(4,3) * t350;
t232 = Icges(4,1) * t278 + Icges(4,4) * t277 - Icges(4,5) * t350;
t231 = Icges(4,4) * t278 + Icges(4,2) * t277 - Icges(4,6) * t350;
t230 = Icges(4,5) * t278 + Icges(4,6) * t277 - Icges(4,3) * t350;
t229 = rSges(3,1) * t273 - rSges(3,2) * t272 + rSges(3,3) * t355;
t228 = rSges(3,1) * t271 - rSges(3,2) * t270 - rSges(3,3) * t354;
t227 = Icges(3,1) * t273 - Icges(3,4) * t272 + Icges(3,5) * t355;
t226 = Icges(3,1) * t271 - Icges(3,4) * t270 - Icges(3,5) * t354;
t225 = Icges(3,4) * t273 - Icges(3,2) * t272 + Icges(3,6) * t355;
t224 = Icges(3,4) * t271 - Icges(3,2) * t270 - Icges(3,6) * t354;
t223 = Icges(3,5) * t273 - Icges(3,6) * t272 + Icges(3,3) * t355;
t222 = Icges(3,5) * t271 - Icges(3,6) * t270 - Icges(3,3) * t354;
t220 = qJD(6) * t264 + t258;
t216 = rSges(5,1) * t265 - rSges(5,2) * t264 - rSges(5,3) * t350;
t215 = Icges(5,1) * t265 - Icges(5,4) * t264 - Icges(5,5) * t350;
t213 = Icges(5,5) * t265 - Icges(5,6) * t264 - Icges(5,3) * t350;
t208 = t247 * t305 + t272 * t304;
t207 = -t247 * t304 + t272 * t305;
t206 = t245 * t305 + t270 * t304;
t205 = -t245 * t304 + t270 * t305;
t201 = rSges(4,1) * t254 + rSges(4,2) * t253 + rSges(4,3) * t272;
t200 = rSges(4,1) * t252 + rSges(4,2) * t251 + rSges(4,3) * t270;
t199 = Icges(4,1) * t254 + Icges(4,4) * t253 + Icges(4,5) * t272;
t198 = Icges(4,1) * t252 + Icges(4,4) * t251 + Icges(4,5) * t270;
t197 = Icges(4,4) * t254 + Icges(4,2) * t253 + Icges(4,6) * t272;
t196 = Icges(4,4) * t252 + Icges(4,2) * t251 + Icges(4,6) * t270;
t195 = Icges(4,5) * t254 + Icges(4,6) * t253 + Icges(4,3) * t272;
t194 = Icges(4,5) * t252 + Icges(4,6) * t251 + Icges(4,3) * t270;
t193 = qJD(6) * t246 + t219;
t192 = qJD(6) * t244 + t218;
t189 = rSges(5,1) * t247 - rSges(5,2) * t246 + rSges(5,3) * t272;
t188 = rSges(5,1) * t245 - rSges(5,2) * t244 + rSges(5,3) * t270;
t187 = Icges(5,1) * t247 - Icges(5,4) * t246 + Icges(5,5) * t272;
t186 = Icges(5,1) * t245 - Icges(5,4) * t244 + Icges(5,5) * t270;
t183 = Icges(5,5) * t247 - Icges(5,6) * t246 + Icges(5,3) * t272;
t182 = Icges(5,5) * t245 - Icges(5,6) * t244 + Icges(5,3) * t270;
t180 = rSges(6,1) * t243 + rSges(6,2) * t242 + rSges(6,3) * t264;
t179 = Icges(6,1) * t243 + Icges(6,4) * t242 + Icges(6,5) * t264;
t178 = Icges(6,4) * t243 + Icges(6,2) * t242 + Icges(6,6) * t264;
t175 = rSges(7,1) * t237 + rSges(7,2) * t236 + rSges(7,3) * t264;
t174 = Icges(7,1) * t237 + Icges(7,4) * t236 + Icges(7,5) * t264;
t173 = Icges(7,4) * t237 + Icges(7,2) * t236 + Icges(7,6) * t264;
t172 = Icges(7,5) * t237 + Icges(7,6) * t236 + Icges(7,3) * t264;
t171 = -pkin(5) * t338 + pkin(10) * t264 + t265 * t360;
t169 = -t228 * t297 + t263 * t285 + t333;
t168 = t229 * t297 - t263 * t286 + t331;
t165 = t228 * t286 - t229 * t285 + t332;
t164 = rSges(6,1) * t212 + rSges(6,2) * t211 + rSges(6,3) * t246;
t163 = rSges(6,1) * t210 + rSges(6,2) * t209 + rSges(6,3) * t244;
t162 = Icges(6,1) * t212 + Icges(6,4) * t211 + Icges(6,5) * t246;
t161 = Icges(6,1) * t210 + Icges(6,4) * t209 + Icges(6,5) * t244;
t160 = Icges(6,4) * t212 + Icges(6,2) * t211 + Icges(6,6) * t246;
t159 = Icges(6,4) * t210 + Icges(6,2) * t209 + Icges(6,6) * t244;
t156 = rSges(7,1) * t208 + rSges(7,2) * t207 + rSges(7,3) * t246;
t155 = rSges(7,1) * t206 + rSges(7,2) * t205 + rSges(7,3) * t244;
t154 = Icges(7,1) * t208 + Icges(7,4) * t207 + Icges(7,5) * t246;
t153 = Icges(7,1) * t206 + Icges(7,4) * t205 + Icges(7,5) * t244;
t152 = Icges(7,4) * t208 + Icges(7,2) * t207 + Icges(7,6) * t246;
t151 = Icges(7,4) * t206 + Icges(7,2) * t205 + Icges(7,6) * t244;
t150 = Icges(7,5) * t208 + Icges(7,6) * t207 + Icges(7,3) * t246;
t149 = Icges(7,5) * t206 + Icges(7,6) * t205 + Icges(7,3) * t244;
t148 = pkin(5) * t356 + pkin(10) * t246 + t247 * t360;
t147 = pkin(5) * t357 + pkin(10) * t244 + t245 * t360;
t146 = -t200 * t274 + t234 * t249 + t330;
t145 = t201 * t274 - t234 * t250 + t328;
t144 = t200 * t250 - t201 * t249 + t329;
t143 = -t188 * t258 + t216 * t218 + t327;
t142 = t189 * t258 - t216 * t219 + t325;
t141 = t188 * t219 - t189 * t218 + t326;
t140 = t180 * t218 + (-t163 - t203) * t258 + t324;
t139 = t164 * t258 + (-t180 - t233) * t219 + t322;
t138 = t163 * t219 + (-t164 - t204) * t218 + t323;
t137 = -t155 * t220 + t171 * t218 + t175 * t192 + (-t147 - t203) * t258 + t324;
t136 = t148 * t258 + t156 * t220 - t175 * t193 + (-t171 - t233) * t219 + t322;
t135 = t155 * t193 - t156 * t192 + (-t148 - t204) * t218 + t147 * t219 + t323;
t1 = m(3) * (t165 ^ 2 + t168 ^ 2 + t169 ^ 2) / 0.2e1 + m(5) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(4) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(6) * (t138 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + m(7) * (t135 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + m(1) * (t282 ^ 2 + t283 ^ 2 + t284 ^ 2) / 0.2e1 + t249 * ((t195 * t270 + t197 * t251 + t199 * t252) * t250 + (t194 * t270 + t196 * t251 + t198 * t252) * t249 + (t230 * t270 + t231 * t251 + t232 * t252) * t274) / 0.2e1 + t250 * ((t195 * t272 + t197 * t253 + t199 * t254) * t250 + (t194 * t272 + t196 * t253 + t198 * t254) * t249 + (t230 * t272 + t231 * t253 + t232 * t254) * t274) / 0.2e1 + t220 * ((t150 * t264 + t152 * t236 + t154 * t237) * t193 + (t149 * t264 + t151 * t236 + t153 * t237) * t192 + (t172 * t264 + t173 * t236 + t174 * t237) * t220) / 0.2e1 + m(2) * (t248 ^ 2 + t256 ^ 2 + t257 ^ 2) / 0.2e1 + t193 * ((t150 * t246 + t152 * t207 + t154 * t208) * t193 + (t149 * t246 + t151 * t207 + t153 * t208) * t192 + (t172 * t246 + t173 * t207 + t174 * t208) * t220) / 0.2e1 + t192 * ((t150 * t244 + t152 * t205 + t154 * t206) * t193 + (t149 * t244 + t151 * t205 + t153 * t206) * t192 + (t172 * t244 + t173 * t205 + t174 * t206) * t220) / 0.2e1 + t274 * ((-t195 * t350 + t197 * t277 + t199 * t278) * t250 + (-t194 * t350 + t196 * t277 + t198 * t278) * t249 + (-t230 * t350 + t231 * t277 + t232 * t278) * t274) / 0.2e1 + t285 * ((-t223 * t354 - t225 * t270 + t227 * t271) * t286 + (-t222 * t354 - t224 * t270 + t226 * t271) * t285 + (-t260 * t354 - t261 * t270 + t262 * t271) * t297) / 0.2e1 + t286 * ((t223 * t355 - t225 * t272 + t227 * t273) * t286 + (t222 * t355 - t224 * t272 + t226 * t273) * t285 + (t260 * t355 - t261 * t272 + t262 * t273) * t297) / 0.2e1 + t297 * ((t222 * t285 + t223 * t286 + t260 * t297) * t315 + ((t225 * t320 + t227 * t318) * t286 + (t224 * t320 + t226 * t318) * t285 + (t261 * t320 + t262 * t318) * t297) * t312) / 0.2e1 + ((t178 * t209 + t179 * t210 + t213 * t270 + t215 * t245 + t244 * t367) * t258 + (t160 * t209 + t162 * t210 + t183 * t270 + t187 * t245 + t244 * t368) * t219 + (t159 * t209 + t161 * t210 + t182 * t270 + t186 * t245 + t369 * t244) * t218) * t218 / 0.2e1 + ((t178 * t211 + t179 * t212 + t213 * t272 + t215 * t247 + t246 * t367) * t258 + (t211 * t160 + t212 * t162 + t183 * t272 + t187 * t247 + t368 * t246) * t219 + (t159 * t211 + t161 * t212 + t182 * t272 + t186 * t247 + t246 * t369) * t218) * t219 / 0.2e1 + ((t178 * t242 + t179 * t243 - t213 * t350 + t215 * t265 + t367 * t264) * t258 + (t160 * t242 + t162 * t243 - t183 * t350 + t187 * t265 + t264 * t368) * t219 + (t159 * t242 + t161 * t243 - t182 * t350 + t186 * t265 + t264 * t369) * t218) * t258 / 0.2e1 + ((-t290 * t311 + t292 * t314 + Icges(1,4)) * V_base(5) + (-t291 * t311 + t293 * t314 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t290 * t314 + t292 * t311 + Icges(1,2)) * V_base(5) + (t291 * t314 + t293 * t311 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t311 + Icges(2,6) * t314 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t314 - Icges(2,6) * t311 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
