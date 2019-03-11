% Calculate kinetic energy for
% S6PRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
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
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPPR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRPPR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPPR1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:58:25
% EndTime: 2019-03-08 20:58:30
% DurationCPUTime: 4.28s
% Computational Cost: add. (3522->432), mult. (5842->600), div. (0->0), fcn. (6932->14), ass. (0->189)
t375 = Icges(5,2) + Icges(6,3);
t374 = Icges(4,3) + Icges(5,3);
t308 = sin(pkin(10));
t311 = cos(pkin(10));
t318 = cos(qJ(2));
t312 = cos(pkin(6));
t316 = sin(qJ(2));
t348 = t312 * t316;
t268 = t308 * t318 + t311 * t348;
t341 = qJ(3) + pkin(11);
t302 = sin(t341);
t332 = cos(t341);
t309 = sin(pkin(6));
t354 = t309 * t311;
t243 = t268 * t332 - t302 * t354;
t347 = t312 * t318;
t267 = t308 * t316 - t311 * t347;
t307 = sin(pkin(12));
t310 = cos(pkin(12));
t209 = -t243 * t307 + t267 * t310;
t357 = t267 * t307;
t210 = t243 * t310 + t357;
t331 = t309 * t332;
t242 = t268 * t302 + t311 * t331;
t373 = -Icges(5,4) * t243 + Icges(6,5) * t210 - Icges(5,6) * t267 + Icges(6,6) * t209 + t242 * t375;
t270 = -t308 * t348 + t311 * t318;
t355 = t308 * t309;
t245 = t270 * t332 + t302 * t355;
t269 = t308 * t347 + t311 * t316;
t211 = -t245 * t307 + t269 * t310;
t356 = t269 * t307;
t212 = t245 * t310 + t356;
t244 = t270 * t302 - t308 * t331;
t372 = -Icges(5,4) * t245 + Icges(6,5) * t212 - Icges(5,6) * t269 + Icges(6,6) * t211 + t244 * t375;
t261 = t312 * t302 + t316 * t331;
t350 = t309 * t318;
t240 = -t261 * t307 - t310 * t350;
t335 = t307 * t350;
t241 = t261 * t310 - t335;
t352 = t309 * t316;
t260 = t302 * t352 - t312 * t332;
t371 = -Icges(5,4) * t261 + Icges(6,5) * t241 + Icges(5,6) * t350 + Icges(6,6) * t240 + t260 * t375;
t315 = sin(qJ(3));
t317 = cos(qJ(3));
t351 = t309 * t317;
t249 = -t268 * t315 - t311 * t351;
t353 = t309 * t315;
t333 = t311 * t353;
t250 = t268 * t317 - t333;
t370 = Icges(4,5) * t250 + Icges(5,5) * t243 + Icges(4,6) * t249 - Icges(5,6) * t242 + t267 * t374;
t251 = -t270 * t315 + t308 * t351;
t334 = t308 * t353;
t252 = t270 * t317 + t334;
t369 = Icges(4,5) * t252 + Icges(5,5) * t245 + Icges(4,6) * t251 - Icges(5,6) * t244 + t269 * t374;
t274 = t312 * t317 - t315 * t352;
t349 = t312 * t315;
t275 = t316 * t351 + t349;
t368 = Icges(4,5) * t275 + Icges(5,5) * t261 + Icges(4,6) * t274 - Icges(5,6) * t260 - t350 * t374;
t362 = pkin(7) * t312;
t361 = pkin(3) * t317;
t360 = pkin(5) * t310;
t358 = Icges(2,4) * t308;
t189 = -pkin(3) * t333 + qJ(4) * t267 + t268 * t361;
t201 = pkin(4) * t243 + qJ(5) * t242;
t345 = -t189 - t201;
t190 = pkin(3) * t334 + qJ(4) * t269 + t270 * t361;
t202 = pkin(4) * t245 + qJ(5) * t244;
t344 = -t190 - t202;
t219 = t261 * pkin(4) + t260 * qJ(5);
t231 = pkin(3) * t349 + (-qJ(4) * t318 + t316 * t361) * t309;
t343 = -t219 - t231;
t342 = qJD(2) * t309;
t340 = V_base(5) * qJ(1) + V_base(1);
t336 = qJD(1) + V_base(3);
t284 = t308 * t342 + V_base(4);
t295 = qJD(2) * t312 + V_base(6);
t248 = qJD(3) * t269 + t284;
t283 = -t311 * t342 + V_base(5);
t247 = qJD(3) * t267 + t283;
t271 = -qJD(3) * t350 + t295;
t277 = pkin(1) * t308 - pkin(7) * t354;
t330 = -t277 * V_base(6) + V_base(5) * t362 + t340;
t278 = pkin(1) * t311 + pkin(7) * t355;
t329 = V_base(4) * t277 - V_base(5) * t278 + t336;
t328 = V_base(6) * t278 + V_base(2) + (-qJ(1) - t362) * V_base(4);
t236 = pkin(2) * t268 + pkin(8) * t267;
t276 = (pkin(2) * t316 - pkin(8) * t318) * t309;
t327 = -t236 * t295 + t283 * t276 + t330;
t237 = pkin(2) * t270 + pkin(8) * t269;
t326 = t284 * t236 - t283 * t237 + t329;
t325 = t295 * t237 - t276 * t284 + t328;
t324 = qJD(4) * t269 + t247 * t231 + t327;
t323 = qJD(4) * t267 + t271 * t190 + t325;
t322 = qJD(5) * t244 + t247 * t219 + t324;
t321 = -qJD(4) * t350 + t248 * t189 + t326;
t320 = qJD(5) * t242 + t271 * t202 + t323;
t319 = qJD(5) * t260 + t248 * t201 + t321;
t306 = pkin(12) + qJ(6);
t304 = Icges(2,4) * t311;
t303 = cos(t306);
t301 = sin(t306);
t292 = rSges(2,1) * t311 - rSges(2,2) * t308;
t291 = rSges(2,1) * t308 + rSges(2,2) * t311;
t290 = Icges(2,1) * t311 - t358;
t289 = Icges(2,1) * t308 + t304;
t288 = -Icges(2,2) * t308 + t304;
t287 = Icges(2,2) * t311 + t358;
t282 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t281 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t280 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t262 = t312 * rSges(3,3) + (rSges(3,1) * t316 + rSges(3,2) * t318) * t309;
t259 = Icges(3,5) * t312 + (Icges(3,1) * t316 + Icges(3,4) * t318) * t309;
t258 = Icges(3,6) * t312 + (Icges(3,4) * t316 + Icges(3,2) * t318) * t309;
t257 = Icges(3,3) * t312 + (Icges(3,5) * t316 + Icges(3,6) * t318) * t309;
t255 = V_base(5) * rSges(2,3) - t291 * V_base(6) + t340;
t254 = t292 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t246 = t291 * V_base(4) - t292 * V_base(5) + t336;
t235 = t261 * t303 - t301 * t350;
t234 = -t261 * t301 - t303 * t350;
t233 = qJD(6) * t260 + t271;
t232 = t275 * rSges(4,1) + t274 * rSges(4,2) - rSges(4,3) * t350;
t230 = Icges(4,1) * t275 + Icges(4,4) * t274 - Icges(4,5) * t350;
t229 = Icges(4,4) * t275 + Icges(4,2) * t274 - Icges(4,6) * t350;
t227 = rSges(3,1) * t270 - rSges(3,2) * t269 + rSges(3,3) * t355;
t226 = rSges(3,1) * t268 - rSges(3,2) * t267 - rSges(3,3) * t354;
t225 = Icges(3,1) * t270 - Icges(3,4) * t269 + Icges(3,5) * t355;
t224 = Icges(3,1) * t268 - Icges(3,4) * t267 - Icges(3,5) * t354;
t223 = Icges(3,4) * t270 - Icges(3,2) * t269 + Icges(3,6) * t355;
t222 = Icges(3,4) * t268 - Icges(3,2) * t267 - Icges(3,6) * t354;
t221 = Icges(3,5) * t270 - Icges(3,6) * t269 + Icges(3,3) * t355;
t220 = Icges(3,5) * t268 - Icges(3,6) * t267 - Icges(3,3) * t354;
t216 = t261 * rSges(5,1) - t260 * rSges(5,2) - rSges(5,3) * t350;
t215 = Icges(5,1) * t261 - Icges(5,4) * t260 - Icges(5,5) * t350;
t208 = t245 * t303 + t269 * t301;
t207 = -t245 * t301 + t269 * t303;
t206 = t243 * t303 + t267 * t301;
t205 = -t243 * t301 + t267 * t303;
t204 = qJD(6) * t244 + t248;
t203 = qJD(6) * t242 + t247;
t198 = rSges(4,1) * t252 + rSges(4,2) * t251 + rSges(4,3) * t269;
t197 = rSges(4,1) * t250 + rSges(4,2) * t249 + rSges(4,3) * t267;
t196 = Icges(4,1) * t252 + Icges(4,4) * t251 + Icges(4,5) * t269;
t195 = Icges(4,1) * t250 + Icges(4,4) * t249 + Icges(4,5) * t267;
t194 = Icges(4,4) * t252 + Icges(4,2) * t251 + Icges(4,6) * t269;
t193 = Icges(4,4) * t250 + Icges(4,2) * t249 + Icges(4,6) * t267;
t187 = rSges(5,1) * t245 - rSges(5,2) * t244 + rSges(5,3) * t269;
t186 = rSges(5,1) * t243 - rSges(5,2) * t242 + rSges(5,3) * t267;
t185 = Icges(5,1) * t245 - Icges(5,4) * t244 + Icges(5,5) * t269;
t184 = Icges(5,1) * t243 - Icges(5,4) * t242 + Icges(5,5) * t267;
t179 = rSges(6,1) * t241 + rSges(6,2) * t240 + rSges(6,3) * t260;
t178 = Icges(6,1) * t241 + Icges(6,4) * t240 + Icges(6,5) * t260;
t177 = Icges(6,4) * t241 + Icges(6,2) * t240 + Icges(6,6) * t260;
t175 = rSges(7,1) * t235 + rSges(7,2) * t234 + rSges(7,3) * t260;
t173 = Icges(7,1) * t235 + Icges(7,4) * t234 + Icges(7,5) * t260;
t172 = Icges(7,4) * t235 + Icges(7,2) * t234 + Icges(7,6) * t260;
t171 = Icges(7,5) * t235 + Icges(7,6) * t234 + Icges(7,3) * t260;
t169 = -pkin(5) * t335 + pkin(9) * t260 + t261 * t360;
t168 = -t226 * t295 + t262 * t283 + t330;
t167 = t227 * t295 - t262 * t284 + t328;
t165 = t226 * t284 - t227 * t283 + t329;
t164 = rSges(6,1) * t212 + rSges(6,2) * t211 + rSges(6,3) * t244;
t163 = rSges(6,1) * t210 + rSges(6,2) * t209 + rSges(6,3) * t242;
t162 = Icges(6,1) * t212 + Icges(6,4) * t211 + Icges(6,5) * t244;
t161 = Icges(6,1) * t210 + Icges(6,4) * t209 + Icges(6,5) * t242;
t160 = Icges(6,4) * t212 + Icges(6,2) * t211 + Icges(6,6) * t244;
t159 = Icges(6,4) * t210 + Icges(6,2) * t209 + Icges(6,6) * t242;
t156 = rSges(7,1) * t208 + rSges(7,2) * t207 + rSges(7,3) * t244;
t155 = rSges(7,1) * t206 + rSges(7,2) * t205 + rSges(7,3) * t242;
t154 = Icges(7,1) * t208 + Icges(7,4) * t207 + Icges(7,5) * t244;
t153 = Icges(7,1) * t206 + Icges(7,4) * t205 + Icges(7,5) * t242;
t152 = Icges(7,4) * t208 + Icges(7,2) * t207 + Icges(7,6) * t244;
t151 = Icges(7,4) * t206 + Icges(7,2) * t205 + Icges(7,6) * t242;
t150 = Icges(7,5) * t208 + Icges(7,6) * t207 + Icges(7,3) * t244;
t149 = Icges(7,5) * t206 + Icges(7,6) * t205 + Icges(7,3) * t242;
t148 = pkin(5) * t356 + pkin(9) * t244 + t245 * t360;
t147 = pkin(5) * t357 + pkin(9) * t242 + t243 * t360;
t146 = -t197 * t271 + t232 * t247 + t327;
t145 = t198 * t271 - t232 * t248 + t325;
t144 = t197 * t248 - t198 * t247 + t326;
t143 = t216 * t247 + (-t186 - t189) * t271 + t324;
t142 = t187 * t271 + (-t216 - t231) * t248 + t323;
t141 = t248 * t186 + (-t187 - t190) * t247 + t321;
t140 = t179 * t247 + (-t163 + t345) * t271 + t322;
t139 = t164 * t271 + (-t179 + t343) * t248 + t320;
t138 = t248 * t163 + (-t164 + t344) * t247 + t319;
t137 = t322 - t155 * t233 + t169 * t247 + t175 * t203 + (-t147 + t345) * t271;
t136 = t148 * t271 + t156 * t233 - t175 * t204 + (-t169 + t343) * t248 + t320;
t135 = t319 + (-t148 + t344) * t247 - t203 * t156 + t204 * t155 + t248 * t147;
t1 = m(7) * (t135 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + m(6) * (t138 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + t295 * ((t220 * t283 + t221 * t284 + t257 * t295) * t312 + ((t223 * t318 + t225 * t316) * t284 + (t222 * t318 + t224 * t316) * t283 + (t258 * t318 + t259 * t316) * t295) * t309) / 0.2e1 + t284 * ((t221 * t355 - t223 * t269 + t225 * t270) * t284 + (t220 * t355 - t222 * t269 + t224 * t270) * t283 + (t257 * t355 - t258 * t269 + t259 * t270) * t295) / 0.2e1 + t283 * ((-t221 * t354 - t223 * t267 + t225 * t268) * t284 + (-t220 * t354 - t222 * t267 + t224 * t268) * t283 + (-t257 * t354 - t258 * t267 + t259 * t268) * t295) / 0.2e1 + m(1) * (t280 ^ 2 + t281 ^ 2 + t282 ^ 2) / 0.2e1 + m(5) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(4) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(3) * (t165 ^ 2 + t167 ^ 2 + t168 ^ 2) / 0.2e1 + t203 * ((t150 * t242 + t152 * t205 + t154 * t206) * t204 + (t242 * t149 + t205 * t151 + t206 * t153) * t203 + (t171 * t242 + t172 * t205 + t173 * t206) * t233) / 0.2e1 + t204 * ((t244 * t150 + t207 * t152 + t208 * t154) * t204 + (t149 * t244 + t151 * t207 + t153 * t208) * t203 + (t171 * t244 + t172 * t207 + t173 * t208) * t233) / 0.2e1 + m(2) * (t246 ^ 2 + t254 ^ 2 + t255 ^ 2) / 0.2e1 + t233 * ((t150 * t260 + t152 * t234 + t154 * t235) * t204 + (t149 * t260 + t151 * t234 + t153 * t235) * t203 + (t171 * t260 + t172 * t234 + t173 * t235) * t233) / 0.2e1 + ((-t287 * t308 + t289 * t311 + Icges(1,4)) * V_base(5) + (-t288 * t308 + t290 * t311 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t287 * t311 + t289 * t308 + Icges(1,2)) * V_base(5) + (t288 * t311 + t290 * t308 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t177 * t209 + t178 * t210 + t215 * t243 + t229 * t249 + t230 * t250 + t242 * t371 + t267 * t368) * t271 + (t160 * t209 + t162 * t210 + t185 * t243 + t194 * t249 + t196 * t250 + t242 * t372 + t267 * t369) * t248 + (t159 * t209 + t161 * t210 + t184 * t243 + t193 * t249 + t195 * t250 + t373 * t242 + t370 * t267) * t247) * t247 / 0.2e1 + ((t177 * t211 + t178 * t212 + t215 * t245 + t229 * t251 + t230 * t252 + t244 * t371 + t269 * t368) * t271 + (t160 * t211 + t162 * t212 + t185 * t245 + t194 * t251 + t196 * t252 + t372 * t244 + t369 * t269) * t248 + (t159 * t211 + t161 * t212 + t184 * t245 + t193 * t251 + t195 * t252 + t244 * t373 + t370 * t269) * t247) * t248 / 0.2e1 + ((t177 * t240 + t178 * t241 + t261 * t215 + t274 * t229 + t275 * t230 + t371 * t260 - t368 * t350) * t271 + (t160 * t240 + t162 * t241 + t261 * t185 + t274 * t194 + t275 * t196 + t260 * t372 - t350 * t369) * t248 + (t159 * t240 + t161 * t241 + t261 * t184 + t274 * t193 + t275 * t195 + t260 * t373 - t370 * t350) * t247) * t271 / 0.2e1 + ((Icges(2,5) * t311 - Icges(2,6) * t308 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t308 + Icges(2,6) * t311 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
