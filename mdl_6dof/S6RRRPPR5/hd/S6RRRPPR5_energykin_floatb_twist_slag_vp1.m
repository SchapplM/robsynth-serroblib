% Calculate kinetic energy for
% S6RRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 15:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPPR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:38:10
% EndTime: 2019-03-09 15:38:15
% DurationCPUTime: 4.41s
% Computational Cost: add. (3582->432), mult. (5842->601), div. (0->0), fcn. (6932->14), ass. (0->190)
t373 = Icges(5,2) + Icges(6,3);
t372 = Icges(4,3) + Icges(5,3);
t311 = cos(pkin(6));
t316 = sin(qJ(1));
t318 = cos(qJ(2));
t347 = t316 * t318;
t315 = sin(qJ(2));
t319 = cos(qJ(1));
t348 = t315 * t319;
t272 = t311 * t348 + t347;
t340 = qJ(3) + pkin(11);
t302 = sin(t340);
t333 = cos(t340);
t309 = sin(pkin(6));
t351 = t309 * t319;
t243 = t272 * t333 - t302 * t351;
t346 = t318 * t319;
t349 = t315 * t316;
t271 = -t311 * t346 + t349;
t308 = sin(pkin(12));
t310 = cos(pkin(12));
t209 = -t243 * t308 + t271 * t310;
t357 = t271 * t308;
t210 = t243 * t310 + t357;
t332 = t309 * t333;
t242 = t272 * t302 + t319 * t332;
t371 = -Icges(5,4) * t243 + Icges(6,5) * t210 - Icges(5,6) * t271 + Icges(6,6) * t209 + t373 * t242;
t274 = -t311 * t349 + t346;
t354 = t309 * t316;
t245 = t274 * t333 + t302 * t354;
t273 = t311 * t347 + t348;
t211 = -t245 * t308 + t273 * t310;
t356 = t273 * t308;
t212 = t245 * t310 + t356;
t244 = t274 * t302 - t316 * t332;
t370 = -Icges(5,4) * t245 + Icges(6,5) * t212 - Icges(5,6) * t273 + Icges(6,6) * t211 + t373 * t244;
t261 = t311 * t302 + t315 * t332;
t352 = t309 * t318;
t240 = -t261 * t308 - t310 * t352;
t336 = t308 * t352;
t241 = t261 * t310 - t336;
t355 = t309 * t315;
t260 = t302 * t355 - t311 * t333;
t369 = -Icges(5,4) * t261 + Icges(6,5) * t241 + Icges(5,6) * t352 + Icges(6,6) * t240 + t373 * t260;
t314 = sin(qJ(3));
t317 = cos(qJ(3));
t248 = -t272 * t314 - t317 * t351;
t334 = t314 * t351;
t249 = t272 * t317 - t334;
t368 = Icges(4,5) * t249 + Icges(5,5) * t243 + Icges(4,6) * t248 - Icges(5,6) * t242 + t372 * t271;
t353 = t309 * t317;
t250 = -t274 * t314 + t316 * t353;
t335 = t314 * t354;
t251 = t274 * t317 + t335;
t367 = Icges(4,5) * t251 + Icges(5,5) * t245 + Icges(4,6) * t250 - Icges(5,6) * t244 + t372 * t273;
t269 = t311 * t317 - t314 * t355;
t350 = t311 * t314;
t270 = t315 * t353 + t350;
t366 = Icges(4,5) * t270 + Icges(5,5) * t261 + Icges(4,6) * t269 - Icges(5,6) * t260 - t372 * t352;
t362 = pkin(8) * t311;
t361 = pkin(3) * t317;
t360 = pkin(5) * t310;
t358 = Icges(2,4) * t316;
t189 = -pkin(3) * t334 + qJ(4) * t271 + t272 * t361;
t201 = pkin(4) * t243 + qJ(5) * t242;
t344 = -t189 - t201;
t190 = pkin(3) * t335 + qJ(4) * t273 + t274 * t361;
t202 = pkin(4) * t245 + qJ(5) * t244;
t343 = -t190 - t202;
t219 = pkin(4) * t261 + qJ(5) * t260;
t223 = pkin(3) * t350 + (-qJ(4) * t318 + t315 * t361) * t309;
t342 = -t219 - t223;
t341 = qJD(2) * t309;
t339 = V_base(5) * pkin(7) + V_base(1);
t284 = t316 * t341 + V_base(4);
t304 = V_base(6) + qJD(1);
t247 = qJD(3) * t273 + t284;
t285 = qJD(2) * t311 + t304;
t283 = -t319 * t341 + V_base(5);
t277 = t316 * pkin(1) - pkin(8) * t351;
t331 = -t277 * t304 + V_base(5) * t362 + t339;
t278 = pkin(1) * t319 + pkin(8) * t354;
t330 = V_base(4) * t277 - t278 * V_base(5) + V_base(3);
t246 = qJD(3) * t271 + t283;
t267 = -qJD(3) * t352 + t285;
t329 = t304 * t278 + V_base(2) + (-pkin(7) - t362) * V_base(4);
t238 = t272 * pkin(2) + t271 * pkin(9);
t276 = (pkin(2) * t315 - pkin(9) * t318) * t309;
t328 = -t238 * t285 + t283 * t276 + t331;
t239 = pkin(2) * t274 + pkin(9) * t273;
t327 = t284 * t238 - t239 * t283 + t330;
t326 = qJD(4) * t273 + t246 * t223 + t328;
t325 = t285 * t239 - t276 * t284 + t329;
t324 = qJD(5) * t244 + t246 * t219 + t326;
t323 = qJD(4) * t271 + t267 * t190 + t325;
t322 = -qJD(4) * t352 + t247 * t189 + t327;
t321 = qJD(5) * t242 + t267 * t202 + t323;
t320 = qJD(5) * t260 + t247 * t201 + t322;
t307 = pkin(12) + qJ(6);
t305 = Icges(2,4) * t319;
t303 = cos(t307);
t301 = sin(t307);
t293 = rSges(2,1) * t319 - t316 * rSges(2,2);
t292 = t316 * rSges(2,1) + rSges(2,2) * t319;
t291 = Icges(2,1) * t319 - t358;
t290 = Icges(2,1) * t316 + t305;
t289 = -Icges(2,2) * t316 + t305;
t288 = Icges(2,2) * t319 + t358;
t281 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t280 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t279 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t262 = rSges(3,3) * t311 + (rSges(3,1) * t315 + rSges(3,2) * t318) * t309;
t259 = Icges(3,5) * t311 + (Icges(3,1) * t315 + Icges(3,4) * t318) * t309;
t258 = Icges(3,6) * t311 + (Icges(3,4) * t315 + Icges(3,2) * t318) * t309;
t257 = Icges(3,3) * t311 + (Icges(3,5) * t315 + Icges(3,6) * t318) * t309;
t255 = V_base(5) * rSges(2,3) - t292 * t304 + t339;
t254 = t293 * t304 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t252 = t292 * V_base(4) - t293 * V_base(5) + V_base(3);
t235 = t261 * t303 - t301 * t352;
t234 = -t261 * t301 - t303 * t352;
t233 = qJD(6) * t260 + t267;
t232 = rSges(3,1) * t274 - rSges(3,2) * t273 + rSges(3,3) * t354;
t231 = t272 * rSges(3,1) - t271 * rSges(3,2) - rSges(3,3) * t351;
t230 = Icges(3,1) * t274 - Icges(3,4) * t273 + Icges(3,5) * t354;
t229 = Icges(3,1) * t272 - Icges(3,4) * t271 - Icges(3,5) * t351;
t228 = Icges(3,4) * t274 - Icges(3,2) * t273 + Icges(3,6) * t354;
t227 = Icges(3,4) * t272 - Icges(3,2) * t271 - Icges(3,6) * t351;
t226 = Icges(3,5) * t274 - Icges(3,6) * t273 + Icges(3,3) * t354;
t225 = Icges(3,5) * t272 - Icges(3,6) * t271 - Icges(3,3) * t351;
t224 = rSges(4,1) * t270 + rSges(4,2) * t269 - rSges(4,3) * t352;
t222 = Icges(4,1) * t270 + Icges(4,4) * t269 - Icges(4,5) * t352;
t221 = Icges(4,4) * t270 + Icges(4,2) * t269 - Icges(4,6) * t352;
t216 = rSges(5,1) * t261 - rSges(5,2) * t260 - rSges(5,3) * t352;
t215 = Icges(5,1) * t261 - Icges(5,4) * t260 - Icges(5,5) * t352;
t208 = t245 * t303 + t273 * t301;
t207 = -t245 * t301 + t273 * t303;
t206 = t243 * t303 + t271 * t301;
t205 = -t243 * t301 + t271 * t303;
t204 = qJD(6) * t244 + t247;
t203 = qJD(6) * t242 + t246;
t198 = rSges(4,1) * t251 + rSges(4,2) * t250 + rSges(4,3) * t273;
t197 = rSges(4,1) * t249 + rSges(4,2) * t248 + rSges(4,3) * t271;
t196 = Icges(4,1) * t251 + Icges(4,4) * t250 + Icges(4,5) * t273;
t195 = Icges(4,1) * t249 + Icges(4,4) * t248 + Icges(4,5) * t271;
t194 = Icges(4,4) * t251 + Icges(4,2) * t250 + Icges(4,6) * t273;
t193 = Icges(4,4) * t249 + Icges(4,2) * t248 + Icges(4,6) * t271;
t188 = rSges(5,1) * t245 - rSges(5,2) * t244 + rSges(5,3) * t273;
t187 = rSges(5,1) * t243 - rSges(5,2) * t242 + rSges(5,3) * t271;
t186 = Icges(5,1) * t245 - Icges(5,4) * t244 + Icges(5,5) * t273;
t185 = Icges(5,1) * t243 - Icges(5,4) * t242 + Icges(5,5) * t271;
t179 = rSges(6,1) * t241 + rSges(6,2) * t240 + rSges(6,3) * t260;
t178 = Icges(6,1) * t241 + Icges(6,4) * t240 + Icges(6,5) * t260;
t177 = Icges(6,4) * t241 + Icges(6,2) * t240 + Icges(6,6) * t260;
t174 = rSges(7,1) * t235 + rSges(7,2) * t234 + rSges(7,3) * t260;
t173 = Icges(7,1) * t235 + Icges(7,4) * t234 + Icges(7,5) * t260;
t172 = Icges(7,4) * t235 + Icges(7,2) * t234 + Icges(7,6) * t260;
t171 = Icges(7,5) * t235 + Icges(7,6) * t234 + Icges(7,3) * t260;
t169 = -pkin(5) * t336 + pkin(10) * t260 + t261 * t360;
t167 = -t231 * t285 + t262 * t283 + t331;
t166 = t232 * t285 - t262 * t284 + t329;
t165 = t231 * t284 - t232 * t283 + t330;
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
t148 = pkin(5) * t356 + pkin(10) * t244 + t245 * t360;
t147 = pkin(5) * t357 + pkin(10) * t242 + t243 * t360;
t146 = -t197 * t267 + t224 * t246 + t328;
t145 = t198 * t267 - t224 * t247 + t325;
t144 = t197 * t247 - t198 * t246 + t327;
t143 = t216 * t246 + (-t187 - t189) * t267 + t326;
t142 = t188 * t267 + (-t216 - t223) * t247 + t323;
t141 = t187 * t247 + (-t188 - t190) * t246 + t322;
t140 = t179 * t246 + (-t163 + t344) * t267 + t324;
t139 = t164 * t267 + (-t179 + t342) * t247 + t321;
t138 = t163 * t247 + (-t164 + t343) * t246 + t320;
t137 = t324 - t155 * t233 + t169 * t246 + t174 * t203 + (-t147 + t344) * t267;
t136 = t148 * t267 + t156 * t233 - t174 * t204 + (-t169 + t342) * t247 + t321;
t135 = t147 * t247 + t155 * t204 - t156 * t203 + (-t148 + t343) * t246 + t320;
t1 = m(7) * (t135 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + m(6) * (t138 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + m(5) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(4) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(3) * (t165 ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + t285 * ((t225 * t283 + t226 * t284 + t257 * t285) * t311 + ((t228 * t318 + t230 * t315) * t284 + (t227 * t318 + t229 * t315) * t283 + (t258 * t318 + t259 * t315) * t285) * t309) / 0.2e1 + t203 * ((t150 * t242 + t152 * t205 + t154 * t206) * t204 + (t242 * t149 + t205 * t151 + t206 * t153) * t203 + (t171 * t242 + t172 * t205 + t173 * t206) * t233) / 0.2e1 + t204 * ((t244 * t150 + t207 * t152 + t208 * t154) * t204 + (t149 * t244 + t151 * t207 + t153 * t208) * t203 + (t171 * t244 + t172 * t207 + t173 * t208) * t233) / 0.2e1 + m(2) * (t252 ^ 2 + t254 ^ 2 + t255 ^ 2) / 0.2e1 + t233 * ((t150 * t260 + t152 * t234 + t154 * t235) * t204 + (t149 * t260 + t151 * t234 + t153 * t235) * t203 + (t171 * t260 + t234 * t172 + t173 * t235) * t233) / 0.2e1 + m(1) * (t279 ^ 2 + t280 ^ 2 + t281 ^ 2) / 0.2e1 + t283 * ((-t226 * t351 - t271 * t228 + t272 * t230) * t284 + (-t225 * t351 - t271 * t227 + t272 * t229) * t283 + (-t257 * t351 - t271 * t258 + t272 * t259) * t285) / 0.2e1 + t284 * ((t226 * t354 - t273 * t228 + t274 * t230) * t284 + (t225 * t354 - t227 * t273 + t229 * t274) * t283 + (t257 * t354 - t258 * t273 + t259 * t274) * t285) / 0.2e1 + ((-t316 * t288 + t290 * t319 + Icges(1,4)) * V_base(5) + (-t316 * t289 + t319 * t291 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t319 * t288 + t316 * t290 + Icges(1,2)) * V_base(5) + (t289 * t319 + t316 * t291 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t177 * t209 + t178 * t210 + t215 * t243 + t221 * t248 + t222 * t249 + t242 * t369 + t271 * t366) * t267 + (t160 * t209 + t162 * t210 + t186 * t243 + t194 * t248 + t196 * t249 + t242 * t370 + t271 * t367) * t247 + (t209 * t159 + t161 * t210 + t243 * t185 + t248 * t193 + t249 * t195 + t371 * t242 + t368 * t271) * t246) * t246 / 0.2e1 + ((t177 * t211 + t178 * t212 + t215 * t245 + t221 * t250 + t222 * t251 + t244 * t369 + t273 * t366) * t267 + (t160 * t211 + t212 * t162 + t245 * t186 + t250 * t194 + t251 * t196 + t370 * t244 + t367 * t273) * t247 + (t159 * t211 + t161 * t212 + t185 * t245 + t193 * t250 + t195 * t251 + t244 * t371 + t368 * t273) * t246) * t247 / 0.2e1 + ((t240 * t177 + t241 * t178 + t261 * t215 + t269 * t221 + t270 * t222 + t369 * t260 - t366 * t352) * t267 + (t160 * t240 + t162 * t241 + t186 * t261 + t194 * t269 + t196 * t270 + t260 * t370 - t352 * t367) * t247 + (t159 * t240 + t161 * t241 + t185 * t261 + t193 * t269 + t195 * t270 + t260 * t371 - t368 * t352) * t246) * t267 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t316 + Icges(2,6) * t319) * V_base(5) + (Icges(2,5) * t319 - Icges(2,6) * t316) * V_base(4) + Icges(2,3) * t304 / 0.2e1) * t304;
T  = t1;
