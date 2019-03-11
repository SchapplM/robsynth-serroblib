% Calculate kinetic energy for
% S6RRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRPR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:28:52
% EndTime: 2019-03-09 10:28:57
% DurationCPUTime: 4.68s
% Computational Cost: add. (3885->434), mult. (8792->610), div. (0->0), fcn. (11044->14), ass. (0->192)
t382 = Icges(5,2) + Icges(6,3);
t322 = cos(pkin(6));
t325 = sin(qJ(2));
t362 = sin(pkin(11));
t363 = cos(pkin(11));
t369 = cos(qJ(2));
t336 = t325 * t363 + t362 * t369;
t276 = t336 * t322;
t285 = -t325 * t362 + t363 * t369;
t326 = sin(qJ(1));
t327 = cos(qJ(1));
t257 = t276 * t327 + t285 * t326;
t324 = sin(qJ(4));
t320 = sin(pkin(6));
t356 = t320 * t327;
t368 = cos(qJ(4));
t235 = t257 * t368 - t324 * t356;
t335 = t322 * t285;
t256 = -t326 * t336 + t327 * t335;
t319 = sin(pkin(12));
t321 = cos(pkin(12));
t202 = -t235 * t319 - t256 * t321;
t360 = t256 * t319;
t203 = t235 * t321 - t360;
t347 = t320 * t368;
t234 = t257 * t324 + t327 * t347;
t381 = -Icges(5,4) * t235 + Icges(6,5) * t203 + Icges(5,6) * t256 + Icges(6,6) * t202 + t234 * t382;
t259 = -t276 * t326 + t285 * t327;
t357 = t320 * t326;
t237 = t259 * t368 + t324 * t357;
t258 = -t326 * t335 - t327 * t336;
t204 = -t237 * t319 - t258 * t321;
t359 = t258 * t319;
t205 = t237 * t321 - t359;
t236 = t259 * t324 - t326 * t347;
t380 = -Icges(5,4) * t237 + Icges(6,5) * t205 + Icges(5,6) * t258 + Icges(6,6) * t204 + t236 * t382;
t275 = t336 * t320;
t263 = t275 * t368 + t322 * t324;
t274 = t285 * t320;
t226 = -t263 * t319 - t274 * t321;
t358 = t274 * t319;
t227 = t263 * t321 - t358;
t262 = t275 * t324 - t322 * t368;
t379 = -Icges(5,4) * t263 + Icges(6,5) * t227 + Icges(5,6) * t274 + Icges(6,6) * t226 + t262 * t382;
t208 = Icges(4,5) * t257 + Icges(4,6) * t256 - Icges(4,3) * t356;
t345 = t327 * t369;
t355 = t326 * t325;
t279 = t322 * t345 - t355;
t346 = t326 * t369;
t354 = t327 * t325;
t280 = t322 * t354 + t346;
t243 = Icges(3,5) * t280 + Icges(3,6) * t279 - Icges(3,3) * t356;
t378 = t208 + t243;
t209 = Icges(4,5) * t259 + Icges(4,6) * t258 + Icges(4,3) * t357;
t281 = -t322 * t346 - t354;
t282 = -t322 * t355 + t345;
t244 = Icges(3,5) * t282 + Icges(3,6) * t281 + Icges(3,3) * t357;
t377 = t209 + t244;
t239 = Icges(4,5) * t275 + Icges(4,6) * t274 + Icges(4,3) * t322;
t270 = Icges(3,3) * t322 + (Icges(3,5) * t325 + Icges(3,6) * t369) * t320;
t376 = t239 + t270;
t367 = pkin(2) * t325;
t366 = pkin(8) * t322;
t365 = pkin(2) * t369;
t364 = pkin(5) * t321;
t361 = Icges(2,4) * t326;
t352 = qJD(2) * t320;
t351 = qJD(3) * t320;
t350 = V_base(5) * pkin(7) + V_base(1);
t293 = t326 * t352 + V_base(4);
t315 = V_base(6) + qJD(1);
t344 = -qJ(3) * t320 + t322 * t367;
t233 = -qJD(4) * t258 + t293;
t294 = qJD(2) * t322 + t315;
t261 = -qJD(4) * t274 + t294;
t292 = -t327 * t352 + V_base(5);
t287 = pkin(1) * t326 - pkin(8) * t356;
t341 = -t287 * t315 + t366 * V_base(5) + t350;
t288 = pkin(1) * t327 + pkin(8) * t357;
t340 = t287 * V_base(4) - t288 * V_base(5) + V_base(3);
t232 = -qJD(4) * t256 + t292;
t286 = qJ(3) * t322 + t320 * t367;
t339 = t286 * t292 + t326 * t351 + t341;
t254 = t326 * t365 + t327 * t344;
t338 = qJD(3) * t322 + t254 * t293 + t340;
t337 = t315 * t288 + V_base(2) + (-pkin(7) - t366) * V_base(4);
t220 = pkin(3) * t257 - pkin(9) * t256;
t251 = pkin(3) * t275 - pkin(9) * t274;
t334 = t292 * t251 + (-t220 - t254) * t294 + t339;
t221 = pkin(3) * t259 - pkin(9) * t258;
t255 = -t326 * t344 + t327 * t365;
t333 = t293 * t220 + (-t221 - t255) * t292 + t338;
t332 = t255 * t294 - t327 * t351 + t337;
t222 = pkin(4) * t263 + qJ(5) * t262;
t331 = qJD(5) * t236 + t222 * t232 + t334;
t194 = pkin(4) * t235 + qJ(5) * t234;
t330 = qJD(5) * t262 + t194 * t233 + t333;
t329 = t294 * t221 + (-t251 - t286) * t293 + t332;
t195 = pkin(4) * t237 + qJ(5) * t236;
t328 = qJD(5) * t234 + t195 * t261 + t329;
t318 = pkin(12) + qJ(6);
t316 = Icges(2,4) * t327;
t314 = cos(t318);
t313 = sin(t318);
t302 = rSges(2,1) * t327 - rSges(2,2) * t326;
t301 = rSges(2,1) * t326 + rSges(2,2) * t327;
t300 = Icges(2,1) * t327 - t361;
t299 = Icges(2,1) * t326 + t316;
t298 = -Icges(2,2) * t326 + t316;
t297 = Icges(2,2) * t327 + t361;
t291 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t290 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t289 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t273 = t322 * rSges(3,3) + (rSges(3,1) * t325 + rSges(3,2) * t369) * t320;
t272 = Icges(3,5) * t322 + (Icges(3,1) * t325 + Icges(3,4) * t369) * t320;
t271 = Icges(3,6) * t322 + (Icges(3,4) * t325 + Icges(3,2) * t369) * t320;
t266 = V_base(5) * rSges(2,3) - t301 * t315 + t350;
t265 = t302 * t315 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t264 = t301 * V_base(4) - t302 * V_base(5) + V_base(3);
t250 = rSges(3,1) * t282 + rSges(3,2) * t281 + rSges(3,3) * t357;
t249 = rSges(3,1) * t280 + rSges(3,2) * t279 - rSges(3,3) * t356;
t248 = Icges(3,1) * t282 + Icges(3,4) * t281 + Icges(3,5) * t357;
t247 = Icges(3,1) * t280 + Icges(3,4) * t279 - Icges(3,5) * t356;
t246 = Icges(3,4) * t282 + Icges(3,2) * t281 + Icges(3,6) * t357;
t245 = Icges(3,4) * t280 + Icges(3,2) * t279 - Icges(3,6) * t356;
t242 = rSges(4,1) * t275 + rSges(4,2) * t274 + rSges(4,3) * t322;
t241 = Icges(4,1) * t275 + Icges(4,4) * t274 + Icges(4,5) * t322;
t240 = Icges(4,4) * t275 + Icges(4,2) * t274 + Icges(4,6) * t322;
t225 = t263 * t314 - t274 * t313;
t224 = -t263 * t313 - t274 * t314;
t223 = qJD(6) * t262 + t261;
t219 = rSges(5,1) * t263 - rSges(5,2) * t262 - rSges(5,3) * t274;
t218 = Icges(5,1) * t263 - Icges(5,4) * t262 - Icges(5,5) * t274;
t216 = Icges(5,5) * t263 - Icges(5,6) * t262 - Icges(5,3) * t274;
t215 = rSges(4,1) * t259 + rSges(4,2) * t258 + rSges(4,3) * t357;
t214 = rSges(4,1) * t257 + rSges(4,2) * t256 - rSges(4,3) * t356;
t213 = Icges(4,1) * t259 + Icges(4,4) * t258 + Icges(4,5) * t357;
t212 = Icges(4,1) * t257 + Icges(4,4) * t256 - Icges(4,5) * t356;
t211 = Icges(4,4) * t259 + Icges(4,2) * t258 + Icges(4,6) * t357;
t210 = Icges(4,4) * t257 + Icges(4,2) * t256 - Icges(4,6) * t356;
t201 = t237 * t314 - t258 * t313;
t200 = -t237 * t313 - t258 * t314;
t199 = t235 * t314 - t256 * t313;
t198 = -t235 * t313 - t256 * t314;
t197 = qJD(6) * t236 + t233;
t196 = qJD(6) * t234 + t232;
t193 = -t249 * t294 + t273 * t292 + t341;
t192 = t250 * t294 - t273 * t293 + t337;
t189 = t249 * t293 - t250 * t292 + t340;
t188 = rSges(6,1) * t227 + rSges(6,2) * t226 + rSges(6,3) * t262;
t187 = Icges(6,1) * t227 + Icges(6,4) * t226 + Icges(6,5) * t262;
t186 = Icges(6,4) * t227 + Icges(6,2) * t226 + Icges(6,6) * t262;
t184 = rSges(5,1) * t237 - rSges(5,2) * t236 - rSges(5,3) * t258;
t183 = rSges(5,1) * t235 - rSges(5,2) * t234 - rSges(5,3) * t256;
t182 = Icges(5,1) * t237 - Icges(5,4) * t236 - Icges(5,5) * t258;
t181 = Icges(5,1) * t235 - Icges(5,4) * t234 - Icges(5,5) * t256;
t178 = Icges(5,5) * t237 - Icges(5,6) * t236 - Icges(5,3) * t258;
t177 = Icges(5,5) * t235 - Icges(5,6) * t234 - Icges(5,3) * t256;
t176 = -pkin(5) * t358 + pkin(10) * t262 + t263 * t364;
t175 = rSges(7,1) * t225 + rSges(7,2) * t224 + rSges(7,3) * t262;
t174 = Icges(7,1) * t225 + Icges(7,4) * t224 + Icges(7,5) * t262;
t173 = Icges(7,4) * t225 + Icges(7,2) * t224 + Icges(7,6) * t262;
t172 = Icges(7,5) * t225 + Icges(7,6) * t224 + Icges(7,3) * t262;
t170 = rSges(6,1) * t205 + rSges(6,2) * t204 + rSges(6,3) * t236;
t169 = rSges(6,1) * t203 + rSges(6,2) * t202 + rSges(6,3) * t234;
t168 = Icges(6,1) * t205 + Icges(6,4) * t204 + Icges(6,5) * t236;
t167 = Icges(6,1) * t203 + Icges(6,4) * t202 + Icges(6,5) * t234;
t166 = Icges(6,4) * t205 + Icges(6,2) * t204 + Icges(6,6) * t236;
t165 = Icges(6,4) * t203 + Icges(6,2) * t202 + Icges(6,6) * t234;
t162 = rSges(7,1) * t201 + rSges(7,2) * t200 + rSges(7,3) * t236;
t161 = rSges(7,1) * t199 + rSges(7,2) * t198 + rSges(7,3) * t234;
t160 = Icges(7,1) * t201 + Icges(7,4) * t200 + Icges(7,5) * t236;
t159 = Icges(7,1) * t199 + Icges(7,4) * t198 + Icges(7,5) * t234;
t158 = Icges(7,4) * t201 + Icges(7,2) * t200 + Icges(7,6) * t236;
t157 = Icges(7,4) * t199 + Icges(7,2) * t198 + Icges(7,6) * t234;
t156 = Icges(7,5) * t201 + Icges(7,6) * t200 + Icges(7,3) * t236;
t155 = Icges(7,5) * t199 + Icges(7,6) * t198 + Icges(7,3) * t234;
t154 = -pkin(5) * t359 + pkin(10) * t236 + t237 * t364;
t153 = -pkin(5) * t360 + pkin(10) * t234 + t235 * t364;
t152 = t242 * t292 + (-t214 - t254) * t294 + t339;
t151 = t294 * t215 + (-t242 - t286) * t293 + t332;
t150 = t214 * t293 + (-t215 - t255) * t292 + t338;
t149 = -t183 * t261 + t219 * t232 + t334;
t148 = t184 * t261 - t219 * t233 + t329;
t147 = t183 * t233 - t184 * t232 + t333;
t146 = t188 * t232 + (-t169 - t194) * t261 + t331;
t145 = t261 * t170 + (-t188 - t222) * t233 + t328;
t144 = t169 * t233 + (-t170 - t195) * t232 + t330;
t143 = -t161 * t223 + t175 * t196 + t176 * t232 + (-t153 - t194) * t261 + t331;
t142 = t261 * t154 + t223 * t162 - t197 * t175 + (-t176 - t222) * t233 + t328;
t141 = t153 * t233 + t161 * t197 - t162 * t196 + (-t154 - t195) * t232 + t330;
t1 = m(7) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(6) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(5) * (t147 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + m(4) * (t150 ^ 2 + t151 ^ 2 + t152 ^ 2) / 0.2e1 + m(3) * (t189 ^ 2 + t192 ^ 2 + t193 ^ 2) / 0.2e1 + t196 * ((t156 * t234 + t158 * t198 + t160 * t199) * t197 + (t155 * t234 + t157 * t198 + t159 * t199) * t196 + (t172 * t234 + t173 * t198 + t174 * t199) * t223) / 0.2e1 + t197 * ((t156 * t236 + t158 * t200 + t160 * t201) * t197 + (t155 * t236 + t157 * t200 + t159 * t201) * t196 + (t172 * t236 + t173 * t200 + t174 * t201) * t223) / 0.2e1 + t223 * ((t156 * t262 + t158 * t224 + t160 * t225) * t197 + (t155 * t262 + t157 * t224 + t159 * t225) * t196 + (t172 * t262 + t173 * t224 + t174 * t225) * t223) / 0.2e1 + m(2) * (t264 ^ 2 + t265 ^ 2 + t266 ^ 2) / 0.2e1 + m(1) * (t289 ^ 2 + t290 ^ 2 + t291 ^ 2) / 0.2e1 + ((t186 * t202 + t187 * t203 - t216 * t256 + t218 * t235 + t234 * t379) * t261 + (t166 * t202 + t168 * t203 - t178 * t256 + t182 * t235 + t234 * t380) * t233 + (t165 * t202 + t167 * t203 - t177 * t256 + t181 * t235 + t381 * t234) * t232) * t232 / 0.2e1 + ((t186 * t204 + t187 * t205 - t216 * t258 + t218 * t237 + t236 * t379) * t261 + (t166 * t204 + t168 * t205 - t178 * t258 + t182 * t237 + t380 * t236) * t233 + (t165 * t204 + t167 * t205 - t177 * t258 + t181 * t237 + t236 * t381) * t232) * t233 / 0.2e1 + ((t186 * t226 + t187 * t227 - t216 * t274 + t218 * t263 + t379 * t262) * t261 + (t166 * t226 + t168 * t227 - t178 * t274 + t182 * t263 + t262 * t380) * t233 + (t165 * t226 + t167 * t227 - t177 * t274 + t181 * t263 + t262 * t381) * t232) * t261 / 0.2e1 + ((t240 * t256 + t241 * t257 + t271 * t279 + t272 * t280 - t356 * t376) * t294 + (t256 * t211 + t257 * t213 + t279 * t246 + t280 * t248 - t356 * t377) * t293 + (t256 * t210 + t257 * t212 + t279 * t245 + t280 * t247 - t356 * t378) * t292) * t292 / 0.2e1 + ((t240 * t258 + t241 * t259 + t271 * t281 + t272 * t282 + t357 * t376) * t294 + (t211 * t258 + t213 * t259 + t246 * t281 + t248 * t282 + t357 * t377) * t293 + (t210 * t258 + t212 * t259 + t245 * t281 + t247 * t282 + t357 * t378) * t292) * t293 / 0.2e1 + ((t243 * t292 + t244 * t293 + t270 * t294) * t322 + ((t246 * t369 + t248 * t325) * t293 + (t245 * t369 + t247 * t325) * t292 + (t271 * t369 + t272 * t325) * t294) * t320 + (t209 * t322 + t211 * t274 + t213 * t275) * t293 + (t208 * t322 + t210 * t274 + t212 * t275) * t292 + (t239 * t322 + t240 * t274 + t241 * t275) * t294) * t294 / 0.2e1 + ((-t297 * t326 + t299 * t327 + Icges(1,4)) * V_base(5) + (-t326 * t298 + t300 * t327 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t297 * t327 + t326 * t299 + Icges(1,2)) * V_base(5) + (t298 * t327 + t300 * t326 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t326 + Icges(2,6) * t327) * V_base(5) + (Icges(2,5) * t327 - Icges(2,6) * t326) * V_base(4) + Icges(2,3) * t315 / 0.2e1) * t315;
T  = t1;
