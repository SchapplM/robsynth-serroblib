% Calculate kinetic energy for
% S6PRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPPR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRPPR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPPR2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:04:02
% EndTime: 2019-03-08 21:04:05
% DurationCPUTime: 3.79s
% Computational Cost: add. (2992->390), mult. (5246->542), div. (0->0), fcn. (6124->12), ass. (0->174)
t365 = Icges(5,1) + Icges(6,2);
t364 = -Icges(5,4) - Icges(6,6);
t363 = Icges(6,4) - Icges(5,5);
t362 = Icges(6,5) - Icges(5,6);
t361 = Icges(5,2) + Icges(6,3);
t360 = Icges(6,1) + Icges(4,3) + Icges(5,3);
t295 = sin(pkin(10));
t297 = cos(pkin(10));
t305 = cos(qJ(2));
t298 = cos(pkin(6));
t302 = sin(qJ(2));
t335 = t298 * t302;
t260 = t295 * t305 + t297 * t335;
t296 = sin(pkin(6));
t329 = qJ(3) + pkin(11);
t321 = cos(t329);
t319 = t296 * t321;
t320 = sin(t329);
t231 = t260 * t320 + t297 * t319;
t318 = t296 * t320;
t232 = t260 * t321 - t297 * t318;
t334 = t298 * t305;
t259 = t295 * t302 - t297 * t334;
t359 = t361 * t231 + t364 * t232 + t362 * t259;
t262 = -t295 * t335 + t297 * t305;
t233 = t262 * t320 - t295 * t319;
t234 = t262 * t321 + t295 * t318;
t261 = t295 * t334 + t297 * t302;
t358 = t361 * t233 + t364 * t234 + t362 * t261;
t357 = t364 * t231 + t365 * t232 - t363 * t259;
t356 = t364 * t233 + t365 * t234 - t363 * t261;
t252 = -t298 * t321 + t302 * t318;
t253 = t298 * t320 + t302 * t319;
t337 = t296 * t305;
t355 = t361 * t252 + t364 * t253 - t362 * t337;
t354 = t364 * t252 + t365 * t253 + t363 * t337;
t301 = sin(qJ(3));
t304 = cos(qJ(3));
t338 = t296 * t304;
t241 = -t260 * t301 - t297 * t338;
t339 = t296 * t301;
t322 = t297 * t339;
t242 = t260 * t304 - t322;
t351 = Icges(4,5) * t242 + Icges(4,6) * t241 + t362 * t231 - t363 * t232 + t360 * t259;
t243 = -t262 * t301 + t295 * t338;
t323 = t295 * t339;
t244 = t262 * t304 + t323;
t350 = Icges(4,5) * t244 + Icges(4,6) * t243 + t362 * t233 - t363 * t234 + t360 * t261;
t266 = t298 * t304 - t302 * t339;
t336 = t298 * t301;
t267 = t302 * t338 + t336;
t349 = Icges(4,5) * t267 + Icges(4,6) * t266 + t362 * t252 - t363 * t253 - t360 * t337;
t345 = pkin(7) * t298;
t344 = pkin(3) * t304;
t342 = Icges(2,4) * t295;
t341 = t295 * t296;
t340 = t296 * t297;
t180 = -pkin(3) * t322 + qJ(4) * t259 + t260 * t344;
t192 = pkin(4) * t232 + qJ(5) * t231;
t333 = -t180 - t192;
t181 = pkin(3) * t323 + qJ(4) * t261 + t262 * t344;
t193 = pkin(4) * t234 + qJ(5) * t233;
t332 = -t181 - t193;
t212 = pkin(4) * t253 + qJ(5) * t252;
t224 = pkin(3) * t336 + (-qJ(4) * t305 + t302 * t344) * t296;
t331 = -t212 - t224;
t330 = qJD(2) * t296;
t328 = V_base(5) * qJ(1) + V_base(1);
t324 = qJD(1) + V_base(3);
t277 = t295 * t330 + V_base(4);
t288 = qJD(2) * t298 + V_base(6);
t240 = qJD(3) * t261 + t277;
t276 = -t297 * t330 + V_base(5);
t239 = qJD(3) * t259 + t276;
t263 = -qJD(3) * t337 + t288;
t269 = pkin(1) * t295 - pkin(7) * t340;
t317 = -t269 * V_base(6) + V_base(5) * t345 + t328;
t270 = pkin(1) * t297 + pkin(7) * t341;
t316 = V_base(4) * t269 - V_base(5) * t270 + t324;
t315 = V_base(6) * t270 + V_base(2) + (-qJ(1) - t345) * V_base(4);
t227 = pkin(2) * t260 + pkin(8) * t259;
t268 = (pkin(2) * t302 - pkin(8) * t305) * t296;
t314 = -t227 * t288 + t276 * t268 + t317;
t228 = pkin(2) * t262 + pkin(8) * t261;
t313 = t277 * t227 - t276 * t228 + t316;
t312 = t288 * t228 - t268 * t277 + t315;
t311 = qJD(4) * t261 + t239 * t224 + t314;
t310 = qJD(4) * t259 + t263 * t181 + t312;
t309 = qJD(5) * t233 + t239 * t212 + t311;
t308 = -qJD(4) * t337 + t240 * t180 + t313;
t307 = qJD(5) * t231 + t263 * t193 + t310;
t306 = qJD(5) * t252 + t240 * t192 + t308;
t303 = cos(qJ(6));
t300 = sin(qJ(6));
t293 = Icges(2,4) * t297;
t285 = rSges(2,1) * t297 - rSges(2,2) * t295;
t284 = rSges(2,1) * t295 + rSges(2,2) * t297;
t283 = Icges(2,1) * t297 - t342;
t282 = Icges(2,1) * t295 + t293;
t281 = -Icges(2,2) * t295 + t293;
t280 = Icges(2,2) * t297 + t342;
t275 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t274 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t273 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t254 = t298 * rSges(3,3) + (rSges(3,1) * t302 + rSges(3,2) * t305) * t296;
t251 = Icges(3,5) * t298 + (Icges(3,1) * t302 + Icges(3,4) * t305) * t296;
t250 = Icges(3,6) * t298 + (Icges(3,4) * t302 + Icges(3,2) * t305) * t296;
t249 = Icges(3,3) * t298 + (Icges(3,5) * t302 + Icges(3,6) * t305) * t296;
t247 = V_base(5) * rSges(2,3) - t284 * V_base(6) + t328;
t246 = t285 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t238 = -pkin(5) * t337 + t253 * pkin(9);
t237 = t284 * V_base(4) - t285 * V_base(5) + t324;
t236 = t252 * t300 - t303 * t337;
t235 = t252 * t303 + t300 * t337;
t226 = qJD(6) * t253 + t263;
t225 = t267 * rSges(4,1) + t266 * rSges(4,2) - rSges(4,3) * t337;
t223 = Icges(4,1) * t267 + Icges(4,4) * t266 - Icges(4,5) * t337;
t222 = Icges(4,4) * t267 + Icges(4,2) * t266 - Icges(4,6) * t337;
t220 = rSges(3,1) * t262 - rSges(3,2) * t261 + rSges(3,3) * t341;
t219 = rSges(3,1) * t260 - rSges(3,2) * t259 - rSges(3,3) * t340;
t218 = Icges(3,1) * t262 - Icges(3,4) * t261 + Icges(3,5) * t341;
t217 = Icges(3,1) * t260 - Icges(3,4) * t259 - Icges(3,5) * t340;
t216 = Icges(3,4) * t262 - Icges(3,2) * t261 + Icges(3,6) * t341;
t215 = Icges(3,4) * t260 - Icges(3,2) * t259 - Icges(3,6) * t340;
t214 = Icges(3,5) * t262 - Icges(3,6) * t261 + Icges(3,3) * t341;
t213 = Icges(3,5) * t260 - Icges(3,6) * t259 - Icges(3,3) * t340;
t209 = t253 * rSges(5,1) - t252 * rSges(5,2) - rSges(5,3) * t337;
t208 = -rSges(6,1) * t337 - t253 * rSges(6,2) + t252 * rSges(6,3);
t201 = pkin(5) * t261 + pkin(9) * t234;
t200 = pkin(5) * t259 + pkin(9) * t232;
t199 = t233 * t300 + t261 * t303;
t198 = t233 * t303 - t261 * t300;
t197 = t231 * t300 + t259 * t303;
t196 = t231 * t303 - t259 * t300;
t195 = qJD(6) * t234 + t240;
t194 = qJD(6) * t232 + t239;
t189 = rSges(4,1) * t244 + rSges(4,2) * t243 + rSges(4,3) * t261;
t188 = rSges(4,1) * t242 + rSges(4,2) * t241 + rSges(4,3) * t259;
t187 = Icges(4,1) * t244 + Icges(4,4) * t243 + Icges(4,5) * t261;
t186 = Icges(4,1) * t242 + Icges(4,4) * t241 + Icges(4,5) * t259;
t185 = Icges(4,4) * t244 + Icges(4,2) * t243 + Icges(4,6) * t261;
t184 = Icges(4,4) * t242 + Icges(4,2) * t241 + Icges(4,6) * t259;
t178 = rSges(5,1) * t234 - rSges(5,2) * t233 + rSges(5,3) * t261;
t177 = rSges(5,1) * t232 - rSges(5,2) * t231 + rSges(5,3) * t259;
t176 = rSges(6,1) * t261 - rSges(6,2) * t234 + rSges(6,3) * t233;
t175 = rSges(6,1) * t259 - rSges(6,2) * t232 + rSges(6,3) * t231;
t162 = rSges(7,1) * t236 + rSges(7,2) * t235 + rSges(7,3) * t253;
t161 = Icges(7,1) * t236 + Icges(7,4) * t235 + Icges(7,5) * t253;
t160 = Icges(7,4) * t236 + Icges(7,2) * t235 + Icges(7,6) * t253;
t159 = Icges(7,5) * t236 + Icges(7,6) * t235 + Icges(7,3) * t253;
t156 = -t219 * t288 + t254 * t276 + t317;
t155 = t220 * t288 - t254 * t277 + t315;
t153 = t219 * t277 - t220 * t276 + t316;
t152 = rSges(7,1) * t199 + rSges(7,2) * t198 + rSges(7,3) * t234;
t151 = rSges(7,1) * t197 + rSges(7,2) * t196 + rSges(7,3) * t232;
t150 = Icges(7,1) * t199 + Icges(7,4) * t198 + Icges(7,5) * t234;
t149 = Icges(7,1) * t197 + Icges(7,4) * t196 + Icges(7,5) * t232;
t148 = Icges(7,4) * t199 + Icges(7,2) * t198 + Icges(7,6) * t234;
t147 = Icges(7,4) * t197 + Icges(7,2) * t196 + Icges(7,6) * t232;
t146 = Icges(7,5) * t199 + Icges(7,6) * t198 + Icges(7,3) * t234;
t145 = Icges(7,5) * t197 + Icges(7,6) * t196 + Icges(7,3) * t232;
t144 = -t188 * t263 + t225 * t239 + t314;
t143 = t189 * t263 - t225 * t240 + t312;
t142 = t188 * t240 - t189 * t239 + t313;
t141 = t209 * t239 + (-t177 - t180) * t263 + t311;
t140 = t178 * t263 + (-t209 - t224) * t240 + t310;
t139 = t240 * t177 + (-t178 - t181) * t239 + t308;
t138 = t208 * t239 + (-t175 + t333) * t263 + t309;
t137 = t176 * t263 + (-t208 + t331) * t240 + t307;
t136 = t240 * t175 + (-t176 + t332) * t239 + t306;
t135 = t309 + (-t200 + t333) * t263 - t151 * t226 + t162 * t194 + t238 * t239;
t134 = t152 * t226 - t162 * t195 + t201 * t263 + (-t238 + t331) * t240 + t307;
t133 = t306 + (-t201 + t332) * t239 - t194 * t152 + t195 * t151 + t240 * t200;
t1 = t194 * ((t146 * t232 + t148 * t196 + t150 * t197) * t195 + (t232 * t145 + t196 * t147 + t197 * t149) * t194 + (t159 * t232 + t160 * t196 + t161 * t197) * t226) / 0.2e1 + m(3) * (t153 ^ 2 + t155 ^ 2 + t156 ^ 2) / 0.2e1 + m(5) * (t139 ^ 2 + t140 ^ 2 + t141 ^ 2) / 0.2e1 + m(4) * (t142 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(6) * (t136 ^ 2 + t137 ^ 2 + t138 ^ 2) / 0.2e1 + m(7) * (t133 ^ 2 + t134 ^ 2 + t135 ^ 2) / 0.2e1 + t276 * ((-t214 * t340 - t216 * t259 + t218 * t260) * t277 + (-t213 * t340 - t215 * t259 + t217 * t260) * t276 + (-t249 * t340 - t250 * t259 + t251 * t260) * t288) / 0.2e1 + t277 * ((t214 * t341 - t216 * t261 + t218 * t262) * t277 + (t213 * t341 - t215 * t261 + t217 * t262) * t276 + (t249 * t341 - t250 * t261 + t251 * t262) * t288) / 0.2e1 + t288 * ((t213 * t276 + t214 * t277 + t249 * t288) * t298 + ((t216 * t305 + t218 * t302) * t277 + (t215 * t305 + t217 * t302) * t276 + (t250 * t305 + t251 * t302) * t288) * t296) / 0.2e1 + t195 * ((t234 * t146 + t198 * t148 + t199 * t150) * t195 + (t145 * t234 + t147 * t198 + t149 * t199) * t194 + (t159 * t234 + t160 * t198 + t161 * t199) * t226) / 0.2e1 + m(2) * (t237 ^ 2 + t246 ^ 2 + t247 ^ 2) / 0.2e1 + t226 * ((t146 * t253 + t148 * t235 + t150 * t236) * t195 + (t145 * t253 + t147 * t235 + t149 * t236) * t194 + (t159 * t253 + t160 * t235 + t161 * t236) * t226) / 0.2e1 + m(1) * (t273 ^ 2 + t274 ^ 2 + t275 ^ 2) / 0.2e1 + ((-t280 * t295 + t282 * t297 + Icges(1,4)) * V_base(5) + (-t281 * t295 + t283 * t297 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t280 * t297 + t282 * t295 + Icges(1,2)) * V_base(5) + (t281 * t297 + t283 * t295 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t222 * t241 + t223 * t242 + t231 * t355 + t232 * t354 + t259 * t349) * t263 + (t185 * t241 + t187 * t242 + t231 * t358 + t232 * t356 + t259 * t350) * t240 + (t184 * t241 + t186 * t242 + t359 * t231 + t357 * t232 + t351 * t259) * t239) * t239 / 0.2e1 + ((t222 * t243 + t223 * t244 + t233 * t355 + t234 * t354 + t261 * t349) * t263 + (t185 * t243 + t187 * t244 + t358 * t233 + t356 * t234 + t350 * t261) * t240 + (t184 * t243 + t186 * t244 + t233 * t359 + t357 * t234 + t351 * t261) * t239) * t240 / 0.2e1 + ((t266 * t222 + t267 * t223 + t355 * t252 + t354 * t253 - t349 * t337) * t263 + (t266 * t185 + t267 * t187 + t252 * t358 + t253 * t356 - t337 * t350) * t240 + (t266 * t184 + t267 * t186 + t252 * t359 + t357 * t253 - t351 * t337) * t239) * t263 / 0.2e1 + ((Icges(2,5) * t297 - Icges(2,6) * t295 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t295 + Icges(2,6) * t297 + Icges(1,6)) * V_base(5) + (Icges(2,3) / 0.2e1 + Icges(1,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
