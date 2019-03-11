% Calculate kinetic energy for
% S6RRRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPPR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:47:40
% EndTime: 2019-03-09 15:47:44
% DurationCPUTime: 3.96s
% Computational Cost: add. (3052->390), mult. (5246->544), div. (0->0), fcn. (6124->12), ass. (0->175)
t363 = Icges(5,1) + Icges(6,2);
t362 = -Icges(5,4) - Icges(6,6);
t361 = Icges(6,4) - Icges(5,5);
t360 = Icges(6,5) - Icges(5,6);
t359 = Icges(5,2) + Icges(6,3);
t358 = Icges(6,1) + Icges(4,3) + Icges(5,3);
t297 = cos(pkin(6));
t302 = sin(qJ(1));
t305 = cos(qJ(2));
t334 = t302 * t305;
t301 = sin(qJ(2));
t306 = cos(qJ(1));
t335 = t301 * t306;
t264 = t297 * t335 + t334;
t296 = sin(pkin(6));
t328 = qJ(3) + pkin(11);
t322 = cos(t328);
t320 = t296 * t322;
t321 = sin(t328);
t233 = t264 * t321 + t306 * t320;
t319 = t296 * t321;
t234 = t264 * t322 - t306 * t319;
t333 = t305 * t306;
t336 = t301 * t302;
t263 = -t297 * t333 + t336;
t357 = t359 * t233 + t362 * t234 + t360 * t263;
t266 = -t297 * t336 + t333;
t235 = t266 * t321 - t302 * t320;
t236 = t266 * t322 + t302 * t319;
t265 = t297 * t334 + t335;
t356 = t359 * t235 + t362 * t236 + t360 * t265;
t355 = t362 * t233 + t363 * t234 - t361 * t263;
t354 = t362 * t235 + t363 * t236 - t361 * t265;
t252 = -t297 * t322 + t301 * t319;
t253 = t297 * t321 + t301 * t320;
t339 = t296 * t305;
t353 = t359 * t252 + t362 * t253 - t360 * t339;
t352 = t362 * t252 + t363 * t253 + t361 * t339;
t300 = sin(qJ(3));
t304 = cos(qJ(3));
t338 = t296 * t306;
t240 = -t264 * t300 - t304 * t338;
t323 = t300 * t338;
t241 = t264 * t304 - t323;
t351 = Icges(4,5) * t241 + Icges(4,6) * t240 + t360 * t233 - t361 * t234 + t358 * t263;
t340 = t296 * t304;
t242 = -t266 * t300 + t302 * t340;
t341 = t296 * t302;
t324 = t300 * t341;
t243 = t266 * t304 + t324;
t350 = Icges(4,5) * t243 + Icges(4,6) * t242 + t360 * t235 - t361 * t236 + t358 * t265;
t261 = -t296 * t300 * t301 + t297 * t304;
t337 = t297 * t300;
t262 = t301 * t340 + t337;
t349 = Icges(4,5) * t262 + Icges(4,6) * t261 + t360 * t252 - t361 * t253 - t358 * t339;
t345 = pkin(8) * t297;
t344 = pkin(3) * t304;
t342 = Icges(2,4) * t302;
t180 = -pkin(3) * t323 + qJ(4) * t263 + t264 * t344;
t192 = pkin(4) * t234 + qJ(5) * t233;
t332 = -t180 - t192;
t181 = pkin(3) * t324 + qJ(4) * t265 + t266 * t344;
t193 = pkin(4) * t236 + qJ(5) * t235;
t331 = -t181 - t193;
t212 = pkin(4) * t253 + qJ(5) * t252;
t216 = pkin(3) * t337 + (-qJ(4) * t305 + t301 * t344) * t296;
t330 = -t212 - t216;
t329 = qJD(2) * t296;
t327 = V_base(5) * pkin(7) + V_base(1);
t277 = t302 * t329 + V_base(4);
t293 = V_base(6) + qJD(1);
t239 = qJD(3) * t265 + t277;
t278 = qJD(2) * t297 + t293;
t276 = -t306 * t329 + V_base(5);
t269 = t302 * pkin(1) - pkin(8) * t338;
t318 = -t269 * t293 + V_base(5) * t345 + t327;
t270 = pkin(1) * t306 + pkin(8) * t341;
t317 = V_base(4) * t269 - t270 * V_base(5) + V_base(3);
t238 = qJD(3) * t263 + t276;
t259 = -qJD(3) * t339 + t278;
t316 = t293 * t270 + V_base(2) + (-pkin(7) - t345) * V_base(4);
t229 = t264 * pkin(2) + t263 * pkin(9);
t268 = (pkin(2) * t301 - pkin(9) * t305) * t296;
t315 = -t229 * t278 + t276 * t268 + t318;
t230 = pkin(2) * t266 + pkin(9) * t265;
t314 = t277 * t229 - t230 * t276 + t317;
t313 = qJD(4) * t265 + t238 * t216 + t315;
t312 = t278 * t230 - t268 * t277 + t316;
t311 = qJD(5) * t235 + t238 * t212 + t313;
t310 = qJD(4) * t263 + t259 * t181 + t312;
t309 = -qJD(4) * t339 + t239 * t180 + t314;
t308 = qJD(5) * t233 + t259 * t193 + t310;
t307 = qJD(5) * t252 + t239 * t192 + t309;
t303 = cos(qJ(6));
t299 = sin(qJ(6));
t294 = Icges(2,4) * t306;
t286 = rSges(2,1) * t306 - t302 * rSges(2,2);
t285 = t302 * rSges(2,1) + rSges(2,2) * t306;
t284 = Icges(2,1) * t306 - t342;
t283 = Icges(2,1) * t302 + t294;
t282 = -Icges(2,2) * t302 + t294;
t281 = Icges(2,2) * t306 + t342;
t273 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t272 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t271 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t254 = rSges(3,3) * t297 + (rSges(3,1) * t301 + rSges(3,2) * t305) * t296;
t251 = Icges(3,5) * t297 + (Icges(3,1) * t301 + Icges(3,4) * t305) * t296;
t250 = Icges(3,6) * t297 + (Icges(3,4) * t301 + Icges(3,2) * t305) * t296;
t249 = Icges(3,3) * t297 + (Icges(3,5) * t301 + Icges(3,6) * t305) * t296;
t247 = V_base(5) * rSges(2,3) - t285 * t293 + t327;
t246 = t286 * t293 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t244 = t285 * V_base(4) - t286 * V_base(5) + V_base(3);
t237 = -pkin(5) * t339 + pkin(10) * t253;
t232 = t252 * t299 - t303 * t339;
t231 = t252 * t303 + t299 * t339;
t226 = qJD(6) * t253 + t259;
t225 = rSges(3,1) * t266 - rSges(3,2) * t265 + rSges(3,3) * t341;
t224 = t264 * rSges(3,1) - t263 * rSges(3,2) - rSges(3,3) * t338;
t223 = Icges(3,1) * t266 - Icges(3,4) * t265 + Icges(3,5) * t341;
t222 = Icges(3,1) * t264 - Icges(3,4) * t263 - Icges(3,5) * t338;
t221 = Icges(3,4) * t266 - Icges(3,2) * t265 + Icges(3,6) * t341;
t220 = Icges(3,4) * t264 - Icges(3,2) * t263 - Icges(3,6) * t338;
t219 = Icges(3,5) * t266 - Icges(3,6) * t265 + Icges(3,3) * t341;
t218 = Icges(3,5) * t264 - Icges(3,6) * t263 - Icges(3,3) * t338;
t217 = rSges(4,1) * t262 + rSges(4,2) * t261 - rSges(4,3) * t339;
t215 = Icges(4,1) * t262 + Icges(4,4) * t261 - Icges(4,5) * t339;
t214 = Icges(4,4) * t262 + Icges(4,2) * t261 - Icges(4,6) * t339;
t209 = pkin(5) * t265 + pkin(10) * t236;
t208 = pkin(5) * t263 + pkin(10) * t234;
t207 = rSges(5,1) * t253 - rSges(5,2) * t252 - rSges(5,3) * t339;
t206 = -rSges(6,1) * t339 - rSges(6,2) * t253 + rSges(6,3) * t252;
t199 = t235 * t299 + t265 * t303;
t198 = t235 * t303 - t265 * t299;
t197 = t233 * t299 + t263 * t303;
t196 = t233 * t303 - t263 * t299;
t195 = qJD(6) * t236 + t239;
t194 = qJD(6) * t234 + t238;
t189 = rSges(4,1) * t243 + rSges(4,2) * t242 + rSges(4,3) * t265;
t188 = rSges(4,1) * t241 + rSges(4,2) * t240 + rSges(4,3) * t263;
t187 = Icges(4,1) * t243 + Icges(4,4) * t242 + Icges(4,5) * t265;
t186 = Icges(4,1) * t241 + Icges(4,4) * t240 + Icges(4,5) * t263;
t185 = Icges(4,4) * t243 + Icges(4,2) * t242 + Icges(4,6) * t265;
t184 = Icges(4,4) * t241 + Icges(4,2) * t240 + Icges(4,6) * t263;
t179 = rSges(5,1) * t236 - rSges(5,2) * t235 + rSges(5,3) * t265;
t178 = rSges(5,1) * t234 - rSges(5,2) * t233 + rSges(5,3) * t263;
t177 = rSges(6,1) * t265 - rSges(6,2) * t236 + rSges(6,3) * t235;
t176 = rSges(6,1) * t263 - rSges(6,2) * t234 + rSges(6,3) * t233;
t162 = rSges(7,1) * t232 + rSges(7,2) * t231 + rSges(7,3) * t253;
t161 = Icges(7,1) * t232 + Icges(7,4) * t231 + Icges(7,5) * t253;
t160 = Icges(7,4) * t232 + Icges(7,2) * t231 + Icges(7,6) * t253;
t159 = Icges(7,5) * t232 + Icges(7,6) * t231 + Icges(7,3) * t253;
t155 = -t224 * t278 + t254 * t276 + t318;
t154 = t225 * t278 - t254 * t277 + t316;
t153 = t224 * t277 - t225 * t276 + t317;
t152 = rSges(7,1) * t199 + rSges(7,2) * t198 + rSges(7,3) * t236;
t151 = rSges(7,1) * t197 + rSges(7,2) * t196 + rSges(7,3) * t234;
t150 = Icges(7,1) * t199 + Icges(7,4) * t198 + Icges(7,5) * t236;
t149 = Icges(7,1) * t197 + Icges(7,4) * t196 + Icges(7,5) * t234;
t148 = Icges(7,4) * t199 + Icges(7,2) * t198 + Icges(7,6) * t236;
t147 = Icges(7,4) * t197 + Icges(7,2) * t196 + Icges(7,6) * t234;
t146 = Icges(7,5) * t199 + Icges(7,6) * t198 + Icges(7,3) * t236;
t145 = Icges(7,5) * t197 + Icges(7,6) * t196 + Icges(7,3) * t234;
t144 = -t188 * t259 + t217 * t238 + t315;
t143 = t189 * t259 - t217 * t239 + t312;
t142 = t188 * t239 - t189 * t238 + t314;
t141 = t207 * t238 + (-t178 - t180) * t259 + t313;
t140 = t179 * t259 + (-t207 - t216) * t239 + t310;
t139 = t178 * t239 + (-t179 - t181) * t238 + t309;
t138 = t206 * t238 + (-t176 + t332) * t259 + t311;
t137 = t177 * t259 + (-t206 + t330) * t239 + t308;
t136 = t176 * t239 + (-t177 + t331) * t238 + t307;
t135 = t311 - t151 * t226 + t162 * t194 + t237 * t238 + (-t208 + t332) * t259;
t134 = t152 * t226 - t162 * t195 + t209 * t259 + (-t237 + t330) * t239 + t308;
t133 = t151 * t195 - t152 * t194 + t208 * t239 + (-t209 + t331) * t238 + t307;
t1 = m(3) * (t153 ^ 2 + t154 ^ 2 + t155 ^ 2) / 0.2e1 + m(4) * (t142 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(6) * (t136 ^ 2 + t137 ^ 2 + t138 ^ 2) / 0.2e1 + m(5) * (t139 ^ 2 + t140 ^ 2 + t141 ^ 2) / 0.2e1 + m(7) * (t133 ^ 2 + t134 ^ 2 + t135 ^ 2) / 0.2e1 + t277 * ((t219 * t341 - t221 * t265 + t223 * t266) * t277 + (t218 * t341 - t220 * t265 + t222 * t266) * t276 + (t249 * t341 - t250 * t265 + t251 * t266) * t278) / 0.2e1 + t276 * ((-t219 * t338 - t263 * t221 + t264 * t223) * t277 + (-t218 * t338 - t263 * t220 + t264 * t222) * t276 + (-t249 * t338 - t263 * t250 + t264 * t251) * t278) / 0.2e1 + t278 * ((t218 * t276 + t219 * t277 + t249 * t278) * t297 + ((t221 * t305 + t223 * t301) * t277 + (t220 * t305 + t222 * t301) * t276 + (t250 * t305 + t251 * t301) * t278) * t296) / 0.2e1 + t194 * ((t146 * t234 + t148 * t196 + t150 * t197) * t195 + (t234 * t145 + t196 * t147 + t197 * t149) * t194 + (t159 * t234 + t160 * t196 + t161 * t197) * t226) / 0.2e1 + t195 * ((t236 * t146 + t198 * t148 + t199 * t150) * t195 + (t145 * t236 + t147 * t198 + t149 * t199) * t194 + (t159 * t236 + t160 * t198 + t161 * t199) * t226) / 0.2e1 + m(2) * (t244 ^ 2 + t246 ^ 2 + t247 ^ 2) / 0.2e1 + t226 * ((t146 * t253 + t148 * t231 + t150 * t232) * t195 + (t145 * t253 + t147 * t231 + t149 * t232) * t194 + (t159 * t253 + t160 * t231 + t161 * t232) * t226) / 0.2e1 + m(1) * (t271 ^ 2 + t272 ^ 2 + t273 ^ 2) / 0.2e1 + ((-t302 * t281 + t283 * t306 + Icges(1,4)) * V_base(5) + (-t302 * t282 + t284 * t306 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t281 * t306 + t302 * t283 + Icges(1,2)) * V_base(5) + (t282 * t306 + t302 * t284 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t214 * t240 + t215 * t241 + t233 * t353 + t234 * t352 + t263 * t349) * t259 + (t185 * t240 + t187 * t241 + t233 * t356 + t234 * t354 + t263 * t350) * t239 + (t184 * t240 + t186 * t241 + t357 * t233 + t355 * t234 + t351 * t263) * t238) * t238 / 0.2e1 + ((t214 * t242 + t215 * t243 + t235 * t353 + t236 * t352 + t265 * t349) * t259 + (t185 * t242 + t187 * t243 + t356 * t235 + t354 * t236 + t350 * t265) * t239 + (t184 * t242 + t186 * t243 + t235 * t357 + t355 * t236 + t351 * t265) * t238) * t239 / 0.2e1 + ((t214 * t261 + t215 * t262 + t353 * t252 + t352 * t253 - t349 * t339) * t259 + (t185 * t261 + t187 * t262 + t252 * t356 + t253 * t354 - t339 * t350) * t239 + (t184 * t261 + t186 * t262 + t252 * t357 + t355 * t253 - t351 * t339) * t238) * t259 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t302 + Icges(2,6) * t306) * V_base(5) + (Icges(2,5) * t306 - Icges(2,6) * t302) * V_base(4) + Icges(2,3) * t293 / 0.2e1) * t293;
T  = t1;
