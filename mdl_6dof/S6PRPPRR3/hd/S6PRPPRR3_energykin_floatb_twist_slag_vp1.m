% Calculate kinetic energy for
% S6PRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPPRR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPPRR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPPRR3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPPRR3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:21:21
% EndTime: 2019-03-08 19:21:25
% DurationCPUTime: 3.95s
% Computational Cost: add. (2676->386), mult. (6399->543), div. (0->0), fcn. (7807->12), ass. (0->167)
t362 = Icges(3,1) + Icges(4,1);
t361 = Icges(3,4) - Icges(4,5);
t360 = Icges(4,4) + Icges(3,5);
t359 = Icges(3,2) + Icges(4,3);
t358 = Icges(4,6) - Icges(3,6);
t357 = -Icges(5,3) - Icges(3,3) - Icges(4,2);
t302 = sin(pkin(10));
t304 = cos(pkin(10));
t308 = sin(qJ(2));
t305 = cos(pkin(6));
t310 = cos(qJ(2));
t334 = t305 * t310;
t269 = t302 * t308 - t304 * t334;
t335 = t305 * t308;
t270 = t302 * t310 + t304 * t335;
t303 = sin(pkin(6));
t337 = t303 * t304;
t354 = t269 * t359 - t270 * t361 - t337 * t358;
t271 = t302 * t334 + t304 * t308;
t272 = -t302 * t335 + t304 * t310;
t338 = t302 * t303;
t353 = t271 * t359 - t272 * t361 + t338 * t358;
t352 = -t361 * t269 + t270 * t362 - t360 * t337;
t351 = -t361 * t271 + t272 * t362 + t360 * t338;
t350 = t358 * t305 + (-t308 * t361 - t310 * t359) * t303;
t349 = t360 * t305 + (t308 * t362 + t361 * t310) * t303;
t301 = sin(pkin(11));
t340 = cos(pkin(11));
t236 = -t269 * t340 + t270 * t301;
t237 = t269 * t301 + t270 * t340;
t348 = Icges(5,5) * t237 - Icges(5,6) * t236 - t269 * t358 - t270 * t360 - t337 * t357;
t238 = -t271 * t340 + t272 * t301;
t239 = t271 * t301 + t272 * t340;
t347 = Icges(5,5) * t239 - Icges(5,6) * t238 - t271 * t358 - t272 * t360 + t338 * t357;
t265 = (t301 * t308 + t310 * t340) * t303;
t266 = (-t301 * t310 + t308 * t340) * t303;
t346 = Icges(5,5) * t266 - Icges(5,6) * t265 + (-t308 * t360 + t310 * t358) * t303 + t357 * t305;
t342 = cos(qJ(5));
t341 = pkin(7) * t305;
t339 = Icges(2,4) * t302;
t307 = sin(qJ(5));
t336 = t303 * t307;
t241 = pkin(2) * t270 + qJ(3) * t269;
t249 = pkin(3) * t270 + qJ(4) * t337;
t333 = -t241 - t249;
t242 = pkin(2) * t272 + qJ(3) * t271;
t250 = pkin(3) * t272 - qJ(4) * t338;
t332 = -t242 - t250;
t275 = (pkin(2) * t308 - qJ(3) * t310) * t303;
t276 = pkin(3) * t303 * t308 - qJ(4) * t305;
t331 = -t275 - t276;
t330 = qJD(2) * t303;
t329 = qJD(4) * t303;
t328 = V_base(5) * qJ(1) + V_base(1);
t324 = qJD(1) + V_base(3);
t323 = t303 * t342;
t283 = t302 * t330 + V_base(4);
t294 = qJD(2) * t305 + V_base(6);
t206 = qJD(5) * t238 + t283;
t245 = qJD(5) * t265 + t294;
t282 = -t304 * t330 + V_base(5);
t205 = qJD(5) * t236 + t282;
t277 = pkin(1) * t302 - pkin(7) * t337;
t322 = -t277 * V_base(6) + V_base(5) * t341 + t328;
t278 = pkin(1) * t304 + pkin(7) * t338;
t321 = V_base(4) * t277 - V_base(5) * t278 + t324;
t320 = V_base(6) * t278 + V_base(2) + (-qJ(1) - t341) * V_base(4);
t319 = qJD(3) * t271 + t282 * t275 + t322;
t318 = qJD(3) * t269 + t294 * t242 + t320;
t317 = -qJD(3) * t303 * t310 + t283 * t241 + t321;
t316 = t294 * t250 + t304 * t329 + t318;
t315 = t282 * t276 - t302 * t329 + t319;
t314 = -qJD(4) * t305 + t283 * t249 + t317;
t195 = pkin(4) * t239 + pkin(8) * t238;
t235 = pkin(4) * t266 + pkin(8) * t265;
t313 = t294 * t195 + (-t235 + t331) * t283 + t316;
t194 = pkin(4) * t237 + pkin(8) * t236;
t312 = t282 * t235 + (-t194 + t333) * t294 + t315;
t311 = t283 * t194 + (-t195 + t332) * t282 + t314;
t309 = cos(qJ(6));
t306 = sin(qJ(6));
t299 = Icges(2,4) * t304;
t291 = rSges(2,1) * t304 - rSges(2,2) * t302;
t290 = rSges(2,1) * t302 + rSges(2,2) * t304;
t289 = Icges(2,1) * t304 - t339;
t288 = Icges(2,1) * t302 + t299;
t287 = -Icges(2,2) * t302 + t299;
t286 = Icges(2,2) * t304 + t339;
t281 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t280 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t279 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t262 = t305 * rSges(3,3) + (rSges(3,1) * t308 + rSges(3,2) * t310) * t303;
t261 = t305 * rSges(4,2) + (rSges(4,1) * t308 - rSges(4,3) * t310) * t303;
t252 = V_base(5) * rSges(2,3) - t290 * V_base(6) + t328;
t251 = t291 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t247 = t266 * t342 - t305 * t307;
t246 = t266 * t307 + t305 * t342;
t244 = t290 * V_base(4) - t291 * V_base(5) + t324;
t232 = rSges(3,1) * t272 - rSges(3,2) * t271 + rSges(3,3) * t338;
t231 = rSges(4,1) * t272 + rSges(4,2) * t338 + rSges(4,3) * t271;
t230 = rSges(3,1) * t270 - rSges(3,2) * t269 - rSges(3,3) * t337;
t229 = rSges(4,1) * t270 - rSges(4,2) * t337 + rSges(4,3) * t269;
t215 = rSges(5,1) * t266 - rSges(5,2) * t265 - rSges(5,3) * t305;
t214 = Icges(5,1) * t266 - Icges(5,4) * t265 - Icges(5,5) * t305;
t213 = Icges(5,4) * t266 - Icges(5,2) * t265 - Icges(5,6) * t305;
t210 = t239 * t342 - t302 * t336;
t209 = t239 * t307 + t302 * t323;
t208 = t237 * t342 + t304 * t336;
t207 = t237 * t307 - t304 * t323;
t203 = t247 * t309 + t265 * t306;
t202 = -t247 * t306 + t265 * t309;
t201 = qJD(6) * t246 + t245;
t200 = pkin(5) * t247 + pkin(9) * t246;
t199 = rSges(6,1) * t247 - rSges(6,2) * t246 + rSges(6,3) * t265;
t198 = Icges(6,1) * t247 - Icges(6,4) * t246 + Icges(6,5) * t265;
t197 = Icges(6,4) * t247 - Icges(6,2) * t246 + Icges(6,6) * t265;
t196 = Icges(6,5) * t247 - Icges(6,6) * t246 + Icges(6,3) * t265;
t193 = rSges(5,1) * t239 - rSges(5,2) * t238 - rSges(5,3) * t338;
t192 = rSges(5,1) * t237 - rSges(5,2) * t236 + rSges(5,3) * t337;
t191 = Icges(5,1) * t239 - Icges(5,4) * t238 - Icges(5,5) * t338;
t190 = Icges(5,1) * t237 - Icges(5,4) * t236 + Icges(5,5) * t337;
t189 = Icges(5,4) * t239 - Icges(5,2) * t238 - Icges(5,6) * t338;
t188 = Icges(5,4) * t237 - Icges(5,2) * t236 + Icges(5,6) * t337;
t183 = t210 * t309 + t238 * t306;
t182 = -t210 * t306 + t238 * t309;
t181 = t208 * t309 + t236 * t306;
t180 = -t208 * t306 + t236 * t309;
t179 = qJD(6) * t209 + t206;
t178 = qJD(6) * t207 + t205;
t177 = pkin(5) * t210 + pkin(9) * t209;
t176 = pkin(5) * t208 + pkin(9) * t207;
t175 = -t230 * t294 + t262 * t282 + t322;
t174 = t232 * t294 - t262 * t283 + t320;
t173 = rSges(7,1) * t203 + rSges(7,2) * t202 + rSges(7,3) * t246;
t172 = Icges(7,1) * t203 + Icges(7,4) * t202 + Icges(7,5) * t246;
t171 = Icges(7,4) * t203 + Icges(7,2) * t202 + Icges(7,6) * t246;
t170 = Icges(7,5) * t203 + Icges(7,6) * t202 + Icges(7,3) * t246;
t169 = t230 * t283 - t232 * t282 + t321;
t168 = rSges(6,1) * t210 - rSges(6,2) * t209 + rSges(6,3) * t238;
t167 = rSges(6,1) * t208 - rSges(6,2) * t207 + rSges(6,3) * t236;
t166 = Icges(6,1) * t210 - Icges(6,4) * t209 + Icges(6,5) * t238;
t165 = Icges(6,1) * t208 - Icges(6,4) * t207 + Icges(6,5) * t236;
t164 = Icges(6,4) * t210 - Icges(6,2) * t209 + Icges(6,6) * t238;
t163 = Icges(6,4) * t208 - Icges(6,2) * t207 + Icges(6,6) * t236;
t162 = Icges(6,5) * t210 - Icges(6,6) * t209 + Icges(6,3) * t238;
t161 = Icges(6,5) * t208 - Icges(6,6) * t207 + Icges(6,3) * t236;
t160 = t261 * t282 + (-t229 - t241) * t294 + t319;
t159 = t231 * t294 + (-t261 - t275) * t283 + t318;
t158 = t283 * t229 + (-t231 - t242) * t282 + t317;
t157 = rSges(7,1) * t183 + rSges(7,2) * t182 + rSges(7,3) * t209;
t156 = rSges(7,1) * t181 + rSges(7,2) * t180 + rSges(7,3) * t207;
t155 = Icges(7,1) * t183 + Icges(7,4) * t182 + Icges(7,5) * t209;
t154 = Icges(7,1) * t181 + Icges(7,4) * t180 + Icges(7,5) * t207;
t153 = Icges(7,4) * t183 + Icges(7,2) * t182 + Icges(7,6) * t209;
t152 = Icges(7,4) * t181 + Icges(7,2) * t180 + Icges(7,6) * t207;
t151 = Icges(7,5) * t183 + Icges(7,6) * t182 + Icges(7,3) * t209;
t150 = Icges(7,5) * t181 + Icges(7,6) * t180 + Icges(7,3) * t207;
t149 = t215 * t282 + (-t192 + t333) * t294 + t315;
t148 = t193 * t294 + (-t215 + t331) * t283 + t316;
t147 = t283 * t192 + (-t193 + t332) * t282 + t314;
t146 = -t167 * t245 + t199 * t205 + t312;
t145 = t168 * t245 - t199 * t206 + t313;
t144 = t206 * t167 - t205 * t168 + t311;
t143 = -t156 * t201 + t173 * t178 - t176 * t245 + t200 * t205 + t312;
t142 = t157 * t201 - t173 * t179 + t177 * t245 - t200 * t206 + t313;
t141 = t179 * t156 - t178 * t157 + t206 * t176 - t205 * t177 + t311;
t1 = m(1) * (t279 ^ 2 + t280 ^ 2 + t281 ^ 2) / 0.2e1 + t245 * ((t162 * t265 - t164 * t246 + t166 * t247) * t206 + (t161 * t265 - t163 * t246 + t165 * t247) * t205 + (t265 * t196 - t246 * t197 + t247 * t198) * t245) / 0.2e1 + m(2) * (t244 ^ 2 + t251 ^ 2 + t252 ^ 2) / 0.2e1 + t206 * ((t238 * t162 - t209 * t164 + t210 * t166) * t206 + (t161 * t238 - t163 * t209 + t165 * t210) * t205 + (t196 * t238 - t197 * t209 + t198 * t210) * t245) / 0.2e1 + t205 * ((t162 * t236 - t164 * t207 + t166 * t208) * t206 + (t236 * t161 - t207 * t163 + t208 * t165) * t205 + (t196 * t236 - t197 * t207 + t198 * t208) * t245) / 0.2e1 + t201 * ((t151 * t246 + t153 * t202 + t155 * t203) * t179 + (t150 * t246 + t152 * t202 + t154 * t203) * t178 + (t246 * t170 + t202 * t171 + t203 * t172) * t201) / 0.2e1 + t178 * ((t151 * t207 + t153 * t180 + t155 * t181) * t179 + (t207 * t150 + t180 * t152 + t181 * t154) * t178 + (t170 * t207 + t171 * t180 + t172 * t181) * t201) / 0.2e1 + t179 * ((t209 * t151 + t182 * t153 + t183 * t155) * t179 + (t150 * t209 + t152 * t182 + t154 * t183) * t178 + (t170 * t209 + t171 * t182 + t172 * t183) * t201) / 0.2e1 + m(3) * (t169 ^ 2 + t174 ^ 2 + t175 ^ 2) / 0.2e1 + m(4) * (t158 ^ 2 + t159 ^ 2 + t160 ^ 2) / 0.2e1 + m(6) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(5) * (t147 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + m(7) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + ((-t286 * t302 + t288 * t304 + Icges(1,4)) * V_base(5) + (-t287 * t302 + t289 * t304 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t286 * t304 + t288 * t302 + Icges(1,2)) * V_base(5) + (t287 * t304 + t289 * t302 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t213 * t236 + t214 * t237 + t269 * t350 + t270 * t349 + t337 * t346) * t294 + (-t189 * t236 + t191 * t237 + t269 * t353 + t270 * t351 + t337 * t347) * t283 + (-t236 * t188 + t237 * t190 + t354 * t269 + t352 * t270 + t348 * t337) * t282) * t282 / 0.2e1 + ((-t213 * t238 + t214 * t239 + t271 * t350 + t272 * t349 - t338 * t346) * t294 + (-t238 * t189 + t239 * t191 + t353 * t271 + t351 * t272 - t347 * t338) * t283 + (-t188 * t238 + t190 * t239 + t271 * t354 + t272 * t352 - t338 * t348) * t282) * t283 / 0.2e1 + ((-t189 * t265 + t191 * t266) * t283 + (-t188 * t265 + t190 * t266) * t282 + (-t265 * t213 + t266 * t214) * t294 + (-t282 * t348 - t283 * t347 - t294 * t346) * t305 + ((t308 * t349 - t310 * t350) * t294 + (t308 * t351 - t310 * t353) * t283 + (t308 * t352 - t310 * t354) * t282) * t303) * t294 / 0.2e1 + ((Icges(2,5) * t302 + Icges(2,6) * t304 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t304 - Icges(2,6) * t302 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
