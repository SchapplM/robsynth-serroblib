% Calculate kinetic energy for
% S6RRRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 23:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR10_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR10_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPR10_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR10_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR10_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR10_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:00:51
% EndTime: 2019-03-09 23:00:56
% DurationCPUTime: 4.57s
% Computational Cost: add. (3157->397), mult. (5456->574), div. (0->0), fcn. (6376->12), ass. (0->178)
t363 = Icges(5,1) + Icges(6,2);
t362 = Icges(6,1) + Icges(5,3);
t361 = -Icges(5,4) - Icges(6,6);
t360 = Icges(6,4) - Icges(5,5);
t359 = Icges(6,5) - Icges(5,6);
t358 = Icges(5,2) + Icges(6,3);
t300 = cos(pkin(6));
t304 = sin(qJ(1));
t307 = cos(qJ(2));
t334 = t304 * t307;
t303 = sin(qJ(2));
t308 = cos(qJ(1));
t335 = t303 * t308;
t267 = t300 * t335 + t334;
t299 = sin(pkin(6));
t332 = qJ(3) + qJ(4);
t325 = cos(t332);
t323 = t299 * t325;
t324 = sin(t332);
t235 = t267 * t324 + t308 * t323;
t322 = t299 * t324;
t236 = t267 * t325 - t308 * t322;
t333 = t307 * t308;
t336 = t303 * t304;
t266 = -t300 * t333 + t336;
t357 = t358 * t235 + t361 * t236 + t359 * t266;
t269 = -t300 * t336 + t333;
t237 = t269 * t324 - t304 * t323;
t238 = t269 * t325 + t304 * t322;
t268 = t300 * t334 + t335;
t356 = t358 * t237 + t361 * t238 + t359 * t268;
t355 = t359 * t235 - t360 * t236 + t362 * t266;
t354 = t359 * t237 - t360 * t238 + t362 * t268;
t353 = t361 * t235 + t363 * t236 - t360 * t266;
t352 = t361 * t237 + t363 * t238 - t360 * t268;
t256 = -t300 * t325 + t303 * t322;
t257 = t300 * t324 + t303 * t323;
t339 = t299 * t307;
t351 = t358 * t256 + t361 * t257 - t359 * t339;
t350 = t361 * t256 + t363 * t257 + t360 * t339;
t349 = t359 * t256 - t360 * t257 - t362 * t339;
t345 = pkin(8) * t300;
t306 = cos(qJ(3));
t344 = pkin(3) * t306;
t342 = Icges(2,4) * t304;
t341 = t299 * t304;
t340 = t299 * t306;
t338 = t299 * t308;
t302 = sin(qJ(3));
t337 = t300 * t302;
t331 = qJD(2) * t299;
t330 = V_base(5) * pkin(7) + V_base(1);
t327 = t302 * t341;
t326 = t302 * t338;
t278 = t304 * t331 + V_base(4);
t296 = V_base(6) + qJD(1);
t241 = qJD(3) * t268 + t278;
t281 = qJD(2) * t300 + t296;
t213 = qJD(4) * t268 + t241;
t277 = -t308 * t331 + V_base(5);
t272 = pkin(1) * t304 - pkin(8) * t338;
t321 = -t272 * t296 + V_base(5) * t345 + t330;
t273 = pkin(1) * t308 + pkin(8) * t341;
t320 = V_base(4) * t272 - t273 * V_base(5) + V_base(3);
t240 = qJD(3) * t266 + t277;
t212 = qJD(4) * t266 + t240;
t319 = t296 * t273 + V_base(2) + (-pkin(7) - t345) * V_base(4);
t250 = (-qJD(3) - qJD(4)) * t339 + t281;
t229 = t267 * pkin(2) + t266 * pkin(9);
t271 = (pkin(2) * t303 - pkin(9) * t307) * t299;
t318 = -t229 * t281 + t277 * t271 + t321;
t230 = t269 * pkin(2) + t268 * pkin(9);
t317 = t278 * t229 - t230 * t277 + t320;
t316 = t281 * t230 - t271 * t278 + t319;
t181 = -pkin(3) * t326 + pkin(10) * t266 + t267 * t344;
t226 = pkin(3) * t337 + (-pkin(10) * t307 + t303 * t344) * t299;
t262 = -qJD(3) * t339 + t281;
t315 = -t181 * t262 + t240 * t226 + t318;
t182 = pkin(3) * t327 + pkin(10) * t268 + t269 * t344;
t314 = t241 * t181 - t182 * t240 + t317;
t218 = pkin(4) * t257 + qJ(5) * t256;
t313 = qJD(5) * t237 + t212 * t218 + t315;
t194 = pkin(4) * t236 + qJ(5) * t235;
t312 = qJD(5) * t256 + t213 * t194 + t314;
t311 = t262 * t182 - t226 * t241 + t316;
t195 = pkin(4) * t238 + qJ(5) * t237;
t310 = qJD(5) * t235 + t250 * t195 + t311;
t305 = cos(qJ(6));
t301 = sin(qJ(6));
t297 = Icges(2,4) * t308;
t289 = rSges(2,1) * t308 - rSges(2,2) * t304;
t288 = rSges(2,1) * t304 + rSges(2,2) * t308;
t287 = Icges(2,1) * t308 - t342;
t286 = Icges(2,1) * t304 + t297;
t285 = -Icges(2,2) * t304 + t297;
t284 = Icges(2,2) * t308 + t342;
t276 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t275 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t274 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t265 = t303 * t340 + t337;
t264 = -t299 * t302 * t303 + t300 * t306;
t255 = rSges(3,3) * t300 + (rSges(3,1) * t303 + rSges(3,2) * t307) * t299;
t254 = Icges(3,5) * t300 + (Icges(3,1) * t303 + Icges(3,4) * t307) * t299;
t253 = Icges(3,6) * t300 + (Icges(3,4) * t303 + Icges(3,2) * t307) * t299;
t252 = Icges(3,3) * t300 + (Icges(3,5) * t303 + Icges(3,6) * t307) * t299;
t249 = V_base(5) * rSges(2,3) - t288 * t296 + t330;
t248 = t289 * t296 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t246 = t288 * V_base(4) - t289 * V_base(5) + V_base(3);
t245 = t269 * t306 + t327;
t244 = -t269 * t302 + t304 * t340;
t243 = t267 * t306 - t326;
t242 = -t267 * t302 - t306 * t338;
t239 = -pkin(5) * t339 + pkin(11) * t257;
t234 = t256 * t301 - t305 * t339;
t233 = t256 * t305 + t301 * t339;
t228 = rSges(3,1) * t269 - rSges(3,2) * t268 + rSges(3,3) * t341;
t227 = rSges(3,1) * t267 - rSges(3,2) * t266 - rSges(3,3) * t338;
t225 = Icges(3,1) * t269 - Icges(3,4) * t268 + Icges(3,5) * t341;
t224 = Icges(3,1) * t267 - Icges(3,4) * t266 - Icges(3,5) * t338;
t223 = Icges(3,4) * t269 - Icges(3,2) * t268 + Icges(3,6) * t341;
t222 = Icges(3,4) * t267 - Icges(3,2) * t266 - Icges(3,6) * t338;
t221 = Icges(3,5) * t269 - Icges(3,6) * t268 + Icges(3,3) * t341;
t220 = Icges(3,5) * t267 - Icges(3,6) * t266 - Icges(3,3) * t338;
t219 = rSges(4,1) * t265 + rSges(4,2) * t264 - rSges(4,3) * t339;
t217 = Icges(4,1) * t265 + Icges(4,4) * t264 - Icges(4,5) * t339;
t216 = Icges(4,4) * t265 + Icges(4,2) * t264 - Icges(4,6) * t339;
t215 = Icges(4,5) * t265 + Icges(4,6) * t264 - Icges(4,3) * t339;
t210 = qJD(6) * t257 + t250;
t209 = rSges(5,1) * t257 - rSges(5,2) * t256 - rSges(5,3) * t339;
t208 = -rSges(6,1) * t339 - rSges(6,2) * t257 + rSges(6,3) * t256;
t201 = pkin(5) * t268 + pkin(11) * t238;
t200 = pkin(5) * t266 + pkin(11) * t236;
t199 = t237 * t301 + t268 * t305;
t198 = t237 * t305 - t268 * t301;
t197 = t235 * t301 + t266 * t305;
t196 = t235 * t305 - t266 * t301;
t192 = rSges(4,1) * t245 + rSges(4,2) * t244 + rSges(4,3) * t268;
t191 = rSges(4,1) * t243 + rSges(4,2) * t242 + rSges(4,3) * t266;
t190 = Icges(4,1) * t245 + Icges(4,4) * t244 + Icges(4,5) * t268;
t189 = Icges(4,1) * t243 + Icges(4,4) * t242 + Icges(4,5) * t266;
t188 = Icges(4,4) * t245 + Icges(4,2) * t244 + Icges(4,6) * t268;
t187 = Icges(4,4) * t243 + Icges(4,2) * t242 + Icges(4,6) * t266;
t186 = Icges(4,5) * t245 + Icges(4,6) * t244 + Icges(4,3) * t268;
t185 = Icges(4,5) * t243 + Icges(4,6) * t242 + Icges(4,3) * t266;
t184 = qJD(6) * t238 + t213;
t183 = qJD(6) * t236 + t212;
t180 = rSges(5,1) * t238 - rSges(5,2) * t237 + rSges(5,3) * t268;
t179 = rSges(5,1) * t236 - rSges(5,2) * t235 + rSges(5,3) * t266;
t178 = rSges(6,1) * t268 - rSges(6,2) * t238 + rSges(6,3) * t237;
t177 = rSges(6,1) * t266 - rSges(6,2) * t236 + rSges(6,3) * t235;
t163 = rSges(7,1) * t234 + rSges(7,2) * t233 + rSges(7,3) * t257;
t162 = Icges(7,1) * t234 + Icges(7,4) * t233 + Icges(7,5) * t257;
t161 = Icges(7,4) * t234 + Icges(7,2) * t233 + Icges(7,6) * t257;
t160 = Icges(7,5) * t234 + Icges(7,6) * t233 + Icges(7,3) * t257;
t155 = -t227 * t281 + t255 * t277 + t321;
t154 = t228 * t281 - t255 * t278 + t319;
t153 = t227 * t278 - t228 * t277 + t320;
t152 = rSges(7,1) * t199 + rSges(7,2) * t198 + rSges(7,3) * t238;
t151 = rSges(7,1) * t197 + rSges(7,2) * t196 + rSges(7,3) * t236;
t150 = Icges(7,1) * t199 + Icges(7,4) * t198 + Icges(7,5) * t238;
t149 = Icges(7,1) * t197 + Icges(7,4) * t196 + Icges(7,5) * t236;
t148 = Icges(7,4) * t199 + Icges(7,2) * t198 + Icges(7,6) * t238;
t147 = Icges(7,4) * t197 + Icges(7,2) * t196 + Icges(7,6) * t236;
t146 = Icges(7,5) * t199 + Icges(7,6) * t198 + Icges(7,3) * t238;
t145 = Icges(7,5) * t197 + Icges(7,6) * t196 + Icges(7,3) * t236;
t144 = -t191 * t262 + t219 * t240 + t318;
t143 = t192 * t262 - t219 * t241 + t316;
t142 = t191 * t241 - t192 * t240 + t317;
t141 = -t179 * t250 + t209 * t212 + t315;
t140 = t180 * t250 - t209 * t213 + t311;
t139 = t179 * t213 - t180 * t212 + t314;
t138 = t208 * t212 + (-t177 - t194) * t250 + t313;
t137 = t178 * t250 + (-t208 - t218) * t213 + t310;
t136 = t177 * t213 + (-t178 - t195) * t212 + t312;
t135 = -t151 * t210 + t163 * t183 + t212 * t239 + (-t194 - t200) * t250 + t313;
t134 = t152 * t210 - t163 * t184 + t201 * t250 + (-t218 - t239) * t213 + t310;
t133 = t151 * t184 - t152 * t183 + t200 * t213 + (-t195 - t201) * t212 + t312;
t1 = m(3) * (t153 ^ 2 + t154 ^ 2 + t155 ^ 2) / 0.2e1 + m(4) * (t142 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(6) * (t136 ^ 2 + t137 ^ 2 + t138 ^ 2) / 0.2e1 + m(5) * (t139 ^ 2 + t140 ^ 2 + t141 ^ 2) / 0.2e1 + m(7) * (t133 ^ 2 + t134 ^ 2 + t135 ^ 2) / 0.2e1 + t278 * ((t221 * t341 - t223 * t268 + t225 * t269) * t278 + (t220 * t341 - t222 * t268 + t224 * t269) * t277 + (t252 * t341 - t253 * t268 + t254 * t269) * t281) / 0.2e1 + t262 * ((-t186 * t339 + t188 * t264 + t190 * t265) * t241 + (-t185 * t339 + t187 * t264 + t189 * t265) * t240 + (-t215 * t339 + t216 * t264 + t217 * t265) * t262) / 0.2e1 + t277 * ((-t221 * t338 - t223 * t266 + t225 * t267) * t278 + (-t220 * t338 - t222 * t266 + t224 * t267) * t277 + (-t252 * t338 - t253 * t266 + t254 * t267) * t281) / 0.2e1 + t281 * ((t220 * t277 + t221 * t278 + t252 * t281) * t300 + ((t223 * t307 + t225 * t303) * t278 + (t222 * t307 + t224 * t303) * t277 + (t253 * t307 + t254 * t303) * t281) * t299) / 0.2e1 + t183 * ((t146 * t236 + t148 * t196 + t150 * t197) * t184 + (t236 * t145 + t196 * t147 + t197 * t149) * t183 + (t160 * t236 + t161 * t196 + t162 * t197) * t210) / 0.2e1 + t184 * ((t238 * t146 + t198 * t148 + t199 * t150) * t184 + (t145 * t238 + t147 * t198 + t149 * t199) * t183 + (t160 * t238 + t161 * t198 + t162 * t199) * t210) / 0.2e1 + m(2) * (t246 ^ 2 + t248 ^ 2 + t249 ^ 2) / 0.2e1 + t210 * ((t146 * t257 + t148 * t233 + t150 * t234) * t184 + (t145 * t257 + t147 * t233 + t149 * t234) * t183 + (t257 * t160 + t233 * t161 + t234 * t162) * t210) / 0.2e1 + t240 * ((t186 * t266 + t188 * t242 + t190 * t243) * t241 + (t185 * t266 + t187 * t242 + t189 * t243) * t240 + (t215 * t266 + t216 * t242 + t217 * t243) * t262) / 0.2e1 + t241 * ((t186 * t268 + t188 * t244 + t190 * t245) * t241 + (t185 * t268 + t187 * t244 + t189 * t245) * t240 + (t215 * t268 + t216 * t244 + t217 * t245) * t262) / 0.2e1 + m(1) * (t274 ^ 2 + t275 ^ 2 + t276 ^ 2) / 0.2e1 + ((t235 * t351 + t236 * t350 + t266 * t349) * t250 + (t235 * t356 + t236 * t352 + t266 * t354) * t213 + (t235 * t357 + t353 * t236 + t355 * t266) * t212) * t212 / 0.2e1 + ((t237 * t351 + t238 * t350 + t268 * t349) * t250 + (t237 * t356 + t238 * t352 + t268 * t354) * t213 + (t237 * t357 + t353 * t238 + t355 * t268) * t212) * t213 / 0.2e1 + ((t256 * t351 + t257 * t350 - t339 * t349) * t250 + (t256 * t356 + t257 * t352 - t339 * t354) * t213 + (t256 * t357 + t353 * t257 - t355 * t339) * t212) * t250 / 0.2e1 + ((-t284 * t304 + t286 * t308 + Icges(1,4)) * V_base(5) + (-t285 * t304 + t287 * t308 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t284 * t308 + t286 * t304 + Icges(1,2)) * V_base(5) + (t285 * t308 + t287 * t304 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t304 + Icges(2,6) * t308) * V_base(5) + (Icges(2,5) * t308 - Icges(2,6) * t304) * V_base(4) + Icges(2,3) * t296 / 0.2e1) * t296;
T  = t1;
