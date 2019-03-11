% Calculate kinetic energy for
% S6RRPRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 11:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR10_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRPR10_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR10_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR10_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR10_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:04:18
% EndTime: 2019-03-09 11:04:22
% DurationCPUTime: 4.44s
% Computational Cost: add. (2997->396), mult. (5136->552), div. (0->0), fcn. (5992->12), ass. (0->178)
t366 = Icges(5,1) + Icges(6,2);
t365 = Icges(6,1) + Icges(5,3);
t364 = -Icges(5,4) - Icges(6,6);
t363 = Icges(6,4) - Icges(5,5);
t362 = Icges(6,5) - Icges(5,6);
t361 = Icges(5,2) + Icges(6,3);
t299 = cos(pkin(6));
t303 = sin(qJ(1));
t305 = cos(qJ(2));
t332 = t303 * t305;
t302 = sin(qJ(2));
t306 = cos(qJ(1));
t333 = t302 * t306;
t264 = t299 * t333 + t332;
t297 = sin(pkin(6));
t328 = pkin(11) + qJ(4);
t322 = cos(t328);
t320 = t297 * t322;
t321 = sin(t328);
t233 = t264 * t321 + t306 * t320;
t319 = t297 * t321;
t234 = t264 * t322 - t306 * t319;
t331 = t305 * t306;
t334 = t302 * t303;
t263 = -t299 * t331 + t334;
t360 = t361 * t233 + t364 * t234 + t362 * t263;
t266 = -t299 * t334 + t331;
t235 = t266 * t321 - t303 * t320;
t236 = t266 * t322 + t303 * t319;
t265 = t299 * t332 + t333;
t359 = t361 * t235 + t364 * t236 + t362 * t265;
t358 = t362 * t233 - t363 * t234 + t365 * t263;
t357 = t362 * t235 - t363 * t236 + t365 * t265;
t356 = t364 * t233 + t366 * t234 - t363 * t263;
t355 = t364 * t235 + t366 * t236 - t363 * t265;
t296 = sin(pkin(11));
t298 = cos(pkin(11));
t335 = t297 * t306;
t238 = -t264 * t296 - t298 * t335;
t323 = t296 * t335;
t239 = t264 * t298 - t323;
t182 = Icges(4,5) * t239 + Icges(4,6) * t238 + Icges(4,3) * t263;
t220 = Icges(3,4) * t264 - Icges(3,2) * t263 - Icges(3,6) * t335;
t354 = t182 - t220;
t337 = t297 * t303;
t240 = -t266 * t296 + t298 * t337;
t324 = t296 * t337;
t241 = t266 * t298 + t324;
t183 = Icges(4,5) * t241 + Icges(4,6) * t240 + Icges(4,3) * t265;
t221 = Icges(3,4) * t266 - Icges(3,2) * t265 + Icges(3,6) * t337;
t353 = t183 - t221;
t252 = -t299 * t322 + t302 * t319;
t253 = t299 * t321 + t302 * t320;
t336 = t297 * t305;
t352 = t361 * t252 + t364 * t253 - t362 * t336;
t351 = t364 * t252 + t366 * t253 + t363 * t336;
t350 = t362 * t252 - t363 * t253 - t365 * t336;
t338 = t297 * t302;
t261 = -t296 * t338 + t298 * t299;
t339 = t296 * t299;
t262 = t298 * t338 + t339;
t212 = Icges(4,5) * t262 + Icges(4,6) * t261 - Icges(4,3) * t336;
t250 = Icges(3,6) * t299 + (Icges(3,4) * t302 + Icges(3,2) * t305) * t297;
t349 = t212 - t250;
t342 = pkin(8) * t299;
t341 = pkin(3) * t298;
t340 = Icges(2,4) * t303;
t329 = qJD(2) * t297;
t327 = V_base(5) * pkin(7) + V_base(1);
t277 = t303 * t329 + V_base(4);
t293 = V_base(6) + qJD(1);
t243 = qJD(4) * t265 + t277;
t278 = qJD(2) * t299 + t293;
t276 = -t306 * t329 + V_base(5);
t269 = t303 * pkin(1) - pkin(8) * t335;
t318 = -t269 * t293 + V_base(5) * t342 + t327;
t270 = pkin(1) * t306 + pkin(8) * t337;
t317 = V_base(4) * t269 - t270 * V_base(5) + V_base(3);
t242 = qJD(4) * t263 + t276;
t259 = -qJD(4) * t336 + t278;
t267 = (pkin(2) * t302 - qJ(3) * t305) * t297;
t316 = qJD(3) * t265 + t276 * t267 + t318;
t315 = t293 * t270 + V_base(2) + (-pkin(7) - t342) * V_base(4);
t228 = pkin(2) * t266 + qJ(3) * t265;
t314 = qJD(3) * t263 + t278 * t228 + t315;
t227 = t264 * pkin(2) + t263 * qJ(3);
t313 = -qJD(3) * t336 + t277 * t227 + t317;
t180 = -pkin(3) * t323 + pkin(9) * t263 + t264 * t341;
t217 = pkin(3) * t339 + (-pkin(9) * t305 + t302 * t341) * t297;
t312 = t276 * t217 + (-t180 - t227) * t278 + t316;
t215 = pkin(4) * t253 + qJ(5) * t252;
t311 = qJD(5) * t235 + t242 * t215 + t312;
t181 = pkin(3) * t324 + pkin(9) * t265 + t266 * t341;
t310 = t278 * t181 + (-t217 - t267) * t277 + t314;
t309 = t277 * t180 + (-t181 - t228) * t276 + t313;
t192 = pkin(4) * t236 + qJ(5) * t235;
t308 = qJD(5) * t233 + t259 * t192 + t310;
t191 = pkin(4) * t234 + qJ(5) * t233;
t307 = qJD(5) * t252 + t243 * t191 + t309;
t304 = cos(qJ(6));
t301 = sin(qJ(6));
t294 = Icges(2,4) * t306;
t286 = rSges(2,1) * t306 - t303 * rSges(2,2);
t285 = t303 * rSges(2,1) + rSges(2,2) * t306;
t284 = Icges(2,1) * t306 - t340;
t283 = Icges(2,1) * t303 + t294;
t282 = -Icges(2,2) * t303 + t294;
t281 = Icges(2,2) * t306 + t340;
t273 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t272 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t271 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t254 = rSges(3,3) * t299 + (rSges(3,1) * t302 + rSges(3,2) * t305) * t297;
t251 = Icges(3,5) * t299 + (Icges(3,1) * t302 + Icges(3,4) * t305) * t297;
t249 = Icges(3,3) * t299 + (Icges(3,5) * t302 + Icges(3,6) * t305) * t297;
t247 = V_base(5) * rSges(2,3) - t285 * t293 + t327;
t246 = t286 * t293 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t244 = t285 * V_base(4) - t286 * V_base(5) + V_base(3);
t237 = -pkin(5) * t336 + pkin(10) * t253;
t232 = t252 * t301 - t304 * t336;
t231 = t252 * t304 + t301 * t336;
t226 = qJD(6) * t253 + t259;
t225 = rSges(3,1) * t266 - rSges(3,2) * t265 + rSges(3,3) * t337;
t224 = t264 * rSges(3,1) - t263 * rSges(3,2) - rSges(3,3) * t335;
t223 = Icges(3,1) * t266 - Icges(3,4) * t265 + Icges(3,5) * t337;
t222 = Icges(3,1) * t264 - Icges(3,4) * t263 - Icges(3,5) * t335;
t219 = Icges(3,5) * t266 - Icges(3,6) * t265 + Icges(3,3) * t337;
t218 = Icges(3,5) * t264 - Icges(3,6) * t263 - Icges(3,3) * t335;
t216 = rSges(4,1) * t262 + rSges(4,2) * t261 - rSges(4,3) * t336;
t214 = Icges(4,1) * t262 + Icges(4,4) * t261 - Icges(4,5) * t336;
t213 = Icges(4,4) * t262 + Icges(4,2) * t261 - Icges(4,6) * t336;
t209 = pkin(5) * t265 + pkin(10) * t236;
t208 = pkin(5) * t263 + pkin(10) * t234;
t207 = rSges(5,1) * t253 - rSges(5,2) * t252 - rSges(5,3) * t336;
t206 = -rSges(6,1) * t336 - rSges(6,2) * t253 + rSges(6,3) * t252;
t199 = t235 * t301 + t265 * t304;
t198 = t235 * t304 - t265 * t301;
t197 = t233 * t301 + t263 * t304;
t196 = t233 * t304 - t263 * t301;
t194 = qJD(6) * t236 + t243;
t193 = qJD(6) * t234 + t242;
t189 = rSges(4,1) * t241 + rSges(4,2) * t240 + rSges(4,3) * t265;
t188 = rSges(4,1) * t239 + rSges(4,2) * t238 + rSges(4,3) * t263;
t187 = Icges(4,1) * t241 + Icges(4,4) * t240 + Icges(4,5) * t265;
t186 = Icges(4,1) * t239 + Icges(4,4) * t238 + Icges(4,5) * t263;
t185 = Icges(4,4) * t241 + Icges(4,2) * t240 + Icges(4,6) * t265;
t184 = Icges(4,4) * t239 + Icges(4,2) * t238 + Icges(4,6) * t263;
t179 = rSges(5,1) * t236 - rSges(5,2) * t235 + rSges(5,3) * t265;
t178 = rSges(5,1) * t234 - rSges(5,2) * t233 + rSges(5,3) * t263;
t177 = rSges(6,1) * t265 - rSges(6,2) * t236 + rSges(6,3) * t235;
t176 = rSges(6,1) * t263 - rSges(6,2) * t234 + rSges(6,3) * t233;
t162 = rSges(7,1) * t232 + rSges(7,2) * t231 + rSges(7,3) * t253;
t161 = Icges(7,1) * t232 + Icges(7,4) * t231 + Icges(7,5) * t253;
t160 = Icges(7,4) * t232 + Icges(7,2) * t231 + Icges(7,6) * t253;
t159 = Icges(7,5) * t232 + Icges(7,6) * t231 + Icges(7,3) * t253;
t155 = -t224 * t278 + t254 * t276 + t318;
t154 = t225 * t278 - t254 * t277 + t315;
t153 = t224 * t277 - t225 * t276 + t317;
t152 = rSges(7,1) * t199 + rSges(7,2) * t198 + rSges(7,3) * t236;
t151 = rSges(7,1) * t197 + rSges(7,2) * t196 + rSges(7,3) * t234;
t150 = Icges(7,1) * t199 + Icges(7,4) * t198 + Icges(7,5) * t236;
t149 = Icges(7,1) * t197 + Icges(7,4) * t196 + Icges(7,5) * t234;
t148 = Icges(7,4) * t199 + Icges(7,2) * t198 + Icges(7,6) * t236;
t147 = Icges(7,4) * t197 + Icges(7,2) * t196 + Icges(7,6) * t234;
t146 = Icges(7,5) * t199 + Icges(7,6) * t198 + Icges(7,3) * t236;
t145 = Icges(7,5) * t197 + Icges(7,6) * t196 + Icges(7,3) * t234;
t144 = t216 * t276 + (-t188 - t227) * t278 + t316;
t143 = t189 * t278 + (-t216 - t267) * t277 + t314;
t142 = t188 * t277 + (-t189 - t228) * t276 + t313;
t141 = -t178 * t259 + t207 * t242 + t312;
t140 = t179 * t259 - t207 * t243 + t310;
t139 = t178 * t243 - t179 * t242 + t309;
t138 = t206 * t242 + (-t176 - t191) * t259 + t311;
t137 = t177 * t259 + (-t206 - t215) * t243 + t308;
t136 = t176 * t243 + (-t177 - t192) * t242 + t307;
t135 = t311 + (-t191 - t208) * t259 - t151 * t226 + t162 * t193 + t237 * t242;
t134 = t152 * t226 - t162 * t194 + t209 * t259 + (-t215 - t237) * t243 + t308;
t133 = t151 * t194 - t152 * t193 + t208 * t243 + (-t192 - t209) * t242 + t307;
t1 = m(7) * (t133 ^ 2 + t134 ^ 2 + t135 ^ 2) / 0.2e1 + m(3) * (t153 ^ 2 + t154 ^ 2 + t155 ^ 2) / 0.2e1 + m(6) * (t136 ^ 2 + t137 ^ 2 + t138 ^ 2) / 0.2e1 + m(5) * (t139 ^ 2 + t140 ^ 2 + t141 ^ 2) / 0.2e1 + m(4) * (t142 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + t193 * ((t146 * t234 + t148 * t196 + t150 * t197) * t194 + (t234 * t145 + t196 * t147 + t197 * t149) * t193 + (t159 * t234 + t160 * t196 + t161 * t197) * t226) / 0.2e1 + m(2) * (t244 ^ 2 + t246 ^ 2 + t247 ^ 2) / 0.2e1 + t226 * ((t146 * t253 + t148 * t231 + t150 * t232) * t194 + (t145 * t253 + t147 * t231 + t149 * t232) * t193 + (t159 * t253 + t160 * t231 + t161 * t232) * t226) / 0.2e1 + m(1) * (t271 ^ 2 + t272 ^ 2 + t273 ^ 2) / 0.2e1 + t194 * ((t236 * t146 + t198 * t148 + t199 * t150) * t194 + (t145 * t236 + t147 * t198 + t149 * t199) * t193 + (t159 * t236 + t160 * t198 + t161 * t199) * t226) / 0.2e1 + ((t233 * t352 + t234 * t351 + t263 * t350) * t259 + (t233 * t359 + t234 * t355 + t263 * t357) * t243 + (t360 * t233 + t356 * t234 + t358 * t263) * t242) * t242 / 0.2e1 + ((t235 * t352 + t236 * t351 + t265 * t350) * t259 + (t359 * t235 + t355 * t236 + t357 * t265) * t243 + (t235 * t360 + t356 * t236 + t358 * t265) * t242) * t243 / 0.2e1 + ((t352 * t252 + t351 * t253 - t350 * t336) * t259 + (t252 * t359 + t253 * t355 - t336 * t357) * t243 + (t252 * t360 + t356 * t253 - t358 * t336) * t242) * t259 / 0.2e1 + ((t213 * t238 + t214 * t239 - t249 * t335 + t264 * t251 + t263 * t349) * t278 + (t185 * t238 + t187 * t239 - t219 * t335 + t264 * t223 + t263 * t353) * t277 + (t184 * t238 + t186 * t239 - t218 * t335 + t264 * t222 + t354 * t263) * t276) * t276 / 0.2e1 + ((t213 * t240 + t214 * t241 + t249 * t337 + t251 * t266 + t265 * t349) * t278 + (t185 * t240 + t187 * t241 + t219 * t337 + t223 * t266 + t353 * t265) * t277 + (t184 * t240 + t186 * t241 + t218 * t337 + t222 * t266 + t265 * t354) * t276) * t277 / 0.2e1 + ((t218 * t276 + t219 * t277 + t249 * t278) * t299 + ((t221 * t305 + t223 * t302) * t277 + (t220 * t305 + t222 * t302) * t276 + (t250 * t305 + t251 * t302) * t278) * t297 + (-t183 * t336 + t185 * t261 + t187 * t262) * t277 + (-t182 * t336 + t184 * t261 + t186 * t262) * t276 + (-t212 * t336 + t213 * t261 + t214 * t262) * t278) * t278 / 0.2e1 + ((-t303 * t281 + t283 * t306 + Icges(1,4)) * V_base(5) + (-t303 * t282 + t284 * t306 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t281 * t306 + t303 * t283 + Icges(1,2)) * V_base(5) + (t282 * t306 + t303 * t284 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t303 + Icges(2,6) * t306) * V_base(5) + (Icges(2,5) * t306 - Icges(2,6) * t303) * V_base(4) + Icges(2,3) * t293 / 0.2e1) * t293;
T  = t1;
