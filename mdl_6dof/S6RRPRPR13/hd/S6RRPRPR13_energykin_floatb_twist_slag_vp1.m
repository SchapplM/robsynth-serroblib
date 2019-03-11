% Calculate kinetic energy for
% S6RRPRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR13_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR13_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRPR13_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR13_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR13_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR13_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:24:03
% EndTime: 2019-03-09 11:24:07
% DurationCPUTime: 4.21s
% Computational Cost: add. (2473->391), mult. (5222->542), div. (0->0), fcn. (6142->12), ass. (0->172)
t364 = Icges(3,1) + Icges(4,2);
t363 = Icges(3,4) + Icges(4,6);
t362 = Icges(3,5) - Icges(4,4);
t361 = Icges(3,2) + Icges(4,3);
t360 = Icges(5,2) + Icges(6,3);
t359 = Icges(3,6) - Icges(4,5);
t358 = Icges(3,3) + Icges(4,1);
t300 = cos(pkin(6));
t304 = sin(qJ(1));
t305 = cos(qJ(2));
t327 = t304 * t305;
t303 = sin(qJ(2));
t306 = cos(qJ(1));
t328 = t303 * t306;
t263 = t300 * t327 + t328;
t302 = sin(qJ(4));
t298 = sin(pkin(6));
t339 = cos(qJ(4));
t319 = t298 * t339;
t233 = t263 * t302 + t304 * t319;
t326 = t305 * t306;
t329 = t303 * t304;
t264 = -t300 * t329 + t326;
t297 = sin(pkin(11));
t299 = cos(pkin(11));
t192 = -t233 * t297 + t264 * t299;
t334 = t264 * t297;
t193 = t233 * t299 + t334;
t332 = t298 * t304;
t232 = -t263 * t339 + t302 * t332;
t357 = -Icges(5,4) * t233 + Icges(6,5) * t193 - Icges(5,6) * t264 + Icges(6,6) * t192 + t360 * t232;
t261 = -t300 * t326 + t329;
t235 = t261 * t302 - t306 * t319;
t262 = t300 * t328 + t327;
t194 = -t235 * t297 + t262 * t299;
t335 = t262 * t297;
t195 = t235 * t299 + t335;
t330 = t298 * t306;
t234 = t261 * t339 + t302 * t330;
t356 = -Icges(5,4) * t235 + Icges(6,5) * t195 - Icges(5,6) * t262 + Icges(6,6) * t194 - t360 * t234;
t331 = t298 * t305;
t260 = t300 * t339 - t302 * t331;
t333 = t298 * t303;
t228 = -t260 * t297 + t299 * t333;
t320 = t297 * t333;
t229 = t260 * t299 + t320;
t259 = t300 * t302 + t305 * t319;
t355 = -Icges(5,4) * t260 + Icges(6,5) * t229 - Icges(5,6) * t333 + Icges(6,6) * t228 + t360 * t259;
t354 = t361 * t263 - t363 * t264 - t359 * t332;
t353 = t361 * t261 - t363 * t262 + t359 * t330;
t352 = -t363 * t263 + t364 * t264 + t362 * t332;
t351 = -t363 * t261 + t364 * t262 - t362 * t330;
t350 = -t359 * t263 + t362 * t264 + t358 * t332;
t349 = -t359 * t261 + t362 * t262 - t358 * t330;
t348 = t358 * t300 + (t362 * t303 + t359 * t305) * t298;
t347 = t359 * t300 + (t363 * t303 + t361 * t305) * t298;
t346 = t362 * t300 + (t364 * t303 + t363 * t305) * t298;
t338 = pkin(8) * t300;
t337 = pkin(5) * t299;
t336 = Icges(2,4) * t304;
t324 = qJD(2) * t298;
t323 = V_base(5) * pkin(7) + V_base(1);
t274 = t304 * t324 + V_base(4);
t293 = V_base(6) + qJD(1);
t231 = qJD(4) * t264 + t274;
t275 = qJD(2) * t300 + t293;
t257 = qJD(4) * t333 + t275;
t273 = -t306 * t324 + V_base(5);
t268 = t304 * pkin(1) - pkin(8) * t330;
t318 = -t268 * t293 + V_base(5) * t338 + t323;
t269 = pkin(1) * t306 + pkin(8) * t332;
t317 = V_base(4) * t268 - t269 * V_base(5) + V_base(3);
t230 = qJD(4) * t262 + t273;
t265 = (pkin(2) * t303 - qJ(3) * t305) * t298;
t316 = qJD(3) * t263 + t273 * t265 + t318;
t315 = t293 * t269 + V_base(2) + (-pkin(7) - t338) * V_base(4);
t223 = pkin(2) * t264 + qJ(3) * t263;
t314 = qJD(3) * t261 + t275 * t223 + t315;
t222 = pkin(2) * t262 + qJ(3) * t261;
t313 = -qJD(3) * t331 + t274 * t222 + t317;
t241 = -pkin(3) * t330 + t262 * pkin(9);
t267 = pkin(3) * t300 + pkin(9) * t333;
t312 = t273 * t267 + (-t222 - t241) * t275 + t316;
t219 = pkin(4) * t260 + qJ(5) * t259;
t311 = qJD(5) * t232 + t230 * t219 + t312;
t240 = pkin(3) * t332 + pkin(9) * t264;
t310 = t275 * t240 + (-t265 - t267) * t274 + t314;
t309 = t274 * t241 + (-t223 - t240) * t273 + t313;
t184 = pkin(4) * t233 + qJ(5) * t232;
t308 = -qJD(5) * t234 + t257 * t184 + t310;
t185 = pkin(4) * t235 - qJ(5) * t234;
t307 = qJD(5) * t259 + t231 * t185 + t309;
t296 = pkin(11) + qJ(6);
t294 = Icges(2,4) * t306;
t292 = cos(t296);
t291 = sin(t296);
t283 = rSges(2,1) * t306 - t304 * rSges(2,2);
t282 = t304 * rSges(2,1) + rSges(2,2) * t306;
t281 = Icges(2,1) * t306 - t336;
t280 = Icges(2,1) * t304 + t294;
t279 = -Icges(2,2) * t304 + t294;
t278 = Icges(2,2) * t306 + t336;
t272 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t271 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t270 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t250 = rSges(4,1) * t300 + (-rSges(4,2) * t303 - rSges(4,3) * t305) * t298;
t249 = rSges(3,3) * t300 + (rSges(3,1) * t303 + rSges(3,2) * t305) * t298;
t239 = V_base(5) * rSges(2,3) - t282 * t293 + t323;
t238 = t283 * t293 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t236 = t282 * V_base(4) - t283 * V_base(5) + V_base(3);
t225 = t260 * t292 + t291 * t333;
t224 = -t260 * t291 + t292 * t333;
t220 = qJD(6) * t259 + t257;
t217 = rSges(3,1) * t264 - rSges(3,2) * t263 + rSges(3,3) * t332;
t216 = t262 * rSges(3,1) - t261 * rSges(3,2) - rSges(3,3) * t330;
t215 = -rSges(4,1) * t330 - t262 * rSges(4,2) + t261 * rSges(4,3);
t214 = rSges(4,1) * t332 - rSges(4,2) * t264 + rSges(4,3) * t263;
t201 = rSges(5,1) * t260 - rSges(5,2) * t259 + rSges(5,3) * t333;
t200 = Icges(5,1) * t260 - Icges(5,4) * t259 + Icges(5,5) * t333;
t198 = Icges(5,5) * t260 - Icges(5,6) * t259 + Icges(5,3) * t333;
t191 = t235 * t292 + t262 * t291;
t190 = -t235 * t291 + t262 * t292;
t189 = t233 * t292 + t264 * t291;
t188 = -t233 * t291 + t264 * t292;
t187 = qJD(6) * t232 + t231;
t186 = -qJD(6) * t234 + t230;
t182 = rSges(5,1) * t235 + rSges(5,2) * t234 + rSges(5,3) * t262;
t181 = rSges(5,1) * t233 - rSges(5,2) * t232 + rSges(5,3) * t264;
t180 = Icges(5,1) * t235 + Icges(5,4) * t234 + Icges(5,5) * t262;
t179 = Icges(5,1) * t233 - Icges(5,4) * t232 + Icges(5,5) * t264;
t176 = Icges(5,5) * t235 + Icges(5,6) * t234 + Icges(5,3) * t262;
t175 = Icges(5,5) * t233 - Icges(5,6) * t232 + Icges(5,3) * t264;
t173 = rSges(6,1) * t229 + rSges(6,2) * t228 + rSges(6,3) * t259;
t172 = Icges(6,1) * t229 + Icges(6,4) * t228 + Icges(6,5) * t259;
t171 = Icges(6,4) * t229 + Icges(6,2) * t228 + Icges(6,6) * t259;
t169 = pkin(5) * t320 + pkin(10) * t259 + t260 * t337;
t168 = rSges(7,1) * t225 + rSges(7,2) * t224 + rSges(7,3) * t259;
t167 = Icges(7,1) * t225 + Icges(7,4) * t224 + Icges(7,5) * t259;
t166 = Icges(7,4) * t225 + Icges(7,2) * t224 + Icges(7,6) * t259;
t165 = Icges(7,5) * t225 + Icges(7,6) * t224 + Icges(7,3) * t259;
t163 = -t216 * t275 + t249 * t273 + t318;
t162 = t217 * t275 - t249 * t274 + t315;
t161 = rSges(6,1) * t195 + rSges(6,2) * t194 - rSges(6,3) * t234;
t160 = rSges(6,1) * t193 + rSges(6,2) * t192 + rSges(6,3) * t232;
t159 = Icges(6,1) * t195 + Icges(6,4) * t194 - Icges(6,5) * t234;
t158 = Icges(6,1) * t193 + Icges(6,4) * t192 + Icges(6,5) * t232;
t157 = Icges(6,4) * t195 + Icges(6,2) * t194 - Icges(6,6) * t234;
t156 = Icges(6,4) * t193 + Icges(6,2) * t192 + Icges(6,6) * t232;
t153 = t216 * t274 - t217 * t273 + t317;
t152 = rSges(7,1) * t191 + rSges(7,2) * t190 - rSges(7,3) * t234;
t151 = rSges(7,1) * t189 + rSges(7,2) * t188 + rSges(7,3) * t232;
t150 = Icges(7,1) * t191 + Icges(7,4) * t190 - Icges(7,5) * t234;
t149 = Icges(7,1) * t189 + Icges(7,4) * t188 + Icges(7,5) * t232;
t148 = Icges(7,4) * t191 + Icges(7,2) * t190 - Icges(7,6) * t234;
t147 = Icges(7,4) * t189 + Icges(7,2) * t188 + Icges(7,6) * t232;
t146 = Icges(7,5) * t191 + Icges(7,6) * t190 - Icges(7,3) * t234;
t145 = Icges(7,5) * t189 + Icges(7,6) * t188 + Icges(7,3) * t232;
t144 = pkin(5) * t335 - pkin(10) * t234 + t235 * t337;
t143 = pkin(5) * t334 + pkin(10) * t232 + t233 * t337;
t142 = t250 * t273 + (-t215 - t222) * t275 + t316;
t141 = t214 * t275 + (-t250 - t265) * t274 + t314;
t140 = t215 * t274 + (-t214 - t223) * t273 + t313;
t139 = -t182 * t257 + t201 * t230 + t312;
t138 = t181 * t257 - t201 * t231 + t310;
t137 = -t181 * t230 + t182 * t231 + t309;
t136 = t173 * t230 + (-t161 - t185) * t257 + t311;
t135 = t160 * t257 + (-t173 - t219) * t231 + t308;
t134 = t161 * t231 + (-t160 - t184) * t230 + t307;
t133 = t311 + (-t144 - t185) * t257 - t152 * t220 + t168 * t186 + t169 * t230;
t132 = t143 * t257 + t151 * t220 - t168 * t187 + (-t169 - t219) * t231 + t308;
t131 = t144 * t231 - t151 * t186 + t152 * t187 + (-t143 - t184) * t230 + t307;
t1 = m(3) * (t153 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(4) * (t140 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + m(6) * (t134 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + m(5) * (t137 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + m(7) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + m(1) * (t270 ^ 2 + t271 ^ 2 + t272 ^ 2) / 0.2e1 + t220 * ((t145 * t259 + t147 * t224 + t149 * t225) * t187 + (t146 * t259 + t148 * t224 + t150 * t225) * t186 + (t259 * t165 + t224 * t166 + t225 * t167) * t220) / 0.2e1 + t186 * ((-t145 * t234 + t147 * t190 + t149 * t191) * t187 + (-t234 * t146 + t190 * t148 + t191 * t150) * t186 + (-t165 * t234 + t166 * t190 + t167 * t191) * t220) / 0.2e1 + m(2) * (t236 ^ 2 + t238 ^ 2 + t239 ^ 2) / 0.2e1 + t187 * ((t232 * t145 + t188 * t147 + t189 * t149) * t187 + (t146 * t232 + t148 * t188 + t150 * t189) * t186 + (t165 * t232 + t166 * t188 + t167 * t189) * t220) / 0.2e1 + ((t171 * t194 + t172 * t195 + t198 * t262 + t200 * t235 - t355 * t234) * t257 + (t156 * t194 + t158 * t195 + t175 * t262 + t179 * t235 - t357 * t234) * t231 + (t157 * t194 + t159 * t195 + t176 * t262 + t180 * t235 - t356 * t234) * t230) * t230 / 0.2e1 + ((t171 * t192 + t172 * t193 + t198 * t264 + t200 * t233 + t232 * t355) * t257 + (t156 * t192 + t158 * t193 + t175 * t264 + t179 * t233 + t357 * t232) * t231 + (t157 * t192 + t159 * t193 + t176 * t264 + t180 * t233 + t232 * t356) * t230) * t231 / 0.2e1 + ((t171 * t228 + t172 * t229 + t198 * t333 + t200 * t260 + t355 * t259) * t257 + (t156 * t228 + t158 * t229 + t175 * t333 + t179 * t260 + t259 * t357) * t231 + (t157 * t228 + t159 * t229 + t176 * t333 + t180 * t260 + t259 * t356) * t230) * t257 / 0.2e1 + ((-t261 * t347 + t262 * t346 - t330 * t348) * t275 + (t261 * t354 + t262 * t352 - t330 * t350) * t274 + (t353 * t261 + t351 * t262 - t349 * t330) * t273) * t273 / 0.2e1 + ((-t263 * t347 + t264 * t346 + t332 * t348) * t275 + (t354 * t263 + t352 * t264 + t350 * t332) * t274 + (t263 * t353 + t264 * t351 + t332 * t349) * t273) * t274 / 0.2e1 + ((t273 * t349 + t274 * t350 + t348 * t275) * t300 + ((t303 * t346 + t305 * t347) * t275 + (t303 * t352 - t305 * t354) * t274 + (t303 * t351 - t305 * t353) * t273) * t298) * t275 / 0.2e1 + ((-t304 * t278 + t280 * t306 + Icges(1,4)) * V_base(5) + (-t304 * t279 + t281 * t306 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t278 * t306 + t304 * t280 + Icges(1,2)) * V_base(5) + (t279 * t306 + t304 * t281 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t304 + Icges(2,6) * t306) * V_base(5) + (Icges(2,5) * t306 - Icges(2,6) * t304) * V_base(4) + Icges(2,3) * t293 / 0.2e1) * t293;
T  = t1;
