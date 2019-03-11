% Calculate kinetic energy for
% S6PRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRPR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRPR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRPR6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:46:34
% EndTime: 2019-03-08 19:46:38
% DurationCPUTime: 4.27s
% Computational Cost: add. (2413->391), mult. (5222->541), div. (0->0), fcn. (6142->12), ass. (0->171)
t366 = Icges(3,1) + Icges(4,2);
t365 = Icges(3,4) + Icges(4,6);
t364 = Icges(3,5) - Icges(4,4);
t363 = Icges(3,2) + Icges(4,3);
t362 = Icges(5,2) + Icges(6,3);
t361 = Icges(3,6) - Icges(4,5);
t360 = Icges(3,3) + Icges(4,1);
t297 = sin(pkin(10));
t300 = cos(pkin(10));
t304 = sin(qJ(2));
t301 = cos(pkin(6));
t305 = cos(qJ(2));
t327 = t301 * t305;
t259 = t297 * t327 + t300 * t304;
t303 = sin(qJ(4));
t298 = sin(pkin(6));
t339 = cos(qJ(4));
t318 = t298 * t339;
t232 = t259 * t303 + t297 * t318;
t328 = t301 * t304;
t260 = -t297 * t328 + t300 * t305;
t296 = sin(pkin(11));
t299 = cos(pkin(11));
t192 = -t232 * t296 + t260 * t299;
t334 = t260 * t296;
t193 = t232 * t299 + t334;
t331 = t298 * t303;
t231 = -t259 * t339 + t297 * t331;
t357 = -Icges(5,4) * t232 + Icges(6,5) * t193 - Icges(5,6) * t260 + Icges(6,6) * t192 + t362 * t231;
t257 = t297 * t304 - t300 * t327;
t234 = t257 * t303 - t300 * t318;
t258 = t297 * t305 + t300 * t328;
t194 = -t234 * t296 + t258 * t299;
t335 = t258 * t296;
t195 = t234 * t299 + t335;
t233 = t257 * t339 + t300 * t331;
t356 = -Icges(5,4) * t234 + Icges(6,5) * t195 - Icges(5,6) * t258 + Icges(6,6) * t194 - t362 * t233;
t329 = t298 * t305;
t265 = t301 * t339 - t303 * t329;
t330 = t298 * t304;
t235 = -t265 * t296 + t299 * t330;
t319 = t296 * t330;
t236 = t265 * t299 + t319;
t264 = t301 * t303 + t305 * t318;
t355 = -Icges(5,4) * t265 + Icges(6,5) * t236 - Icges(5,6) * t330 + Icges(6,6) * t235 + t362 * t264;
t333 = t297 * t298;
t354 = t363 * t259 - t365 * t260 - t361 * t333;
t332 = t298 * t300;
t353 = t363 * t257 - t365 * t258 + t361 * t332;
t352 = -t365 * t259 + t366 * t260 + t364 * t333;
t351 = -t365 * t257 + t366 * t258 - t364 * t332;
t350 = -t361 * t259 + t364 * t260 + t360 * t333;
t349 = -t361 * t257 + t364 * t258 - t360 * t332;
t348 = t360 * t301 + (t364 * t304 + t361 * t305) * t298;
t347 = t361 * t301 + (t365 * t304 + t363 * t305) * t298;
t346 = t364 * t301 + (t366 * t304 + t365 * t305) * t298;
t338 = pkin(7) * t301;
t337 = pkin(5) * t299;
t336 = Icges(2,4) * t297;
t325 = qJD(2) * t298;
t324 = V_base(5) * qJ(1) + V_base(1);
t320 = qJD(1) + V_base(3);
t274 = t297 * t325 + V_base(4);
t285 = qJD(2) * t301 + V_base(6);
t230 = qJD(4) * t260 + t274;
t261 = qJD(4) * t330 + t285;
t273 = -t300 * t325 + V_base(5);
t229 = qJD(4) * t258 + t273;
t267 = pkin(1) * t297 - pkin(7) * t332;
t317 = -t267 * V_base(6) + V_base(5) * t338 + t324;
t268 = pkin(1) * t300 + pkin(7) * t333;
t316 = V_base(4) * t267 - V_base(5) * t268 + t320;
t315 = V_base(6) * t268 + V_base(2) + (-qJ(1) - t338) * V_base(4);
t266 = (pkin(2) * t304 - qJ(3) * t305) * t298;
t314 = qJD(3) * t259 + t273 * t266 + t317;
t220 = pkin(2) * t260 + qJ(3) * t259;
t313 = qJD(3) * t257 + t285 * t220 + t315;
t219 = pkin(2) * t258 + qJ(3) * t257;
t312 = -qJD(3) * t329 + t274 * t219 + t316;
t239 = -pkin(3) * t332 + pkin(8) * t258;
t269 = pkin(3) * t301 + pkin(8) * t330;
t311 = t273 * t269 + (-t219 - t239) * t285 + t314;
t238 = pkin(3) * t333 + pkin(8) * t260;
t310 = t285 * t238 + (-t266 - t269) * t274 + t313;
t221 = pkin(4) * t265 + qJ(5) * t264;
t309 = qJD(5) * t231 + t229 * t221 + t311;
t308 = t274 * t239 + (-t220 - t238) * t273 + t312;
t184 = pkin(4) * t232 + qJ(5) * t231;
t307 = -qJD(5) * t233 + t261 * t184 + t310;
t185 = pkin(4) * t234 - qJ(5) * t233;
t306 = qJD(5) * t264 + t230 * t185 + t308;
t295 = pkin(11) + qJ(6);
t293 = Icges(2,4) * t300;
t292 = cos(t295);
t291 = sin(t295);
t282 = rSges(2,1) * t300 - rSges(2,2) * t297;
t281 = rSges(2,1) * t297 + rSges(2,2) * t300;
t280 = Icges(2,1) * t300 - t336;
t279 = Icges(2,1) * t297 + t293;
t278 = -Icges(2,2) * t297 + t293;
t277 = Icges(2,2) * t300 + t336;
t272 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t271 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t270 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t250 = t301 * rSges(4,1) + (-rSges(4,2) * t304 - rSges(4,3) * t305) * t298;
t249 = t301 * rSges(3,3) + (rSges(3,1) * t304 + rSges(3,2) * t305) * t298;
t241 = V_base(5) * rSges(2,3) - t281 * V_base(6) + t324;
t240 = t282 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t228 = t281 * V_base(4) - t282 * V_base(5) + t320;
t225 = t265 * t292 + t291 * t330;
t224 = -t265 * t291 + t292 * t330;
t223 = qJD(6) * t264 + t261;
t217 = rSges(5,1) * t265 - rSges(5,2) * t264 + rSges(5,3) * t330;
t216 = Icges(5,1) * t265 - Icges(5,4) * t264 + Icges(5,5) * t330;
t214 = Icges(5,5) * t265 - Icges(5,6) * t264 + Icges(5,3) * t330;
t213 = rSges(3,1) * t260 - rSges(3,2) * t259 + rSges(3,3) * t333;
t212 = rSges(3,1) * t258 - rSges(3,2) * t257 - rSges(3,3) * t332;
t211 = -rSges(4,1) * t332 - rSges(4,2) * t258 + rSges(4,3) * t257;
t210 = rSges(4,1) * t333 - rSges(4,2) * t260 + rSges(4,3) * t259;
t191 = t234 * t292 + t258 * t291;
t190 = -t234 * t291 + t258 * t292;
t189 = t232 * t292 + t260 * t291;
t188 = -t232 * t291 + t260 * t292;
t187 = qJD(6) * t231 + t230;
t186 = -qJD(6) * t233 + t229;
t181 = rSges(6,1) * t236 + rSges(6,2) * t235 + rSges(6,3) * t264;
t180 = rSges(5,1) * t234 + rSges(5,2) * t233 + rSges(5,3) * t258;
t179 = rSges(5,1) * t232 - rSges(5,2) * t231 + rSges(5,3) * t260;
t178 = Icges(6,1) * t236 + Icges(6,4) * t235 + Icges(6,5) * t264;
t177 = Icges(6,4) * t236 + Icges(6,2) * t235 + Icges(6,6) * t264;
t175 = Icges(5,1) * t234 + Icges(5,4) * t233 + Icges(5,5) * t258;
t174 = Icges(5,1) * t232 - Icges(5,4) * t231 + Icges(5,5) * t260;
t171 = Icges(5,5) * t234 + Icges(5,6) * t233 + Icges(5,3) * t258;
t170 = Icges(5,5) * t232 - Icges(5,6) * t231 + Icges(5,3) * t260;
t169 = pkin(5) * t319 + pkin(9) * t264 + t265 * t337;
t168 = rSges(7,1) * t225 + rSges(7,2) * t224 + rSges(7,3) * t264;
t167 = Icges(7,1) * t225 + Icges(7,4) * t224 + Icges(7,5) * t264;
t166 = Icges(7,4) * t225 + Icges(7,2) * t224 + Icges(7,6) * t264;
t165 = Icges(7,5) * t225 + Icges(7,6) * t224 + Icges(7,3) * t264;
t163 = -t212 * t285 + t249 * t273 + t317;
t162 = t213 * t285 - t249 * t274 + t315;
t161 = rSges(6,1) * t195 + rSges(6,2) * t194 - rSges(6,3) * t233;
t160 = rSges(6,1) * t193 + rSges(6,2) * t192 + rSges(6,3) * t231;
t159 = Icges(6,1) * t195 + Icges(6,4) * t194 - Icges(6,5) * t233;
t158 = Icges(6,1) * t193 + Icges(6,4) * t192 + Icges(6,5) * t231;
t157 = Icges(6,4) * t195 + Icges(6,2) * t194 - Icges(6,6) * t233;
t156 = Icges(6,4) * t193 + Icges(6,2) * t192 + Icges(6,6) * t231;
t153 = t212 * t274 - t213 * t273 + t316;
t152 = rSges(7,1) * t191 + rSges(7,2) * t190 - rSges(7,3) * t233;
t151 = rSges(7,1) * t189 + rSges(7,2) * t188 + rSges(7,3) * t231;
t150 = Icges(7,1) * t191 + Icges(7,4) * t190 - Icges(7,5) * t233;
t149 = Icges(7,1) * t189 + Icges(7,4) * t188 + Icges(7,5) * t231;
t148 = Icges(7,4) * t191 + Icges(7,2) * t190 - Icges(7,6) * t233;
t147 = Icges(7,4) * t189 + Icges(7,2) * t188 + Icges(7,6) * t231;
t146 = Icges(7,5) * t191 + Icges(7,6) * t190 - Icges(7,3) * t233;
t145 = Icges(7,5) * t189 + Icges(7,6) * t188 + Icges(7,3) * t231;
t144 = pkin(5) * t335 - pkin(9) * t233 + t234 * t337;
t143 = pkin(5) * t334 + pkin(9) * t231 + t232 * t337;
t142 = t250 * t273 + (-t211 - t219) * t285 + t314;
t141 = t210 * t285 + (-t250 - t266) * t274 + t313;
t140 = t274 * t211 + (-t210 - t220) * t273 + t312;
t139 = -t180 * t261 + t217 * t229 + t311;
t138 = t179 * t261 - t217 * t230 + t310;
t137 = -t229 * t179 + t230 * t180 + t308;
t136 = t181 * t229 + (-t161 - t185) * t261 + t309;
t135 = t160 * t261 + (-t181 - t221) * t230 + t307;
t134 = t230 * t161 + (-t160 - t184) * t229 + t306;
t133 = (-t144 - t185) * t261 - t152 * t223 + t168 * t186 + t169 * t229 + t309;
t132 = t143 * t261 + t151 * t223 - t168 * t187 + (-t169 - t221) * t230 + t307;
t131 = -t186 * t151 + t187 * t152 + (-t143 - t184) * t229 + t306 + t230 * t144;
t1 = m(3) * (t153 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(5) * (t137 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + m(4) * (t140 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + m(6) * (t134 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + m(7) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + t187 * ((t231 * t145 + t188 * t147 + t189 * t149) * t187 + (t146 * t231 + t148 * t188 + t150 * t189) * t186 + (t165 * t231 + t166 * t188 + t167 * t189) * t223) / 0.2e1 + t186 * ((-t145 * t233 + t147 * t190 + t149 * t191) * t187 + (-t233 * t146 + t190 * t148 + t191 * t150) * t186 + (-t165 * t233 + t166 * t190 + t167 * t191) * t223) / 0.2e1 + m(2) * (t228 ^ 2 + t240 ^ 2 + t241 ^ 2) / 0.2e1 + t223 * ((t145 * t264 + t147 * t224 + t149 * t225) * t187 + (t146 * t264 + t148 * t224 + t150 * t225) * t186 + (t264 * t165 + t224 * t166 + t225 * t167) * t223) / 0.2e1 + m(1) * (t270 ^ 2 + t271 ^ 2 + t272 ^ 2) / 0.2e1 + ((t177 * t194 + t178 * t195 + t214 * t258 + t216 * t234 - t233 * t355) * t261 + (t156 * t194 + t158 * t195 + t170 * t258 + t174 * t234 - t233 * t357) * t230 + (t194 * t157 + t195 * t159 + t258 * t171 + t234 * t175 - t233 * t356) * t229) * t229 / 0.2e1 + ((t177 * t192 + t178 * t193 + t214 * t260 + t216 * t232 + t231 * t355) * t261 + (t192 * t156 + t193 * t158 + t260 * t170 + t232 * t174 + t357 * t231) * t230 + (t157 * t192 + t159 * t193 + t171 * t260 + t175 * t232 + t231 * t356) * t229) * t230 / 0.2e1 + ((t235 * t177 + t236 * t178 + t214 * t330 + t265 * t216 + t264 * t355) * t261 + (t156 * t235 + t158 * t236 + t170 * t330 + t174 * t265 + t264 * t357) * t230 + (t157 * t235 + t159 * t236 + t171 * t330 + t175 * t265 + t264 * t356) * t229) * t261 / 0.2e1 + ((-t257 * t347 + t258 * t346 - t332 * t348) * t285 + (t257 * t354 + t258 * t352 - t332 * t350) * t274 + (t257 * t353 + t258 * t351 - t332 * t349) * t273) * t273 / 0.2e1 + ((-t259 * t347 + t260 * t346 + t333 * t348) * t285 + (t259 * t354 + t260 * t352 + t333 * t350) * t274 + (t259 * t353 + t260 * t351 + t333 * t349) * t273) * t274 / 0.2e1 + ((t273 * t349 + t274 * t350 + t285 * t348) * t301 + ((t304 * t346 + t305 * t347) * t285 + (t304 * t352 - t305 * t354) * t274 + (t304 * t351 - t353 * t305) * t273) * t298) * t285 / 0.2e1 + ((-t277 * t297 + t279 * t300 + Icges(1,4)) * V_base(5) + (-t297 * t278 + t300 * t280 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t300 * t277 + t297 * t279 + Icges(1,2)) * V_base(5) + (t278 * t300 + t280 * t297 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t297 + Icges(2,6) * t300 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t300 - Icges(2,6) * t297 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
