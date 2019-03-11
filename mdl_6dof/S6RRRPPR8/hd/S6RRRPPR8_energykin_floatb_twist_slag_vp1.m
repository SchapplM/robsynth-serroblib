% Calculate kinetic energy for
% S6RRRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
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
% Datum: 2019-03-09 16:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR8_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPPR8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR8_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR8_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR8_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:03:30
% EndTime: 2019-03-09 16:03:34
% DurationCPUTime: 4.65s
% Computational Cost: add. (2275->344), mult. (5148->485), div. (0->0), fcn. (6035->10), ass. (0->153)
t347 = Icges(4,1) + Icges(5,1) + Icges(6,2);
t346 = Icges(6,1) + Icges(4,2) + Icges(5,3);
t345 = -Icges(4,4) - Icges(6,4) + Icges(5,5);
t344 = Icges(5,4) + Icges(4,5) + Icges(6,6);
t343 = Icges(6,5) + Icges(4,6) - Icges(5,6);
t342 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t294 = sin(qJ(2));
t295 = sin(qJ(1));
t297 = cos(qJ(2));
t298 = cos(qJ(1));
t327 = cos(pkin(6));
t311 = t298 * t327;
t260 = t294 * t311 + t295 * t297;
t292 = sin(pkin(6));
t329 = cos(qJ(3));
t315 = t292 * t329;
t328 = sin(qJ(3));
t237 = t260 * t328 + t298 * t315;
t314 = t292 * t328;
t238 = t260 * t329 - t298 * t314;
t259 = t294 * t295 - t297 * t311;
t341 = t237 * t343 - t238 * t344 - t259 * t342;
t312 = t295 * t327;
t262 = -t294 * t312 + t298 * t297;
t239 = t262 * t328 - t295 * t315;
t240 = t262 * t329 + t295 * t314;
t261 = t298 * t294 + t297 * t312;
t340 = t239 * t343 - t240 * t344 - t261 * t342;
t339 = t237 * t346 + t238 * t345 - t259 * t343;
t338 = t239 * t346 + t240 * t345 - t261 * t343;
t337 = t345 * t237 + t238 * t347 + t344 * t259;
t336 = t345 * t239 + t240 * t347 + t344 * t261;
t257 = t294 * t314 - t327 * t329;
t258 = t294 * t315 + t327 * t328;
t324 = t292 * t297;
t335 = t257 * t343 - t258 * t344 + t324 * t342;
t334 = t257 * t346 + t258 * t345 + t324 * t343;
t333 = t345 * t257 + t258 * t347 - t344 * t324;
t326 = Icges(2,4) * t295;
t325 = t292 * t295;
t323 = t292 * t298;
t191 = pkin(3) * t238 + qJ(4) * t237;
t202 = pkin(4) * t238 - qJ(5) * t259;
t322 = -t191 - t202;
t192 = pkin(3) * t240 + qJ(4) * t239;
t203 = pkin(4) * t240 - qJ(5) * t261;
t321 = -t192 - t203;
t226 = pkin(3) * t258 + qJ(4) * t257;
t245 = pkin(4) * t258 + qJ(5) * t324;
t320 = -t226 - t245;
t319 = qJD(2) * t292;
t318 = V_base(5) * pkin(7) + V_base(1);
t313 = t327 * pkin(8);
t271 = t295 * t319 + V_base(4);
t289 = V_base(6) + qJD(1);
t236 = qJD(3) * t261 + t271;
t272 = qJD(2) * t327 + t289;
t270 = -t298 * t319 + V_base(5);
t265 = t295 * pkin(1) - pkin(8) * t323;
t310 = -t265 * t289 + V_base(5) * t313 + t318;
t266 = pkin(1) * t298 + pkin(8) * t325;
t309 = V_base(4) * t265 - t266 * V_base(5) + V_base(3);
t235 = qJD(3) * t259 + t270;
t255 = -qJD(3) * t324 + t272;
t229 = pkin(2) * t260 + pkin(9) * t259;
t264 = (pkin(2) * t294 - pkin(9) * t297) * t292;
t308 = -t229 * t272 + t270 * t264 + t310;
t230 = pkin(2) * t262 + pkin(9) * t261;
t307 = t271 * t229 - t230 * t270 + t309;
t306 = t289 * t266 + V_base(2) + (-t313 - pkin(7)) * V_base(4);
t305 = qJD(4) * t239 + t235 * t226 + t308;
t304 = qJD(4) * t257 + t236 * t191 + t307;
t303 = -qJD(5) * t261 + t235 * t245 + t305;
t302 = t272 * t230 - t271 * t264 + t306;
t301 = qJD(5) * t324 + t236 * t202 + t304;
t300 = qJD(4) * t237 + t255 * t192 + t302;
t299 = -qJD(5) * t259 + t255 * t203 + t300;
t296 = cos(qJ(6));
t293 = sin(qJ(6));
t290 = Icges(2,4) * t298;
t280 = rSges(2,1) * t298 - t295 * rSges(2,2);
t279 = t295 * rSges(2,1) + rSges(2,2) * t298;
t278 = Icges(2,1) * t298 - t326;
t277 = Icges(2,1) * t295 + t290;
t276 = -Icges(2,2) * t295 + t290;
t275 = Icges(2,2) * t298 + t326;
t269 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t268 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t267 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t249 = t327 * rSges(3,3) + (rSges(3,1) * t294 + rSges(3,2) * t297) * t292;
t248 = Icges(3,5) * t327 + (Icges(3,1) * t294 + Icges(3,4) * t297) * t292;
t247 = Icges(3,6) * t327 + (Icges(3,4) * t294 + Icges(3,2) * t297) * t292;
t246 = Icges(3,3) * t327 + (Icges(3,5) * t294 + Icges(3,6) * t297) * t292;
t244 = V_base(5) * rSges(2,3) - t279 * t289 + t318;
t243 = t280 * t289 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t241 = t279 * V_base(4) - t280 * V_base(5) + V_base(3);
t234 = t257 * t296 + t293 * t324;
t233 = -t257 * t293 + t296 * t324;
t228 = qJD(6) * t258 + t255;
t227 = pkin(5) * t257 + pkin(10) * t258;
t225 = rSges(3,1) * t262 - rSges(3,2) * t261 + rSges(3,3) * t325;
t224 = t260 * rSges(3,1) - t259 * rSges(3,2) - rSges(3,3) * t323;
t223 = Icges(3,1) * t262 - Icges(3,4) * t261 + Icges(3,5) * t325;
t222 = Icges(3,1) * t260 - Icges(3,4) * t259 - Icges(3,5) * t323;
t221 = Icges(3,4) * t262 - Icges(3,2) * t261 + Icges(3,6) * t325;
t220 = Icges(3,4) * t260 - Icges(3,2) * t259 - Icges(3,6) * t323;
t219 = Icges(3,5) * t262 - Icges(3,6) * t261 + Icges(3,3) * t325;
t218 = Icges(3,5) * t260 - Icges(3,6) * t259 - Icges(3,3) * t323;
t217 = rSges(4,1) * t258 - rSges(4,2) * t257 - rSges(4,3) * t324;
t216 = rSges(5,1) * t258 - rSges(5,2) * t324 + rSges(5,3) * t257;
t215 = rSges(6,1) * t257 - rSges(6,2) * t258 + rSges(6,3) * t324;
t201 = t239 * t296 - t261 * t293;
t200 = -t239 * t293 - t261 * t296;
t199 = t237 * t296 - t259 * t293;
t198 = -t237 * t293 - t259 * t296;
t196 = qJD(6) * t240 + t236;
t195 = qJD(6) * t238 + t235;
t194 = pkin(5) * t239 + pkin(10) * t240;
t193 = pkin(5) * t237 + pkin(10) * t238;
t187 = rSges(4,1) * t240 - rSges(4,2) * t239 + rSges(4,3) * t261;
t186 = rSges(5,1) * t240 + rSges(5,2) * t261 + rSges(5,3) * t239;
t185 = rSges(6,1) * t239 - rSges(6,2) * t240 - rSges(6,3) * t261;
t184 = rSges(4,1) * t238 - rSges(4,2) * t237 + rSges(4,3) * t259;
t183 = rSges(5,1) * t238 + rSges(5,2) * t259 + rSges(5,3) * t237;
t182 = rSges(6,1) * t237 - rSges(6,2) * t238 - rSges(6,3) * t259;
t162 = rSges(7,1) * t234 + rSges(7,2) * t233 + rSges(7,3) * t258;
t161 = Icges(7,1) * t234 + Icges(7,4) * t233 + Icges(7,5) * t258;
t160 = Icges(7,4) * t234 + Icges(7,2) * t233 + Icges(7,6) * t258;
t159 = Icges(7,5) * t234 + Icges(7,6) * t233 + Icges(7,3) * t258;
t157 = -t224 * t272 + t249 * t270 + t310;
t156 = t272 * t225 - t271 * t249 + t306;
t155 = rSges(7,1) * t201 + rSges(7,2) * t200 + rSges(7,3) * t240;
t154 = rSges(7,1) * t199 + rSges(7,2) * t198 + rSges(7,3) * t238;
t153 = Icges(7,1) * t201 + Icges(7,4) * t200 + Icges(7,5) * t240;
t152 = Icges(7,1) * t199 + Icges(7,4) * t198 + Icges(7,5) * t238;
t151 = Icges(7,4) * t201 + Icges(7,2) * t200 + Icges(7,6) * t240;
t150 = Icges(7,4) * t199 + Icges(7,2) * t198 + Icges(7,6) * t238;
t149 = Icges(7,5) * t201 + Icges(7,6) * t200 + Icges(7,3) * t240;
t148 = Icges(7,5) * t199 + Icges(7,6) * t198 + Icges(7,3) * t238;
t147 = t224 * t271 - t225 * t270 + t309;
t146 = -t184 * t255 + t217 * t235 + t308;
t145 = t255 * t187 - t236 * t217 + t302;
t144 = t184 * t236 - t187 * t235 + t307;
t143 = t216 * t235 + (-t183 - t191) * t255 + t305;
t142 = t255 * t186 + (-t216 - t226) * t236 + t300;
t141 = t183 * t236 + (-t186 - t192) * t235 + t304;
t140 = t215 * t235 + (-t182 + t322) * t255 + t303;
t139 = t255 * t185 + (-t215 + t320) * t236 + t299;
t138 = t182 * t236 + (-t185 + t321) * t235 + t301;
t137 = t303 - t154 * t228 + t162 * t195 + t227 * t235 + (-t193 + t322) * t255;
t136 = t228 * t155 - t196 * t162 + t255 * t194 + (-t227 + t320) * t236 + t299;
t135 = t154 * t196 - t155 * t195 + t193 * t236 + (-t194 + t321) * t235 + t301;
t1 = t272 * (((t221 * t297 + t223 * t294) * t271 + (t220 * t297 + t222 * t294) * t270 + (t247 * t297 + t248 * t294) * t272) * t292 + (t218 * t270 + t219 * t271 + t246 * t272) * t327) / 0.2e1 + t271 * ((t219 * t325 - t261 * t221 + t262 * t223) * t271 + (t218 * t325 - t220 * t261 + t222 * t262) * t270 + (t246 * t325 - t247 * t261 + t248 * t262) * t272) / 0.2e1 + m(3) * (t147 ^ 2 + t156 ^ 2 + t157 ^ 2) / 0.2e1 + m(4) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(6) * (t138 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + m(5) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(7) * (t135 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + t270 * ((-t219 * t323 - t259 * t221 + t260 * t223) * t271 + (-t218 * t323 - t259 * t220 + t260 * t222) * t270 + (-t246 * t323 - t259 * t247 + t260 * t248) * t272) / 0.2e1 + t195 * ((t149 * t238 + t151 * t198 + t153 * t199) * t196 + (t238 * t148 + t198 * t150 + t199 * t152) * t195 + (t159 * t238 + t160 * t198 + t161 * t199) * t228) / 0.2e1 + t196 * ((t240 * t149 + t200 * t151 + t201 * t153) * t196 + (t148 * t240 + t150 * t200 + t152 * t201) * t195 + (t159 * t240 + t160 * t200 + t161 * t201) * t228) / 0.2e1 + m(2) * (t241 ^ 2 + t243 ^ 2 + t244 ^ 2) / 0.2e1 + t228 * ((t149 * t258 + t151 * t233 + t153 * t234) * t196 + (t148 * t258 + t150 * t233 + t152 * t234) * t195 + (t258 * t159 + t233 * t160 + t234 * t161) * t228) / 0.2e1 + m(1) * (t267 ^ 2 + t268 ^ 2 + t269 ^ 2) / 0.2e1 + ((-t295 * t275 + t277 * t298 + Icges(1,4)) * V_base(5) + (-t295 * t276 + t298 * t278 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t298 * t275 + t295 * t277 + Icges(1,2)) * V_base(5) + (t276 * t298 + t295 * t278 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t237 * t334 + t238 * t333 - t259 * t335) * t255 + (t237 * t338 + t238 * t336 - t259 * t340) * t236 + (t339 * t237 + t337 * t238 - t341 * t259) * t235) * t235 / 0.2e1 + ((t239 * t334 + t240 * t333 - t261 * t335) * t255 + (t338 * t239 + t336 * t240 - t340 * t261) * t236 + (t339 * t239 + t337 * t240 - t341 * t261) * t235) * t236 / 0.2e1 + ((t334 * t257 + t333 * t258 + t335 * t324) * t255 + (t257 * t338 + t336 * t258 + t340 * t324) * t236 + (t339 * t257 + t337 * t258 + t324 * t341) * t235) * t255 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t295 + Icges(2,6) * t298) * V_base(5) + (Icges(2,5) * t298 - Icges(2,6) * t295) * V_base(4) + Icges(2,3) * t289 / 0.2e1) * t289;
T  = t1;
