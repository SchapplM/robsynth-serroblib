% Calculate kinetic energy for
% S6PRPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRRR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRR3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:31:36
% EndTime: 2019-03-08 20:31:40
% DurationCPUTime: 4.31s
% Computational Cost: add. (3577->443), mult. (5384->639), div. (0->0), fcn. (6256->14), ass. (0->192)
t313 = sin(pkin(11));
t316 = cos(pkin(11));
t322 = cos(qJ(2));
t317 = cos(pkin(6));
t320 = sin(qJ(2));
t351 = t317 * t320;
t272 = t313 * t322 + t316 * t351;
t312 = sin(pkin(12));
t315 = cos(pkin(12));
t314 = sin(pkin(6));
t354 = t314 * t316;
t245 = -t272 * t312 - t315 * t354;
t339 = t312 * t354;
t246 = t272 * t315 - t339;
t350 = t317 * t322;
t271 = t313 * t320 - t316 * t350;
t189 = Icges(4,5) * t246 + Icges(4,6) * t245 + Icges(4,3) * t271;
t224 = Icges(3,4) * t272 - Icges(3,2) * t271 - Icges(3,6) * t354;
t367 = t189 - t224;
t274 = -t313 * t351 + t316 * t322;
t355 = t313 * t314;
t247 = -t274 * t312 + t315 * t355;
t340 = t312 * t355;
t248 = t274 * t315 + t340;
t273 = t313 * t350 + t316 * t320;
t190 = Icges(4,5) * t248 + Icges(4,6) * t247 + Icges(4,3) * t273;
t225 = Icges(3,4) * t274 - Icges(3,2) * t273 + Icges(3,6) * t355;
t366 = t190 - t225;
t353 = t314 * t320;
t269 = -t312 * t353 + t315 * t317;
t356 = t312 * t317;
t270 = t315 * t353 + t356;
t352 = t314 * t322;
t219 = Icges(4,5) * t270 + Icges(4,6) * t269 - Icges(4,3) * t352;
t258 = Icges(3,6) * t317 + (Icges(3,4) * t320 + Icges(3,2) * t322) * t314;
t365 = t219 - t258;
t359 = pkin(7) * t317;
t358 = t315 * pkin(3);
t357 = Icges(2,4) * t313;
t311 = pkin(12) + qJ(4);
t306 = cos(t311);
t348 = pkin(4) * t306;
t346 = qJD(2) * t314;
t345 = V_base(5) * qJ(1) + V_base(1);
t341 = qJD(1) + V_base(3);
t288 = t313 * t346 + V_base(4);
t299 = qJD(2) * t317 + V_base(6);
t338 = qJ(5) + t311;
t305 = sin(t311);
t337 = pkin(4) * t305;
t250 = qJD(4) * t273 + t288;
t336 = cos(t338);
t217 = qJD(5) * t273 + t250;
t287 = -t316 * t346 + V_base(5);
t335 = t314 * t336;
t249 = qJD(4) * t271 + t287;
t279 = pkin(1) * t313 - pkin(7) * t354;
t334 = -t279 * V_base(6) + V_base(5) * t359 + t345;
t280 = pkin(1) * t316 + pkin(7) * t355;
t333 = V_base(4) * t279 - V_base(5) * t280 + t341;
t216 = qJD(5) * t271 + t249;
t254 = (-qJD(4) - qJD(5)) * t352 + t299;
t332 = V_base(6) * t280 + V_base(2) + (-qJ(1) - t359) * V_base(4);
t278 = (pkin(2) * t320 - qJ(3) * t322) * t314;
t331 = qJD(3) * t273 + t287 * t278 + t334;
t233 = pkin(2) * t274 + qJ(3) * t273;
t330 = qJD(3) * t271 + t299 * t233 + t332;
t232 = pkin(2) * t272 + qJ(3) * t271;
t329 = -qJD(3) * t352 + t288 * t232 + t333;
t185 = -pkin(3) * t339 + pkin(8) * t271 + t272 * t358;
t229 = pkin(3) * t356 + (-pkin(8) * t322 + t320 * t358) * t314;
t328 = t287 * t229 + (-t185 - t232) * t299 + t331;
t186 = pkin(3) * t340 + pkin(8) * t273 + t274 * t358;
t327 = t299 * t186 + (-t229 - t278) * t288 + t330;
t326 = t288 * t185 + (-t186 - t233) * t287 + t329;
t158 = pkin(9) * t271 + t272 * t348 - t337 * t354;
t199 = t337 * t317 + (-pkin(9) * t322 + t320 * t348) * t314;
t275 = -qJD(4) * t352 + t299;
t325 = -t158 * t275 + t249 * t199 + t328;
t159 = pkin(9) * t273 + t274 * t348 + t337 * t355;
t324 = t275 * t159 - t199 * t250 + t327;
t323 = t250 * t158 - t249 * t159 + t326;
t321 = cos(qJ(6));
t319 = sin(qJ(6));
t307 = Icges(2,4) * t316;
t302 = sin(t338);
t296 = rSges(2,1) * t316 - rSges(2,2) * t313;
t295 = rSges(2,1) * t313 + rSges(2,2) * t316;
t294 = Icges(2,1) * t316 - t357;
t293 = Icges(2,1) * t313 + t307;
t292 = -Icges(2,2) * t313 + t307;
t291 = Icges(2,2) * t316 + t357;
t286 = -rSges(1,1) * V_base(5) + V_base(4) * rSges(1,2) + V_base(3);
t285 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t284 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t262 = t317 * rSges(3,3) + (rSges(3,1) * t320 + rSges(3,2) * t322) * t314;
t261 = t305 * t317 + t306 * t353;
t260 = -t305 * t353 + t306 * t317;
t259 = Icges(3,5) * t317 + (Icges(3,1) * t320 + Icges(3,4) * t322) * t314;
t257 = Icges(3,3) * t317 + (Icges(3,5) * t320 + Icges(3,6) * t322) * t314;
t256 = t317 * t302 + t320 * t335;
t255 = t302 * t353 - t317 * t336;
t253 = V_base(5) * rSges(2,3) - t295 * V_base(6) + t345;
t252 = t296 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t244 = t295 * V_base(4) - t296 * V_base(5) + t341;
t243 = t274 * t306 + t305 * t355;
t242 = -t274 * t305 + t306 * t355;
t241 = t272 * t306 - t305 * t354;
t240 = -t272 * t305 - t306 * t354;
t239 = t256 * t321 - t319 * t352;
t238 = -t256 * t319 - t321 * t352;
t237 = t274 * t336 + t302 * t355;
t236 = t274 * t302 - t313 * t335;
t235 = t272 * t336 - t302 * t354;
t234 = t272 * t302 + t316 * t335;
t231 = rSges(3,1) * t274 - rSges(3,2) * t273 + rSges(3,3) * t355;
t230 = rSges(3,1) * t272 - rSges(3,2) * t271 - rSges(3,3) * t354;
t228 = t270 * rSges(4,1) + t269 * rSges(4,2) - rSges(4,3) * t352;
t227 = Icges(3,1) * t274 - Icges(3,4) * t273 + Icges(3,5) * t355;
t226 = Icges(3,1) * t272 - Icges(3,4) * t271 - Icges(3,5) * t354;
t223 = Icges(3,5) * t274 - Icges(3,6) * t273 + Icges(3,3) * t355;
t222 = Icges(3,5) * t272 - Icges(3,6) * t271 - Icges(3,3) * t354;
t221 = Icges(4,1) * t270 + Icges(4,4) * t269 - Icges(4,5) * t352;
t220 = Icges(4,4) * t270 + Icges(4,2) * t269 - Icges(4,6) * t352;
t214 = qJD(6) * t255 + t254;
t213 = pkin(5) * t256 + pkin(10) * t255;
t212 = t261 * rSges(5,1) + t260 * rSges(5,2) - rSges(5,3) * t352;
t211 = Icges(5,1) * t261 + Icges(5,4) * t260 - Icges(5,5) * t352;
t210 = Icges(5,4) * t261 + Icges(5,2) * t260 - Icges(5,6) * t352;
t209 = Icges(5,5) * t261 + Icges(5,6) * t260 - Icges(5,3) * t352;
t207 = t256 * rSges(6,1) - t255 * rSges(6,2) - rSges(6,3) * t352;
t206 = Icges(6,1) * t256 - Icges(6,4) * t255 - Icges(6,5) * t352;
t205 = Icges(6,4) * t256 - Icges(6,2) * t255 - Icges(6,6) * t352;
t204 = Icges(6,5) * t256 - Icges(6,6) * t255 - Icges(6,3) * t352;
t203 = t237 * t321 + t273 * t319;
t202 = -t237 * t319 + t273 * t321;
t201 = t235 * t321 + t271 * t319;
t200 = -t235 * t319 + t271 * t321;
t198 = pkin(5) * t237 + pkin(10) * t236;
t197 = pkin(5) * t235 + pkin(10) * t234;
t196 = rSges(4,1) * t248 + rSges(4,2) * t247 + rSges(4,3) * t273;
t195 = rSges(4,1) * t246 + rSges(4,2) * t245 + rSges(4,3) * t271;
t194 = Icges(4,1) * t248 + Icges(4,4) * t247 + Icges(4,5) * t273;
t193 = Icges(4,1) * t246 + Icges(4,4) * t245 + Icges(4,5) * t271;
t192 = Icges(4,4) * t248 + Icges(4,2) * t247 + Icges(4,6) * t273;
t191 = Icges(4,4) * t246 + Icges(4,2) * t245 + Icges(4,6) * t271;
t188 = qJD(6) * t236 + t217;
t187 = qJD(6) * t234 + t216;
t184 = rSges(5,1) * t243 + rSges(5,2) * t242 + rSges(5,3) * t273;
t183 = rSges(5,1) * t241 + rSges(5,2) * t240 + rSges(5,3) * t271;
t182 = Icges(5,1) * t243 + Icges(5,4) * t242 + Icges(5,5) * t273;
t181 = Icges(5,1) * t241 + Icges(5,4) * t240 + Icges(5,5) * t271;
t180 = Icges(5,4) * t243 + Icges(5,2) * t242 + Icges(5,6) * t273;
t179 = Icges(5,4) * t241 + Icges(5,2) * t240 + Icges(5,6) * t271;
t178 = Icges(5,5) * t243 + Icges(5,6) * t242 + Icges(5,3) * t273;
t177 = Icges(5,5) * t241 + Icges(5,6) * t240 + Icges(5,3) * t271;
t176 = rSges(6,1) * t237 - rSges(6,2) * t236 + rSges(6,3) * t273;
t175 = rSges(6,1) * t235 - rSges(6,2) * t234 + rSges(6,3) * t271;
t174 = Icges(6,1) * t237 - Icges(6,4) * t236 + Icges(6,5) * t273;
t173 = Icges(6,1) * t235 - Icges(6,4) * t234 + Icges(6,5) * t271;
t172 = Icges(6,4) * t237 - Icges(6,2) * t236 + Icges(6,6) * t273;
t171 = Icges(6,4) * t235 - Icges(6,2) * t234 + Icges(6,6) * t271;
t170 = Icges(6,5) * t237 - Icges(6,6) * t236 + Icges(6,3) * t273;
t169 = Icges(6,5) * t235 - Icges(6,6) * t234 + Icges(6,3) * t271;
t167 = rSges(7,1) * t239 + rSges(7,2) * t238 + rSges(7,3) * t255;
t166 = Icges(7,1) * t239 + Icges(7,4) * t238 + Icges(7,5) * t255;
t165 = Icges(7,4) * t239 + Icges(7,2) * t238 + Icges(7,6) * t255;
t164 = Icges(7,5) * t239 + Icges(7,6) * t238 + Icges(7,3) * t255;
t161 = -t230 * t299 + t262 * t287 + t334;
t160 = t231 * t299 - t262 * t288 + t332;
t156 = t230 * t288 - t231 * t287 + t333;
t154 = rSges(7,1) * t203 + rSges(7,2) * t202 + rSges(7,3) * t236;
t153 = rSges(7,1) * t201 + rSges(7,2) * t200 + rSges(7,3) * t234;
t152 = Icges(7,1) * t203 + Icges(7,4) * t202 + Icges(7,5) * t236;
t151 = Icges(7,1) * t201 + Icges(7,4) * t200 + Icges(7,5) * t234;
t150 = Icges(7,4) * t203 + Icges(7,2) * t202 + Icges(7,6) * t236;
t149 = Icges(7,4) * t201 + Icges(7,2) * t200 + Icges(7,6) * t234;
t148 = Icges(7,5) * t203 + Icges(7,6) * t202 + Icges(7,3) * t236;
t147 = Icges(7,5) * t201 + Icges(7,6) * t200 + Icges(7,3) * t234;
t146 = t228 * t287 + (-t195 - t232) * t299 + t331;
t145 = t196 * t299 + (-t228 - t278) * t288 + t330;
t144 = t288 * t195 + (-t196 - t233) * t287 + t329;
t143 = -t183 * t275 + t212 * t249 + t328;
t142 = t184 * t275 - t212 * t250 + t327;
t141 = t250 * t183 - t249 * t184 + t326;
t140 = -t175 * t254 + t207 * t216 + t325;
t139 = t176 * t254 - t207 * t217 + t324;
t138 = t217 * t175 - t216 * t176 + t323;
t137 = -t153 * t214 + t167 * t187 - t197 * t254 + t213 * t216 + t325;
t136 = t154 * t214 - t167 * t188 + t198 * t254 - t213 * t217 + t324;
t135 = t188 * t153 - t187 * t154 + t217 * t197 - t216 * t198 + t323;
t1 = m(2) * (t244 ^ 2 + t252 ^ 2 + t253 ^ 2) / 0.2e1 + t187 * ((t148 * t234 + t150 * t200 + t152 * t201) * t188 + (t234 * t147 + t200 * t149 + t201 * t151) * t187 + (t164 * t234 + t165 * t200 + t166 * t201) * t214) / 0.2e1 + t188 * ((t236 * t148 + t202 * t150 + t203 * t152) * t188 + (t147 * t236 + t149 * t202 + t151 * t203) * t187 + (t164 * t236 + t165 * t202 + t166 * t203) * t214) / 0.2e1 + t275 * ((-t178 * t352 + t260 * t180 + t261 * t182) * t250 + (-t177 * t352 + t260 * t179 + t261 * t181) * t249 + (-t209 * t352 + t260 * t210 + t261 * t211) * t275) / 0.2e1 + t254 * ((-t170 * t352 - t255 * t172 + t256 * t174) * t217 + (-t169 * t352 - t255 * t171 + t256 * t173) * t216 + (-t204 * t352 - t255 * t205 + t256 * t206) * t254) / 0.2e1 + m(3) * (t156 ^ 2 + t160 ^ 2 + t161 ^ 2) / 0.2e1 + m(5) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(4) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(6) * (t138 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + m(7) * (t135 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + t214 * ((t148 * t255 + t150 * t238 + t152 * t239) * t188 + (t147 * t255 + t149 * t238 + t151 * t239) * t187 + (t255 * t164 + t238 * t165 + t239 * t166) * t214) / 0.2e1 + t216 * ((t170 * t271 - t172 * t234 + t174 * t235) * t217 + (t271 * t169 - t234 * t171 + t235 * t173) * t216 + (t204 * t271 - t205 * t234 + t206 * t235) * t254) / 0.2e1 + t217 * ((t273 * t170 - t236 * t172 + t237 * t174) * t217 + (t169 * t273 - t171 * t236 + t173 * t237) * t216 + (t204 * t273 - t205 * t236 + t206 * t237) * t254) / 0.2e1 + t249 * ((t178 * t271 + t180 * t240 + t182 * t241) * t250 + (t177 * t271 + t179 * t240 + t181 * t241) * t249 + (t209 * t271 + t210 * t240 + t211 * t241) * t275) / 0.2e1 + t250 * ((t178 * t273 + t180 * t242 + t182 * t243) * t250 + (t177 * t273 + t179 * t242 + t181 * t243) * t249 + (t209 * t273 + t210 * t242 + t211 * t243) * t275) / 0.2e1 + m(1) * (t284 ^ 2 + t285 ^ 2 + t286 ^ 2) / 0.2e1 + ((t220 * t245 + t221 * t246 - t257 * t354 + t259 * t272 + t271 * t365) * t299 + (t192 * t245 + t194 * t246 - t223 * t354 + t227 * t272 + t271 * t366) * t288 + (t191 * t245 + t193 * t246 - t222 * t354 + t226 * t272 + t367 * t271) * t287) * t287 / 0.2e1 + ((t220 * t247 + t221 * t248 + t257 * t355 + t259 * t274 + t273 * t365) * t299 + (t192 * t247 + t194 * t248 + t223 * t355 + t227 * t274 + t273 * t366) * t288 + (t191 * t247 + t193 * t248 + t222 * t355 + t226 * t274 + t273 * t367) * t287) * t288 / 0.2e1 + ((t222 * t287 + t223 * t288 + t257 * t299) * t317 + ((t225 * t322 + t227 * t320) * t288 + (t224 * t322 + t226 * t320) * t287 + (t258 * t322 + t259 * t320) * t299) * t314 + (-t190 * t352 + t269 * t192 + t270 * t194) * t288 + (-t189 * t352 + t269 * t191 + t270 * t193) * t287 + (-t219 * t352 + t269 * t220 + t270 * t221) * t299) * t299 / 0.2e1 + ((-t291 * t313 + t293 * t316 + Icges(1,4)) * V_base(5) + (-t292 * t313 + t294 * t316 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t291 * t316 + t293 * t313 + Icges(1,2)) * V_base(5) + (t292 * t316 + t294 * t313 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t316 - Icges(2,6) * t313 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t313 + Icges(2,6) * t316 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
