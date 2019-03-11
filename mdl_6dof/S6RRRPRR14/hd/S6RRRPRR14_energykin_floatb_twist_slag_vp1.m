% Calculate kinetic energy for
% S6RRRPRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 20:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR14_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR14_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRR14_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR14_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR14_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR14_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:12:17
% EndTime: 2019-03-09 20:12:21
% DurationCPUTime: 4.30s
% Computational Cost: add. (2775->396), mult. (5944->571), div. (0->0), fcn. (7097->12), ass. (0->177)
t361 = Icges(4,1) + Icges(5,2);
t360 = Icges(5,1) + Icges(4,3);
t359 = -Icges(4,4) - Icges(5,6);
t358 = Icges(5,4) - Icges(4,5);
t357 = Icges(5,5) - Icges(4,6);
t356 = Icges(4,2) + Icges(5,3);
t302 = cos(pkin(6));
t305 = sin(qJ(1));
t307 = cos(qJ(2));
t329 = t305 * t307;
t304 = sin(qJ(2));
t308 = cos(qJ(1));
t330 = t304 * t308;
t266 = t302 * t330 + t329;
t301 = sin(pkin(6));
t343 = cos(qJ(3));
t323 = t301 * t343;
t342 = sin(qJ(3));
t244 = t266 * t342 + t308 * t323;
t322 = t301 * t342;
t245 = t266 * t343 - t308 * t322;
t328 = t307 * t308;
t331 = t304 * t305;
t265 = -t302 * t328 + t331;
t355 = t244 * t356 + t245 * t359 + t265 * t357;
t268 = -t302 * t331 + t328;
t246 = t268 * t342 - t305 * t323;
t247 = t268 * t343 + t305 * t322;
t267 = t302 * t329 + t330;
t354 = t246 * t356 + t247 * t359 + t267 * t357;
t353 = t244 * t357 - t245 * t358 + t265 * t360;
t352 = t246 * t357 - t247 * t358 + t267 * t360;
t351 = t244 * t359 + t245 * t361 - t265 * t358;
t350 = t246 * t359 + t247 * t361 - t267 * t358;
t263 = -t302 * t343 + t304 * t322;
t264 = t302 * t342 + t304 * t323;
t333 = t301 * t307;
t349 = t263 * t356 + t264 * t359 - t333 * t357;
t348 = t263 * t359 + t264 * t361 + t333 * t358;
t347 = t263 * t357 - t264 * t358 - t333 * t360;
t341 = pkin(8) * t302;
t306 = cos(qJ(5));
t340 = pkin(5) * t306;
t338 = Icges(2,4) * t305;
t303 = sin(qJ(5));
t337 = t244 * t303;
t336 = t246 * t303;
t335 = t263 * t303;
t334 = t301 * t305;
t332 = t301 * t308;
t327 = qJD(2) * t301;
t326 = V_base(5) * pkin(7) + V_base(1);
t277 = t305 * t327 + V_base(4);
t295 = V_base(6) + qJD(1);
t243 = qJD(3) * t267 + t277;
t278 = qJD(2) * t302 + t295;
t199 = qJD(5) * t247 + t243;
t276 = -t308 * t327 + V_base(5);
t271 = pkin(1) * t305 - pkin(8) * t332;
t321 = -t271 * t295 + t341 * V_base(5) + t326;
t272 = pkin(1) * t308 + pkin(8) * t334;
t320 = t271 * V_base(4) - t272 * V_base(5) + V_base(3);
t242 = qJD(3) * t265 + t276;
t198 = qJD(5) * t245 + t242;
t261 = -qJD(3) * t333 + t278;
t319 = t295 * t272 + V_base(2) + (-pkin(7) - t341) * V_base(4);
t231 = qJD(5) * t264 + t261;
t232 = pkin(2) * t266 + pkin(9) * t265;
t270 = (pkin(2) * t304 - pkin(9) * t307) * t301;
t318 = -t232 * t278 + t270 * t276 + t321;
t233 = pkin(2) * t268 + pkin(9) * t267;
t317 = t232 * t277 - t233 * t276 + t320;
t230 = pkin(3) * t264 + qJ(4) * t263;
t316 = qJD(4) * t246 + t230 * t242 + t318;
t196 = pkin(3) * t245 + qJ(4) * t244;
t315 = qJD(4) * t263 + t196 * t243 + t317;
t314 = t233 * t278 - t270 * t277 + t319;
t197 = pkin(3) * t247 + qJ(4) * t246;
t313 = qJD(4) * t244 + t197 * t261 + t314;
t210 = pkin(4) * t265 + pkin(10) * t245;
t252 = -pkin(4) * t333 + pkin(10) * t264;
t312 = t242 * t252 + (-t196 - t210) * t261 + t316;
t211 = pkin(4) * t267 + pkin(10) * t247;
t311 = t243 * t210 + (-t197 - t211) * t242 + t315;
t310 = t261 * t211 + (-t230 - t252) * t243 + t313;
t300 = qJ(5) + qJ(6);
t298 = Icges(2,4) * t308;
t297 = cos(t300);
t296 = sin(t300);
t286 = rSges(2,1) * t308 - rSges(2,2) * t305;
t285 = rSges(2,1) * t305 + rSges(2,2) * t308;
t284 = Icges(2,1) * t308 - t338;
t283 = Icges(2,1) * t305 + t298;
t282 = -Icges(2,2) * t305 + t298;
t281 = Icges(2,2) * t308 + t338;
t275 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t274 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t273 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t256 = rSges(3,3) * t302 + (rSges(3,1) * t304 + rSges(3,2) * t307) * t301;
t255 = Icges(3,5) * t302 + (Icges(3,1) * t304 + Icges(3,4) * t307) * t301;
t254 = Icges(3,6) * t302 + (Icges(3,4) * t304 + Icges(3,2) * t307) * t301;
t253 = Icges(3,3) * t302 + (Icges(3,5) * t304 + Icges(3,6) * t307) * t301;
t251 = V_base(5) * rSges(2,3) - t285 * t295 + t326;
t250 = t286 * t295 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t248 = t285 * V_base(4) - t286 * V_base(5) + V_base(3);
t241 = -t306 * t333 + t335;
t240 = t263 * t306 + t303 * t333;
t235 = t263 * t296 - t297 * t333;
t234 = t263 * t297 + t296 * t333;
t229 = rSges(3,1) * t268 - rSges(3,2) * t267 + rSges(3,3) * t334;
t228 = rSges(3,1) * t266 - rSges(3,2) * t265 - rSges(3,3) * t332;
t227 = Icges(3,1) * t268 - Icges(3,4) * t267 + Icges(3,5) * t334;
t226 = Icges(3,1) * t266 - Icges(3,4) * t265 - Icges(3,5) * t332;
t225 = Icges(3,4) * t268 - Icges(3,2) * t267 + Icges(3,6) * t334;
t224 = Icges(3,4) * t266 - Icges(3,2) * t265 - Icges(3,6) * t332;
t223 = Icges(3,5) * t268 - Icges(3,6) * t267 + Icges(3,3) * t334;
t222 = Icges(3,5) * t266 - Icges(3,6) * t265 - Icges(3,3) * t332;
t221 = rSges(4,1) * t264 - rSges(4,2) * t263 - rSges(4,3) * t333;
t220 = -rSges(5,1) * t333 - rSges(5,2) * t264 + rSges(5,3) * t263;
t209 = t267 * t306 + t336;
t208 = t246 * t306 - t267 * t303;
t207 = t265 * t306 + t337;
t206 = t244 * t306 - t265 * t303;
t205 = qJD(6) * t264 + t231;
t203 = t246 * t296 + t267 * t297;
t202 = t246 * t297 - t267 * t296;
t201 = t244 * t296 + t265 * t297;
t200 = t244 * t297 - t265 * t296;
t193 = pkin(5) * t335 + pkin(11) * t264 - t333 * t340;
t191 = rSges(4,1) * t247 - rSges(4,2) * t246 + rSges(4,3) * t267;
t190 = rSges(4,1) * t245 - rSges(4,2) * t244 + rSges(4,3) * t265;
t189 = rSges(5,1) * t267 - rSges(5,2) * t247 + rSges(5,3) * t246;
t188 = rSges(5,1) * t265 - rSges(5,2) * t245 + rSges(5,3) * t244;
t174 = rSges(6,1) * t241 + rSges(6,2) * t240 + rSges(6,3) * t264;
t173 = Icges(6,1) * t241 + Icges(6,4) * t240 + Icges(6,5) * t264;
t172 = Icges(6,4) * t241 + Icges(6,2) * t240 + Icges(6,6) * t264;
t171 = Icges(6,5) * t241 + Icges(6,6) * t240 + Icges(6,3) * t264;
t170 = rSges(7,1) * t235 + rSges(7,2) * t234 + rSges(7,3) * t264;
t169 = Icges(7,1) * t235 + Icges(7,4) * t234 + Icges(7,5) * t264;
t168 = Icges(7,4) * t235 + Icges(7,2) * t234 + Icges(7,6) * t264;
t167 = Icges(7,5) * t235 + Icges(7,6) * t234 + Icges(7,3) * t264;
t166 = qJD(6) * t247 + t199;
t165 = qJD(6) * t245 + t198;
t163 = -t228 * t278 + t256 * t276 + t321;
t162 = t229 * t278 - t256 * t277 + t319;
t161 = pkin(5) * t336 + pkin(11) * t247 + t267 * t340;
t160 = pkin(5) * t337 + pkin(11) * t245 + t265 * t340;
t159 = rSges(6,1) * t209 + rSges(6,2) * t208 + rSges(6,3) * t247;
t158 = rSges(6,1) * t207 + rSges(6,2) * t206 + rSges(6,3) * t245;
t157 = Icges(6,1) * t209 + Icges(6,4) * t208 + Icges(6,5) * t247;
t156 = Icges(6,1) * t207 + Icges(6,4) * t206 + Icges(6,5) * t245;
t155 = Icges(6,4) * t209 + Icges(6,2) * t208 + Icges(6,6) * t247;
t154 = Icges(6,4) * t207 + Icges(6,2) * t206 + Icges(6,6) * t245;
t153 = Icges(6,5) * t209 + Icges(6,6) * t208 + Icges(6,3) * t247;
t152 = Icges(6,5) * t207 + Icges(6,6) * t206 + Icges(6,3) * t245;
t151 = t228 * t277 - t229 * t276 + t320;
t150 = rSges(7,1) * t203 + rSges(7,2) * t202 + rSges(7,3) * t247;
t149 = rSges(7,1) * t201 + rSges(7,2) * t200 + rSges(7,3) * t245;
t148 = Icges(7,1) * t203 + Icges(7,4) * t202 + Icges(7,5) * t247;
t147 = Icges(7,1) * t201 + Icges(7,4) * t200 + Icges(7,5) * t245;
t146 = Icges(7,4) * t203 + Icges(7,2) * t202 + Icges(7,6) * t247;
t145 = Icges(7,4) * t201 + Icges(7,2) * t200 + Icges(7,6) * t245;
t144 = Icges(7,5) * t203 + Icges(7,6) * t202 + Icges(7,3) * t247;
t143 = Icges(7,5) * t201 + Icges(7,6) * t200 + Icges(7,3) * t245;
t142 = -t190 * t261 + t221 * t242 + t318;
t141 = t191 * t261 - t221 * t243 + t314;
t140 = t190 * t243 - t191 * t242 + t317;
t139 = t220 * t242 + (-t188 - t196) * t261 + t316;
t138 = t189 * t261 + (-t220 - t230) * t243 + t313;
t137 = t188 * t243 + (-t189 - t197) * t242 + t315;
t136 = -t158 * t231 + t174 * t198 + t312;
t135 = t159 * t231 - t174 * t199 + t310;
t134 = t158 * t199 - t159 * t198 + t311;
t133 = -t149 * t205 - t160 * t231 + t165 * t170 + t193 * t198 + t312;
t132 = t150 * t205 + t161 * t231 - t166 * t170 - t193 * t199 + t310;
t131 = t149 * t166 - t150 * t165 + t160 * t199 - t161 * t198 + t311;
t1 = m(3) * (t151 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(1) * (t273 ^ 2 + t274 ^ 2 + t275 ^ 2) / 0.2e1 + t231 * ((t153 * t264 + t155 * t240 + t157 * t241) * t199 + (t152 * t264 + t154 * t240 + t156 * t241) * t198 + (t171 * t264 + t172 * t240 + t173 * t241) * t231) / 0.2e1 + t205 * ((t144 * t264 + t146 * t234 + t148 * t235) * t166 + (t143 * t264 + t145 * t234 + t147 * t235) * t165 + (t167 * t264 + t168 * t234 + t169 * t235) * t205) / 0.2e1 + t199 * ((t153 * t247 + t155 * t208 + t157 * t209) * t199 + (t152 * t247 + t154 * t208 + t156 * t209) * t198 + (t171 * t247 + t172 * t208 + t173 * t209) * t231) / 0.2e1 + t166 * ((t144 * t247 + t146 * t202 + t148 * t203) * t166 + (t143 * t247 + t145 * t202 + t147 * t203) * t165 + (t167 * t247 + t168 * t202 + t169 * t203) * t205) / 0.2e1 + m(2) * (t248 ^ 2 + t250 ^ 2 + t251 ^ 2) / 0.2e1 + t198 * ((t153 * t245 + t155 * t206 + t157 * t207) * t199 + (t152 * t245 + t154 * t206 + t156 * t207) * t198 + (t171 * t245 + t172 * t206 + t173 * t207) * t231) / 0.2e1 + t165 * ((t144 * t245 + t146 * t200 + t148 * t201) * t166 + (t143 * t245 + t145 * t200 + t147 * t201) * t165 + (t167 * t245 + t168 * t200 + t169 * t201) * t205) / 0.2e1 + t278 * ((t222 * t276 + t223 * t277 + t253 * t278) * t302 + ((t225 * t307 + t227 * t304) * t277 + (t224 * t307 + t226 * t304) * t276 + (t254 * t307 + t255 * t304) * t278) * t301) / 0.2e1 + m(7) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + m(6) * (t134 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + m(5) * (t137 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + m(4) * (t140 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + t276 * ((-t223 * t332 - t225 * t265 + t227 * t266) * t277 + (-t222 * t332 - t224 * t265 + t226 * t266) * t276 + (-t253 * t332 - t254 * t265 + t255 * t266) * t278) / 0.2e1 + t277 * ((t223 * t334 - t225 * t267 + t227 * t268) * t277 + (t222 * t334 - t224 * t267 + t226 * t268) * t276 + (t253 * t334 - t254 * t267 + t255 * t268) * t278) / 0.2e1 + ((t244 * t349 + t245 * t348 + t265 * t347) * t261 + (t244 * t354 + t245 * t350 + t265 * t352) * t243 + (t355 * t244 + t351 * t245 + t353 * t265) * t242) * t242 / 0.2e1 + ((t246 * t349 + t247 * t348 + t267 * t347) * t261 + (t354 * t246 + t350 * t247 + t352 * t267) * t243 + (t246 * t355 + t247 * t351 + t267 * t353) * t242) * t243 / 0.2e1 + ((t349 * t263 + t348 * t264 - t347 * t333) * t261 + (t263 * t354 + t264 * t350 - t333 * t352) * t243 + (t263 * t355 + t264 * t351 - t333 * t353) * t242) * t261 / 0.2e1 + ((-t281 * t305 + t283 * t308 + Icges(1,4)) * V_base(5) + (-t282 * t305 + t284 * t308 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t281 * t308 + t283 * t305 + Icges(1,2)) * V_base(5) + (t282 * t308 + t284 * t305 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t305 + Icges(2,6) * t308) * V_base(5) + (Icges(2,5) * t308 - Icges(2,6) * t305) * V_base(4) + Icges(2,3) * t295 / 0.2e1) * t295;
T  = t1;
