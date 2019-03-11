% Calculate kinetic energy for
% S6RRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR8_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRR8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR8_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR8_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR8_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:44:00
% EndTime: 2019-03-09 18:44:05
% DurationCPUTime: 4.94s
% Computational Cost: add. (3750->437), mult. (6058->631), div. (0->0), fcn. (7208->14), ass. (0->192)
t369 = Icges(4,3) + Icges(5,3);
t312 = cos(pkin(6));
t317 = sin(qJ(1));
t320 = cos(qJ(2));
t346 = t317 * t320;
t316 = sin(qJ(2));
t321 = cos(qJ(1));
t347 = t316 * t321;
t275 = t312 * t347 + t346;
t343 = qJ(3) + pkin(12);
t304 = sin(t343);
t311 = sin(pkin(6));
t336 = cos(t343);
t335 = t311 * t336;
t245 = t275 * t304 + t321 * t335;
t350 = t311 * t321;
t246 = t275 * t336 - t304 * t350;
t315 = sin(qJ(3));
t319 = cos(qJ(3));
t251 = -t275 * t315 - t319 * t350;
t337 = t315 * t350;
t252 = t275 * t319 - t337;
t345 = t320 * t321;
t348 = t316 * t317;
t274 = -t312 * t345 + t348;
t368 = Icges(4,5) * t252 + Icges(5,5) * t246 + Icges(4,6) * t251 - Icges(5,6) * t245 + t274 * t369;
t277 = -t312 * t348 + t345;
t247 = t277 * t304 - t317 * t335;
t353 = t311 * t317;
t248 = t277 * t336 + t304 * t353;
t352 = t311 * t319;
t253 = -t277 * t315 + t317 * t352;
t338 = t315 * t353;
t254 = t277 * t319 + t338;
t276 = t312 * t346 + t347;
t367 = Icges(4,5) * t254 + Icges(5,5) * t248 + Icges(4,6) * t253 - Icges(5,6) * t247 + t276 * t369;
t354 = t311 * t316;
t263 = t304 * t354 - t312 * t336;
t264 = t304 * t312 + t316 * t335;
t272 = t312 * t319 - t315 * t354;
t349 = t312 * t315;
t273 = t316 * t352 + t349;
t351 = t311 * t320;
t366 = Icges(4,5) * t273 + Icges(5,5) * t264 + Icges(4,6) * t272 - Icges(5,6) * t263 - t351 * t369;
t362 = pkin(8) * t312;
t361 = pkin(3) * t319;
t318 = cos(qJ(5));
t360 = pkin(5) * t318;
t357 = Icges(2,4) * t317;
t314 = sin(qJ(5));
t356 = t274 * t314;
t355 = t276 * t314;
t344 = qJD(2) * t311;
t342 = V_base(5) * pkin(7) + V_base(1);
t339 = t314 * t351;
t287 = t317 * t344 + V_base(4);
t305 = V_base(6) + qJD(1);
t250 = qJD(3) * t276 + t287;
t288 = qJD(2) * t312 + t305;
t206 = qJD(5) * t247 + t250;
t286 = -t321 * t344 + V_base(5);
t280 = pkin(1) * t317 - pkin(8) * t350;
t334 = -t280 * t305 + t362 * V_base(5) + t342;
t281 = pkin(1) * t321 + pkin(8) * t353;
t333 = t280 * V_base(4) - t281 * V_base(5) + V_base(3);
t249 = qJD(3) * t274 + t286;
t205 = qJD(5) * t245 + t249;
t270 = -qJD(3) * t351 + t288;
t332 = t305 * t281 + V_base(2) + (-pkin(7) - t362) * V_base(4);
t236 = qJD(5) * t263 + t270;
t241 = pkin(2) * t275 + pkin(9) * t274;
t279 = (pkin(2) * t316 - pkin(9) * t320) * t311;
t331 = -t241 * t288 + t279 * t286 + t334;
t242 = pkin(2) * t277 + pkin(9) * t276;
t330 = t241 * t287 - t242 * t286 + t333;
t226 = pkin(3) * t349 + (-qJ(4) * t320 + t316 * t361) * t311;
t329 = qJD(4) * t276 + t226 * t249 + t331;
t328 = t242 * t288 - t279 * t287 + t332;
t192 = pkin(3) * t338 + qJ(4) * t276 + t277 * t361;
t327 = qJD(4) * t274 + t192 * t270 + t328;
t191 = -pkin(3) * t337 + qJ(4) * t274 + t275 * t361;
t326 = -qJD(4) * t351 + t191 * t250 + t330;
t203 = pkin(4) * t246 + pkin(10) * t245;
t222 = pkin(4) * t264 + pkin(10) * t263;
t325 = t249 * t222 + (-t191 - t203) * t270 + t329;
t204 = pkin(4) * t248 + pkin(10) * t247;
t324 = t270 * t204 + (-t222 - t226) * t250 + t327;
t323 = t250 * t203 + (-t192 - t204) * t249 + t326;
t310 = qJ(5) + qJ(6);
t308 = Icges(2,4) * t321;
t307 = cos(t310);
t306 = sin(t310);
t296 = rSges(2,1) * t321 - rSges(2,2) * t317;
t295 = rSges(2,1) * t317 + rSges(2,2) * t321;
t294 = Icges(2,1) * t321 - t357;
t293 = Icges(2,1) * t317 + t308;
t292 = -Icges(2,2) * t317 + t308;
t291 = Icges(2,2) * t321 + t357;
t284 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t283 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t282 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t265 = rSges(3,3) * t312 + (rSges(3,1) * t316 + rSges(3,2) * t320) * t311;
t262 = Icges(3,5) * t312 + (Icges(3,1) * t316 + Icges(3,4) * t320) * t311;
t261 = Icges(3,6) * t312 + (Icges(3,4) * t316 + Icges(3,2) * t320) * t311;
t260 = Icges(3,3) * t312 + (Icges(3,5) * t316 + Icges(3,6) * t320) * t311;
t258 = V_base(5) * rSges(2,3) - t295 * t305 + t342;
t257 = t296 * t305 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t255 = t295 * V_base(4) - t296 * V_base(5) + V_base(3);
t244 = t264 * t318 - t339;
t243 = -t264 * t314 - t318 * t351;
t238 = t264 * t307 - t306 * t351;
t237 = -t264 * t306 - t307 * t351;
t235 = rSges(3,1) * t277 - rSges(3,2) * t276 + rSges(3,3) * t353;
t234 = rSges(3,1) * t275 - rSges(3,2) * t274 - rSges(3,3) * t350;
t233 = Icges(3,1) * t277 - Icges(3,4) * t276 + Icges(3,5) * t353;
t232 = Icges(3,1) * t275 - Icges(3,4) * t274 - Icges(3,5) * t350;
t231 = Icges(3,4) * t277 - Icges(3,2) * t276 + Icges(3,6) * t353;
t230 = Icges(3,4) * t275 - Icges(3,2) * t274 - Icges(3,6) * t350;
t229 = Icges(3,5) * t277 - Icges(3,6) * t276 + Icges(3,3) * t353;
t228 = Icges(3,5) * t275 - Icges(3,6) * t274 - Icges(3,3) * t350;
t227 = rSges(4,1) * t273 + rSges(4,2) * t272 - rSges(4,3) * t351;
t225 = Icges(4,1) * t273 + Icges(4,4) * t272 - Icges(4,5) * t351;
t224 = Icges(4,4) * t273 + Icges(4,2) * t272 - Icges(4,6) * t351;
t219 = rSges(5,1) * t264 - rSges(5,2) * t263 - rSges(5,3) * t351;
t218 = Icges(5,1) * t264 - Icges(5,4) * t263 - Icges(5,5) * t351;
t217 = Icges(5,4) * t264 - Icges(5,2) * t263 - Icges(5,6) * t351;
t215 = t248 * t318 + t355;
t214 = -t248 * t314 + t276 * t318;
t213 = t246 * t318 + t356;
t212 = -t246 * t314 + t274 * t318;
t211 = t248 * t307 + t276 * t306;
t210 = -t248 * t306 + t276 * t307;
t209 = t246 * t307 + t274 * t306;
t208 = -t246 * t306 + t274 * t307;
t207 = qJD(6) * t263 + t236;
t200 = rSges(4,1) * t254 + rSges(4,2) * t253 + rSges(4,3) * t276;
t199 = rSges(4,1) * t252 + rSges(4,2) * t251 + rSges(4,3) * t274;
t198 = Icges(4,1) * t254 + Icges(4,4) * t253 + Icges(4,5) * t276;
t197 = Icges(4,1) * t252 + Icges(4,4) * t251 + Icges(4,5) * t274;
t196 = Icges(4,4) * t254 + Icges(4,2) * t253 + Icges(4,6) * t276;
t195 = Icges(4,4) * t252 + Icges(4,2) * t251 + Icges(4,6) * t274;
t190 = rSges(5,1) * t248 - rSges(5,2) * t247 + rSges(5,3) * t276;
t189 = rSges(5,1) * t246 - rSges(5,2) * t245 + rSges(5,3) * t274;
t187 = Icges(5,1) * t248 - Icges(5,4) * t247 + Icges(5,5) * t276;
t186 = Icges(5,1) * t246 - Icges(5,4) * t245 + Icges(5,5) * t274;
t185 = Icges(5,4) * t248 - Icges(5,2) * t247 + Icges(5,6) * t276;
t184 = Icges(5,4) * t246 - Icges(5,2) * t245 + Icges(5,6) * t274;
t181 = rSges(6,1) * t244 + rSges(6,2) * t243 + rSges(6,3) * t263;
t180 = Icges(6,1) * t244 + Icges(6,4) * t243 + Icges(6,5) * t263;
t179 = Icges(6,4) * t244 + Icges(6,2) * t243 + Icges(6,6) * t263;
t178 = Icges(6,5) * t244 + Icges(6,6) * t243 + Icges(6,3) * t263;
t177 = qJD(6) * t247 + t206;
t176 = qJD(6) * t245 + t205;
t174 = rSges(7,1) * t238 + rSges(7,2) * t237 + rSges(7,3) * t263;
t173 = Icges(7,1) * t238 + Icges(7,4) * t237 + Icges(7,5) * t263;
t172 = Icges(7,4) * t238 + Icges(7,2) * t237 + Icges(7,6) * t263;
t171 = Icges(7,5) * t238 + Icges(7,6) * t237 + Icges(7,3) * t263;
t169 = -pkin(5) * t339 + pkin(11) * t263 + t264 * t360;
t167 = -t234 * t288 + t265 * t286 + t334;
t166 = t235 * t288 - t265 * t287 + t332;
t165 = t234 * t287 - t235 * t286 + t333;
t164 = rSges(6,1) * t215 + rSges(6,2) * t214 + rSges(6,3) * t247;
t163 = rSges(6,1) * t213 + rSges(6,2) * t212 + rSges(6,3) * t245;
t162 = Icges(6,1) * t215 + Icges(6,4) * t214 + Icges(6,5) * t247;
t161 = Icges(6,1) * t213 + Icges(6,4) * t212 + Icges(6,5) * t245;
t160 = Icges(6,4) * t215 + Icges(6,2) * t214 + Icges(6,6) * t247;
t159 = Icges(6,4) * t213 + Icges(6,2) * t212 + Icges(6,6) * t245;
t158 = Icges(6,5) * t215 + Icges(6,6) * t214 + Icges(6,3) * t247;
t157 = Icges(6,5) * t213 + Icges(6,6) * t212 + Icges(6,3) * t245;
t156 = rSges(7,1) * t211 + rSges(7,2) * t210 + rSges(7,3) * t247;
t155 = rSges(7,1) * t209 + rSges(7,2) * t208 + rSges(7,3) * t245;
t154 = Icges(7,1) * t211 + Icges(7,4) * t210 + Icges(7,5) * t247;
t153 = Icges(7,1) * t209 + Icges(7,4) * t208 + Icges(7,5) * t245;
t152 = Icges(7,4) * t211 + Icges(7,2) * t210 + Icges(7,6) * t247;
t151 = Icges(7,4) * t209 + Icges(7,2) * t208 + Icges(7,6) * t245;
t150 = Icges(7,5) * t211 + Icges(7,6) * t210 + Icges(7,3) * t247;
t149 = Icges(7,5) * t209 + Icges(7,6) * t208 + Icges(7,3) * t245;
t148 = pkin(5) * t355 + pkin(11) * t247 + t248 * t360;
t147 = pkin(5) * t356 + pkin(11) * t245 + t246 * t360;
t146 = -t199 * t270 + t227 * t249 + t331;
t145 = t200 * t270 - t227 * t250 + t328;
t144 = t199 * t250 - t200 * t249 + t330;
t143 = t219 * t249 + (-t189 - t191) * t270 + t329;
t142 = t190 * t270 + (-t219 - t226) * t250 + t327;
t141 = t189 * t250 + (-t190 - t192) * t249 + t326;
t140 = -t163 * t236 + t181 * t205 + t325;
t139 = t164 * t236 - t181 * t206 + t324;
t138 = t163 * t206 - t164 * t205 + t323;
t137 = -t147 * t236 - t155 * t207 + t169 * t205 + t174 * t176 + t325;
t136 = t148 * t236 + t156 * t207 - t169 * t206 - t174 * t177 + t324;
t135 = t147 * t206 - t148 * t205 + t155 * t177 - t156 * t176 + t323;
t1 = m(3) * (t165 ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + m(4) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(5) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(6) * (t138 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + m(7) * (t135 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + t288 * ((t228 * t286 + t229 * t287 + t260 * t288) * t312 + ((t231 * t320 + t233 * t316) * t287 + (t230 * t320 + t232 * t316) * t286 + (t261 * t320 + t262 * t316) * t288) * t311) / 0.2e1 + t287 * ((t229 * t353 - t231 * t276 + t233 * t277) * t287 + (t228 * t353 - t230 * t276 + t232 * t277) * t286 + (t260 * t353 - t261 * t276 + t262 * t277) * t288) / 0.2e1 + t286 * ((-t229 * t350 - t231 * t274 + t233 * t275) * t287 + (-t228 * t350 - t230 * t274 + t232 * t275) * t286 + (-t260 * t350 - t261 * t274 + t262 * t275) * t288) / 0.2e1 + t176 * ((t150 * t245 + t152 * t208 + t154 * t209) * t177 + (t245 * t149 + t208 * t151 + t209 * t153) * t176 + (t171 * t245 + t172 * t208 + t173 * t209) * t207) / 0.2e1 + t205 * ((t158 * t245 + t160 * t212 + t162 * t213) * t206 + (t245 * t157 + t212 * t159 + t213 * t161) * t205 + (t178 * t245 + t179 * t212 + t180 * t213) * t236) / 0.2e1 + t177 * ((t247 * t150 + t210 * t152 + t211 * t154) * t177 + (t149 * t247 + t151 * t210 + t153 * t211) * t176 + (t171 * t247 + t172 * t210 + t173 * t211) * t207) / 0.2e1 + t206 * ((t247 * t158 + t214 * t160 + t215 * t162) * t206 + (t157 * t247 + t159 * t214 + t161 * t215) * t205 + (t178 * t247 + t179 * t214 + t180 * t215) * t236) / 0.2e1 + m(2) * (t255 ^ 2 + t257 ^ 2 + t258 ^ 2) / 0.2e1 + t207 * ((t150 * t263 + t152 * t237 + t154 * t238) * t177 + (t149 * t263 + t151 * t237 + t153 * t238) * t176 + (t263 * t171 + t237 * t172 + t238 * t173) * t207) / 0.2e1 + t236 * ((t158 * t263 + t160 * t243 + t162 * t244) * t206 + (t157 * t263 + t159 * t243 + t161 * t244) * t205 + (t178 * t263 + t179 * t243 + t180 * t244) * t236) / 0.2e1 + m(1) * (t282 ^ 2 + t283 ^ 2 + t284 ^ 2) / 0.2e1 + ((-t217 * t245 + t218 * t246 + t224 * t251 + t225 * t252 + t274 * t366) * t270 + (-t185 * t245 + t187 * t246 + t196 * t251 + t198 * t252 + t274 * t367) * t250 + (-t184 * t245 + t186 * t246 + t195 * t251 + t197 * t252 + t274 * t368) * t249) * t249 / 0.2e1 + ((-t217 * t247 + t218 * t248 + t224 * t253 + t225 * t254 + t276 * t366) * t270 + (-t185 * t247 + t187 * t248 + t196 * t253 + t198 * t254 + t276 * t367) * t250 + (-t184 * t247 + t186 * t248 + t195 * t253 + t197 * t254 + t276 * t368) * t249) * t250 / 0.2e1 + ((-t217 * t263 + t218 * t264 + t224 * t272 + t225 * t273 - t351 * t366) * t270 + (-t185 * t263 + t187 * t264 + t196 * t272 + t198 * t273 - t351 * t367) * t250 + (-t184 * t263 + t186 * t264 + t195 * t272 + t197 * t273 - t351 * t368) * t249) * t270 / 0.2e1 + ((-t291 * t317 + t293 * t321 + Icges(1,4)) * V_base(5) + (-t292 * t317 + t294 * t321 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t291 * t321 + t293 * t317 + Icges(1,2)) * V_base(5) + (t292 * t321 + t294 * t317 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t317 + Icges(2,6) * t321) * V_base(5) + (Icges(2,5) * t321 - Icges(2,6) * t317) * V_base(4) + Icges(2,3) * t305 / 0.2e1) * t305;
T  = t1;
