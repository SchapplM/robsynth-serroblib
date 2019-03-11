% Calculate kinetic energy for
% S6RRRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR9_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPR9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR9_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR9_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR9_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:47:41
% EndTime: 2019-03-09 22:47:46
% DurationCPUTime: 4.75s
% Computational Cost: add. (3687->438), mult. (6052->631), div. (0->0), fcn. (7184->14), ass. (0->192)
t369 = Icges(5,2) + Icges(6,3);
t314 = cos(pkin(6));
t318 = sin(qJ(1));
t320 = cos(qJ(2));
t347 = t318 * t320;
t317 = sin(qJ(2));
t321 = cos(qJ(1));
t348 = t317 * t321;
t275 = t314 * t348 + t347;
t345 = qJ(3) + qJ(4);
t307 = sin(t345);
t336 = cos(t345);
t312 = sin(pkin(6));
t351 = t312 * t321;
t245 = t275 * t336 - t307 * t351;
t346 = t320 * t321;
t349 = t317 * t318;
t274 = -t314 * t346 + t349;
t311 = sin(pkin(12));
t313 = cos(pkin(12));
t209 = -t245 * t311 + t274 * t313;
t357 = t274 * t311;
t210 = t245 * t313 + t357;
t335 = t312 * t336;
t244 = t275 * t307 + t321 * t335;
t368 = -Icges(5,4) * t245 + Icges(6,5) * t210 - Icges(5,6) * t274 + Icges(6,6) * t209 + t369 * t244;
t277 = -t314 * t349 + t346;
t354 = t312 * t318;
t247 = t277 * t336 + t307 * t354;
t276 = t314 * t347 + t348;
t211 = -t247 * t311 + t276 * t313;
t356 = t276 * t311;
t212 = t247 * t313 + t356;
t246 = t277 * t307 - t318 * t335;
t367 = -Icges(5,4) * t247 + Icges(6,5) * t212 - Icges(5,6) * t276 + Icges(6,6) * t211 + t369 * t246;
t265 = t314 * t307 + t317 * t335;
t352 = t312 * t320;
t242 = -t265 * t311 - t313 * t352;
t339 = t311 * t352;
t243 = t265 * t313 - t339;
t355 = t312 * t317;
t264 = t307 * t355 - t314 * t336;
t366 = -Icges(5,4) * t265 + Icges(6,5) * t243 + Icges(5,6) * t352 + Icges(6,6) * t242 + t369 * t264;
t362 = pkin(8) * t314;
t319 = cos(qJ(3));
t361 = pkin(3) * t319;
t360 = pkin(5) * t313;
t358 = Icges(2,4) * t318;
t353 = t312 * t319;
t316 = sin(qJ(3));
t350 = t314 * t316;
t343 = qJD(2) * t312;
t342 = V_base(5) * pkin(7) + V_base(1);
t338 = t316 * t354;
t337 = t316 * t351;
t286 = t318 * t343 + V_base(4);
t306 = V_base(6) + qJD(1);
t249 = qJD(3) * t276 + t286;
t288 = qJD(2) * t314 + t306;
t220 = qJD(4) * t276 + t249;
t285 = -t321 * t343 + V_base(5);
t280 = pkin(1) * t318 - pkin(8) * t351;
t334 = -t280 * t306 + V_base(5) * t362 + t342;
t281 = pkin(1) * t321 + pkin(8) * t354;
t333 = V_base(4) * t280 - t281 * V_base(5) + V_base(3);
t248 = qJD(3) * t274 + t285;
t219 = qJD(4) * t274 + t248;
t332 = t306 * t281 + V_base(2) + (-pkin(7) - t362) * V_base(4);
t258 = (-qJD(3) - qJD(4)) * t352 + t288;
t238 = t275 * pkin(2) + t274 * pkin(9);
t279 = (pkin(2) * t317 - pkin(9) * t320) * t312;
t331 = -t238 * t288 + t285 * t279 + t334;
t239 = t277 * pkin(2) + t276 * pkin(9);
t330 = t286 * t238 - t239 * t285 + t333;
t329 = t288 * t239 - t279 * t286 + t332;
t190 = -pkin(3) * t337 + pkin(10) * t274 + t275 * t361;
t233 = pkin(3) * t350 + (-pkin(10) * t320 + t317 * t361) * t312;
t270 = -qJD(3) * t352 + t288;
t328 = -t190 * t270 + t248 * t233 + t331;
t191 = pkin(3) * t338 + pkin(10) * t276 + t277 * t361;
t327 = t249 * t190 - t191 * t248 + t330;
t225 = pkin(4) * t265 + qJ(5) * t264;
t326 = qJD(5) * t246 + t219 * t225 + t328;
t203 = pkin(4) * t245 + qJ(5) * t244;
t325 = qJD(5) * t264 + t220 * t203 + t327;
t324 = t270 * t191 - t233 * t249 + t329;
t204 = pkin(4) * t247 + qJ(5) * t246;
t323 = qJD(5) * t244 + t258 * t204 + t324;
t310 = pkin(12) + qJ(6);
t308 = Icges(2,4) * t321;
t305 = cos(t310);
t304 = sin(t310);
t296 = rSges(2,1) * t321 - rSges(2,2) * t318;
t295 = rSges(2,1) * t318 + rSges(2,2) * t321;
t294 = Icges(2,1) * t321 - t358;
t293 = Icges(2,1) * t318 + t308;
t292 = -Icges(2,2) * t318 + t308;
t291 = Icges(2,2) * t321 + t358;
t284 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t283 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t282 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t273 = t317 * t353 + t350;
t272 = t314 * t319 - t316 * t355;
t263 = rSges(3,3) * t314 + (rSges(3,1) * t317 + rSges(3,2) * t320) * t312;
t262 = Icges(3,5) * t314 + (Icges(3,1) * t317 + Icges(3,4) * t320) * t312;
t261 = Icges(3,6) * t314 + (Icges(3,4) * t317 + Icges(3,2) * t320) * t312;
t260 = Icges(3,3) * t314 + (Icges(3,5) * t317 + Icges(3,6) * t320) * t312;
t257 = V_base(5) * rSges(2,3) - t295 * t306 + t342;
t256 = t296 * t306 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t254 = t295 * V_base(4) - t296 * V_base(5) + V_base(3);
t253 = t277 * t319 + t338;
t252 = -t277 * t316 + t318 * t353;
t251 = t275 * t319 - t337;
t250 = -t275 * t316 - t319 * t351;
t237 = t265 * t305 - t304 * t352;
t236 = -t265 * t304 - t305 * t352;
t235 = rSges(3,1) * t277 - rSges(3,2) * t276 + rSges(3,3) * t354;
t234 = rSges(3,1) * t275 - rSges(3,2) * t274 - rSges(3,3) * t351;
t232 = Icges(3,1) * t277 - Icges(3,4) * t276 + Icges(3,5) * t354;
t231 = Icges(3,1) * t275 - Icges(3,4) * t274 - Icges(3,5) * t351;
t230 = Icges(3,4) * t277 - Icges(3,2) * t276 + Icges(3,6) * t354;
t229 = Icges(3,4) * t275 - Icges(3,2) * t274 - Icges(3,6) * t351;
t228 = Icges(3,5) * t277 - Icges(3,6) * t276 + Icges(3,3) * t354;
t227 = Icges(3,5) * t275 - Icges(3,6) * t274 - Icges(3,3) * t351;
t226 = rSges(4,1) * t273 + rSges(4,2) * t272 - rSges(4,3) * t352;
t224 = Icges(4,1) * t273 + Icges(4,4) * t272 - Icges(4,5) * t352;
t223 = Icges(4,4) * t273 + Icges(4,2) * t272 - Icges(4,6) * t352;
t222 = Icges(4,5) * t273 + Icges(4,6) * t272 - Icges(4,3) * t352;
t217 = qJD(6) * t264 + t258;
t216 = rSges(5,1) * t265 - rSges(5,2) * t264 - rSges(5,3) * t352;
t215 = Icges(5,1) * t265 - Icges(5,4) * t264 - Icges(5,5) * t352;
t213 = Icges(5,5) * t265 - Icges(5,6) * t264 - Icges(5,3) * t352;
t208 = t247 * t305 + t276 * t304;
t207 = -t247 * t304 + t276 * t305;
t206 = t245 * t305 + t274 * t304;
t205 = -t245 * t304 + t274 * t305;
t201 = rSges(4,1) * t253 + rSges(4,2) * t252 + rSges(4,3) * t276;
t200 = rSges(4,1) * t251 + rSges(4,2) * t250 + rSges(4,3) * t274;
t199 = Icges(4,1) * t253 + Icges(4,4) * t252 + Icges(4,5) * t276;
t198 = Icges(4,1) * t251 + Icges(4,4) * t250 + Icges(4,5) * t274;
t197 = Icges(4,4) * t253 + Icges(4,2) * t252 + Icges(4,6) * t276;
t196 = Icges(4,4) * t251 + Icges(4,2) * t250 + Icges(4,6) * t274;
t195 = Icges(4,5) * t253 + Icges(4,6) * t252 + Icges(4,3) * t276;
t194 = Icges(4,5) * t251 + Icges(4,6) * t250 + Icges(4,3) * t274;
t193 = qJD(6) * t246 + t220;
t192 = qJD(6) * t244 + t219;
t189 = rSges(5,1) * t247 - rSges(5,2) * t246 + rSges(5,3) * t276;
t188 = rSges(5,1) * t245 - rSges(5,2) * t244 + rSges(5,3) * t274;
t187 = Icges(5,1) * t247 - Icges(5,4) * t246 + Icges(5,5) * t276;
t186 = Icges(5,1) * t245 - Icges(5,4) * t244 + Icges(5,5) * t274;
t183 = Icges(5,5) * t247 - Icges(5,6) * t246 + Icges(5,3) * t276;
t182 = Icges(5,5) * t245 - Icges(5,6) * t244 + Icges(5,3) * t274;
t180 = rSges(6,1) * t243 + rSges(6,2) * t242 + rSges(6,3) * t264;
t179 = Icges(6,1) * t243 + Icges(6,4) * t242 + Icges(6,5) * t264;
t178 = Icges(6,4) * t243 + Icges(6,2) * t242 + Icges(6,6) * t264;
t175 = rSges(7,1) * t237 + rSges(7,2) * t236 + rSges(7,3) * t264;
t174 = Icges(7,1) * t237 + Icges(7,4) * t236 + Icges(7,5) * t264;
t173 = Icges(7,4) * t237 + Icges(7,2) * t236 + Icges(7,6) * t264;
t172 = Icges(7,5) * t237 + Icges(7,6) * t236 + Icges(7,3) * t264;
t170 = -pkin(5) * t339 + pkin(11) * t264 + t265 * t360;
t167 = -t234 * t288 + t263 * t285 + t334;
t166 = t235 * t288 - t263 * t286 + t332;
t165 = t234 * t286 - t235 * t285 + t333;
t164 = rSges(6,1) * t212 + rSges(6,2) * t211 + rSges(6,3) * t246;
t163 = rSges(6,1) * t210 + rSges(6,2) * t209 + rSges(6,3) * t244;
t162 = Icges(6,1) * t212 + Icges(6,4) * t211 + Icges(6,5) * t246;
t161 = Icges(6,1) * t210 + Icges(6,4) * t209 + Icges(6,5) * t244;
t160 = Icges(6,4) * t212 + Icges(6,2) * t211 + Icges(6,6) * t246;
t159 = Icges(6,4) * t210 + Icges(6,2) * t209 + Icges(6,6) * t244;
t156 = rSges(7,1) * t208 + rSges(7,2) * t207 + rSges(7,3) * t246;
t155 = rSges(7,1) * t206 + rSges(7,2) * t205 + rSges(7,3) * t244;
t154 = Icges(7,1) * t208 + Icges(7,4) * t207 + Icges(7,5) * t246;
t153 = Icges(7,1) * t206 + Icges(7,4) * t205 + Icges(7,5) * t244;
t152 = Icges(7,4) * t208 + Icges(7,2) * t207 + Icges(7,6) * t246;
t151 = Icges(7,4) * t206 + Icges(7,2) * t205 + Icges(7,6) * t244;
t150 = Icges(7,5) * t208 + Icges(7,6) * t207 + Icges(7,3) * t246;
t149 = Icges(7,5) * t206 + Icges(7,6) * t205 + Icges(7,3) * t244;
t148 = pkin(5) * t356 + pkin(11) * t246 + t247 * t360;
t147 = pkin(5) * t357 + pkin(11) * t244 + t245 * t360;
t146 = -t200 * t270 + t226 * t248 + t331;
t145 = t201 * t270 - t226 * t249 + t329;
t144 = t200 * t249 - t201 * t248 + t330;
t143 = -t188 * t258 + t216 * t219 + t328;
t142 = t189 * t258 - t216 * t220 + t324;
t141 = t188 * t220 - t189 * t219 + t327;
t140 = t180 * t219 + (-t163 - t203) * t258 + t326;
t139 = t164 * t258 + (-t180 - t225) * t220 + t323;
t138 = t163 * t220 + (-t164 - t204) * t219 + t325;
t137 = t326 + (-t147 - t203) * t258 - t155 * t217 + t170 * t219 + t175 * t192;
t136 = t148 * t258 + t156 * t217 - t175 * t193 + (-t170 - t225) * t220 + t323;
t135 = t147 * t220 + t155 * t193 - t156 * t192 + (-t148 - t204) * t219 + t325;
t1 = m(6) * (t138 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + m(7) * (t135 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + m(1) * (t282 ^ 2 + t283 ^ 2 + t284 ^ 2) / 0.2e1 + t249 * ((t195 * t276 + t252 * t197 + t253 * t199) * t249 + (t194 * t276 + t196 * t252 + t198 * t253) * t248 + (t222 * t276 + t223 * t252 + t224 * t253) * t270) / 0.2e1 + t248 * ((t195 * t274 + t197 * t250 + t199 * t251) * t249 + (t274 * t194 + t250 * t196 + t251 * t198) * t248 + (t222 * t274 + t223 * t250 + t224 * t251) * t270) / 0.2e1 + t217 * ((t150 * t264 + t152 * t236 + t154 * t237) * t193 + (t149 * t264 + t151 * t236 + t153 * t237) * t192 + (t264 * t172 + t236 * t173 + t174 * t237) * t217) / 0.2e1 + m(2) * (t254 ^ 2 + t256 ^ 2 + t257 ^ 2) / 0.2e1 + t193 * ((t150 * t246 + t207 * t152 + t154 * t208) * t193 + (t149 * t246 + t151 * t207 + t153 * t208) * t192 + (t172 * t246 + t173 * t207 + t174 * t208) * t217) / 0.2e1 + t192 * ((t150 * t244 + t152 * t205 + t154 * t206) * t193 + (t244 * t149 + t205 * t151 + t206 * t153) * t192 + (t172 * t244 + t173 * t205 + t174 * t206) * t217) / 0.2e1 + m(3) * (t165 ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + t286 * ((t228 * t354 - t276 * t230 + t277 * t232) * t286 + (t227 * t354 - t229 * t276 + t231 * t277) * t285 + (t260 * t354 - t261 * t276 + t262 * t277) * t288) / 0.2e1 + t270 * ((-t195 * t352 + t197 * t272 + t199 * t273) * t249 + (-t194 * t352 + t196 * t272 + t198 * t273) * t248 + (-t222 * t352 + t272 * t223 + t273 * t224) * t270) / 0.2e1 + t285 * ((-t228 * t351 - t230 * t274 + t232 * t275) * t286 + (-t227 * t351 - t274 * t229 + t275 * t231) * t285 + (-t260 * t351 - t261 * t274 + t262 * t275) * t288) / 0.2e1 + t288 * ((t227 * t285 + t228 * t286 + t260 * t288) * t314 + ((t230 * t320 + t232 * t317) * t286 + (t229 * t320 + t231 * t317) * t285 + (t261 * t320 + t262 * t317) * t288) * t312) / 0.2e1 + m(4) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(5) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + ((t178 * t209 + t179 * t210 + t213 * t274 + t215 * t245 + t244 * t366) * t258 + (t160 * t209 + t162 * t210 + t183 * t274 + t187 * t245 + t244 * t367) * t220 + (t209 * t159 + t210 * t161 + t274 * t182 + t245 * t186 + t368 * t244) * t219) * t219 / 0.2e1 + ((t178 * t211 + t179 * t212 + t213 * t276 + t215 * t247 + t246 * t366) * t258 + (t211 * t160 + t212 * t162 + t276 * t183 + t247 * t187 + t367 * t246) * t220 + (t159 * t211 + t161 * t212 + t182 * t276 + t186 * t247 + t246 * t368) * t219) * t220 / 0.2e1 + ((t242 * t178 + t243 * t179 - t213 * t352 + t265 * t215 + t366 * t264) * t258 + (t160 * t242 + t162 * t243 - t183 * t352 + t187 * t265 + t264 * t367) * t220 + (t159 * t242 + t161 * t243 - t182 * t352 + t186 * t265 + t264 * t368) * t219) * t258 / 0.2e1 + ((-t291 * t318 + t293 * t321 + Icges(1,4)) * V_base(5) + (-t318 * t292 + t321 * t294 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t321 * t291 + t318 * t293 + Icges(1,2)) * V_base(5) + (t292 * t321 + t294 * t318 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t318 + Icges(2,6) * t321) * V_base(5) + (Icges(2,5) * t321 - Icges(2,6) * t318) * V_base(4) + Icges(2,3) * t306 / 0.2e1) * t306;
T  = t1;
