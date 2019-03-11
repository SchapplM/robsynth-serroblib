% Calculate kinetic energy for
% S6PRRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRPRR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRR2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:55:48
% EndTime: 2019-03-08 21:55:53
% DurationCPUTime: 5.03s
% Computational Cost: add. (3690->437), mult. (6058->630), div. (0->0), fcn. (7208->14), ass. (0->191)
t371 = Icges(4,3) + Icges(5,3);
t310 = sin(pkin(11));
t312 = cos(pkin(11));
t320 = cos(qJ(2));
t313 = cos(pkin(6));
t317 = sin(qJ(2));
t347 = t313 * t317;
t271 = t310 * t320 + t312 * t347;
t344 = qJ(3) + pkin(12);
t304 = sin(t344);
t311 = sin(pkin(6));
t335 = cos(t344);
t334 = t311 * t335;
t243 = t271 * t304 + t312 * t334;
t353 = t311 * t312;
t244 = t271 * t335 - t304 * t353;
t316 = sin(qJ(3));
t319 = cos(qJ(3));
t350 = t311 * t319;
t252 = -t271 * t316 - t312 * t350;
t352 = t311 * t316;
t337 = t312 * t352;
t253 = t271 * t319 - t337;
t346 = t313 * t320;
t270 = t310 * t317 - t312 * t346;
t370 = Icges(4,5) * t253 + Icges(5,5) * t244 + Icges(4,6) * t252 - Icges(5,6) * t243 + t371 * t270;
t273 = -t310 * t347 + t312 * t320;
t245 = t273 * t304 - t310 * t334;
t354 = t310 * t311;
t246 = t273 * t335 + t304 * t354;
t254 = -t273 * t316 + t310 * t350;
t338 = t310 * t352;
t255 = t273 * t319 + t338;
t272 = t310 * t346 + t312 * t317;
t369 = Icges(4,5) * t255 + Icges(5,5) * t246 + Icges(4,6) * t254 - Icges(5,6) * t245 + t371 * t272;
t351 = t311 * t317;
t263 = t304 * t351 - t313 * t335;
t264 = t313 * t304 + t317 * t334;
t277 = t313 * t319 - t316 * t351;
t348 = t313 * t316;
t278 = t317 * t350 + t348;
t349 = t311 * t320;
t368 = Icges(4,5) * t278 + Icges(5,5) * t264 + Icges(4,6) * t277 - Icges(5,6) * t263 - t371 * t349;
t362 = pkin(7) * t313;
t361 = pkin(3) * t319;
t318 = cos(qJ(5));
t360 = pkin(5) * t318;
t357 = Icges(2,4) * t310;
t315 = sin(qJ(5));
t356 = t270 * t315;
t355 = t272 * t315;
t345 = qJD(2) * t311;
t343 = V_base(5) * qJ(1) + V_base(1);
t339 = qJD(1) + V_base(3);
t336 = t315 * t349;
t287 = t310 * t345 + V_base(4);
t298 = qJD(2) * t313 + V_base(6);
t251 = qJD(3) * t272 + t287;
t206 = qJD(5) * t245 + t251;
t286 = -t312 * t345 + V_base(5);
t250 = qJD(3) * t270 + t286;
t274 = -qJD(3) * t349 + t298;
t280 = pkin(1) * t310 - pkin(7) * t353;
t333 = -t280 * V_base(6) + V_base(5) * t362 + t343;
t281 = pkin(1) * t312 + pkin(7) * t354;
t332 = V_base(4) * t280 - t281 * V_base(5) + t339;
t205 = qJD(5) * t243 + t250;
t236 = qJD(5) * t263 + t274;
t331 = V_base(6) * t281 + V_base(2) + (-qJ(1) - t362) * V_base(4);
t239 = pkin(2) * t271 + pkin(8) * t270;
t279 = (pkin(2) * t317 - pkin(8) * t320) * t311;
t330 = -t239 * t298 + t286 * t279 + t333;
t240 = pkin(2) * t273 + pkin(8) * t272;
t329 = t287 * t239 - t240 * t286 + t332;
t328 = t298 * t240 - t279 * t287 + t331;
t234 = pkin(3) * t348 + (-qJ(4) * t320 + t317 * t361) * t311;
t327 = qJD(4) * t272 + t250 * t234 + t330;
t191 = pkin(3) * t338 + qJ(4) * t272 + t273 * t361;
t326 = qJD(4) * t270 + t274 * t191 + t328;
t190 = -pkin(3) * t337 + qJ(4) * t270 + t271 * t361;
t325 = -qJD(4) * t349 + t251 * t190 + t329;
t203 = t244 * pkin(4) + t243 * pkin(9);
t228 = t264 * pkin(4) + t263 * pkin(9);
t324 = t250 * t228 + (-t190 - t203) * t274 + t327;
t204 = t246 * pkin(4) + t245 * pkin(9);
t323 = t274 * t204 + (-t228 - t234) * t251 + t326;
t322 = t251 * t203 + (-t191 - t204) * t250 + t325;
t309 = qJ(5) + qJ(6);
t307 = cos(t309);
t306 = sin(t309);
t305 = Icges(2,4) * t312;
t295 = rSges(2,1) * t312 - rSges(2,2) * t310;
t294 = rSges(2,1) * t310 + rSges(2,2) * t312;
t293 = Icges(2,1) * t312 - t357;
t292 = Icges(2,1) * t310 + t305;
t291 = -Icges(2,2) * t310 + t305;
t290 = Icges(2,2) * t312 + t357;
t285 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t284 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t283 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t265 = rSges(3,3) * t313 + (rSges(3,1) * t317 + rSges(3,2) * t320) * t311;
t262 = Icges(3,5) * t313 + (Icges(3,1) * t317 + Icges(3,4) * t320) * t311;
t261 = Icges(3,6) * t313 + (Icges(3,4) * t317 + Icges(3,2) * t320) * t311;
t260 = Icges(3,3) * t313 + (Icges(3,5) * t317 + Icges(3,6) * t320) * t311;
t258 = V_base(5) * rSges(2,3) - t294 * V_base(6) + t343;
t257 = t295 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t249 = t294 * V_base(4) - t295 * V_base(5) + t339;
t248 = t264 * t318 - t336;
t247 = -t264 * t315 - t318 * t349;
t238 = t264 * t307 - t306 * t349;
t237 = -t264 * t306 - t307 * t349;
t235 = rSges(4,1) * t278 + rSges(4,2) * t277 - rSges(4,3) * t349;
t233 = Icges(4,1) * t278 + Icges(4,4) * t277 - Icges(4,5) * t349;
t232 = Icges(4,4) * t278 + Icges(4,2) * t277 - Icges(4,6) * t349;
t230 = rSges(3,1) * t273 - rSges(3,2) * t272 + rSges(3,3) * t354;
t229 = rSges(3,1) * t271 - rSges(3,2) * t270 - rSges(3,3) * t353;
t227 = Icges(3,1) * t273 - Icges(3,4) * t272 + Icges(3,5) * t354;
t226 = Icges(3,1) * t271 - Icges(3,4) * t270 - Icges(3,5) * t353;
t225 = Icges(3,4) * t273 - Icges(3,2) * t272 + Icges(3,6) * t354;
t224 = Icges(3,4) * t271 - Icges(3,2) * t270 - Icges(3,6) * t353;
t223 = Icges(3,5) * t273 - Icges(3,6) * t272 + Icges(3,3) * t354;
t222 = Icges(3,5) * t271 - Icges(3,6) * t270 - Icges(3,3) * t353;
t219 = rSges(5,1) * t264 - rSges(5,2) * t263 - rSges(5,3) * t349;
t218 = Icges(5,1) * t264 - Icges(5,4) * t263 - Icges(5,5) * t349;
t217 = Icges(5,4) * t264 - Icges(5,2) * t263 - Icges(5,6) * t349;
t215 = t246 * t318 + t355;
t214 = -t246 * t315 + t272 * t318;
t213 = t244 * t318 + t356;
t212 = -t244 * t315 + t270 * t318;
t211 = qJD(6) * t263 + t236;
t210 = t246 * t307 + t272 * t306;
t209 = -t246 * t306 + t272 * t307;
t208 = t244 * t307 + t270 * t306;
t207 = -t244 * t306 + t270 * t307;
t200 = rSges(4,1) * t255 + rSges(4,2) * t254 + rSges(4,3) * t272;
t199 = rSges(4,1) * t253 + rSges(4,2) * t252 + rSges(4,3) * t270;
t198 = Icges(4,1) * t255 + Icges(4,4) * t254 + Icges(4,5) * t272;
t197 = Icges(4,1) * t253 + Icges(4,4) * t252 + Icges(4,5) * t270;
t196 = Icges(4,4) * t255 + Icges(4,2) * t254 + Icges(4,6) * t272;
t195 = Icges(4,4) * t253 + Icges(4,2) * t252 + Icges(4,6) * t270;
t189 = rSges(5,1) * t246 - rSges(5,2) * t245 + rSges(5,3) * t272;
t188 = rSges(5,1) * t244 - rSges(5,2) * t243 + rSges(5,3) * t270;
t187 = Icges(5,1) * t246 - Icges(5,4) * t245 + Icges(5,5) * t272;
t186 = Icges(5,1) * t244 - Icges(5,4) * t243 + Icges(5,5) * t270;
t185 = Icges(5,4) * t246 - Icges(5,2) * t245 + Icges(5,6) * t272;
t184 = Icges(5,4) * t244 - Icges(5,2) * t243 + Icges(5,6) * t270;
t181 = rSges(6,1) * t248 + rSges(6,2) * t247 + rSges(6,3) * t263;
t180 = Icges(6,1) * t248 + Icges(6,4) * t247 + Icges(6,5) * t263;
t179 = Icges(6,4) * t248 + Icges(6,2) * t247 + Icges(6,6) * t263;
t178 = Icges(6,5) * t248 + Icges(6,6) * t247 + Icges(6,3) * t263;
t177 = rSges(7,1) * t238 + rSges(7,2) * t237 + rSges(7,3) * t263;
t176 = Icges(7,1) * t238 + Icges(7,4) * t237 + Icges(7,5) * t263;
t175 = Icges(7,4) * t238 + Icges(7,2) * t237 + Icges(7,6) * t263;
t174 = Icges(7,5) * t238 + Icges(7,6) * t237 + Icges(7,3) * t263;
t173 = qJD(6) * t245 + t206;
t172 = qJD(6) * t243 + t205;
t170 = -pkin(5) * t336 + pkin(10) * t263 + t264 * t360;
t168 = -t229 * t298 + t265 * t286 + t333;
t167 = t230 * t298 - t265 * t287 + t331;
t165 = t229 * t287 - t230 * t286 + t332;
t164 = rSges(6,1) * t215 + rSges(6,2) * t214 + rSges(6,3) * t245;
t163 = rSges(6,1) * t213 + rSges(6,2) * t212 + rSges(6,3) * t243;
t162 = Icges(6,1) * t215 + Icges(6,4) * t214 + Icges(6,5) * t245;
t161 = Icges(6,1) * t213 + Icges(6,4) * t212 + Icges(6,5) * t243;
t160 = Icges(6,4) * t215 + Icges(6,2) * t214 + Icges(6,6) * t245;
t159 = Icges(6,4) * t213 + Icges(6,2) * t212 + Icges(6,6) * t243;
t158 = Icges(6,5) * t215 + Icges(6,6) * t214 + Icges(6,3) * t245;
t157 = Icges(6,5) * t213 + Icges(6,6) * t212 + Icges(6,3) * t243;
t156 = rSges(7,1) * t210 + rSges(7,2) * t209 + rSges(7,3) * t245;
t155 = rSges(7,1) * t208 + rSges(7,2) * t207 + rSges(7,3) * t243;
t154 = Icges(7,1) * t210 + Icges(7,4) * t209 + Icges(7,5) * t245;
t153 = Icges(7,1) * t208 + Icges(7,4) * t207 + Icges(7,5) * t243;
t152 = Icges(7,4) * t210 + Icges(7,2) * t209 + Icges(7,6) * t245;
t151 = Icges(7,4) * t208 + Icges(7,2) * t207 + Icges(7,6) * t243;
t150 = Icges(7,5) * t210 + Icges(7,6) * t209 + Icges(7,3) * t245;
t149 = Icges(7,5) * t208 + Icges(7,6) * t207 + Icges(7,3) * t243;
t148 = pkin(5) * t355 + pkin(10) * t245 + t246 * t360;
t147 = pkin(5) * t356 + pkin(10) * t243 + t244 * t360;
t146 = -t199 * t274 + t235 * t250 + t330;
t145 = t200 * t274 - t235 * t251 + t328;
t144 = t199 * t251 - t200 * t250 + t329;
t143 = t219 * t250 + (-t188 - t190) * t274 + t327;
t142 = t189 * t274 + (-t219 - t234) * t251 + t326;
t141 = t188 * t251 + (-t189 - t191) * t250 + t325;
t140 = -t163 * t236 + t181 * t205 + t324;
t139 = t164 * t236 - t181 * t206 + t323;
t138 = t163 * t206 - t164 * t205 + t322;
t137 = -t147 * t236 - t155 * t211 + t170 * t205 + t172 * t177 + t324;
t136 = t148 * t236 + t156 * t211 - t170 * t206 - t173 * t177 + t323;
t135 = t147 * t206 - t148 * t205 + t155 * t173 - t156 * t172 + t322;
t1 = m(4) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(5) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(6) * (t138 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + m(7) * (t135 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + t298 * ((t222 * t286 + t223 * t287 + t260 * t298) * t313 + ((t225 * t320 + t227 * t317) * t287 + (t224 * t320 + t226 * t317) * t286 + (t261 * t320 + t262 * t317) * t298) * t311) / 0.2e1 + t205 * ((t158 * t243 + t160 * t212 + t162 * t213) * t206 + (t243 * t157 + t212 * t159 + t213 * t161) * t205 + (t178 * t243 + t179 * t212 + t180 * t213) * t236) / 0.2e1 + t172 * ((t150 * t243 + t152 * t207 + t154 * t208) * t173 + (t243 * t149 + t207 * t151 + t208 * t153) * t172 + (t174 * t243 + t175 * t207 + t176 * t208) * t211) / 0.2e1 + t206 * ((t245 * t158 + t214 * t160 + t215 * t162) * t206 + (t157 * t245 + t159 * t214 + t161 * t215) * t205 + (t178 * t245 + t179 * t214 + t180 * t215) * t236) / 0.2e1 + t173 * ((t245 * t150 + t209 * t152 + t210 * t154) * t173 + (t149 * t245 + t151 * t209 + t153 * t210) * t172 + (t174 * t245 + t175 * t209 + t176 * t210) * t211) / 0.2e1 + m(2) * (t249 ^ 2 + t257 ^ 2 + t258 ^ 2) / 0.2e1 + t211 * ((t150 * t263 + t152 * t237 + t154 * t238) * t173 + (t149 * t263 + t151 * t237 + t153 * t238) * t172 + (t263 * t174 + t237 * t175 + t238 * t176) * t211) / 0.2e1 + t236 * ((t158 * t263 + t160 * t247 + t162 * t248) * t206 + (t157 * t263 + t159 * t247 + t161 * t248) * t205 + (t178 * t263 + t179 * t247 + t180 * t248) * t236) / 0.2e1 + m(1) * (t283 ^ 2 + t284 ^ 2 + t285 ^ 2) / 0.2e1 + t286 * ((-t223 * t353 - t225 * t270 + t227 * t271) * t287 + (-t222 * t353 - t224 * t270 + t226 * t271) * t286 + (-t260 * t353 - t261 * t270 + t262 * t271) * t298) / 0.2e1 + t287 * ((t223 * t354 - t225 * t272 + t227 * t273) * t287 + (t222 * t354 - t224 * t272 + t226 * t273) * t286 + (t260 * t354 - t261 * t272 + t262 * t273) * t298) / 0.2e1 + m(3) * (t165 ^ 2 + t167 ^ 2 + t168 ^ 2) / 0.2e1 + ((-t217 * t243 + t218 * t244 + t232 * t252 + t233 * t253 + t270 * t368) * t274 + (-t185 * t243 + t187 * t244 + t196 * t252 + t198 * t253 + t270 * t369) * t251 + (-t184 * t243 + t186 * t244 + t195 * t252 + t197 * t253 + t370 * t270) * t250) * t250 / 0.2e1 + ((-t217 * t245 + t218 * t246 + t232 * t254 + t233 * t255 + t272 * t368) * t274 + (-t185 * t245 + t187 * t246 + t196 * t254 + t198 * t255 + t369 * t272) * t251 + (-t184 * t245 + t186 * t246 + t195 * t254 + t197 * t255 + t272 * t370) * t250) * t251 / 0.2e1 + ((-t217 * t263 + t218 * t264 + t232 * t277 + t233 * t278 - t349 * t368) * t274 + (-t185 * t263 + t187 * t264 + t196 * t277 + t198 * t278 - t349 * t369) * t251 + (-t184 * t263 + t186 * t264 + t195 * t277 + t197 * t278 - t349 * t370) * t250) * t274 / 0.2e1 + ((-t290 * t310 + t292 * t312 + Icges(1,4)) * V_base(5) + (-t291 * t310 + t293 * t312 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t290 * t312 + t292 * t310 + Icges(1,2)) * V_base(5) + (t291 * t312 + t293 * t310 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t312 - Icges(2,6) * t310 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t310 + Icges(2,6) * t312 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
