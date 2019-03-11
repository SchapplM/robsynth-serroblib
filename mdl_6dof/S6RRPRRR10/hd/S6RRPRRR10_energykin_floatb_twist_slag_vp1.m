% Calculate kinetic energy for
% S6RRPRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR10_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR10_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRR10_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR10_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR10_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR10_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:18:14
% EndTime: 2019-03-09 14:18:18
% DurationCPUTime: 4.18s
% Computational Cost: add. (3695->442), mult. (5948->639), div. (0->0), fcn. (7076->14), ass. (0->196)
t314 = cos(pkin(6));
t318 = sin(qJ(1));
t320 = cos(qJ(2));
t347 = t318 * t320;
t317 = sin(qJ(2));
t321 = cos(qJ(1));
t348 = t317 * t321;
t275 = t314 * t348 + t347;
t311 = sin(pkin(12));
t313 = cos(pkin(12));
t312 = sin(pkin(6));
t350 = t312 * t321;
t249 = -t275 * t311 - t313 * t350;
t338 = t311 * t350;
t250 = t275 * t313 - t338;
t346 = t320 * t321;
t349 = t317 * t318;
t274 = -t314 * t346 + t349;
t193 = Icges(4,5) * t250 + Icges(4,6) * t249 + Icges(4,3) * t274;
t230 = Icges(3,4) * t275 - Icges(3,2) * t274 - Icges(3,6) * t350;
t367 = t193 - t230;
t277 = -t314 * t349 + t346;
t352 = t312 * t318;
t251 = -t277 * t311 + t313 * t352;
t339 = t311 * t352;
t252 = t277 * t313 + t339;
t276 = t314 * t347 + t348;
t194 = Icges(4,5) * t252 + Icges(4,6) * t251 + Icges(4,3) * t276;
t231 = Icges(3,4) * t277 - Icges(3,2) * t276 + Icges(3,6) * t352;
t366 = t194 - t231;
t353 = t312 * t317;
t272 = -t311 * t353 + t313 * t314;
t354 = t311 * t314;
t273 = t313 * t353 + t354;
t351 = t312 * t320;
t222 = Icges(4,5) * t273 + Icges(4,6) * t272 - Icges(4,3) * t351;
t261 = Icges(3,6) * t314 + (Icges(3,4) * t317 + Icges(3,2) * t320) * t312;
t365 = t222 - t261;
t361 = pkin(8) * t314;
t360 = pkin(3) * t313;
t319 = cos(qJ(5));
t359 = pkin(5) * t319;
t357 = Icges(2,4) * t318;
t316 = sin(qJ(5));
t356 = t274 * t316;
t355 = t276 * t316;
t344 = qJD(2) * t312;
t343 = pkin(12) + qJ(4);
t342 = V_base(5) * pkin(7) + V_base(1);
t337 = t316 * t351;
t287 = t318 * t344 + V_base(4);
t305 = V_base(6) + qJD(1);
t336 = cos(t343);
t254 = qJD(4) * t276 + t287;
t288 = qJD(2) * t314 + t305;
t335 = t312 * t336;
t304 = sin(t343);
t247 = t277 * t304 - t318 * t335;
t205 = qJD(5) * t247 + t254;
t286 = -t321 * t344 + V_base(5);
t280 = pkin(1) * t318 - pkin(8) * t350;
t334 = -t280 * t305 + V_base(5) * t361 + t342;
t281 = pkin(1) * t321 + pkin(8) * t352;
t333 = V_base(4) * t280 - t281 * V_base(5) + V_base(3);
t253 = qJD(4) * t274 + t286;
t245 = t275 * t304 + t321 * t335;
t204 = qJD(5) * t245 + t253;
t270 = -qJD(4) * t351 + t288;
t278 = (pkin(2) * t317 - qJ(3) * t320) * t312;
t332 = qJD(3) * t276 + t286 * t278 + t334;
t331 = t305 * t281 + V_base(2) + (-pkin(7) - t361) * V_base(4);
t263 = t304 * t353 - t314 * t336;
t236 = qJD(5) * t263 + t270;
t240 = pkin(2) * t277 + qJ(3) * t276;
t330 = qJD(3) * t274 + t288 * t240 + t331;
t239 = pkin(2) * t275 + qJ(3) * t274;
t329 = -qJD(3) * t351 + t287 * t239 + t333;
t191 = -pkin(3) * t338 + pkin(9) * t274 + t275 * t360;
t227 = pkin(3) * t354 + (-pkin(9) * t320 + t317 * t360) * t312;
t328 = t286 * t227 + (-t191 - t239) * t288 + t332;
t192 = pkin(3) * t339 + pkin(9) * t276 + t277 * t360;
t327 = t288 * t192 + (-t227 - t278) * t287 + t330;
t326 = t287 * t191 + (-t192 - t240) * t286 + t329;
t246 = t275 * t336 - t304 * t350;
t202 = t246 * pkin(4) + t245 * pkin(10);
t264 = t314 * t304 + t317 * t335;
t226 = t264 * pkin(4) + t263 * pkin(10);
t325 = -t202 * t270 + t253 * t226 + t328;
t248 = t277 * t336 + t304 * t352;
t203 = t248 * pkin(4) + t247 * pkin(10);
t324 = t270 * t203 - t226 * t254 + t327;
t323 = t254 * t202 - t203 * t253 + t326;
t310 = qJ(5) + qJ(6);
t308 = Icges(2,4) * t321;
t307 = cos(t310);
t306 = sin(t310);
t296 = rSges(2,1) * t321 - rSges(2,2) * t318;
t295 = rSges(2,1) * t318 + rSges(2,2) * t321;
t294 = Icges(2,1) * t321 - t357;
t293 = Icges(2,1) * t318 + t308;
t292 = -Icges(2,2) * t318 + t308;
t291 = Icges(2,2) * t321 + t357;
t284 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t283 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t282 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t265 = rSges(3,3) * t314 + (rSges(3,1) * t317 + rSges(3,2) * t320) * t312;
t262 = Icges(3,5) * t314 + (Icges(3,1) * t317 + Icges(3,4) * t320) * t312;
t260 = Icges(3,3) * t314 + (Icges(3,5) * t317 + Icges(3,6) * t320) * t312;
t258 = V_base(5) * rSges(2,3) - t295 * t305 + t342;
t257 = t296 * t305 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t255 = t295 * V_base(4) - t296 * V_base(5) + V_base(3);
t244 = t264 * t319 - t337;
t243 = -t264 * t316 - t319 * t351;
t238 = t264 * t307 - t306 * t351;
t237 = -t264 * t306 - t307 * t351;
t235 = rSges(3,1) * t277 - rSges(3,2) * t276 + rSges(3,3) * t352;
t234 = rSges(3,1) * t275 - rSges(3,2) * t274 - rSges(3,3) * t350;
t233 = Icges(3,1) * t277 - Icges(3,4) * t276 + Icges(3,5) * t352;
t232 = Icges(3,1) * t275 - Icges(3,4) * t274 - Icges(3,5) * t350;
t229 = Icges(3,5) * t277 - Icges(3,6) * t276 + Icges(3,3) * t352;
t228 = Icges(3,5) * t275 - Icges(3,6) * t274 - Icges(3,3) * t350;
t225 = rSges(4,1) * t273 + rSges(4,2) * t272 - rSges(4,3) * t351;
t224 = Icges(4,1) * t273 + Icges(4,4) * t272 - Icges(4,5) * t351;
t223 = Icges(4,4) * t273 + Icges(4,2) * t272 - Icges(4,6) * t351;
t219 = rSges(5,1) * t264 - rSges(5,2) * t263 - rSges(5,3) * t351;
t218 = Icges(5,1) * t264 - Icges(5,4) * t263 - Icges(5,5) * t351;
t217 = Icges(5,4) * t264 - Icges(5,2) * t263 - Icges(5,6) * t351;
t216 = Icges(5,5) * t264 - Icges(5,6) * t263 - Icges(5,3) * t351;
t215 = t248 * t319 + t355;
t214 = -t248 * t316 + t276 * t319;
t213 = t246 * t319 + t356;
t212 = -t246 * t316 + t274 * t319;
t210 = t248 * t307 + t276 * t306;
t209 = -t248 * t306 + t276 * t307;
t208 = t246 * t307 + t274 * t306;
t207 = -t246 * t306 + t274 * t307;
t206 = qJD(6) * t263 + t236;
t200 = rSges(4,1) * t252 + rSges(4,2) * t251 + rSges(4,3) * t276;
t199 = rSges(4,1) * t250 + rSges(4,2) * t249 + rSges(4,3) * t274;
t198 = Icges(4,1) * t252 + Icges(4,4) * t251 + Icges(4,5) * t276;
t197 = Icges(4,1) * t250 + Icges(4,4) * t249 + Icges(4,5) * t274;
t196 = Icges(4,4) * t252 + Icges(4,2) * t251 + Icges(4,6) * t276;
t195 = Icges(4,4) * t250 + Icges(4,2) * t249 + Icges(4,6) * t274;
t190 = rSges(5,1) * t248 - rSges(5,2) * t247 + rSges(5,3) * t276;
t189 = rSges(5,1) * t246 - rSges(5,2) * t245 + rSges(5,3) * t274;
t187 = Icges(5,1) * t248 - Icges(5,4) * t247 + Icges(5,5) * t276;
t186 = Icges(5,1) * t246 - Icges(5,4) * t245 + Icges(5,5) * t274;
t185 = Icges(5,4) * t248 - Icges(5,2) * t247 + Icges(5,6) * t276;
t184 = Icges(5,4) * t246 - Icges(5,2) * t245 + Icges(5,6) * t274;
t183 = Icges(5,5) * t248 - Icges(5,6) * t247 + Icges(5,3) * t276;
t182 = Icges(5,5) * t246 - Icges(5,6) * t245 + Icges(5,3) * t274;
t181 = rSges(6,1) * t244 + rSges(6,2) * t243 + rSges(6,3) * t263;
t180 = Icges(6,1) * t244 + Icges(6,4) * t243 + Icges(6,5) * t263;
t179 = Icges(6,4) * t244 + Icges(6,2) * t243 + Icges(6,6) * t263;
t178 = Icges(6,5) * t244 + Icges(6,6) * t243 + Icges(6,3) * t263;
t175 = qJD(6) * t247 + t205;
t174 = qJD(6) * t245 + t204;
t172 = rSges(7,1) * t238 + rSges(7,2) * t237 + rSges(7,3) * t263;
t171 = Icges(7,1) * t238 + Icges(7,4) * t237 + Icges(7,5) * t263;
t170 = Icges(7,4) * t238 + Icges(7,2) * t237 + Icges(7,6) * t263;
t169 = Icges(7,5) * t238 + Icges(7,6) * t237 + Icges(7,3) * t263;
t168 = -pkin(5) * t337 + pkin(11) * t263 + t264 * t359;
t167 = -t234 * t288 + t265 * t286 + t334;
t166 = t235 * t288 - t265 * t287 + t331;
t165 = t234 * t287 - t235 * t286 + t333;
t164 = rSges(6,1) * t215 + rSges(6,2) * t214 + rSges(6,3) * t247;
t163 = rSges(6,1) * t213 + rSges(6,2) * t212 + rSges(6,3) * t245;
t162 = Icges(6,1) * t215 + Icges(6,4) * t214 + Icges(6,5) * t247;
t161 = Icges(6,1) * t213 + Icges(6,4) * t212 + Icges(6,5) * t245;
t160 = Icges(6,4) * t215 + Icges(6,2) * t214 + Icges(6,6) * t247;
t159 = Icges(6,4) * t213 + Icges(6,2) * t212 + Icges(6,6) * t245;
t158 = Icges(6,5) * t215 + Icges(6,6) * t214 + Icges(6,3) * t247;
t157 = Icges(6,5) * t213 + Icges(6,6) * t212 + Icges(6,3) * t245;
t156 = rSges(7,1) * t210 + rSges(7,2) * t209 + rSges(7,3) * t247;
t155 = rSges(7,1) * t208 + rSges(7,2) * t207 + rSges(7,3) * t245;
t154 = Icges(7,1) * t210 + Icges(7,4) * t209 + Icges(7,5) * t247;
t153 = Icges(7,1) * t208 + Icges(7,4) * t207 + Icges(7,5) * t245;
t152 = Icges(7,4) * t210 + Icges(7,2) * t209 + Icges(7,6) * t247;
t151 = Icges(7,4) * t208 + Icges(7,2) * t207 + Icges(7,6) * t245;
t150 = Icges(7,5) * t210 + Icges(7,6) * t209 + Icges(7,3) * t247;
t149 = Icges(7,5) * t208 + Icges(7,6) * t207 + Icges(7,3) * t245;
t148 = pkin(5) * t355 + pkin(11) * t247 + t248 * t359;
t147 = pkin(5) * t356 + pkin(11) * t245 + t246 * t359;
t146 = t225 * t286 + (-t199 - t239) * t288 + t332;
t145 = t200 * t288 + (-t225 - t278) * t287 + t330;
t144 = t199 * t287 + (-t200 - t240) * t286 + t329;
t143 = -t189 * t270 + t219 * t253 + t328;
t142 = t190 * t270 - t219 * t254 + t327;
t141 = t189 * t254 - t190 * t253 + t326;
t140 = -t163 * t236 + t181 * t204 + t325;
t139 = t164 * t236 - t181 * t205 + t324;
t138 = t163 * t205 - t164 * t204 + t323;
t137 = -t147 * t236 - t155 * t206 + t168 * t204 + t172 * t174 + t325;
t136 = t148 * t236 + t156 * t206 - t168 * t205 - t172 * t175 + t324;
t135 = t147 * t205 - t148 * t204 + t155 * t175 - t156 * t174 + t323;
t1 = m(7) * (t135 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + m(6) * (t138 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + m(5) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(4) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(3) * (t165 ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + t270 * ((-t183 * t351 - t185 * t263 + t187 * t264) * t254 + (-t182 * t351 - t184 * t263 + t186 * t264) * t253 + (-t216 * t351 - t217 * t263 + t218 * t264) * t270) / 0.2e1 + t174 * ((t150 * t245 + t152 * t207 + t154 * t208) * t175 + (t245 * t149 + t207 * t151 + t208 * t153) * t174 + (t169 * t245 + t170 * t207 + t171 * t208) * t206) / 0.2e1 + t204 * ((t158 * t245 + t160 * t212 + t162 * t213) * t205 + (t245 * t157 + t212 * t159 + t213 * t161) * t204 + (t178 * t245 + t179 * t212 + t180 * t213) * t236) / 0.2e1 + t175 * ((t247 * t150 + t209 * t152 + t210 * t154) * t175 + (t149 * t247 + t151 * t209 + t153 * t210) * t174 + (t169 * t247 + t170 * t209 + t171 * t210) * t206) / 0.2e1 + t205 * ((t247 * t158 + t214 * t160 + t215 * t162) * t205 + (t157 * t247 + t159 * t214 + t161 * t215) * t204 + (t178 * t247 + t179 * t214 + t180 * t215) * t236) / 0.2e1 + m(2) * (t255 ^ 2 + t257 ^ 2 + t258 ^ 2) / 0.2e1 + t206 * ((t150 * t263 + t152 * t237 + t154 * t238) * t175 + (t149 * t263 + t151 * t237 + t153 * t238) * t174 + (t263 * t169 + t237 * t170 + t238 * t171) * t206) / 0.2e1 + t236 * ((t158 * t263 + t160 * t243 + t162 * t244) * t205 + (t157 * t263 + t159 * t243 + t161 * t244) * t204 + (t178 * t263 + t179 * t243 + t180 * t244) * t236) / 0.2e1 + t253 * ((t183 * t274 - t185 * t245 + t187 * t246) * t254 + (t182 * t274 - t184 * t245 + t186 * t246) * t253 + (t216 * t274 - t217 * t245 + t218 * t246) * t270) / 0.2e1 + t254 * ((t183 * t276 - t185 * t247 + t187 * t248) * t254 + (t182 * t276 - t184 * t247 + t186 * t248) * t253 + (t216 * t276 - t217 * t247 + t218 * t248) * t270) / 0.2e1 + m(1) * (t282 ^ 2 + t283 ^ 2 + t284 ^ 2) / 0.2e1 + ((t223 * t249 + t224 * t250 - t260 * t350 + t262 * t275 + t274 * t365) * t288 + (t196 * t249 + t198 * t250 - t229 * t350 + t233 * t275 + t274 * t366) * t287 + (t195 * t249 + t197 * t250 - t228 * t350 + t232 * t275 + t274 * t367) * t286) * t286 / 0.2e1 + ((t223 * t251 + t224 * t252 + t260 * t352 + t262 * t277 + t276 * t365) * t288 + (t196 * t251 + t198 * t252 + t229 * t352 + t233 * t277 + t276 * t366) * t287 + (t195 * t251 + t197 * t252 + t228 * t352 + t232 * t277 + t276 * t367) * t286) * t287 / 0.2e1 + ((t228 * t286 + t229 * t287 + t260 * t288) * t314 + ((t231 * t320 + t233 * t317) * t287 + (t230 * t320 + t232 * t317) * t286 + (t261 * t320 + t262 * t317) * t288) * t312 + (-t194 * t351 + t196 * t272 + t198 * t273) * t287 + (-t193 * t351 + t195 * t272 + t197 * t273) * t286 + (-t222 * t351 + t223 * t272 + t224 * t273) * t288) * t288 / 0.2e1 + ((-t291 * t318 + t293 * t321 + Icges(1,4)) * V_base(5) + (-t292 * t318 + t294 * t321 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t291 * t321 + t293 * t318 + Icges(1,2)) * V_base(5) + (t292 * t321 + t294 * t318 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t318 + Icges(2,6) * t321) * V_base(5) + (Icges(2,5) * t321 - Icges(2,6) * t318) * V_base(4) + Icges(2,3) * t305 / 0.2e1) * t305;
T  = t1;
