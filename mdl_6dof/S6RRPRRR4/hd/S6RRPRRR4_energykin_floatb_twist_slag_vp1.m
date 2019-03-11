% Calculate kinetic energy for
% S6RRPRRR4
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
% Datum: 2019-03-09 13:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:27:40
% EndTime: 2019-03-09 13:27:46
% DurationCPUTime: 5.82s
% Computational Cost: add. (3997->441), mult. (8194->646), div. (0->0), fcn. (10180->14), ass. (0->196)
t337 = sin(qJ(2));
t378 = sin(pkin(12));
t379 = cos(pkin(12));
t384 = cos(qJ(2));
t301 = -t337 * t379 - t378 * t384;
t338 = sin(qJ(1));
t341 = cos(qJ(1));
t380 = cos(pkin(6));
t356 = t380 * t378;
t357 = t380 * t379;
t348 = -t337 * t356 + t357 * t384;
t271 = t301 * t338 + t341 * t348;
t293 = t337 * t357 + t356 * t384;
t302 = -t337 * t378 + t379 * t384;
t272 = t293 * t341 + t302 * t338;
t334 = sin(pkin(6));
t375 = t334 * t341;
t220 = Icges(4,5) * t272 + Icges(4,6) * t271 - Icges(4,3) * t375;
t360 = t380 * t384;
t296 = -t337 * t338 + t341 * t360;
t363 = t337 * t380;
t297 = t338 * t384 + t341 * t363;
t256 = Icges(3,5) * t297 + Icges(3,6) * t296 - Icges(3,3) * t375;
t390 = t220 + t256;
t273 = t301 * t341 - t338 * t348;
t274 = -t293 * t338 + t302 * t341;
t376 = t334 * t338;
t221 = Icges(4,5) * t274 + Icges(4,6) * t273 + Icges(4,3) * t376;
t298 = -t337 * t341 - t338 * t360;
t299 = -t338 * t363 + t341 * t384;
t257 = Icges(3,5) * t299 + Icges(3,6) * t298 + Icges(3,3) * t376;
t389 = t221 + t257;
t291 = t302 * t334;
t292 = t301 * t334;
t252 = -Icges(4,5) * t292 + Icges(4,6) * t291 + Icges(4,3) * t380;
t287 = Icges(3,3) * t380 + (Icges(3,5) * t337 + Icges(3,6) * t384) * t334;
t388 = t252 + t287;
t383 = pkin(2) * t384;
t340 = cos(qJ(4));
t382 = pkin(4) * t340;
t377 = Icges(2,4) * t338;
t374 = qJ(4) + qJ(5);
t373 = qJD(2) * t334;
t372 = qJD(3) * t334;
t371 = V_base(5) * pkin(7) + V_base(1);
t336 = sin(qJ(4));
t368 = t336 * t376;
t367 = t336 * t375;
t366 = t380 * pkin(8);
t310 = t338 * t373 + V_base(4);
t365 = cos(t374);
t330 = V_base(6) + qJD(1);
t364 = pkin(2) * t363 - qJ(3) * t334;
t361 = t380 * t336;
t245 = -qJD(4) * t273 + t310;
t312 = qJD(2) * t380 + t330;
t358 = t334 * t365;
t214 = -qJD(5) * t273 + t245;
t277 = -qJD(4) * t291 + t312;
t309 = -t341 * t373 + V_base(5);
t304 = pkin(1) * t338 - pkin(8) * t375;
t355 = -t304 * t330 + t366 * V_base(5) + t371;
t250 = -qJD(5) * t291 + t277;
t305 = pkin(1) * t341 + pkin(8) * t376;
t354 = t304 * V_base(4) - t305 * V_base(5) + V_base(3);
t244 = -qJD(4) * t271 + t309;
t213 = -qJD(5) * t271 + t244;
t303 = pkin(2) * t334 * t337 + qJ(3) * t380;
t353 = t303 * t309 + t338 * t372 + t355;
t269 = t338 * t383 + t341 * t364;
t352 = qJD(3) * t380 + t269 * t310 + t354;
t351 = t330 * t305 + V_base(2) + (-t366 - pkin(7)) * V_base(4);
t233 = pkin(3) * t272 - pkin(9) * t271;
t264 = -pkin(3) * t292 - pkin(9) * t291;
t350 = t309 * t264 + (-t233 - t269) * t312 + t353;
t234 = pkin(3) * t274 - pkin(9) * t273;
t270 = -t338 * t364 + t341 * t383;
t349 = t310 * t233 + (-t234 - t270) * t309 + t352;
t347 = t270 * t312 - t341 * t372 + t351;
t177 = -pkin(4) * t367 - pkin(10) * t271 + t272 * t382;
t211 = pkin(4) * t361 - pkin(10) * t291 - t292 * t382;
t346 = -t177 * t277 + t211 * t244 + t350;
t178 = pkin(4) * t368 - pkin(10) * t273 + t274 * t382;
t345 = t177 * t245 - t178 * t244 + t349;
t344 = t312 * t234 + (-t264 - t303) * t310 + t347;
t343 = t178 * t277 - t211 * t245 + t344;
t339 = cos(qJ(6));
t335 = sin(qJ(6));
t332 = Icges(2,4) * t341;
t331 = sin(t374);
t320 = rSges(2,1) * t341 - rSges(2,2) * t338;
t319 = rSges(2,1) * t338 + rSges(2,2) * t341;
t318 = Icges(2,1) * t341 - t377;
t317 = Icges(2,1) * t338 + t332;
t316 = -Icges(2,2) * t338 + t332;
t315 = Icges(2,2) * t341 + t377;
t308 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t307 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t306 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t290 = t380 * rSges(3,3) + (rSges(3,1) * t337 + rSges(3,2) * t384) * t334;
t289 = Icges(3,5) * t380 + (Icges(3,1) * t337 + Icges(3,4) * t384) * t334;
t288 = Icges(3,6) * t380 + (Icges(3,4) * t337 + Icges(3,2) * t384) * t334;
t282 = V_base(5) * rSges(2,3) - t319 * t330 + t371;
t281 = t320 * t330 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t280 = t319 * V_base(4) - t320 * V_base(5) + V_base(3);
t279 = -t292 * t340 + t361;
t278 = t292 * t336 + t340 * t380;
t276 = -t292 * t365 + t331 * t380;
t275 = -t292 * t331 - t365 * t380;
t263 = rSges(3,1) * t299 + rSges(3,2) * t298 + rSges(3,3) * t376;
t262 = rSges(3,1) * t297 + rSges(3,2) * t296 - rSges(3,3) * t375;
t261 = Icges(3,1) * t299 + Icges(3,4) * t298 + Icges(3,5) * t376;
t260 = Icges(3,1) * t297 + Icges(3,4) * t296 - Icges(3,5) * t375;
t259 = Icges(3,4) * t299 + Icges(3,2) * t298 + Icges(3,6) * t376;
t258 = Icges(3,4) * t297 + Icges(3,2) * t296 - Icges(3,6) * t375;
t255 = -rSges(4,1) * t292 + rSges(4,2) * t291 + rSges(4,3) * t380;
t254 = -Icges(4,1) * t292 + Icges(4,4) * t291 + Icges(4,5) * t380;
t253 = -Icges(4,4) * t292 + Icges(4,2) * t291 + Icges(4,6) * t380;
t249 = t274 * t340 + t368;
t248 = -t274 * t336 + t340 * t376;
t247 = t272 * t340 - t367;
t246 = -t272 * t336 - t340 * t375;
t241 = t274 * t365 + t331 * t376;
t240 = t274 * t331 - t338 * t358;
t239 = t272 * t365 - t331 * t375;
t238 = t272 * t331 + t341 * t358;
t237 = t276 * t339 - t291 * t335;
t236 = -t276 * t335 - t291 * t339;
t235 = pkin(5) * t276 + pkin(11) * t275;
t232 = rSges(5,1) * t279 + rSges(5,2) * t278 - rSges(5,3) * t291;
t231 = Icges(5,1) * t279 + Icges(5,4) * t278 - Icges(5,5) * t291;
t230 = Icges(5,4) * t279 + Icges(5,2) * t278 - Icges(5,6) * t291;
t229 = Icges(5,5) * t279 + Icges(5,6) * t278 - Icges(5,3) * t291;
t228 = qJD(6) * t275 + t250;
t227 = rSges(4,1) * t274 + rSges(4,2) * t273 + rSges(4,3) * t376;
t226 = rSges(4,1) * t272 + rSges(4,2) * t271 - rSges(4,3) * t375;
t225 = Icges(4,1) * t274 + Icges(4,4) * t273 + Icges(4,5) * t376;
t224 = Icges(4,1) * t272 + Icges(4,4) * t271 - Icges(4,5) * t375;
t223 = Icges(4,4) * t274 + Icges(4,2) * t273 + Icges(4,6) * t376;
t222 = Icges(4,4) * t272 + Icges(4,2) * t271 - Icges(4,6) * t375;
t219 = rSges(6,1) * t276 - rSges(6,2) * t275 - rSges(6,3) * t291;
t218 = Icges(6,1) * t276 - Icges(6,4) * t275 - Icges(6,5) * t291;
t217 = Icges(6,4) * t276 - Icges(6,2) * t275 - Icges(6,6) * t291;
t216 = Icges(6,5) * t276 - Icges(6,6) * t275 - Icges(6,3) * t291;
t210 = t241 * t339 - t273 * t335;
t209 = -t241 * t335 - t273 * t339;
t208 = t239 * t339 - t271 * t335;
t207 = -t239 * t335 - t271 * t339;
t206 = -t262 * t312 + t290 * t309 + t355;
t205 = t263 * t312 - t290 * t310 + t351;
t204 = pkin(5) * t241 + pkin(11) * t240;
t203 = pkin(5) * t239 + pkin(11) * t238;
t202 = t262 * t310 - t263 * t309 + t354;
t200 = rSges(5,1) * t249 + rSges(5,2) * t248 - rSges(5,3) * t273;
t199 = rSges(5,1) * t247 + rSges(5,2) * t246 - rSges(5,3) * t271;
t198 = Icges(5,1) * t249 + Icges(5,4) * t248 - Icges(5,5) * t273;
t197 = Icges(5,1) * t247 + Icges(5,4) * t246 - Icges(5,5) * t271;
t196 = Icges(5,4) * t249 + Icges(5,2) * t248 - Icges(5,6) * t273;
t195 = Icges(5,4) * t247 + Icges(5,2) * t246 - Icges(5,6) * t271;
t194 = Icges(5,5) * t249 + Icges(5,6) * t248 - Icges(5,3) * t273;
t193 = Icges(5,5) * t247 + Icges(5,6) * t246 - Icges(5,3) * t271;
t192 = qJD(6) * t240 + t214;
t191 = qJD(6) * t238 + t213;
t190 = rSges(6,1) * t241 - rSges(6,2) * t240 - rSges(6,3) * t273;
t189 = rSges(6,1) * t239 - rSges(6,2) * t238 - rSges(6,3) * t271;
t188 = rSges(7,1) * t237 + rSges(7,2) * t236 + rSges(7,3) * t275;
t187 = Icges(6,1) * t241 - Icges(6,4) * t240 - Icges(6,5) * t273;
t186 = Icges(6,1) * t239 - Icges(6,4) * t238 - Icges(6,5) * t271;
t185 = Icges(6,4) * t241 - Icges(6,2) * t240 - Icges(6,6) * t273;
t184 = Icges(6,4) * t239 - Icges(6,2) * t238 - Icges(6,6) * t271;
t183 = Icges(6,5) * t241 - Icges(6,6) * t240 - Icges(6,3) * t273;
t182 = Icges(6,5) * t239 - Icges(6,6) * t238 - Icges(6,3) * t271;
t181 = Icges(7,1) * t237 + Icges(7,4) * t236 + Icges(7,5) * t275;
t180 = Icges(7,4) * t237 + Icges(7,2) * t236 + Icges(7,6) * t275;
t179 = Icges(7,5) * t237 + Icges(7,6) * t236 + Icges(7,3) * t275;
t174 = rSges(7,1) * t210 + rSges(7,2) * t209 + rSges(7,3) * t240;
t173 = rSges(7,1) * t208 + rSges(7,2) * t207 + rSges(7,3) * t238;
t172 = Icges(7,1) * t210 + Icges(7,4) * t209 + Icges(7,5) * t240;
t171 = Icges(7,1) * t208 + Icges(7,4) * t207 + Icges(7,5) * t238;
t170 = Icges(7,4) * t210 + Icges(7,2) * t209 + Icges(7,6) * t240;
t169 = Icges(7,4) * t208 + Icges(7,2) * t207 + Icges(7,6) * t238;
t168 = Icges(7,5) * t210 + Icges(7,6) * t209 + Icges(7,3) * t240;
t167 = Icges(7,5) * t208 + Icges(7,6) * t207 + Icges(7,3) * t238;
t166 = t255 * t309 + (-t226 - t269) * t312 + t353;
t165 = t312 * t227 + (-t255 - t303) * t310 + t347;
t164 = t226 * t310 + (-t227 - t270) * t309 + t352;
t163 = -t199 * t277 + t232 * t244 + t350;
t162 = t200 * t277 - t232 * t245 + t344;
t161 = t199 * t245 - t200 * t244 + t349;
t160 = -t189 * t250 + t213 * t219 + t346;
t159 = t190 * t250 - t214 * t219 + t343;
t158 = t189 * t214 - t190 * t213 + t345;
t157 = -t173 * t228 + t188 * t191 - t203 * t250 + t213 * t235 + t346;
t156 = t174 * t228 - t188 * t192 + t204 * t250 - t214 * t235 + t343;
t155 = t173 * t192 - t174 * t191 + t203 * t214 - t204 * t213 + t345;
t1 = m(7) * (t155 ^ 2 + t156 ^ 2 + t157 ^ 2) / 0.2e1 + m(6) * (t158 ^ 2 + t159 ^ 2 + t160 ^ 2) / 0.2e1 + m(5) * (t161 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(4) * (t164 ^ 2 + t165 ^ 2 + t166 ^ 2) / 0.2e1 + m(3) * (t202 ^ 2 + t205 ^ 2 + t206 ^ 2) / 0.2e1 + t191 * ((t168 * t238 + t170 * t207 + t172 * t208) * t192 + (t238 * t167 + t207 * t169 + t208 * t171) * t191 + (t179 * t238 + t180 * t207 + t181 * t208) * t228) / 0.2e1 + t192 * ((t240 * t168 + t209 * t170 + t210 * t172) * t192 + (t167 * t240 + t169 * t209 + t171 * t210) * t191 + (t179 * t240 + t180 * t209 + t181 * t210) * t228) / 0.2e1 + t213 * ((-t183 * t271 - t185 * t238 + t187 * t239) * t214 + (-t271 * t182 - t238 * t184 + t239 * t186) * t213 + (-t216 * t271 - t217 * t238 + t218 * t239) * t250) / 0.2e1 + t214 * ((-t273 * t183 - t240 * t185 + t241 * t187) * t214 + (-t182 * t273 - t184 * t240 + t186 * t241) * t213 + (-t216 * t273 - t217 * t240 + t218 * t241) * t250) / 0.2e1 + t228 * ((t168 * t275 + t170 * t236 + t172 * t237) * t192 + (t167 * t275 + t169 * t236 + t171 * t237) * t191 + (t275 * t179 + t236 * t180 + t237 * t181) * t228) / 0.2e1 + t244 * ((-t194 * t271 + t196 * t246 + t198 * t247) * t245 + (-t271 * t193 + t246 * t195 + t247 * t197) * t244 + (-t229 * t271 + t230 * t246 + t231 * t247) * t277) / 0.2e1 + t245 * ((-t273 * t194 + t248 * t196 + t249 * t198) * t245 + (-t193 * t273 + t195 * t248 + t197 * t249) * t244 + (-t229 * t273 + t230 * t248 + t231 * t249) * t277) / 0.2e1 + m(2) * (t280 ^ 2 + t281 ^ 2 + t282 ^ 2) / 0.2e1 + t250 * ((-t183 * t291 - t185 * t275 + t187 * t276) * t214 + (-t182 * t291 - t184 * t275 + t186 * t276) * t213 + (-t291 * t216 - t275 * t217 + t276 * t218) * t250) / 0.2e1 + t277 * ((-t194 * t291 + t196 * t278 + t198 * t279) * t245 + (-t193 * t291 + t195 * t278 + t197 * t279) * t244 + (-t291 * t229 + t278 * t230 + t279 * t231) * t277) / 0.2e1 + m(1) * (t306 ^ 2 + t307 ^ 2 + t308 ^ 2) / 0.2e1 + ((t253 * t271 + t254 * t272 + t288 * t296 + t289 * t297 - t375 * t388) * t312 + (t223 * t271 + t225 * t272 + t259 * t296 + t261 * t297 - t375 * t389) * t310 + (t271 * t222 + t272 * t224 + t296 * t258 + t297 * t260 - t390 * t375) * t309) * t309 / 0.2e1 + ((t253 * t273 + t254 * t274 + t288 * t298 + t289 * t299 + t376 * t388) * t312 + (t273 * t223 + t274 * t225 + t298 * t259 + t299 * t261 + t389 * t376) * t310 + (t222 * t273 + t224 * t274 + t258 * t298 + t260 * t299 + t376 * t390) * t309) * t310 / 0.2e1 + ((t221 * t380 + t291 * t223 - t292 * t225) * t310 + (t220 * t380 + t291 * t222 - t292 * t224) * t309 + (t252 * t380 + t291 * t253 - t292 * t254) * t312 + ((t259 * t384 + t261 * t337) * t310 + (t258 * t384 + t260 * t337) * t309 + (t288 * t384 + t289 * t337) * t312) * t334 + (t256 * t309 + t257 * t310 + t287 * t312) * t380) * t312 / 0.2e1 + ((-t315 * t338 + t317 * t341 + Icges(1,4)) * V_base(5) + (-t316 * t338 + t318 * t341 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t315 * t341 + t317 * t338 + Icges(1,2)) * V_base(5) + (t316 * t341 + t318 * t338 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t338 + Icges(2,6) * t341) * V_base(5) + (Icges(2,5) * t341 - Icges(2,6) * t338) * V_base(4) + Icges(2,3) * t330 / 0.2e1) * t330;
T  = t1;
