% Calculate kinetic energy for
% S6PRRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRPP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPP2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPP2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:48:39
% EndTime: 2019-03-08 22:48:44
% DurationCPUTime: 4.61s
% Computational Cost: add. (2785->334), mult. (6592->472), div. (0->0), fcn. (8001->10), ass. (0->149)
t356 = Icges(5,1) + Icges(6,1) + Icges(7,1);
t355 = -Icges(5,4) + Icges(7,4) + Icges(6,5);
t354 = Icges(6,4) + Icges(5,5) - Icges(7,5);
t353 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t352 = Icges(6,2) + Icges(5,3) + Icges(7,3);
t351 = -Icges(5,6) + Icges(6,6) - Icges(7,6);
t350 = rSges(7,1) + pkin(5);
t349 = -rSges(7,3) - qJ(6);
t295 = sin(pkin(10));
t297 = cos(pkin(10));
t301 = cos(qJ(2));
t300 = sin(qJ(2));
t332 = cos(pkin(6));
t315 = t300 * t332;
t262 = t295 * t301 + t297 * t315;
t296 = sin(pkin(6));
t299 = sin(qJ(3));
t328 = t296 * t299;
t334 = cos(qJ(3));
t245 = t262 * t334 - t297 * t328;
t314 = t301 * t332;
t261 = t295 * t300 - t297 * t314;
t298 = sin(qJ(4));
t333 = cos(qJ(4));
t217 = t245 * t298 - t261 * t333;
t218 = t245 * t333 + t261 * t298;
t317 = t296 * t334;
t244 = t262 * t299 + t297 * t317;
t346 = t217 * t351 + t218 * t354 + t244 * t352;
t264 = -t295 * t315 + t297 * t301;
t247 = t264 * t334 + t295 * t328;
t263 = t295 * t314 + t297 * t300;
t219 = t247 * t298 - t263 * t333;
t220 = t247 * t333 + t263 * t298;
t246 = t264 * t299 - t295 * t317;
t345 = t219 * t351 + t220 * t354 + t246 * t352;
t344 = t217 * t353 + t218 * t355 + t244 * t351;
t343 = t219 * t353 + t220 * t355 + t246 * t351;
t342 = t217 * t355 + t218 * t356 + t244 * t354;
t341 = t219 * t355 + t220 * t356 + t246 * t354;
t269 = t299 * t332 + t300 * t317;
t327 = t296 * t301;
t248 = t269 * t298 + t327 * t333;
t249 = t269 * t333 - t298 * t327;
t268 = t300 * t328 - t332 * t334;
t340 = t248 * t351 + t249 * t354 + t268 * t352;
t339 = t248 * t353 + t249 * t355 + t268 * t351;
t338 = t248 * t355 + t249 * t356 + t268 * t354;
t331 = Icges(2,4) * t295;
t330 = t295 * t296;
t329 = t296 * t297;
t326 = rSges(7,2) * t217 + t218 * t350 + t244 * t349;
t325 = rSges(7,2) * t219 + t220 * t350 + t246 * t349;
t324 = rSges(7,2) * t248 + t249 * t350 + t268 * t349;
t323 = qJD(2) * t296;
t322 = V_base(5) * qJ(1) + V_base(1);
t318 = qJD(1) + V_base(3);
t316 = t332 * pkin(7);
t277 = t295 * t323 + V_base(4);
t288 = qJD(2) * t332 + V_base(6);
t243 = qJD(3) * t263 + t277;
t276 = -t297 * t323 + V_base(5);
t242 = qJD(3) * t261 + t276;
t265 = -qJD(3) * t327 + t288;
t271 = pkin(1) * t295 - pkin(7) * t329;
t313 = -t271 * V_base(6) + t316 * V_base(5) + t322;
t272 = pkin(1) * t297 + pkin(7) * t330;
t312 = t271 * V_base(4) - t272 * V_base(5) + t318;
t311 = V_base(6) * t272 + V_base(2) + (-t316 - qJ(1)) * V_base(4);
t236 = pkin(2) * t262 + pkin(8) * t261;
t270 = (pkin(2) * t300 - pkin(8) * t301) * t296;
t310 = -t236 * t288 + t270 * t276 + t313;
t237 = pkin(2) * t264 + pkin(8) * t263;
t309 = t236 * t277 - t237 * t276 + t312;
t308 = t237 * t288 - t270 * t277 + t311;
t210 = pkin(3) * t245 + pkin(9) * t244;
t238 = pkin(3) * t269 + pkin(9) * t268;
t307 = -t210 * t265 + t238 * t242 + t310;
t211 = pkin(3) * t247 + pkin(9) * t246;
t306 = t210 * t243 - t211 * t242 + t309;
t212 = pkin(4) * t249 + qJ(5) * t248;
t213 = qJD(4) * t244 + t242;
t305 = qJD(5) * t219 + t212 * t213 + t307;
t183 = pkin(4) * t218 + qJ(5) * t217;
t214 = qJD(4) * t246 + t243;
t304 = qJD(5) * t248 + t183 * t214 + t306;
t303 = t211 * t265 - t238 * t243 + t308;
t184 = pkin(4) * t220 + qJ(5) * t219;
t239 = qJD(4) * t268 + t265;
t302 = qJD(5) * t217 + t184 * t239 + t303;
t293 = Icges(2,4) * t297;
t285 = rSges(2,1) * t297 - rSges(2,2) * t295;
t284 = rSges(2,1) * t295 + rSges(2,2) * t297;
t283 = Icges(2,1) * t297 - t331;
t282 = Icges(2,1) * t295 + t293;
t281 = -Icges(2,2) * t295 + t293;
t280 = Icges(2,2) * t297 + t331;
t275 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t274 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t273 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t256 = t332 * rSges(3,3) + (rSges(3,1) * t300 + rSges(3,2) * t301) * t296;
t255 = Icges(3,5) * t332 + (Icges(3,1) * t300 + Icges(3,4) * t301) * t296;
t254 = Icges(3,6) * t332 + (Icges(3,4) * t300 + Icges(3,2) * t301) * t296;
t253 = Icges(3,3) * t332 + (Icges(3,5) * t300 + Icges(3,6) * t301) * t296;
t252 = V_base(5) * rSges(2,3) - t284 * V_base(6) + t322;
t251 = t285 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t240 = t284 * V_base(4) - t285 * V_base(5) + t318;
t235 = rSges(4,1) * t269 - rSges(4,2) * t268 - rSges(4,3) * t327;
t234 = Icges(4,1) * t269 - Icges(4,4) * t268 - Icges(4,5) * t327;
t233 = Icges(4,4) * t269 - Icges(4,2) * t268 - Icges(4,6) * t327;
t232 = Icges(4,5) * t269 - Icges(4,6) * t268 - Icges(4,3) * t327;
t231 = rSges(3,1) * t264 - rSges(3,2) * t263 + rSges(3,3) * t330;
t230 = rSges(3,1) * t262 - rSges(3,2) * t261 - rSges(3,3) * t329;
t229 = Icges(3,1) * t264 - Icges(3,4) * t263 + Icges(3,5) * t330;
t228 = Icges(3,1) * t262 - Icges(3,4) * t261 - Icges(3,5) * t329;
t227 = Icges(3,4) * t264 - Icges(3,2) * t263 + Icges(3,6) * t330;
t226 = Icges(3,4) * t262 - Icges(3,2) * t261 - Icges(3,6) * t329;
t225 = Icges(3,5) * t264 - Icges(3,6) * t263 + Icges(3,3) * t330;
t224 = Icges(3,5) * t262 - Icges(3,6) * t261 - Icges(3,3) * t329;
t208 = rSges(5,1) * t249 - rSges(5,2) * t248 + rSges(5,3) * t268;
t207 = rSges(6,1) * t249 + rSges(6,2) * t268 + rSges(6,3) * t248;
t195 = rSges(4,1) * t247 - rSges(4,2) * t246 + rSges(4,3) * t263;
t194 = rSges(4,1) * t245 - rSges(4,2) * t244 + rSges(4,3) * t261;
t193 = Icges(4,1) * t247 - Icges(4,4) * t246 + Icges(4,5) * t263;
t192 = Icges(4,1) * t245 - Icges(4,4) * t244 + Icges(4,5) * t261;
t191 = Icges(4,4) * t247 - Icges(4,2) * t246 + Icges(4,6) * t263;
t190 = Icges(4,4) * t245 - Icges(4,2) * t244 + Icges(4,6) * t261;
t189 = Icges(4,5) * t247 - Icges(4,6) * t246 + Icges(4,3) * t263;
t188 = Icges(4,5) * t245 - Icges(4,6) * t244 + Icges(4,3) * t261;
t182 = -t230 * t288 + t256 * t276 + t313;
t181 = t231 * t288 - t256 * t277 + t311;
t179 = rSges(5,1) * t220 - rSges(5,2) * t219 + rSges(5,3) * t246;
t178 = rSges(6,1) * t220 + rSges(6,2) * t246 + rSges(6,3) * t219;
t176 = rSges(5,1) * t218 - rSges(5,2) * t217 + rSges(5,3) * t244;
t175 = rSges(6,1) * t218 + rSges(6,2) * t244 + rSges(6,3) * t217;
t154 = t230 * t277 - t231 * t276 + t312;
t152 = -t194 * t265 + t235 * t242 + t310;
t151 = t195 * t265 - t235 * t243 + t308;
t150 = t194 * t243 - t195 * t242 + t309;
t149 = -t176 * t239 + t208 * t213 + t307;
t148 = t179 * t239 - t208 * t214 + t303;
t147 = t176 * t214 - t179 * t213 + t306;
t146 = t207 * t213 + (-t175 - t183) * t239 + t305;
t145 = t239 * t178 + (-t207 - t212) * t214 + t302;
t144 = t175 * t214 + (-t178 - t184) * t213 + t304;
t143 = -qJD(6) * t246 + t324 * t213 + (-t183 - t326) * t239 + t305;
t142 = -qJD(6) * t244 + t325 * t239 + (-t212 - t324) * t214 + t302;
t141 = -qJD(6) * t268 + t326 * t214 + (-t184 - t325) * t213 + t304;
t1 = m(1) * (t273 ^ 2 + t274 ^ 2 + t275 ^ 2) / 0.2e1 + t242 * ((t189 * t261 - t191 * t244 + t193 * t245) * t243 + (t261 * t188 - t244 * t190 + t245 * t192) * t242 + (t232 * t261 - t233 * t244 + t234 * t245) * t265) / 0.2e1 + t243 * ((t263 * t189 - t246 * t191 + t247 * t193) * t243 + (t188 * t263 - t190 * t246 + t192 * t247) * t242 + (t232 * t263 - t233 * t246 + t234 * t247) * t265) / 0.2e1 + m(2) * (t240 ^ 2 + t251 ^ 2 + t252 ^ 2) / 0.2e1 + m(3) * (t154 ^ 2 + t181 ^ 2 + t182 ^ 2) / 0.2e1 + m(4) * (t150 ^ 2 + t151 ^ 2 + t152 ^ 2) / 0.2e1 + m(5) * (t147 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + m(6) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(7) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + t288 * (((t227 * t301 + t229 * t300) * t277 + (t226 * t301 + t228 * t300) * t276 + (t254 * t301 + t255 * t300) * t288) * t296 + (t224 * t276 + t225 * t277 + t253 * t288) * t332) / 0.2e1 + t277 * ((t225 * t330 - t263 * t227 + t264 * t229) * t277 + (t224 * t330 - t226 * t263 + t228 * t264) * t276 + (t253 * t330 - t254 * t263 + t255 * t264) * t288) / 0.2e1 + t265 * ((-t189 * t327 - t191 * t268 + t193 * t269) * t243 + (-t188 * t327 - t190 * t268 + t192 * t269) * t242 + (-t232 * t327 - t268 * t233 + t269 * t234) * t265) / 0.2e1 + t276 * ((-t225 * t329 - t227 * t261 + t229 * t262) * t277 + (-t224 * t329 - t261 * t226 + t262 * t228) * t276 + (-t253 * t329 - t254 * t261 + t255 * t262) * t288) / 0.2e1 + ((-t280 * t295 + t282 * t297 + Icges(1,4)) * V_base(5) + (-t295 * t281 + t297 * t283 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t297 * t280 + t295 * t282 + Icges(1,2)) * V_base(5) + (t281 * t297 + t283 * t295 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t217 * t339 + t218 * t338 + t244 * t340) * t239 + (t217 * t343 + t218 * t341 + t244 * t345) * t214 + (t344 * t217 + t342 * t218 + t346 * t244) * t213) * t213 / 0.2e1 + ((t219 * t339 + t220 * t338 + t246 * t340) * t239 + (t343 * t219 + t341 * t220 + t345 * t246) * t214 + (t219 * t344 + t220 * t342 + t246 * t346) * t213) * t214 / 0.2e1 + ((t339 * t248 + t338 * t249 + t340 * t268) * t239 + (t248 * t343 + t249 * t341 + t268 * t345) * t214 + (t248 * t344 + t249 * t342 + t268 * t346) * t213) * t239 / 0.2e1 + ((Icges(2,5) * t295 + Icges(2,6) * t297 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t297 - Icges(2,6) * t295 + Icges(1,5)) * V_base(4) + (Icges(2,3) / 0.2e1 + Icges(1,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
