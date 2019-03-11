% Calculate kinetic energy for
% S6PRPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRRP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRP1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:55:08
% EndTime: 2019-03-08 19:55:13
% DurationCPUTime: 4.56s
% Computational Cost: add. (3748->388), mult. (9009->549), div. (0->0), fcn. (11331->12), ass. (0->178)
t384 = Icges(6,1) + Icges(7,1);
t383 = Icges(6,4) + Icges(7,4);
t382 = Icges(6,5) + Icges(7,5);
t381 = Icges(6,2) + Icges(7,2);
t380 = Icges(6,6) + Icges(7,6);
t379 = Icges(6,3) + Icges(7,3);
t378 = rSges(7,3) + qJ(6);
t307 = cos(pkin(6));
t311 = sin(qJ(2));
t313 = cos(qJ(2));
t351 = sin(pkin(11));
t352 = cos(pkin(11));
t325 = t311 * t352 + t313 * t351;
t265 = t325 * t307;
t274 = -t311 * t351 + t313 * t352;
t304 = sin(pkin(10));
t306 = cos(pkin(10));
t247 = t265 * t306 + t274 * t304;
t305 = sin(pkin(6));
t310 = sin(qJ(4));
t344 = t305 * t310;
t357 = cos(qJ(4));
t225 = t247 * t357 - t306 * t344;
t321 = t307 * t274;
t246 = -t304 * t325 + t306 * t321;
t309 = sin(qJ(5));
t312 = cos(qJ(5));
t196 = -t225 * t309 - t246 * t312;
t349 = t246 * t309;
t197 = t225 * t312 - t349;
t331 = t305 * t357;
t224 = t247 * t310 + t306 * t331;
t376 = t196 * t380 + t197 * t382 + t224 * t379;
t249 = -t265 * t304 + t274 * t306;
t227 = t249 * t357 + t304 * t344;
t248 = -t304 * t321 - t306 * t325;
t198 = -t227 * t309 - t248 * t312;
t348 = t248 * t309;
t199 = t227 * t312 - t348;
t226 = t249 * t310 - t304 * t331;
t375 = t198 * t380 + t199 * t382 + t226 * t379;
t374 = t196 * t381 + t197 * t383 + t224 * t380;
t373 = t198 * t381 + t199 * t383 + t226 * t380;
t372 = t383 * t196 + t197 * t384 + t382 * t224;
t371 = t383 * t198 + t199 * t384 + t382 * t226;
t264 = t325 * t305;
t253 = t264 * t357 + t307 * t310;
t263 = t274 * t305;
t218 = -t253 * t309 - t263 * t312;
t347 = t263 * t309;
t219 = t253 * t312 - t347;
t252 = t264 * t310 - t307 * t357;
t370 = t218 * t380 + t219 * t382 + t252 * t379;
t369 = t218 * t381 + t219 * t383 + t252 * t380;
t368 = t383 * t218 + t219 * t384 + t382 * t252;
t345 = t305 * t306;
t202 = Icges(4,5) * t247 + Icges(4,6) * t246 - Icges(4,3) * t345;
t342 = t307 * t313;
t267 = -t304 * t311 + t306 * t342;
t343 = t307 * t311;
t268 = t304 * t313 + t306 * t343;
t233 = Icges(3,5) * t268 + Icges(3,6) * t267 - Icges(3,3) * t345;
t367 = t202 + t233;
t346 = t304 * t305;
t203 = Icges(4,5) * t249 + Icges(4,6) * t248 + Icges(4,3) * t346;
t269 = -t304 * t342 - t306 * t311;
t270 = -t304 * t343 + t306 * t313;
t234 = Icges(3,5) * t270 + Icges(3,6) * t269 + Icges(3,3) * t346;
t366 = t203 + t234;
t229 = Icges(4,5) * t264 + Icges(4,6) * t263 + Icges(4,3) * t307;
t259 = Icges(3,3) * t307 + (Icges(3,5) * t311 + Icges(3,6) * t313) * t305;
t365 = t229 + t259;
t356 = pkin(7) * t307;
t355 = pkin(2) * t313;
t354 = pkin(5) * t312;
t350 = Icges(2,4) * t304;
t341 = rSges(7,1) * t197 + rSges(7,2) * t196 - pkin(5) * t349 + t378 * t224 + t225 * t354;
t340 = rSges(7,1) * t199 + rSges(7,2) * t198 - pkin(5) * t348 + t378 * t226 + t227 * t354;
t339 = rSges(7,1) * t219 + rSges(7,2) * t218 - pkin(5) * t347 + t378 * t252 + t253 * t354;
t338 = qJD(2) * t305;
t337 = qJD(3) * t305;
t336 = V_base(5) * qJ(1) + V_base(1);
t332 = qJD(1) + V_base(3);
t282 = t304 * t338 + V_base(4);
t293 = qJD(2) * t307 + V_base(6);
t330 = pkin(2) * t343 - qJ(3) * t305;
t223 = -qJD(4) * t248 + t282;
t251 = -qJD(4) * t263 + t293;
t281 = -t306 * t338 + V_base(5);
t222 = -qJD(4) * t246 + t281;
t276 = pkin(1) * t304 - pkin(7) * t345;
t327 = -t276 * V_base(6) + V_base(5) * t356 + t336;
t277 = pkin(1) * t306 + pkin(7) * t346;
t326 = V_base(4) * t276 - t277 * V_base(5) + t332;
t324 = V_base(6) * t277 + V_base(2) + (-qJ(1) - t356) * V_base(4);
t275 = pkin(2) * t305 * t311 + qJ(3) * t307;
t323 = t281 * t275 + t304 * t337 + t327;
t241 = t304 * t355 + t306 * t330;
t322 = qJD(3) * t307 + t282 * t241 + t326;
t242 = -t304 * t330 + t306 * t355;
t320 = t293 * t242 - t306 * t337 + t324;
t210 = pkin(3) * t247 - pkin(8) * t246;
t245 = pkin(3) * t264 - pkin(8) * t263;
t319 = t281 * t245 + (-t210 - t241) * t293 + t323;
t211 = pkin(3) * t249 - pkin(8) * t248;
t318 = t282 * t210 + (-t211 - t242) * t281 + t322;
t192 = pkin(4) * t225 + pkin(9) * t224;
t216 = pkin(4) * t253 + pkin(9) * t252;
t317 = -t192 * t251 + t222 * t216 + t319;
t193 = pkin(4) * t227 + pkin(9) * t226;
t316 = t223 * t192 - t193 * t222 + t318;
t315 = t293 * t211 + (-t245 - t275) * t282 + t320;
t314 = t251 * t193 - t216 * t223 + t315;
t302 = Icges(2,4) * t306;
t290 = rSges(2,1) * t306 - rSges(2,2) * t304;
t289 = rSges(2,1) * t304 + rSges(2,2) * t306;
t288 = Icges(2,1) * t306 - t350;
t287 = Icges(2,1) * t304 + t302;
t286 = -Icges(2,2) * t304 + t302;
t285 = Icges(2,2) * t306 + t350;
t280 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t279 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t278 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t262 = t307 * rSges(3,3) + (rSges(3,1) * t311 + rSges(3,2) * t313) * t305;
t261 = Icges(3,5) * t307 + (Icges(3,1) * t311 + Icges(3,4) * t313) * t305;
t260 = Icges(3,6) * t307 + (Icges(3,4) * t311 + Icges(3,2) * t313) * t305;
t255 = V_base(5) * rSges(2,3) - t289 * V_base(6) + t336;
t254 = t290 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t250 = t289 * V_base(4) - t290 * V_base(5) + t332;
t240 = rSges(3,1) * t270 + rSges(3,2) * t269 + rSges(3,3) * t346;
t239 = rSges(3,1) * t268 + rSges(3,2) * t267 - rSges(3,3) * t345;
t238 = Icges(3,1) * t270 + Icges(3,4) * t269 + Icges(3,5) * t346;
t237 = Icges(3,1) * t268 + Icges(3,4) * t267 - Icges(3,5) * t345;
t236 = Icges(3,4) * t270 + Icges(3,2) * t269 + Icges(3,6) * t346;
t235 = Icges(3,4) * t268 + Icges(3,2) * t267 - Icges(3,6) * t345;
t232 = rSges(4,1) * t264 + rSges(4,2) * t263 + rSges(4,3) * t307;
t231 = Icges(4,1) * t264 + Icges(4,4) * t263 + Icges(4,5) * t307;
t230 = Icges(4,4) * t264 + Icges(4,2) * t263 + Icges(4,6) * t307;
t217 = qJD(5) * t252 + t251;
t215 = rSges(5,1) * t253 - rSges(5,2) * t252 - rSges(5,3) * t263;
t214 = Icges(5,1) * t253 - Icges(5,4) * t252 - Icges(5,5) * t263;
t213 = Icges(5,4) * t253 - Icges(5,2) * t252 - Icges(5,6) * t263;
t212 = Icges(5,5) * t253 - Icges(5,6) * t252 - Icges(5,3) * t263;
t209 = rSges(4,1) * t249 + rSges(4,2) * t248 + rSges(4,3) * t346;
t208 = rSges(4,1) * t247 + rSges(4,2) * t246 - rSges(4,3) * t345;
t207 = Icges(4,1) * t249 + Icges(4,4) * t248 + Icges(4,5) * t346;
t206 = Icges(4,1) * t247 + Icges(4,4) * t246 - Icges(4,5) * t345;
t205 = Icges(4,4) * t249 + Icges(4,2) * t248 + Icges(4,6) * t346;
t204 = Icges(4,4) * t247 + Icges(4,2) * t246 - Icges(4,6) * t345;
t195 = qJD(5) * t226 + t223;
t194 = qJD(5) * t224 + t222;
t191 = -t239 * t293 + t262 * t281 + t327;
t190 = t240 * t293 - t262 * t282 + t324;
t187 = rSges(6,1) * t219 + rSges(6,2) * t218 + rSges(6,3) * t252;
t179 = t239 * t282 - t240 * t281 + t326;
t178 = rSges(5,1) * t227 - rSges(5,2) * t226 - rSges(5,3) * t248;
t177 = rSges(5,1) * t225 - rSges(5,2) * t224 - rSges(5,3) * t246;
t176 = Icges(5,1) * t227 - Icges(5,4) * t226 - Icges(5,5) * t248;
t175 = Icges(5,1) * t225 - Icges(5,4) * t224 - Icges(5,5) * t246;
t174 = Icges(5,4) * t227 - Icges(5,2) * t226 - Icges(5,6) * t248;
t173 = Icges(5,4) * t225 - Icges(5,2) * t224 - Icges(5,6) * t246;
t172 = Icges(5,5) * t227 - Icges(5,6) * t226 - Icges(5,3) * t248;
t171 = Icges(5,5) * t225 - Icges(5,6) * t224 - Icges(5,3) * t246;
t168 = rSges(6,1) * t199 + rSges(6,2) * t198 + rSges(6,3) * t226;
t166 = rSges(6,1) * t197 + rSges(6,2) * t196 + rSges(6,3) * t224;
t150 = t232 * t281 + (-t208 - t241) * t293 + t323;
t149 = t209 * t293 + (-t232 - t275) * t282 + t320;
t148 = t208 * t282 + (-t209 - t242) * t281 + t322;
t147 = -t177 * t251 + t215 * t222 + t319;
t146 = t178 * t251 - t215 * t223 + t315;
t145 = t177 * t223 - t178 * t222 + t318;
t144 = -t166 * t217 + t187 * t194 + t317;
t143 = t168 * t217 - t187 * t195 + t314;
t142 = t166 * t195 - t168 * t194 + t316;
t141 = qJD(6) * t226 + t194 * t339 - t217 * t341 + t317;
t140 = qJD(6) * t224 - t195 * t339 + t217 * t340 + t314;
t139 = qJD(6) * t252 - t194 * t340 + t195 * t341 + t316;
t1 = m(4) * (t148 ^ 2 + t149 ^ 2 + t150 ^ 2) / 0.2e1 + m(5) * (t145 ^ 2 + t146 ^ 2 + t147 ^ 2) / 0.2e1 + m(6) * (t142 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(7) * (t139 ^ 2 + t140 ^ 2 + t141 ^ 2) / 0.2e1 + m(3) * (t179 ^ 2 + t190 ^ 2 + t191 ^ 2) / 0.2e1 + m(1) * (t278 ^ 2 + t279 ^ 2 + t280 ^ 2) / 0.2e1 + t251 * ((-t172 * t263 - t174 * t252 + t176 * t253) * t223 + (-t171 * t263 - t173 * t252 + t175 * t253) * t222 + (-t212 * t263 - t213 * t252 + t214 * t253) * t251) / 0.2e1 + m(2) * (t250 ^ 2 + t254 ^ 2 + t255 ^ 2) / 0.2e1 + t222 * ((-t172 * t246 - t174 * t224 + t176 * t225) * t223 + (-t171 * t246 - t173 * t224 + t175 * t225) * t222 + (-t212 * t246 - t213 * t224 + t214 * t225) * t251) / 0.2e1 + t223 * ((-t172 * t248 - t174 * t226 + t176 * t227) * t223 + (-t171 * t248 - t173 * t226 + t175 * t227) * t222 + (-t212 * t248 - t213 * t226 + t214 * t227) * t251) / 0.2e1 + ((t196 * t369 + t197 * t368 + t224 * t370) * t217 + (t196 * t373 + t197 * t371 + t224 * t375) * t195 + (t374 * t196 + t372 * t197 + t376 * t224) * t194) * t194 / 0.2e1 + ((t198 * t369 + t199 * t368 + t226 * t370) * t217 + (t373 * t198 + t371 * t199 + t375 * t226) * t195 + (t198 * t374 + t199 * t372 + t226 * t376) * t194) * t195 / 0.2e1 + ((t369 * t218 + t368 * t219 + t370 * t252) * t217 + (t218 * t373 + t219 * t371 + t252 * t375) * t195 + (t218 * t374 + t372 * t219 + t376 * t252) * t194) * t217 / 0.2e1 + ((t230 * t246 + t231 * t247 + t260 * t267 + t261 * t268 - t345 * t365) * t293 + (t205 * t246 + t207 * t247 + t236 * t267 + t238 * t268 - t345 * t366) * t282 + (t204 * t246 + t206 * t247 + t235 * t267 + t237 * t268 - t367 * t345) * t281) * t281 / 0.2e1 + ((t230 * t248 + t231 * t249 + t260 * t269 + t261 * t270 + t346 * t365) * t293 + (t205 * t248 + t207 * t249 + t236 * t269 + t238 * t270 + t366 * t346) * t282 + (t204 * t248 + t206 * t249 + t235 * t269 + t237 * t270 + t346 * t367) * t281) * t282 / 0.2e1 + ((t233 * t281 + t234 * t282 + t259 * t293) * t307 + ((t236 * t313 + t238 * t311) * t282 + (t235 * t313 + t237 * t311) * t281 + (t260 * t313 + t261 * t311) * t293) * t305 + (t203 * t307 + t205 * t263 + t207 * t264) * t282 + (t202 * t307 + t204 * t263 + t206 * t264) * t281 + (t229 * t307 + t230 * t263 + t231 * t264) * t293) * t293 / 0.2e1 + ((-t285 * t304 + t287 * t306 + Icges(1,4)) * V_base(5) + (-t286 * t304 + t288 * t306 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t285 * t306 + t287 * t304 + Icges(1,2)) * V_base(5) + (t286 * t306 + t288 * t304 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t306 - Icges(2,6) * t304 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t304 + Icges(2,6) * t306 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
