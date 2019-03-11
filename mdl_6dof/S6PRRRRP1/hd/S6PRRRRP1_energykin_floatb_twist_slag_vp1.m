% Calculate kinetic energy for
% S6PRRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRRP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRP1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:55:38
% EndTime: 2019-03-08 23:55:43
% DurationCPUTime: 4.59s
% Computational Cost: add. (3557->392), mult. (6178->570), div. (0->0), fcn. (7345->12), ass. (0->179)
t380 = Icges(6,1) + Icges(7,1);
t379 = Icges(6,4) + Icges(7,4);
t378 = Icges(6,5) + Icges(7,5);
t377 = Icges(6,2) + Icges(7,2);
t376 = Icges(6,6) + Icges(7,6);
t375 = Icges(6,3) + Icges(7,3);
t374 = rSges(7,3) + qJ(6);
t304 = sin(pkin(11));
t306 = cos(pkin(11));
t314 = cos(qJ(2));
t307 = cos(pkin(6));
t311 = sin(qJ(2));
t344 = t307 * t311;
t268 = t304 * t314 + t306 * t344;
t342 = qJ(3) + qJ(4);
t302 = sin(t342);
t329 = cos(t342);
t305 = sin(pkin(6));
t350 = t305 * t306;
t241 = t268 * t329 - t302 * t350;
t343 = t307 * t314;
t267 = t304 * t311 - t306 * t343;
t309 = sin(qJ(5));
t312 = cos(qJ(5));
t211 = -t241 * t309 + t267 * t312;
t353 = t267 * t309;
t212 = t241 * t312 + t353;
t328 = t305 * t329;
t240 = t268 * t302 + t306 * t328;
t372 = t211 * t376 + t212 * t378 + t240 * t375;
t270 = -t304 * t344 + t306 * t314;
t351 = t304 * t305;
t243 = t270 * t329 + t302 * t351;
t269 = t304 * t343 + t306 * t311;
t213 = -t243 * t309 + t269 * t312;
t352 = t269 * t309;
t214 = t243 * t312 + t352;
t242 = t270 * t302 - t304 * t328;
t371 = t213 * t376 + t214 * t378 + t242 * t375;
t370 = t211 * t377 + t212 * t379 + t240 * t376;
t369 = t213 * t377 + t214 * t379 + t242 * t376;
t368 = t379 * t211 + t212 * t380 + t378 * t240;
t367 = t379 * t213 + t214 * t380 + t378 * t242;
t262 = t307 * t302 + t311 * t328;
t346 = t305 * t314;
t244 = -t262 * t309 - t312 * t346;
t330 = t309 * t346;
t245 = t262 * t312 - t330;
t348 = t305 * t311;
t261 = t302 * t348 - t307 * t329;
t366 = t244 * t376 + t245 * t378 + t261 * t375;
t365 = t244 * t377 + t245 * t379 + t261 * t376;
t364 = t379 * t244 + t245 * t380 + t378 * t261;
t359 = pkin(7) * t307;
t313 = cos(qJ(3));
t358 = pkin(3) * t313;
t357 = pkin(5) * t312;
t354 = Icges(2,4) * t304;
t310 = sin(qJ(3));
t349 = t305 * t310;
t347 = t305 * t313;
t345 = t307 * t310;
t341 = rSges(7,1) * t212 + rSges(7,2) * t211 + pkin(5) * t353 + t240 * t374 + t241 * t357;
t340 = rSges(7,1) * t214 + rSges(7,2) * t213 + pkin(5) * t352 + t242 * t374 + t243 * t357;
t339 = rSges(7,1) * t245 + rSges(7,2) * t244 - pkin(5) * t330 + t261 * t374 + t262 * t357;
t338 = qJD(2) * t305;
t337 = V_base(5) * qJ(1) + V_base(1);
t333 = qJD(1) + V_base(3);
t332 = t304 * t349;
t331 = t306 * t349;
t283 = t304 * t338 + V_base(4);
t294 = qJD(2) * t307 + V_base(6);
t248 = qJD(3) * t269 + t283;
t221 = qJD(4) * t269 + t248;
t282 = -t306 * t338 + V_base(5);
t247 = qJD(3) * t267 + t282;
t277 = pkin(1) * t304 - pkin(7) * t350;
t327 = -t277 * V_base(6) + V_base(5) * t359 + t337;
t278 = pkin(1) * t306 + pkin(7) * t351;
t326 = V_base(4) * t277 - t278 * V_base(5) + t333;
t220 = qJD(4) * t267 + t247;
t256 = (-qJD(3) - qJD(4)) * t346 + t294;
t325 = V_base(6) * t278 + V_base(2) + (-qJ(1) - t359) * V_base(4);
t238 = t268 * pkin(2) + t267 * pkin(8);
t276 = (pkin(2) * t311 - pkin(8) * t314) * t305;
t324 = -t238 * t294 + t282 * t276 + t327;
t239 = t270 * pkin(2) + t269 * pkin(8);
t323 = t283 * t238 - t239 * t282 + t326;
t322 = t294 * t239 - t276 * t283 + t325;
t196 = -pkin(3) * t331 + pkin(9) * t267 + t268 * t358;
t237 = pkin(3) * t345 + (-pkin(9) * t314 + t311 * t358) * t305;
t271 = -qJD(3) * t346 + t294;
t321 = -t196 * t271 + t247 * t237 + t324;
t197 = pkin(3) * t332 + pkin(9) * t269 + t270 * t358;
t320 = t248 * t196 - t197 * t247 + t323;
t319 = t271 * t197 - t237 * t248 + t322;
t209 = pkin(4) * t241 + pkin(10) * t240;
t236 = pkin(4) * t262 + pkin(10) * t261;
t318 = -t209 * t256 + t220 * t236 + t321;
t210 = pkin(4) * t243 + pkin(10) * t242;
t317 = t221 * t209 - t210 * t220 + t320;
t316 = t256 * t210 - t221 * t236 + t319;
t301 = Icges(2,4) * t306;
t292 = rSges(2,1) * t306 - rSges(2,2) * t304;
t291 = rSges(2,1) * t304 + rSges(2,2) * t306;
t290 = Icges(2,1) * t306 - t354;
t289 = Icges(2,1) * t304 + t301;
t288 = -Icges(2,2) * t304 + t301;
t287 = Icges(2,2) * t306 + t354;
t281 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t280 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t279 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t275 = t311 * t347 + t345;
t274 = t307 * t313 - t310 * t348;
t260 = rSges(3,3) * t307 + (rSges(3,1) * t311 + rSges(3,2) * t314) * t305;
t259 = Icges(3,5) * t307 + (Icges(3,1) * t311 + Icges(3,4) * t314) * t305;
t258 = Icges(3,6) * t307 + (Icges(3,4) * t311 + Icges(3,2) * t314) * t305;
t257 = Icges(3,3) * t307 + (Icges(3,5) * t311 + Icges(3,6) * t314) * t305;
t255 = V_base(5) * rSges(2,3) - t291 * V_base(6) + t337;
t254 = t292 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t252 = t270 * t313 + t332;
t251 = -t270 * t310 + t304 * t347;
t250 = t268 * t313 - t331;
t249 = -t268 * t310 - t306 * t347;
t246 = t291 * V_base(4) - t292 * V_base(5) + t333;
t235 = rSges(4,1) * t275 + rSges(4,2) * t274 - rSges(4,3) * t346;
t234 = Icges(4,1) * t275 + Icges(4,4) * t274 - Icges(4,5) * t346;
t233 = Icges(4,4) * t275 + Icges(4,2) * t274 - Icges(4,6) * t346;
t232 = Icges(4,5) * t275 + Icges(4,6) * t274 - Icges(4,3) * t346;
t231 = rSges(3,1) * t270 - rSges(3,2) * t269 + rSges(3,3) * t351;
t230 = rSges(3,1) * t268 - rSges(3,2) * t267 - rSges(3,3) * t350;
t229 = Icges(3,1) * t270 - Icges(3,4) * t269 + Icges(3,5) * t351;
t228 = Icges(3,1) * t268 - Icges(3,4) * t267 - Icges(3,5) * t350;
t227 = Icges(3,4) * t270 - Icges(3,2) * t269 + Icges(3,6) * t351;
t226 = Icges(3,4) * t268 - Icges(3,2) * t267 - Icges(3,6) * t350;
t225 = Icges(3,5) * t270 - Icges(3,6) * t269 + Icges(3,3) * t351;
t224 = Icges(3,5) * t268 - Icges(3,6) * t267 - Icges(3,3) * t350;
t222 = qJD(5) * t261 + t256;
t218 = rSges(5,1) * t262 - rSges(5,2) * t261 - rSges(5,3) * t346;
t217 = Icges(5,1) * t262 - Icges(5,4) * t261 - Icges(5,5) * t346;
t216 = Icges(5,4) * t262 - Icges(5,2) * t261 - Icges(5,6) * t346;
t215 = Icges(5,5) * t262 - Icges(5,6) * t261 - Icges(5,3) * t346;
t207 = rSges(4,1) * t252 + rSges(4,2) * t251 + rSges(4,3) * t269;
t206 = rSges(4,1) * t250 + rSges(4,2) * t249 + rSges(4,3) * t267;
t205 = Icges(4,1) * t252 + Icges(4,4) * t251 + Icges(4,5) * t269;
t204 = Icges(4,1) * t250 + Icges(4,4) * t249 + Icges(4,5) * t267;
t203 = Icges(4,4) * t252 + Icges(4,2) * t251 + Icges(4,6) * t269;
t202 = Icges(4,4) * t250 + Icges(4,2) * t249 + Icges(4,6) * t267;
t201 = Icges(4,5) * t252 + Icges(4,6) * t251 + Icges(4,3) * t269;
t200 = Icges(4,5) * t250 + Icges(4,6) * t249 + Icges(4,3) * t267;
t199 = qJD(5) * t242 + t221;
t198 = qJD(5) * t240 + t220;
t195 = rSges(5,1) * t243 - rSges(5,2) * t242 + rSges(5,3) * t269;
t194 = rSges(5,1) * t241 - rSges(5,2) * t240 + rSges(5,3) * t267;
t193 = Icges(5,1) * t243 - Icges(5,4) * t242 + Icges(5,5) * t269;
t192 = Icges(5,1) * t241 - Icges(5,4) * t240 + Icges(5,5) * t267;
t191 = Icges(5,4) * t243 - Icges(5,2) * t242 + Icges(5,6) * t269;
t190 = Icges(5,4) * t241 - Icges(5,2) * t240 + Icges(5,6) * t267;
t189 = Icges(5,5) * t243 - Icges(5,6) * t242 + Icges(5,3) * t269;
t188 = Icges(5,5) * t241 - Icges(5,6) * t240 + Icges(5,3) * t267;
t187 = rSges(6,1) * t245 + rSges(6,2) * t244 + rSges(6,3) * t261;
t175 = -t230 * t294 + t260 * t282 + t327;
t174 = t231 * t294 - t260 * t283 + t325;
t171 = t230 * t283 - t231 * t282 + t326;
t170 = rSges(6,1) * t214 + rSges(6,2) * t213 + rSges(6,3) * t242;
t168 = rSges(6,1) * t212 + rSges(6,2) * t211 + rSges(6,3) * t240;
t152 = -t206 * t271 + t235 * t247 + t324;
t151 = t207 * t271 - t235 * t248 + t322;
t150 = t206 * t248 - t207 * t247 + t323;
t149 = -t194 * t256 + t218 * t220 + t321;
t148 = t195 * t256 - t218 * t221 + t319;
t147 = t194 * t221 - t195 * t220 + t320;
t146 = -t168 * t222 + t187 * t198 + t318;
t145 = t170 * t222 - t187 * t199 + t316;
t144 = t168 * t199 - t170 * t198 + t317;
t143 = qJD(6) * t242 + t198 * t339 - t222 * t341 + t318;
t142 = qJD(6) * t240 - t199 * t339 + t222 * t340 + t316;
t141 = qJD(6) * t261 - t198 * t340 + t199 * t341 + t317;
t1 = m(6) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(5) * (t147 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + m(4) * (t150 ^ 2 + t151 ^ 2 + t152 ^ 2) / 0.2e1 + m(3) * (t171 ^ 2 + t174 ^ 2 + t175 ^ 2) / 0.2e1 + m(7) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + t294 * ((t224 * t282 + t225 * t283 + t257 * t294) * t307 + ((t227 * t314 + t229 * t311) * t283 + (t226 * t314 + t228 * t311) * t282 + (t258 * t314 + t259 * t311) * t294) * t305) / 0.2e1 + m(2) * (t246 ^ 2 + t254 ^ 2 + t255 ^ 2) / 0.2e1 + t220 * ((t189 * t267 - t191 * t240 + t193 * t241) * t221 + (t267 * t188 - t240 * t190 + t241 * t192) * t220 + (t215 * t267 - t216 * t240 + t217 * t241) * t256) / 0.2e1 + t221 * ((t269 * t189 - t242 * t191 + t243 * t193) * t221 + (t188 * t269 - t190 * t242 + t192 * t243) * t220 + (t215 * t269 - t216 * t242 + t217 * t243) * t256) / 0.2e1 + t247 * ((t201 * t267 + t203 * t249 + t205 * t250) * t248 + (t200 * t267 + t202 * t249 + t204 * t250) * t247 + (t232 * t267 + t233 * t249 + t234 * t250) * t271) / 0.2e1 + t248 * ((t201 * t269 + t203 * t251 + t205 * t252) * t248 + (t200 * t269 + t202 * t251 + t204 * t252) * t247 + (t232 * t269 + t233 * t251 + t234 * t252) * t271) / 0.2e1 + m(1) * (t279 ^ 2 + t280 ^ 2 + t281 ^ 2) / 0.2e1 + t271 * ((-t201 * t346 + t203 * t274 + t205 * t275) * t248 + (-t200 * t346 + t202 * t274 + t204 * t275) * t247 + (-t232 * t346 + t233 * t274 + t234 * t275) * t271) / 0.2e1 + t256 * ((-t189 * t346 - t191 * t261 + t193 * t262) * t221 + (-t188 * t346 - t190 * t261 + t192 * t262) * t220 + (-t215 * t346 - t216 * t261 + t217 * t262) * t256) / 0.2e1 + t282 * ((-t225 * t350 - t227 * t267 + t229 * t268) * t283 + (-t224 * t350 - t226 * t267 + t228 * t268) * t282 + (-t257 * t350 - t258 * t267 + t259 * t268) * t294) / 0.2e1 + t283 * ((t225 * t351 - t227 * t269 + t229 * t270) * t283 + (t224 * t351 - t226 * t269 + t228 * t270) * t282 + (t257 * t351 - t258 * t269 + t259 * t270) * t294) / 0.2e1 + ((t211 * t365 + t212 * t364 + t240 * t366) * t222 + (t211 * t369 + t212 * t367 + t240 * t371) * t199 + (t370 * t211 + t368 * t212 + t372 * t240) * t198) * t198 / 0.2e1 + ((t213 * t365 + t214 * t364 + t242 * t366) * t222 + (t369 * t213 + t367 * t214 + t371 * t242) * t199 + (t213 * t370 + t214 * t368 + t242 * t372) * t198) * t199 / 0.2e1 + ((t365 * t244 + t364 * t245 + t366 * t261) * t222 + (t244 * t369 + t245 * t367 + t261 * t371) * t199 + (t244 * t370 + t245 * t368 + t261 * t372) * t198) * t222 / 0.2e1 + ((-t287 * t304 + t289 * t306 + Icges(1,4)) * V_base(5) + (-t288 * t304 + t290 * t306 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t287 * t306 + t289 * t304 + Icges(1,2)) * V_base(5) + (t288 * t306 + t290 * t304 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t306 - Icges(2,6) * t304 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t304 + Icges(2,6) * t306 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
