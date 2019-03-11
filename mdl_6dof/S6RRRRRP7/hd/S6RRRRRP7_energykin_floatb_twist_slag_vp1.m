% Calculate kinetic energy for
% S6RRRRRP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP7_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRP7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP7_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP7_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP7_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:34:23
% EndTime: 2019-03-10 01:34:28
% DurationCPUTime: 4.57s
% Computational Cost: add. (3617->392), mult. (6178->571), div. (0->0), fcn. (7345->12), ass. (0->180)
t378 = Icges(6,1) + Icges(7,1);
t377 = Icges(6,4) + Icges(7,4);
t376 = Icges(6,5) + Icges(7,5);
t375 = Icges(6,2) + Icges(7,2);
t374 = Icges(6,6) + Icges(7,6);
t373 = Icges(6,3) + Icges(7,3);
t372 = rSges(7,3) + qJ(6);
t306 = cos(pkin(6));
t311 = sin(qJ(1));
t314 = cos(qJ(2));
t343 = t311 * t314;
t310 = sin(qJ(2));
t315 = cos(qJ(1));
t344 = t310 * t315;
t272 = t306 * t344 + t343;
t341 = qJ(3) + qJ(4);
t302 = sin(t341);
t330 = cos(t341);
t305 = sin(pkin(6));
t347 = t305 * t315;
t243 = t272 * t330 - t302 * t347;
t342 = t314 * t315;
t345 = t310 * t311;
t271 = -t306 * t342 + t345;
t308 = sin(qJ(5));
t312 = cos(qJ(5));
t211 = -t243 * t308 + t271 * t312;
t353 = t271 * t308;
t212 = t243 * t312 + t353;
t329 = t305 * t330;
t242 = t272 * t302 + t315 * t329;
t371 = t211 * t374 + t212 * t376 + t242 * t373;
t274 = -t306 * t345 + t342;
t350 = t305 * t311;
t245 = t274 * t330 + t302 * t350;
t273 = t306 * t343 + t344;
t213 = -t245 * t308 + t273 * t312;
t352 = t273 * t308;
t214 = t245 * t312 + t352;
t244 = t274 * t302 - t311 * t329;
t370 = t213 * t374 + t214 * t376 + t244 * t373;
t369 = t211 * t375 + t212 * t377 + t242 * t374;
t368 = t213 * t375 + t214 * t377 + t244 * t374;
t367 = t377 * t211 + t212 * t378 + t376 * t242;
t366 = t377 * t213 + t214 * t378 + t376 * t244;
t262 = t306 * t302 + t310 * t329;
t348 = t305 * t314;
t240 = -t262 * t308 - t312 * t348;
t333 = t308 * t348;
t241 = t262 * t312 - t333;
t351 = t305 * t310;
t261 = t302 * t351 - t306 * t330;
t365 = t240 * t374 + t241 * t376 + t261 * t373;
t364 = t240 * t375 + t241 * t377 + t261 * t374;
t363 = t377 * t240 + t241 * t378 + t376 * t261;
t359 = pkin(8) * t306;
t313 = cos(qJ(3));
t358 = pkin(3) * t313;
t357 = pkin(5) * t312;
t354 = Icges(2,4) * t311;
t349 = t305 * t313;
t309 = sin(qJ(3));
t346 = t306 * t309;
t340 = rSges(7,1) * t212 + rSges(7,2) * t211 + pkin(5) * t353 + t242 * t372 + t243 * t357;
t339 = rSges(7,1) * t214 + rSges(7,2) * t213 + pkin(5) * t352 + t244 * t372 + t245 * t357;
t338 = rSges(7,1) * t241 + rSges(7,2) * t240 - pkin(5) * t333 + t261 * t372 + t262 * t357;
t337 = qJD(2) * t305;
t336 = V_base(5) * pkin(7) + V_base(1);
t332 = t309 * t350;
t331 = t309 * t347;
t283 = t311 * t337 + V_base(4);
t301 = V_base(6) + qJD(1);
t247 = qJD(3) * t273 + t283;
t285 = qJD(2) * t306 + t301;
t222 = qJD(4) * t273 + t247;
t282 = -t315 * t337 + V_base(5);
t277 = pkin(1) * t311 - pkin(8) * t347;
t328 = -t277 * t301 + V_base(5) * t359 + t336;
t278 = pkin(1) * t315 + pkin(8) * t350;
t327 = V_base(4) * t277 - t278 * V_base(5) + V_base(3);
t246 = qJD(3) * t271 + t282;
t221 = qJD(4) * t271 + t246;
t326 = t301 * t278 + V_base(2) + (-pkin(7) - t359) * V_base(4);
t256 = (-qJD(3) - qJD(4)) * t348 + t285;
t238 = t272 * pkin(2) + t271 * pkin(9);
t276 = (pkin(2) * t310 - pkin(9) * t314) * t305;
t325 = -t238 * t285 + t282 * t276 + t328;
t239 = t274 * pkin(2) + t273 * pkin(9);
t324 = t283 * t238 - t239 * t282 + t327;
t323 = t285 * t239 - t276 * t283 + t326;
t196 = -pkin(3) * t331 + pkin(10) * t271 + t272 * t358;
t235 = pkin(3) * t346 + (-pkin(10) * t314 + t310 * t358) * t305;
t267 = -qJD(3) * t348 + t285;
t322 = -t196 * t267 + t246 * t235 + t325;
t197 = pkin(3) * t332 + pkin(10) * t273 + t274 * t358;
t321 = t247 * t196 - t197 * t246 + t324;
t320 = t267 * t197 - t235 * t247 + t323;
t209 = pkin(4) * t243 + pkin(11) * t242;
t228 = pkin(4) * t262 + pkin(11) * t261;
t319 = -t209 * t256 + t221 * t228 + t322;
t210 = pkin(4) * t245 + pkin(11) * t244;
t318 = t222 * t209 - t210 * t221 + t321;
t317 = t256 * t210 - t222 * t228 + t320;
t303 = Icges(2,4) * t315;
t293 = rSges(2,1) * t315 - rSges(2,2) * t311;
t292 = rSges(2,1) * t311 + rSges(2,2) * t315;
t291 = Icges(2,1) * t315 - t354;
t290 = Icges(2,1) * t311 + t303;
t289 = -Icges(2,2) * t311 + t303;
t288 = Icges(2,2) * t315 + t354;
t281 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t280 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t279 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t270 = t310 * t349 + t346;
t269 = t306 * t313 - t309 * t351;
t260 = rSges(3,3) * t306 + (rSges(3,1) * t310 + rSges(3,2) * t314) * t305;
t259 = Icges(3,5) * t306 + (Icges(3,1) * t310 + Icges(3,4) * t314) * t305;
t258 = Icges(3,6) * t306 + (Icges(3,4) * t310 + Icges(3,2) * t314) * t305;
t257 = Icges(3,3) * t306 + (Icges(3,5) * t310 + Icges(3,6) * t314) * t305;
t255 = V_base(5) * rSges(2,3) - t292 * t301 + t336;
t254 = t293 * t301 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t252 = t292 * V_base(4) - t293 * V_base(5) + V_base(3);
t251 = t274 * t313 + t332;
t250 = -t274 * t309 + t311 * t349;
t249 = t272 * t313 - t331;
t248 = -t272 * t309 - t313 * t347;
t237 = rSges(3,1) * t274 - rSges(3,2) * t273 + rSges(3,3) * t350;
t236 = rSges(3,1) * t272 - rSges(3,2) * t271 - rSges(3,3) * t347;
t234 = Icges(3,1) * t274 - Icges(3,4) * t273 + Icges(3,5) * t350;
t233 = Icges(3,1) * t272 - Icges(3,4) * t271 - Icges(3,5) * t347;
t232 = Icges(3,4) * t274 - Icges(3,2) * t273 + Icges(3,6) * t350;
t231 = Icges(3,4) * t272 - Icges(3,2) * t271 - Icges(3,6) * t347;
t230 = Icges(3,5) * t274 - Icges(3,6) * t273 + Icges(3,3) * t350;
t229 = Icges(3,5) * t272 - Icges(3,6) * t271 - Icges(3,3) * t347;
t227 = rSges(4,1) * t270 + rSges(4,2) * t269 - rSges(4,3) * t348;
t226 = Icges(4,1) * t270 + Icges(4,4) * t269 - Icges(4,5) * t348;
t225 = Icges(4,4) * t270 + Icges(4,2) * t269 - Icges(4,6) * t348;
t224 = Icges(4,5) * t270 + Icges(4,6) * t269 - Icges(4,3) * t348;
t219 = qJD(5) * t261 + t256;
t218 = rSges(5,1) * t262 - rSges(5,2) * t261 - rSges(5,3) * t348;
t217 = Icges(5,1) * t262 - Icges(5,4) * t261 - Icges(5,5) * t348;
t216 = Icges(5,4) * t262 - Icges(5,2) * t261 - Icges(5,6) * t348;
t215 = Icges(5,5) * t262 - Icges(5,6) * t261 - Icges(5,3) * t348;
t207 = rSges(4,1) * t251 + rSges(4,2) * t250 + rSges(4,3) * t273;
t206 = rSges(4,1) * t249 + rSges(4,2) * t248 + rSges(4,3) * t271;
t205 = Icges(4,1) * t251 + Icges(4,4) * t250 + Icges(4,5) * t273;
t204 = Icges(4,1) * t249 + Icges(4,4) * t248 + Icges(4,5) * t271;
t203 = Icges(4,4) * t251 + Icges(4,2) * t250 + Icges(4,6) * t273;
t202 = Icges(4,4) * t249 + Icges(4,2) * t248 + Icges(4,6) * t271;
t201 = Icges(4,5) * t251 + Icges(4,6) * t250 + Icges(4,3) * t273;
t200 = Icges(4,5) * t249 + Icges(4,6) * t248 + Icges(4,3) * t271;
t199 = qJD(5) * t244 + t222;
t198 = qJD(5) * t242 + t221;
t195 = rSges(5,1) * t245 - rSges(5,2) * t244 + rSges(5,3) * t273;
t194 = rSges(5,1) * t243 - rSges(5,2) * t242 + rSges(5,3) * t271;
t193 = Icges(5,1) * t245 - Icges(5,4) * t244 + Icges(5,5) * t273;
t192 = Icges(5,1) * t243 - Icges(5,4) * t242 + Icges(5,5) * t271;
t191 = Icges(5,4) * t245 - Icges(5,2) * t244 + Icges(5,6) * t273;
t190 = Icges(5,4) * t243 - Icges(5,2) * t242 + Icges(5,6) * t271;
t189 = Icges(5,5) * t245 - Icges(5,6) * t244 + Icges(5,3) * t273;
t188 = Icges(5,5) * t243 - Icges(5,6) * t242 + Icges(5,3) * t271;
t186 = rSges(6,1) * t241 + rSges(6,2) * t240 + rSges(6,3) * t261;
t173 = -t236 * t285 + t260 * t282 + t328;
t172 = t237 * t285 - t260 * t283 + t326;
t171 = t236 * t283 - t237 * t282 + t327;
t170 = rSges(6,1) * t214 + rSges(6,2) * t213 + rSges(6,3) * t244;
t168 = rSges(6,1) * t212 + rSges(6,2) * t211 + rSges(6,3) * t242;
t152 = -t206 * t267 + t227 * t246 + t325;
t151 = t207 * t267 - t227 * t247 + t323;
t150 = t206 * t247 - t207 * t246 + t324;
t149 = -t194 * t256 + t218 * t221 + t322;
t148 = t195 * t256 - t218 * t222 + t320;
t147 = t194 * t222 - t195 * t221 + t321;
t146 = -t168 * t219 + t186 * t198 + t319;
t145 = t170 * t219 - t186 * t199 + t317;
t144 = t168 * t199 - t170 * t198 + t318;
t143 = qJD(6) * t244 + t198 * t338 - t219 * t340 + t319;
t142 = qJD(6) * t242 - t199 * t338 + t219 * t339 + t317;
t141 = qJD(6) * t261 - t198 * t339 + t199 * t340 + t318;
t1 = m(7) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(6) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(5) * (t147 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + m(4) * (t150 ^ 2 + t151 ^ 2 + t152 ^ 2) / 0.2e1 + m(3) * (t171 ^ 2 + t172 ^ 2 + t173 ^ 2) / 0.2e1 + t285 * ((t229 * t282 + t230 * t283 + t257 * t285) * t306 + ((t232 * t314 + t234 * t310) * t283 + (t231 * t314 + t233 * t310) * t282 + (t258 * t314 + t259 * t310) * t285) * t305) / 0.2e1 + t282 * ((-t230 * t347 - t232 * t271 + t234 * t272) * t283 + (-t229 * t347 - t231 * t271 + t233 * t272) * t282 + (-t257 * t347 - t258 * t271 + t259 * t272) * t285) / 0.2e1 + t267 * ((-t201 * t348 + t203 * t269 + t205 * t270) * t247 + (-t200 * t348 + t202 * t269 + t204 * t270) * t246 + (-t224 * t348 + t225 * t269 + t226 * t270) * t267) / 0.2e1 + t256 * ((-t189 * t348 - t191 * t261 + t193 * t262) * t222 + (-t188 * t348 - t190 * t261 + t192 * t262) * t221 + (-t215 * t348 - t216 * t261 + t217 * t262) * t256) / 0.2e1 + t283 * ((t230 * t350 - t232 * t273 + t234 * t274) * t283 + (t229 * t350 - t231 * t273 + t233 * t274) * t282 + (t257 * t350 - t258 * t273 + t259 * t274) * t285) / 0.2e1 + m(2) * (t252 ^ 2 + t254 ^ 2 + t255 ^ 2) / 0.2e1 + t221 * ((t189 * t271 - t191 * t242 + t193 * t243) * t222 + (t271 * t188 - t242 * t190 + t243 * t192) * t221 + (t215 * t271 - t216 * t242 + t217 * t243) * t256) / 0.2e1 + t246 * ((t201 * t271 + t203 * t248 + t205 * t249) * t247 + (t200 * t271 + t202 * t248 + t204 * t249) * t246 + (t224 * t271 + t225 * t248 + t226 * t249) * t267) / 0.2e1 + t222 * ((t273 * t189 - t244 * t191 + t245 * t193) * t222 + (t188 * t273 - t190 * t244 + t192 * t245) * t221 + (t215 * t273 - t216 * t244 + t217 * t245) * t256) / 0.2e1 + t247 * ((t201 * t273 + t203 * t250 + t205 * t251) * t247 + (t200 * t273 + t202 * t250 + t204 * t251) * t246 + (t224 * t273 + t225 * t250 + t226 * t251) * t267) / 0.2e1 + m(1) * (t279 ^ 2 + t280 ^ 2 + t281 ^ 2) / 0.2e1 + ((t211 * t364 + t212 * t363 + t242 * t365) * t219 + (t211 * t368 + t212 * t366 + t242 * t370) * t199 + (t369 * t211 + t367 * t212 + t371 * t242) * t198) * t198 / 0.2e1 + ((t213 * t364 + t214 * t363 + t244 * t365) * t219 + (t368 * t213 + t366 * t214 + t370 * t244) * t199 + (t369 * t213 + t367 * t214 + t244 * t371) * t198) * t199 / 0.2e1 + ((t364 * t240 + t363 * t241 + t365 * t261) * t219 + (t240 * t368 + t241 * t366 + t261 * t370) * t199 + (t369 * t240 + t367 * t241 + t261 * t371) * t198) * t219 / 0.2e1 + ((-t288 * t311 + t290 * t315 + Icges(1,4)) * V_base(5) + (-t289 * t311 + t291 * t315 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t288 * t315 + t290 * t311 + Icges(1,2)) * V_base(5) + (t289 * t315 + t291 * t311 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t311 + Icges(2,6) * t315) * V_base(5) + (Icges(2,5) * t315 - Icges(2,6) * t311) * V_base(4) + Icges(2,3) * t301 / 0.2e1) * t301;
T  = t1;
