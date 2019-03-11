% Calculate kinetic energy for
% S6PRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
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
% Datum: 2019-03-08 22:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRPP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPP1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPP1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:42:45
% EndTime: 2019-03-08 22:42:49
% DurationCPUTime: 3.75s
% Computational Cost: add. (3224->380), mult. (6709->530), div. (0->0), fcn. (8131->12), ass. (0->168)
t382 = Icges(6,1) + Icges(7,1);
t381 = -Icges(6,4) + Icges(7,5);
t380 = Icges(7,4) + Icges(6,5);
t379 = Icges(6,2) + Icges(7,3);
t378 = -Icges(6,6) + Icges(7,6);
t377 = Icges(7,2) + Icges(5,3) + Icges(6,3);
t376 = rSges(7,1) + pkin(5);
t375 = rSges(7,3) + qJ(6);
t313 = sin(pkin(10));
t315 = cos(pkin(10));
t322 = cos(qJ(2));
t316 = cos(pkin(6));
t320 = sin(qJ(2));
t349 = t316 * t320;
t278 = t313 * t322 + t315 * t349;
t314 = sin(pkin(6));
t319 = sin(qJ(3));
t351 = t314 * t319;
t360 = cos(qJ(3));
t260 = t278 * t360 - t315 * t351;
t348 = t316 * t322;
t277 = t313 * t320 - t315 * t348;
t343 = qJ(4) + pkin(11);
t310 = sin(t343);
t335 = cos(t343);
t226 = t260 * t310 - t277 * t335;
t227 = t260 * t335 + t277 * t310;
t336 = t314 * t360;
t259 = t278 * t319 + t315 * t336;
t373 = t226 * t379 + t227 * t381 + t259 * t378;
t280 = -t313 * t349 + t315 * t322;
t262 = t280 * t360 + t313 * t351;
t279 = t313 * t348 + t315 * t320;
t228 = t262 * t310 - t279 * t335;
t229 = t262 * t335 + t279 * t310;
t261 = t280 * t319 - t313 * t336;
t372 = t228 * t379 + t229 * t381 + t261 * t378;
t371 = t226 * t381 + t227 * t382 + t259 * t380;
t370 = t228 * t381 + t229 * t382 + t261 * t380;
t285 = t316 * t319 + t320 * t336;
t350 = t314 * t322;
t252 = t285 * t310 + t335 * t350;
t253 = t285 * t335 - t310 * t350;
t284 = -t316 * t360 + t320 * t351;
t369 = t252 * t379 + t253 * t381 + t284 * t378;
t368 = t252 * t381 + t253 * t382 + t284 * t380;
t318 = sin(qJ(4));
t321 = cos(qJ(4));
t230 = -t260 * t318 + t277 * t321;
t355 = t277 * t318;
t231 = t260 * t321 + t355;
t366 = Icges(5,5) * t231 + Icges(5,6) * t230 + t226 * t378 + t227 * t380 + t259 * t377;
t232 = -t262 * t318 + t279 * t321;
t354 = t279 * t318;
t233 = t262 * t321 + t354;
t365 = Icges(5,5) * t233 + Icges(5,6) * t232 + t228 * t378 + t229 * t380 + t261 * t377;
t263 = -t285 * t318 - t321 * t350;
t337 = t318 * t350;
t264 = t285 * t321 - t337;
t364 = Icges(5,5) * t264 + Icges(5,6) * t263 + t252 * t378 + t253 * t380 + t284 * t377;
t359 = pkin(7) * t316;
t358 = pkin(4) * t321;
t356 = Icges(2,4) * t313;
t353 = t313 * t314;
t352 = t314 * t315;
t347 = rSges(7,2) * t259 + t226 * t375 + t227 * t376;
t346 = rSges(7,2) * t261 + t228 * t375 + t229 * t376;
t345 = rSges(7,2) * t284 + t252 * t375 + t253 * t376;
t344 = qJD(2) * t314;
t342 = V_base(5) * qJ(1) + V_base(1);
t338 = qJD(1) + V_base(3);
t293 = t313 * t344 + V_base(4);
t304 = qJD(2) * t316 + V_base(6);
t258 = qJD(3) * t279 + t293;
t292 = -t315 * t344 + V_base(5);
t257 = qJD(3) * t277 + t292;
t281 = -qJD(3) * t350 + t304;
t287 = pkin(1) * t313 - pkin(7) * t352;
t334 = -t287 * V_base(6) + t359 * V_base(5) + t342;
t288 = pkin(1) * t315 + pkin(7) * t353;
t333 = t287 * V_base(4) - t288 * V_base(5) + t338;
t332 = V_base(6) * t288 + V_base(2) + (-qJ(1) - t359) * V_base(4);
t248 = pkin(2) * t278 + pkin(8) * t277;
t286 = (pkin(2) * t320 - pkin(8) * t322) * t314;
t331 = -t248 * t304 + t286 * t292 + t334;
t249 = pkin(2) * t280 + pkin(8) * t279;
t330 = t248 * t293 - t249 * t292 + t333;
t329 = t249 * t304 - t286 * t293 + t332;
t222 = pkin(3) * t260 + pkin(9) * t259;
t250 = pkin(3) * t285 + t284 * pkin(9);
t328 = -t222 * t281 + t250 * t257 + t331;
t223 = pkin(3) * t262 + pkin(9) * t261;
t327 = t222 * t258 - t223 * t257 + t330;
t326 = t223 * t281 - t250 * t258 + t329;
t206 = -pkin(4) * t337 + qJ(5) * t284 + t285 * t358;
t224 = qJD(4) * t259 + t257;
t325 = qJD(5) * t261 + t206 * t224 + t328;
t165 = pkin(4) * t355 + qJ(5) * t259 + t260 * t358;
t225 = qJD(4) * t261 + t258;
t324 = qJD(5) * t284 + t165 * t225 + t327;
t166 = pkin(4) * t354 + qJ(5) * t261 + t262 * t358;
t251 = qJD(4) * t284 + t281;
t323 = qJD(5) * t259 + t166 * t251 + t326;
t311 = Icges(2,4) * t315;
t301 = rSges(2,1) * t315 - rSges(2,2) * t313;
t300 = rSges(2,1) * t313 + rSges(2,2) * t315;
t299 = Icges(2,1) * t315 - t356;
t298 = Icges(2,1) * t313 + t311;
t297 = -Icges(2,2) * t313 + t311;
t296 = Icges(2,2) * t315 + t356;
t291 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t290 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t289 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t273 = t316 * rSges(3,3) + (rSges(3,1) * t320 + rSges(3,2) * t322) * t314;
t272 = Icges(3,5) * t316 + (Icges(3,1) * t320 + Icges(3,4) * t322) * t314;
t271 = Icges(3,6) * t316 + (Icges(3,4) * t320 + Icges(3,2) * t322) * t314;
t270 = Icges(3,3) * t316 + (Icges(3,5) * t320 + Icges(3,6) * t322) * t314;
t267 = V_base(5) * rSges(2,3) - t300 * V_base(6) + t342;
t266 = t301 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t256 = t300 * V_base(4) - t301 * V_base(5) + t338;
t247 = rSges(4,1) * t285 - rSges(4,2) * t284 - rSges(4,3) * t350;
t246 = Icges(4,1) * t285 - Icges(4,4) * t284 - Icges(4,5) * t350;
t245 = Icges(4,4) * t285 - Icges(4,2) * t284 - Icges(4,6) * t350;
t244 = Icges(4,5) * t285 - Icges(4,6) * t284 - Icges(4,3) * t350;
t243 = rSges(3,1) * t280 - rSges(3,2) * t279 + rSges(3,3) * t353;
t242 = rSges(3,1) * t278 - rSges(3,2) * t277 - rSges(3,3) * t352;
t241 = Icges(3,1) * t280 - Icges(3,4) * t279 + Icges(3,5) * t353;
t240 = Icges(3,1) * t278 - Icges(3,4) * t277 - Icges(3,5) * t352;
t239 = Icges(3,4) * t280 - Icges(3,2) * t279 + Icges(3,6) * t353;
t238 = Icges(3,4) * t278 - Icges(3,2) * t277 - Icges(3,6) * t352;
t237 = Icges(3,5) * t280 - Icges(3,6) * t279 + Icges(3,3) * t353;
t236 = Icges(3,5) * t278 - Icges(3,6) * t277 - Icges(3,3) * t352;
t219 = rSges(5,1) * t264 + rSges(5,2) * t263 + rSges(5,3) * t284;
t217 = Icges(5,1) * t264 + Icges(5,4) * t263 + Icges(5,5) * t284;
t216 = Icges(5,4) * t264 + Icges(5,2) * t263 + Icges(5,6) * t284;
t214 = rSges(4,1) * t262 - rSges(4,2) * t261 + rSges(4,3) * t279;
t213 = rSges(4,1) * t260 - rSges(4,2) * t259 + rSges(4,3) * t277;
t212 = Icges(4,1) * t262 - Icges(4,4) * t261 + Icges(4,5) * t279;
t211 = Icges(4,1) * t260 - Icges(4,4) * t259 + Icges(4,5) * t277;
t210 = Icges(4,4) * t262 - Icges(4,2) * t261 + Icges(4,6) * t279;
t209 = Icges(4,4) * t260 - Icges(4,2) * t259 + Icges(4,6) * t277;
t208 = Icges(4,5) * t262 - Icges(4,6) * t261 + Icges(4,3) * t279;
t207 = Icges(4,5) * t260 - Icges(4,6) * t259 + Icges(4,3) * t277;
t205 = rSges(6,1) * t253 - rSges(6,2) * t252 + rSges(6,3) * t284;
t196 = -t242 * t304 + t273 * t292 + t334;
t195 = t243 * t304 - t273 * t293 + t332;
t192 = rSges(5,1) * t233 + rSges(5,2) * t232 + rSges(5,3) * t261;
t191 = rSges(5,1) * t231 + rSges(5,2) * t230 + rSges(5,3) * t259;
t190 = Icges(5,1) * t233 + Icges(5,4) * t232 + Icges(5,5) * t261;
t189 = Icges(5,1) * t231 + Icges(5,4) * t230 + Icges(5,5) * t259;
t188 = Icges(5,4) * t233 + Icges(5,2) * t232 + Icges(5,6) * t261;
t187 = Icges(5,4) * t231 + Icges(5,2) * t230 + Icges(5,6) * t259;
t184 = t242 * t293 - t243 * t292 + t333;
t182 = rSges(6,1) * t229 - rSges(6,2) * t228 + rSges(6,3) * t261;
t180 = rSges(6,1) * t227 - rSges(6,2) * t226 + rSges(6,3) * t259;
t162 = -t213 * t281 + t247 * t257 + t331;
t161 = t214 * t281 - t247 * t258 + t329;
t160 = t213 * t258 - t214 * t257 + t330;
t159 = -t191 * t251 + t219 * t224 + t328;
t158 = t192 * t251 - t219 * t225 + t326;
t157 = t191 * t225 - t192 * t224 + t327;
t156 = t205 * t224 + (-t165 - t180) * t251 + t325;
t155 = t182 * t251 + (-t205 - t206) * t225 + t323;
t154 = t180 * t225 + (-t166 - t182) * t224 + t324;
t153 = qJD(6) * t228 + t345 * t224 + (-t165 - t347) * t251 + t325;
t152 = qJD(6) * t226 + t346 * t251 + (-t206 - t345) * t225 + t323;
t151 = qJD(6) * t252 + t347 * t225 + (-t166 - t346) * t224 + t324;
t1 = m(3) * (t184 ^ 2 + t195 ^ 2 + t196 ^ 2) / 0.2e1 + m(1) * (t289 ^ 2 + t290 ^ 2 + t291 ^ 2) / 0.2e1 + t258 * ((t208 * t279 - t210 * t261 + t212 * t262) * t258 + (t207 * t279 - t209 * t261 + t211 * t262) * t257 + (t244 * t279 - t245 * t261 + t246 * t262) * t281) / 0.2e1 + t257 * ((t208 * t277 - t210 * t259 + t212 * t260) * t258 + (t207 * t277 - t209 * t259 + t211 * t260) * t257 + (t244 * t277 - t245 * t259 + t246 * t260) * t281) / 0.2e1 + m(2) * (t256 ^ 2 + t266 ^ 2 + t267 ^ 2) / 0.2e1 + t304 * ((t236 * t292 + t237 * t293 + t270 * t304) * t316 + ((t239 * t322 + t241 * t320) * t293 + (t238 * t322 + t240 * t320) * t292 + (t271 * t322 + t272 * t320) * t304) * t314) / 0.2e1 + t292 * ((-t237 * t352 - t239 * t277 + t241 * t278) * t293 + (-t236 * t352 - t238 * t277 + t240 * t278) * t292 + (-t270 * t352 - t271 * t277 + t272 * t278) * t304) / 0.2e1 + t293 * ((t237 * t353 - t239 * t279 + t241 * t280) * t293 + (t236 * t353 - t238 * t279 + t240 * t280) * t292 + (t270 * t353 - t271 * t279 + t272 * t280) * t304) / 0.2e1 + t281 * ((-t208 * t350 - t210 * t284 + t212 * t285) * t258 + (-t207 * t350 - t209 * t284 + t211 * t285) * t257 + (-t244 * t350 - t284 * t245 + t285 * t246) * t281) / 0.2e1 + m(5) * (t157 ^ 2 + t158 ^ 2 + t159 ^ 2) / 0.2e1 + m(4) * (t160 ^ 2 + t161 ^ 2 + t162 ^ 2) / 0.2e1 + m(6) * (t154 ^ 2 + t155 ^ 2 + t156 ^ 2) / 0.2e1 + m(7) * (t151 ^ 2 + t152 ^ 2 + t153 ^ 2) / 0.2e1 + ((-t296 * t313 + t298 * t315 + Icges(1,4)) * V_base(5) + (-t297 * t313 + t299 * t315 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t296 * t315 + t298 * t313 + Icges(1,2)) * V_base(5) + (t297 * t315 + t299 * t313 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t216 * t230 + t217 * t231 + t226 * t369 + t227 * t368 + t259 * t364) * t251 + (t188 * t230 + t190 * t231 + t226 * t372 + t227 * t370 + t259 * t365) * t225 + (t187 * t230 + t189 * t231 + t373 * t226 + t371 * t227 + t366 * t259) * t224) * t224 / 0.2e1 + ((t216 * t232 + t217 * t233 + t228 * t369 + t229 * t368 + t261 * t364) * t251 + (t188 * t232 + t190 * t233 + t372 * t228 + t370 * t229 + t365 * t261) * t225 + (t187 * t232 + t189 * t233 + t228 * t373 + t229 * t371 + t261 * t366) * t224) * t225 / 0.2e1 + ((t263 * t216 + t264 * t217 + t369 * t252 + t368 * t253 + t364 * t284) * t251 + (t188 * t263 + t190 * t264 + t252 * t372 + t253 * t370 + t284 * t365) * t225 + (t187 * t263 + t189 * t264 + t252 * t373 + t253 * t371 + t284 * t366) * t224) * t251 / 0.2e1 + ((Icges(2,5) * t313 + Icges(2,6) * t315 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t315 - Icges(2,6) * t313 + Icges(1,5)) * V_base(4) + (Icges(2,3) / 0.2e1 + Icges(1,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
