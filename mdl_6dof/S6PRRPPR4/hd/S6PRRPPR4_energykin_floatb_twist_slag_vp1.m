% Calculate kinetic energy for
% S6PRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPPR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRPPR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPPR4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:12:25
% EndTime: 2019-03-08 21:12:30
% DurationCPUTime: 5.08s
% Computational Cost: add. (3021->386), mult. (7256->535), div. (0->0), fcn. (8949->12), ass. (0->167)
t365 = Icges(5,1) + Icges(6,1);
t364 = -Icges(5,4) + Icges(6,5);
t363 = Icges(6,4) + Icges(5,5);
t362 = Icges(5,2) + Icges(6,3);
t361 = -Icges(5,6) + Icges(6,6);
t360 = Icges(4,2) + Icges(6,2) + Icges(5,3);
t305 = sin(pkin(10));
t307 = cos(pkin(10));
t312 = cos(qJ(2));
t310 = sin(qJ(2));
t344 = cos(pkin(6));
t326 = t310 * t344;
t271 = t305 * t312 + t307 * t326;
t306 = sin(pkin(6));
t309 = sin(qJ(3));
t339 = t306 * t309;
t345 = cos(qJ(3));
t255 = t271 * t345 - t307 * t339;
t325 = t312 * t344;
t270 = t305 * t310 - t307 * t325;
t304 = sin(pkin(11));
t343 = cos(pkin(11));
t223 = t255 * t304 - t270 * t343;
t224 = t255 * t343 + t270 * t304;
t328 = t306 * t345;
t254 = t271 * t309 + t307 * t328;
t358 = t362 * t223 + t364 * t224 + t361 * t254;
t273 = -t305 * t326 + t307 * t312;
t257 = t273 * t345 + t305 * t339;
t272 = t305 * t325 + t307 * t310;
t225 = t257 * t304 - t272 * t343;
t226 = t257 * t343 + t272 * t304;
t256 = t273 * t309 - t305 * t328;
t357 = t362 * t225 + t364 * t226 + t361 * t256;
t356 = t364 * t223 + t365 * t224 + t363 * t254;
t355 = t364 * t225 + t365 * t226 + t363 * t256;
t278 = t309 * t344 + t310 * t328;
t338 = t306 * t312;
t252 = t278 * t304 + t338 * t343;
t253 = t278 * t343 - t304 * t338;
t277 = t310 * t339 - t344 * t345;
t354 = t362 * t252 + t364 * t253 + t361 * t277;
t353 = t364 * t252 + t365 * t253 + t363 * t277;
t351 = -Icges(4,4) * t255 - Icges(4,6) * t270 + t361 * t223 + t363 * t224 + t360 * t254;
t350 = -Icges(4,4) * t257 - Icges(4,6) * t272 + t361 * t225 + t363 * t226 + t360 * t256;
t349 = -Icges(4,4) * t278 + Icges(4,6) * t338 + t361 * t252 + t363 * t253 + t360 * t277;
t342 = Icges(2,4) * t305;
t341 = t305 * t306;
t340 = t306 * t307;
t190 = pkin(4) * t224 + qJ(5) * t223;
t217 = pkin(3) * t255 + qJ(4) * t254;
t337 = -t190 - t217;
t191 = pkin(4) * t226 + qJ(5) * t225;
t218 = pkin(3) * t257 + qJ(4) * t256;
t336 = -t191 - t218;
t216 = pkin(4) * t253 + qJ(5) * t252;
t244 = pkin(3) * t278 + qJ(4) * t277;
t335 = -t216 - t244;
t334 = qJD(2) * t306;
t333 = V_base(5) * qJ(1) + V_base(1);
t329 = qJD(1) + V_base(3);
t327 = t344 * pkin(7);
t286 = t305 * t334 + V_base(4);
t297 = qJD(2) * t344 + V_base(6);
t251 = qJD(3) * t272 + t286;
t285 = -t307 * t334 + V_base(5);
t250 = qJD(3) * t270 + t285;
t274 = -qJD(3) * t338 + t297;
t280 = pkin(1) * t305 - pkin(7) * t340;
t324 = -t280 * V_base(6) + V_base(5) * t327 + t333;
t281 = pkin(1) * t307 + pkin(7) * t341;
t323 = V_base(4) * t280 - t281 * V_base(5) + t329;
t322 = V_base(6) * t281 + V_base(2) + (-t327 - qJ(1)) * V_base(4);
t242 = pkin(2) * t271 + pkin(8) * t270;
t279 = (pkin(2) * t310 - pkin(8) * t312) * t306;
t321 = -t242 * t297 + t285 * t279 + t324;
t243 = pkin(2) * t273 + pkin(8) * t272;
t320 = t286 * t242 - t243 * t285 + t323;
t319 = qJD(4) * t256 + t250 * t244 + t321;
t318 = qJD(4) * t277 + t251 * t217 + t320;
t317 = t297 * t243 - t286 * t279 + t322;
t316 = qJD(5) * t225 + t250 * t216 + t319;
t315 = qJD(5) * t252 + t251 * t190 + t318;
t314 = qJD(4) * t254 + t274 * t218 + t317;
t313 = qJD(5) * t223 + t274 * t191 + t314;
t311 = cos(qJ(6));
t308 = sin(qJ(6));
t302 = Icges(2,4) * t307;
t294 = rSges(2,1) * t307 - rSges(2,2) * t305;
t293 = rSges(2,1) * t305 + rSges(2,2) * t307;
t292 = Icges(2,1) * t307 - t342;
t291 = Icges(2,1) * t305 + t302;
t290 = -Icges(2,2) * t305 + t302;
t289 = Icges(2,2) * t307 + t342;
t284 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t283 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t282 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t264 = t344 * rSges(3,3) + (rSges(3,1) * t310 + rSges(3,2) * t312) * t306;
t263 = Icges(3,5) * t344 + (Icges(3,1) * t310 + Icges(3,4) * t312) * t306;
t262 = Icges(3,6) * t344 + (Icges(3,4) * t310 + Icges(3,2) * t312) * t306;
t261 = Icges(3,3) * t344 + (Icges(3,5) * t310 + Icges(3,6) * t312) * t306;
t260 = V_base(5) * rSges(2,3) - t293 * V_base(6) + t333;
t259 = t294 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t249 = t293 * V_base(4) - t294 * V_base(5) + t329;
t245 = -qJD(6) * t277 + t274;
t241 = t278 * rSges(4,1) - t277 * rSges(4,2) - rSges(4,3) * t338;
t240 = Icges(4,1) * t278 - Icges(4,4) * t277 - Icges(4,5) * t338;
t238 = Icges(4,5) * t278 - Icges(4,6) * t277 - Icges(4,3) * t338;
t237 = rSges(3,1) * t273 - rSges(3,2) * t272 + rSges(3,3) * t341;
t236 = rSges(3,1) * t271 - rSges(3,2) * t270 - rSges(3,3) * t340;
t235 = Icges(3,1) * t273 - Icges(3,4) * t272 + Icges(3,5) * t341;
t234 = Icges(3,1) * t271 - Icges(3,4) * t270 - Icges(3,5) * t340;
t233 = Icges(3,4) * t273 - Icges(3,2) * t272 + Icges(3,6) * t341;
t232 = Icges(3,4) * t271 - Icges(3,2) * t270 - Icges(3,6) * t340;
t231 = Icges(3,5) * t273 - Icges(3,6) * t272 + Icges(3,3) * t341;
t230 = Icges(3,5) * t271 - Icges(3,6) * t270 - Icges(3,3) * t340;
t227 = pkin(5) * t253 - pkin(9) * t277;
t220 = -qJD(6) * t256 + t251;
t219 = -qJD(6) * t254 + t250;
t214 = t252 * t308 + t253 * t311;
t213 = t252 * t311 - t253 * t308;
t211 = rSges(5,1) * t253 - rSges(5,2) * t252 + rSges(5,3) * t277;
t210 = rSges(6,1) * t253 + rSges(6,2) * t277 + rSges(6,3) * t252;
t209 = rSges(4,1) * t257 - rSges(4,2) * t256 + rSges(4,3) * t272;
t208 = rSges(4,1) * t255 - rSges(4,2) * t254 + rSges(4,3) * t270;
t201 = Icges(4,1) * t257 - Icges(4,4) * t256 + Icges(4,5) * t272;
t200 = Icges(4,1) * t255 - Icges(4,4) * t254 + Icges(4,5) * t270;
t197 = Icges(4,5) * t257 - Icges(4,6) * t256 + Icges(4,3) * t272;
t196 = Icges(4,5) * t255 - Icges(4,6) * t254 + Icges(4,3) * t270;
t195 = pkin(5) * t226 - pkin(9) * t256;
t194 = pkin(5) * t224 - pkin(9) * t254;
t189 = t225 * t308 + t226 * t311;
t188 = t225 * t311 - t226 * t308;
t187 = t223 * t308 + t224 * t311;
t186 = t223 * t311 - t224 * t308;
t185 = -t236 * t297 + t264 * t285 + t324;
t184 = t297 * t237 - t286 * t264 + t322;
t181 = rSges(5,1) * t226 - rSges(5,2) * t225 + rSges(5,3) * t256;
t180 = rSges(6,1) * t226 + rSges(6,2) * t256 + rSges(6,3) * t225;
t179 = rSges(5,1) * t224 - rSges(5,2) * t223 + rSges(5,3) * t254;
t178 = rSges(6,1) * t224 + rSges(6,2) * t254 + rSges(6,3) * t223;
t165 = t236 * t286 - t237 * t285 + t323;
t164 = rSges(7,1) * t214 + rSges(7,2) * t213 - rSges(7,3) * t277;
t163 = Icges(7,1) * t214 + Icges(7,4) * t213 - Icges(7,5) * t277;
t162 = Icges(7,4) * t214 + Icges(7,2) * t213 - Icges(7,6) * t277;
t161 = Icges(7,5) * t214 + Icges(7,6) * t213 - Icges(7,3) * t277;
t160 = rSges(7,1) * t189 + rSges(7,2) * t188 - rSges(7,3) * t256;
t159 = rSges(7,1) * t187 + rSges(7,2) * t186 - rSges(7,3) * t254;
t158 = Icges(7,1) * t189 + Icges(7,4) * t188 - Icges(7,5) * t256;
t157 = Icges(7,1) * t187 + Icges(7,4) * t186 - Icges(7,5) * t254;
t156 = Icges(7,4) * t189 + Icges(7,2) * t188 - Icges(7,6) * t256;
t155 = Icges(7,4) * t187 + Icges(7,2) * t186 - Icges(7,6) * t254;
t154 = Icges(7,5) * t189 + Icges(7,6) * t188 - Icges(7,3) * t256;
t153 = Icges(7,5) * t187 + Icges(7,6) * t186 - Icges(7,3) * t254;
t152 = -t208 * t274 + t241 * t250 + t321;
t151 = t274 * t209 - t251 * t241 + t317;
t150 = t208 * t251 - t209 * t250 + t320;
t149 = t211 * t250 + (-t179 - t217) * t274 + t319;
t148 = t274 * t181 + (-t211 - t244) * t251 + t314;
t147 = t179 * t251 + (-t181 - t218) * t250 + t318;
t146 = t210 * t250 + (-t178 + t337) * t274 + t316;
t145 = t274 * t180 + (-t210 + t335) * t251 + t313;
t144 = t178 * t251 + (-t180 + t336) * t250 + t315;
t143 = t316 + (-t194 + t337) * t274 - t159 * t245 + t164 * t219 + t227 * t250;
t142 = t245 * t160 - t220 * t164 + t274 * t195 + (-t227 + t335) * t251 + t313;
t141 = t315 + (-t195 + t336) * t250 + t159 * t220 - t160 * t219 + t194 * t251;
t1 = t285 * ((-t231 * t340 - t233 * t270 + t235 * t271) * t286 + (-t230 * t340 - t270 * t232 + t271 * t234) * t285 + (-t261 * t340 - t262 * t270 + t263 * t271) * t297) / 0.2e1 + t286 * ((t231 * t341 - t272 * t233 + t273 * t235) * t286 + (t230 * t341 - t232 * t272 + t234 * t273) * t285 + (t261 * t341 - t262 * t272 + t263 * t273) * t297) / 0.2e1 + m(7) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(5) * (t147 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + m(6) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(4) * (t150 ^ 2 + t151 ^ 2 + t152 ^ 2) / 0.2e1 + m(3) * (t165 ^ 2 + t184 ^ 2 + t185 ^ 2) / 0.2e1 + t219 * ((-t154 * t254 + t156 * t186 + t158 * t187) * t220 + (-t153 * t254 + t155 * t186 + t157 * t187) * t219 + (-t161 * t254 + t162 * t186 + t163 * t187) * t245) / 0.2e1 + t220 * ((-t256 * t154 + t188 * t156 + t158 * t189) * t220 + (-t153 * t256 + t155 * t188 + t189 * t157) * t219 + (-t161 * t256 + t162 * t188 + t163 * t189) * t245) / 0.2e1 + m(2) * (t249 ^ 2 + t259 ^ 2 + t260 ^ 2) / 0.2e1 + t245 * ((-t154 * t277 + t156 * t213 + t158 * t214) * t220 + (-t153 * t277 + t155 * t213 + t157 * t214) * t219 + (-t277 * t161 + t213 * t162 + t163 * t214) * t245) / 0.2e1 + m(1) * (t282 ^ 2 + t283 ^ 2 + t284 ^ 2) / 0.2e1 + t297 * (((t233 * t312 + t235 * t310) * t286 + (t232 * t312 + t234 * t310) * t285 + (t262 * t312 + t263 * t310) * t297) * t306 + (t230 * t285 + t231 * t286 + t261 * t297) * t344) / 0.2e1 + ((-t289 * t305 + t291 * t307 + Icges(1,4)) * V_base(5) + (-t305 * t290 + t307 * t292 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t307 * t289 + t305 * t291 + Icges(1,2)) * V_base(5) + (t290 * t307 + t292 * t305 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t223 * t354 + t224 * t353 + t238 * t270 + t240 * t255 + t254 * t349) * t274 + (t197 * t270 + t201 * t255 + t223 * t357 + t224 * t355 + t254 * t350) * t251 + (t270 * t196 + t255 * t200 + t358 * t223 + t356 * t224 + t351 * t254) * t250) * t250 / 0.2e1 + ((t225 * t354 + t226 * t353 + t238 * t272 + t240 * t257 + t256 * t349) * t274 + (t272 * t197 + t257 * t201 + t357 * t225 + t355 * t226 + t350 * t256) * t251 + (t196 * t272 + t200 * t257 + t225 * t358 + t226 * t356 + t256 * t351) * t250) * t251 / 0.2e1 + ((-t238 * t338 + t278 * t240 + t354 * t252 + t353 * t253 + t349 * t277) * t274 + (-t197 * t338 + t278 * t201 + t252 * t357 + t253 * t355 + t277 * t350) * t251 + (-t196 * t338 + t278 * t200 + t252 * t358 + t253 * t356 + t277 * t351) * t250) * t274 / 0.2e1 + ((Icges(2,5) * t307 - Icges(2,6) * t305 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t305 + Icges(2,6) * t307 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
