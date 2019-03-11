% Calculate kinetic energy for
% S6RRRPPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
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
% Datum: 2019-03-09 15:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPP1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPPP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPP1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPP1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPP1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:14:23
% EndTime: 2019-03-09 15:14:26
% DurationCPUTime: 3.01s
% Computational Cost: add. (2313->319), mult. (5550->443), div. (0->0), fcn. (6561->10), ass. (0->155)
t372 = Icges(5,1) + Icges(6,2) + Icges(7,3);
t371 = Icges(6,1) + Icges(7,1) + Icges(5,3);
t370 = -Icges(5,4) - Icges(6,6) + Icges(7,6);
t369 = -Icges(6,4) + Icges(5,5) + Icges(7,5);
t368 = Icges(7,4) + Icges(6,5) - Icges(5,6);
t367 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t366 = rSges(7,1) + pkin(5);
t365 = rSges(7,3) + qJ(6);
t298 = sin(qJ(3));
t301 = cos(qJ(3));
t303 = cos(qJ(1));
t300 = sin(qJ(1));
t302 = cos(qJ(2));
t340 = t300 * t302;
t261 = -t298 * t340 - t301 * t303;
t262 = -t298 * t303 + t301 * t340;
t299 = sin(qJ(2));
t349 = sin(pkin(6));
t350 = cos(pkin(10));
t324 = t349 * t350;
t318 = t299 * t324;
t351 = cos(pkin(6));
t326 = t351 * t350;
t348 = sin(pkin(10));
t213 = -t261 * t326 + t262 * t348 - t300 * t318;
t323 = t349 * t348;
t317 = t299 * t323;
t325 = t351 * t348;
t214 = t261 * t325 + t262 * t350 + t300 * t317;
t328 = t299 * t351;
t234 = -t261 * t349 + t300 * t328;
t364 = t370 * t213 + t372 * t214 + t369 * t234;
t339 = t302 * t303;
t263 = -t298 * t339 + t300 * t301;
t264 = t300 * t298 + t301 * t339;
t215 = -t263 * t326 + t264 * t348 - t303 * t318;
t216 = t263 * t325 + t264 * t350 + t303 * t317;
t235 = -t263 * t349 + t303 * t328;
t363 = t370 * t215 + t372 * t216 + t369 * t235;
t362 = t367 * t213 + t370 * t214 + t368 * t234;
t361 = t367 * t215 + t370 * t216 + t368 * t235;
t360 = t368 * t213 + t369 * t214 + t371 * t234;
t359 = t368 * t215 + t369 * t216 + t371 * t235;
t230 = t302 * t324 + (t298 * t326 + t301 * t348) * t299;
t342 = t299 * t301;
t344 = t298 * t299;
t231 = -t302 * t323 - t325 * t344 + t342 * t350;
t258 = -t302 * t351 + t344 * t349;
t358 = t370 * t230 + t372 * t231 + t369 * t258;
t357 = t367 * t230 + t370 * t231 + t368 * t258;
t356 = t368 * t230 + t369 * t231 + t371 * t258;
t347 = Icges(2,4) * t300;
t346 = Icges(3,4) * t299;
t345 = Icges(3,4) * t302;
t343 = t299 * t300;
t341 = t299 * t303;
t338 = rSges(7,2) * t213 + t365 * t214 + t366 * t234;
t337 = rSges(7,2) * t215 + t365 * t216 + t366 * t235;
t188 = pkin(4) * t214 + qJ(5) * t213;
t218 = t262 * pkin(3) + qJ(4) * t234;
t336 = -t188 - t218;
t189 = pkin(4) * t216 + qJ(5) * t215;
t219 = t264 * pkin(3) + qJ(4) * t235;
t335 = -t189 - t219;
t334 = rSges(7,2) * t230 + t365 * t231 + t366 * t258;
t206 = pkin(4) * t231 + qJ(5) * t230;
t239 = pkin(3) * t342 + qJ(4) * t258;
t333 = -t206 - t239;
t332 = qJD(3) * t299;
t331 = V_base(5) * pkin(7) + V_base(1);
t293 = qJD(2) * t300 + V_base(4);
t294 = V_base(6) + qJD(1);
t327 = pkin(2) * t302 + pkin(9) * t299;
t292 = -qJD(2) * t303 + V_base(5);
t322 = rSges(3,1) * t302 - rSges(3,2) * t299;
t321 = Icges(3,1) * t302 - t346;
t320 = -Icges(3,2) * t299 + t345;
t319 = Icges(3,5) * t302 - Icges(3,6) * t299;
t291 = pkin(1) * t303 + t300 * pkin(8);
t316 = -V_base(4) * pkin(7) + t294 * t291 + V_base(2);
t290 = t300 * pkin(1) - pkin(8) * t303;
t315 = V_base(4) * t290 - t291 * V_base(5) + V_base(3);
t266 = t327 * t300;
t289 = pkin(2) * t299 - pkin(9) * t302;
t314 = t292 * t289 + (-t266 - t290) * t294 + t331;
t313 = (-Icges(3,3) * t303 + t300 * t319) * t292 + (Icges(3,3) * t300 + t303 * t319) * t293 + (Icges(3,5) * t299 + Icges(3,6) * t302) * t294;
t267 = t327 * t303;
t312 = t294 * t267 - t289 * t293 + t316;
t259 = t300 * t332 + t292;
t311 = qJD(4) * t235 + t259 * t239 + t314;
t310 = t293 * t266 - t267 * t292 + t315;
t284 = -qJD(3) * t302 + t294;
t309 = qJD(4) * t234 + t284 * t219 + t312;
t308 = qJD(5) * t215 + t259 * t206 + t311;
t260 = t303 * t332 + t293;
t307 = qJD(4) * t258 + t260 * t218 + t310;
t306 = qJD(5) * t213 + t284 * t189 + t309;
t305 = qJD(5) * t230 + t260 * t188 + t307;
t248 = -Icges(3,6) * t303 + t300 * t320;
t249 = Icges(3,6) * t300 + t303 * t320;
t251 = -Icges(3,5) * t303 + t300 * t321;
t252 = Icges(3,5) * t300 + t303 * t321;
t278 = Icges(3,2) * t302 + t346;
t281 = Icges(3,1) * t299 + t345;
t304 = (-t249 * t299 + t252 * t302) * t293 + (-t248 * t299 + t251 * t302) * t292 + (-t278 * t299 + t281 * t302) * t294;
t296 = Icges(2,4) * t303;
t287 = rSges(2,1) * t303 - t300 * rSges(2,2);
t286 = t300 * rSges(2,1) + rSges(2,2) * t303;
t285 = rSges(3,1) * t299 + rSges(3,2) * t302;
t283 = Icges(2,1) * t303 - t347;
t282 = Icges(2,1) * t300 + t296;
t280 = -Icges(2,2) * t300 + t296;
t279 = Icges(2,2) * t303 + t347;
t274 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t273 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t272 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t255 = t300 * rSges(3,3) + t303 * t322;
t254 = -rSges(3,3) * t303 + t300 * t322;
t253 = -rSges(4,3) * t302 + (rSges(4,1) * t301 - rSges(4,2) * t298) * t299;
t250 = -Icges(4,5) * t302 + (Icges(4,1) * t301 - Icges(4,4) * t298) * t299;
t247 = -Icges(4,6) * t302 + (Icges(4,4) * t301 - Icges(4,2) * t298) * t299;
t244 = -Icges(4,3) * t302 + (Icges(4,5) * t301 - Icges(4,6) * t298) * t299;
t238 = V_base(5) * rSges(2,3) - t286 * t294 + t331;
t237 = t287 * t294 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t236 = t286 * V_base(4) - t287 * V_base(5) + V_base(3);
t227 = t264 * rSges(4,1) + t263 * rSges(4,2) + rSges(4,3) * t341;
t226 = rSges(4,1) * t262 + rSges(4,2) * t261 + rSges(4,3) * t343;
t225 = Icges(4,1) * t264 + Icges(4,4) * t263 + Icges(4,5) * t341;
t224 = Icges(4,1) * t262 + Icges(4,4) * t261 + Icges(4,5) * t343;
t223 = Icges(4,4) * t264 + Icges(4,2) * t263 + Icges(4,6) * t341;
t222 = Icges(4,4) * t262 + Icges(4,2) * t261 + Icges(4,6) * t343;
t221 = Icges(4,5) * t264 + Icges(4,6) * t263 + Icges(4,3) * t341;
t220 = Icges(4,5) * t262 + Icges(4,6) * t261 + Icges(4,3) * t343;
t208 = t285 * t292 + (-t254 - t290) * t294 + t331;
t207 = t255 * t294 - t285 * t293 + t316;
t204 = t254 * t293 - t255 * t292 + t315;
t203 = rSges(6,1) * t258 - rSges(6,2) * t231 + rSges(6,3) * t230;
t201 = rSges(5,1) * t231 - rSges(5,2) * t230 + rSges(5,3) * t258;
t185 = rSges(6,1) * t235 - rSges(6,2) * t216 + rSges(6,3) * t215;
t183 = rSges(6,1) * t234 - rSges(6,2) * t214 + rSges(6,3) * t213;
t181 = rSges(5,1) * t216 - rSges(5,2) * t215 + rSges(5,3) * t235;
t180 = rSges(5,1) * t214 - rSges(5,2) * t213 + rSges(5,3) * t234;
t161 = -t226 * t284 + t253 * t259 + t314;
t160 = t227 * t284 - t253 * t260 + t312;
t159 = t226 * t260 - t227 * t259 + t310;
t158 = t201 * t259 + (-t180 - t218) * t284 + t311;
t157 = t181 * t284 + (-t201 - t239) * t260 + t309;
t156 = t180 * t260 + (-t181 - t219) * t259 + t307;
t155 = t203 * t259 + (-t183 + t336) * t284 + t308;
t154 = t185 * t284 + (-t203 + t333) * t260 + t306;
t153 = t183 * t260 + (-t185 + t335) * t259 + t305;
t152 = qJD(6) * t216 + t334 * t259 + (t336 - t338) * t284 + t308;
t151 = qJD(6) * t214 + t337 * t284 + (t333 - t334) * t260 + t306;
t150 = qJD(6) * t231 + t338 * t260 + (t335 - t337) * t259 + t305;
t1 = m(4) * (t159 ^ 2 + t160 ^ 2 + t161 ^ 2) / 0.2e1 + m(5) * (t156 ^ 2 + t157 ^ 2 + t158 ^ 2) / 0.2e1 + m(7) * (t150 ^ 2 + t151 ^ 2 + t152 ^ 2) / 0.2e1 + m(6) * (t153 ^ 2 + t154 ^ 2 + t155 ^ 2) / 0.2e1 + m(3) * (t204 ^ 2 + t207 ^ 2 + t208 ^ 2) / 0.2e1 + m(2) * (t236 ^ 2 + t237 ^ 2 + t238 ^ 2) / 0.2e1 + m(1) * (t272 ^ 2 + t273 ^ 2 + t274 ^ 2) / 0.2e1 + t293 * (t313 * t300 + t304 * t303) / 0.2e1 + t292 * (t304 * t300 - t313 * t303) / 0.2e1 + ((t249 * t302 + t252 * t299) * t293 + (t248 * t302 + t251 * t299) * t292 + (t278 * t302 + t281 * t299 + Icges(2,3)) * t294) * t294 / 0.2e1 + ((-t300 * t279 + t282 * t303 + Icges(1,4)) * V_base(5) + (-t300 * t280 + t283 * t303 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t279 * t303 + t300 * t282 + Icges(1,2)) * V_base(5) + (t280 * t303 + t300 * t283 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t213 * t357 + t214 * t358 + t234 * t356 + t244 * t343 + t247 * t261 + t250 * t262) * t284 + (t213 * t361 + t214 * t363 + t221 * t343 + t223 * t261 + t225 * t262 + t234 * t359) * t260 + (t362 * t213 + t364 * t214 + t220 * t343 + t222 * t261 + t224 * t262 + t360 * t234) * t259) * t259 / 0.2e1 + ((t215 * t357 + t216 * t358 + t235 * t356 + t244 * t341 + t263 * t247 + t264 * t250) * t284 + (t361 * t215 + t363 * t216 + t221 * t341 + t263 * t223 + t264 * t225 + t359 * t235) * t260 + (t362 * t215 + t216 * t364 + t220 * t341 + t263 * t222 + t264 * t224 + t360 * t235) * t259) * t260 / 0.2e1 + ((-t244 * t302 + (-t247 * t298 + t250 * t301) * t299 + t356 * t258 + t358 * t231 + t357 * t230) * t284 + (-t221 * t302 + (-t223 * t298 + t225 * t301) * t299 + t359 * t258 + t363 * t231 + t361 * t230) * t260 + (-t220 * t302 + (-t222 * t298 + t224 * t301) * t299 + t360 * t258 + t364 * t231 + t362 * t230) * t259) * t284 / 0.2e1 + V_base(4) * t294 * (Icges(2,5) * t303 - Icges(2,6) * t300) + V_base(5) * t294 * (Icges(2,5) * t300 + Icges(2,6) * t303) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
