% Calculate kinetic energy for
% S6RPRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR12_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR12_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPR12_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR12_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR12_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR12_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:47:58
% EndTime: 2019-03-09 05:48:03
% DurationCPUTime: 4.96s
% Computational Cost: add. (4786->406), mult. (12529->575), div. (0->0), fcn. (15975->14), ass. (0->179)
t412 = Icges(5,1) + Icges(6,2);
t411 = Icges(6,1) + Icges(5,3);
t410 = -Icges(5,4) - Icges(6,6);
t409 = -Icges(6,4) + Icges(5,5);
t408 = Icges(6,5) - Icges(5,6);
t407 = Icges(5,2) + Icges(6,3);
t352 = cos(pkin(12));
t350 = sin(pkin(12));
t357 = sin(qJ(1));
t383 = t357 * t350;
t353 = cos(pkin(6));
t359 = cos(qJ(1));
t384 = t353 * t359;
t323 = t352 * t384 - t383;
t382 = t357 * t352;
t324 = t350 * t384 + t382;
t356 = sin(qJ(3));
t351 = sin(pkin(6));
t390 = sin(pkin(7));
t375 = t351 * t390;
t391 = cos(pkin(7));
t393 = cos(qJ(3));
t282 = t324 * t393 + (t323 * t391 - t359 * t375) * t356;
t376 = t351 * t391;
t305 = -t323 * t390 - t359 * t376;
t355 = sin(qJ(4));
t392 = cos(qJ(4));
t264 = t282 * t355 - t305 * t392;
t265 = t282 * t392 + t305 * t355;
t372 = t393 * t390;
t370 = t351 * t372;
t373 = t391 * t393;
t281 = -t323 * t373 + t324 * t356 + t359 * t370;
t406 = t264 * t407 + t265 * t410 + t281 * t408;
t325 = -t350 * t359 - t353 * t382;
t326 = t352 * t359 - t353 * t383;
t284 = t326 * t393 + (t325 * t391 + t357 * t375) * t356;
t306 = -t325 * t390 + t357 * t376;
t266 = t284 * t355 - t306 * t392;
t267 = t284 * t392 + t306 * t355;
t283 = -t325 * t373 + t326 * t356 - t357 * t370;
t405 = t266 * t407 + t267 * t410 + t283 * t408;
t404 = t264 * t408 + t265 * t409 + t281 * t411;
t403 = t266 * t408 + t267 * t409 + t283 * t411;
t402 = t410 * t264 + t265 * t412 + t409 * t281;
t401 = t410 * t266 + t267 * t412 + t409 * t283;
t304 = t353 * t390 * t356 + (t352 * t356 * t391 + t350 * t393) * t351;
t322 = -t352 * t375 + t353 * t391;
t278 = t304 * t355 - t322 * t392;
t279 = t304 * t392 + t322 * t355;
t387 = t350 * t351;
t303 = -t351 * t352 * t373 - t353 * t372 + t356 * t387;
t400 = t278 * t407 + t279 * t410 + t303 * t408;
t399 = t278 * t408 + t279 * t409 + t303 * t411;
t398 = t410 * t278 + t279 * t412 + t409 * t303;
t389 = Icges(2,4) * t357;
t388 = qJ(2) * t353;
t386 = t351 * t357;
t385 = t351 * t359;
t381 = qJD(2) * t351;
t380 = V_base(5) * pkin(8) + V_base(1);
t297 = qJD(3) * t306 + V_base(4);
t296 = qJD(3) * t305 + V_base(5);
t347 = V_base(6) + qJD(1);
t377 = -pkin(8) - t388;
t328 = t357 * pkin(1) - qJ(2) * t385;
t374 = qJD(2) * t353 + V_base(4) * t328 + V_base(3);
t263 = qJD(4) * t283 + t297;
t262 = qJD(4) * t281 + t296;
t312 = qJD(3) * t322 + t347;
t371 = t357 * t381 + V_base(5) * t388 + t380;
t275 = qJD(4) * t303 + t312;
t329 = pkin(1) * t359 + qJ(2) * t386;
t369 = t347 * t329 - t359 * t381 + V_base(2);
t286 = t324 * pkin(2) + pkin(9) * t305;
t309 = pkin(2) * t387 + pkin(9) * t322;
t368 = V_base(5) * t309 + (-t286 - t328) * t347 + t371;
t287 = t326 * pkin(2) + pkin(9) * t306;
t367 = V_base(4) * t286 + (-t287 - t329) * V_base(5) + t374;
t255 = pkin(3) * t282 + pkin(10) * t281;
t273 = pkin(3) * t304 + pkin(10) * t303;
t366 = -t255 * t312 + t296 * t273 + t368;
t256 = pkin(3) * t284 + pkin(10) * t283;
t365 = t297 * t255 - t256 * t296 + t367;
t364 = t347 * t287 + (-t309 + t377) * V_base(4) + t369;
t251 = pkin(4) * t279 + qJ(5) * t278;
t363 = qJD(5) * t266 + t262 * t251 + t366;
t222 = pkin(4) * t265 + qJ(5) * t264;
t362 = qJD(5) * t278 + t263 * t222 + t365;
t361 = t312 * t256 - t297 * t273 + t364;
t223 = pkin(4) * t267 + qJ(5) * t266;
t360 = qJD(5) * t264 + t275 * t223 + t361;
t358 = cos(qJ(6));
t354 = sin(qJ(6));
t348 = Icges(2,4) * t359;
t342 = rSges(2,1) * t359 - t357 * rSges(2,2);
t341 = t357 * rSges(2,1) + rSges(2,2) * t359;
t340 = Icges(2,1) * t359 - t389;
t339 = Icges(2,1) * t357 + t348;
t338 = -Icges(2,2) * t357 + t348;
t337 = Icges(2,2) * t359 + t389;
t336 = Icges(2,5) * t359 - Icges(2,6) * t357;
t335 = Icges(2,5) * t357 + Icges(2,6) * t359;
t334 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t333 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t332 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t318 = rSges(3,3) * t353 + (rSges(3,1) * t350 + rSges(3,2) * t352) * t351;
t317 = Icges(3,5) * t353 + (Icges(3,1) * t350 + Icges(3,4) * t352) * t351;
t316 = Icges(3,6) * t353 + (Icges(3,4) * t350 + Icges(3,2) * t352) * t351;
t315 = Icges(3,3) * t353 + (Icges(3,5) * t350 + Icges(3,6) * t352) * t351;
t311 = V_base(5) * rSges(2,3) - t341 * t347 + t380;
t310 = t342 * t347 + V_base(2) + (-rSges(2,3) - pkin(8)) * V_base(4);
t308 = t341 * V_base(4) - t342 * V_base(5) + V_base(3);
t295 = rSges(3,1) * t326 + rSges(3,2) * t325 + rSges(3,3) * t386;
t294 = t324 * rSges(3,1) + t323 * rSges(3,2) - rSges(3,3) * t385;
t293 = Icges(3,1) * t326 + Icges(3,4) * t325 + Icges(3,5) * t386;
t292 = Icges(3,1) * t324 + Icges(3,4) * t323 - Icges(3,5) * t385;
t291 = Icges(3,4) * t326 + Icges(3,2) * t325 + Icges(3,6) * t386;
t290 = Icges(3,4) * t324 + Icges(3,2) * t323 - Icges(3,6) * t385;
t289 = Icges(3,5) * t326 + Icges(3,6) * t325 + Icges(3,3) * t386;
t288 = Icges(3,5) * t324 + Icges(3,6) * t323 - Icges(3,3) * t385;
t272 = rSges(4,1) * t304 - rSges(4,2) * t303 + rSges(4,3) * t322;
t271 = Icges(4,1) * t304 - Icges(4,4) * t303 + Icges(4,5) * t322;
t270 = Icges(4,4) * t304 - Icges(4,2) * t303 + Icges(4,6) * t322;
t269 = Icges(4,5) * t304 - Icges(4,6) * t303 + Icges(4,3) * t322;
t268 = pkin(5) * t303 + pkin(11) * t279;
t259 = t278 * t354 + t303 * t358;
t258 = t278 * t358 - t303 * t354;
t254 = t318 * V_base(5) + (-t294 - t328) * t347 + t371;
t253 = t347 * t295 + (-t318 + t377) * V_base(4) + t369;
t252 = qJD(6) * t279 + t275;
t250 = t294 * V_base(4) + (-t295 - t329) * V_base(5) + t374;
t248 = rSges(4,1) * t284 - rSges(4,2) * t283 + rSges(4,3) * t306;
t247 = rSges(4,1) * t282 - rSges(4,2) * t281 + rSges(4,3) * t305;
t246 = Icges(4,1) * t284 - Icges(4,4) * t283 + Icges(4,5) * t306;
t245 = Icges(4,1) * t282 - Icges(4,4) * t281 + Icges(4,5) * t305;
t244 = Icges(4,4) * t284 - Icges(4,2) * t283 + Icges(4,6) * t306;
t243 = Icges(4,4) * t282 - Icges(4,2) * t281 + Icges(4,6) * t305;
t242 = Icges(4,5) * t284 - Icges(4,6) * t283 + Icges(4,3) * t306;
t241 = Icges(4,5) * t282 - Icges(4,6) * t281 + Icges(4,3) * t305;
t240 = rSges(5,1) * t279 - rSges(5,2) * t278 + rSges(5,3) * t303;
t239 = rSges(6,1) * t303 - rSges(6,2) * t279 + rSges(6,3) * t278;
t231 = pkin(5) * t283 + pkin(11) * t267;
t230 = pkin(5) * t281 + pkin(11) * t265;
t229 = t266 * t354 + t283 * t358;
t228 = t266 * t358 - t283 * t354;
t227 = t264 * t354 + t281 * t358;
t226 = t264 * t358 - t281 * t354;
t225 = qJD(6) * t267 + t263;
t224 = qJD(6) * t265 + t262;
t220 = rSges(5,1) * t267 - rSges(5,2) * t266 + rSges(5,3) * t283;
t219 = rSges(5,1) * t265 - rSges(5,2) * t264 + rSges(5,3) * t281;
t218 = rSges(6,1) * t283 - rSges(6,2) * t267 + rSges(6,3) * t266;
t217 = rSges(6,1) * t281 - rSges(6,2) * t265 + rSges(6,3) * t264;
t203 = rSges(7,1) * t259 + rSges(7,2) * t258 + rSges(7,3) * t279;
t202 = Icges(7,1) * t259 + Icges(7,4) * t258 + Icges(7,5) * t279;
t201 = Icges(7,4) * t259 + Icges(7,2) * t258 + Icges(7,6) * t279;
t200 = Icges(7,5) * t259 + Icges(7,6) * t258 + Icges(7,3) * t279;
t198 = rSges(7,1) * t229 + rSges(7,2) * t228 + rSges(7,3) * t267;
t197 = rSges(7,1) * t227 + rSges(7,2) * t226 + rSges(7,3) * t265;
t196 = Icges(7,1) * t229 + Icges(7,4) * t228 + Icges(7,5) * t267;
t195 = Icges(7,1) * t227 + Icges(7,4) * t226 + Icges(7,5) * t265;
t194 = Icges(7,4) * t229 + Icges(7,2) * t228 + Icges(7,6) * t267;
t193 = Icges(7,4) * t227 + Icges(7,2) * t226 + Icges(7,6) * t265;
t192 = Icges(7,5) * t229 + Icges(7,6) * t228 + Icges(7,3) * t267;
t191 = Icges(7,5) * t227 + Icges(7,6) * t226 + Icges(7,3) * t265;
t190 = -t247 * t312 + t272 * t296 + t368;
t189 = t312 * t248 - t297 * t272 + t364;
t188 = t247 * t297 - t248 * t296 + t367;
t187 = -t219 * t275 + t240 * t262 + t366;
t186 = t275 * t220 - t263 * t240 + t361;
t185 = t219 * t263 - t220 * t262 + t365;
t184 = t239 * t262 + (-t217 - t222) * t275 + t363;
t183 = t275 * t218 + (-t239 - t251) * t263 + t360;
t182 = t217 * t263 + (-t218 - t223) * t262 + t362;
t181 = t363 - t197 * t252 + t203 * t224 + t262 * t268 + (-t222 - t230) * t275;
t180 = t252 * t198 - t225 * t203 + t275 * t231 + (-t251 - t268) * t263 + t360;
t179 = t197 * t225 - t198 * t224 + t230 * t263 + (-t223 - t231) * t262 + t362;
t1 = m(7) * (t179 ^ 2 + t180 ^ 2 + t181 ^ 2) / 0.2e1 + m(4) * (t188 ^ 2 + t189 ^ 2 + t190 ^ 2) / 0.2e1 + m(5) * (t185 ^ 2 + t186 ^ 2 + t187 ^ 2) / 0.2e1 + m(6) * (t182 ^ 2 + t183 ^ 2 + t184 ^ 2) / 0.2e1 + m(3) * (t250 ^ 2 + t253 ^ 2 + t254 ^ 2) / 0.2e1 + t224 * ((t192 * t265 + t194 * t226 + t196 * t227) * t225 + (t265 * t191 + t226 * t193 + t227 * t195) * t224 + (t200 * t265 + t201 * t226 + t202 * t227) * t252) / 0.2e1 + t225 * ((t267 * t192 + t228 * t194 + t229 * t196) * t225 + (t191 * t267 + t193 * t228 + t195 * t229) * t224 + (t200 * t267 + t201 * t228 + t202 * t229) * t252) / 0.2e1 + t252 * ((t192 * t279 + t194 * t258 + t196 * t259) * t225 + (t191 * t279 + t193 * t258 + t195 * t259) * t224 + (t200 * t279 + t201 * t258 + t202 * t259) * t252) / 0.2e1 + m(2) * (t308 ^ 2 + t310 ^ 2 + t311 ^ 2) / 0.2e1 + t296 * ((t242 * t305 - t244 * t281 + t246 * t282) * t297 + (t241 * t305 - t243 * t281 + t245 * t282) * t296 + (t269 * t305 - t270 * t281 + t271 * t282) * t312) / 0.2e1 + t297 * ((t242 * t306 - t244 * t283 + t246 * t284) * t297 + (t241 * t306 - t243 * t283 + t245 * t284) * t296 + (t269 * t306 - t270 * t283 + t271 * t284) * t312) / 0.2e1 + t312 * ((t242 * t322 - t244 * t303 + t246 * t304) * t297 + (t241 * t322 - t243 * t303 + t245 * t304) * t296 + (t269 * t322 - t270 * t303 + t271 * t304) * t312) / 0.2e1 + m(1) * (t332 ^ 2 + t333 ^ 2 + t334 ^ 2) / 0.2e1 + ((t264 * t400 + t265 * t398 + t281 * t399) * t275 + (t264 * t405 + t265 * t401 + t281 * t403) * t263 + (t406 * t264 + t402 * t265 + t404 * t281) * t262) * t262 / 0.2e1 + ((t266 * t400 + t267 * t398 + t283 * t399) * t275 + (t405 * t266 + t401 * t267 + t403 * t283) * t263 + (t266 * t406 + t402 * t267 + t404 * t283) * t262) * t263 / 0.2e1 + ((t278 * t400 + t279 * t398 + t303 * t399) * t275 + (t278 * t405 + t279 * t401 + t303 * t403) * t263 + (t278 * t406 + t402 * t279 + t404 * t303) * t262) * t275 / 0.2e1 + ((t288 * V_base(5) + t289 * V_base(4) + t315 * t347) * t353 + ((t291 * t352 + t293 * t350) * V_base(4) + (t290 * t352 + t292 * t350) * V_base(5) + (t316 * t352 + t317 * t350) * t347) * t351 + Icges(2,3) * t347 + t335 * V_base(5) + t336 * V_base(4)) * t347 / 0.2e1 + ((t315 * t386 + t316 * t325 + t317 * t326 + t336) * t347 + (t288 * t386 + t290 * t325 + t292 * t326 - t357 * t337 + t339 * t359 + Icges(1,4)) * V_base(5) + (t289 * t386 + t291 * t325 + t293 * t326 - t357 * t338 + t340 * t359 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t315 * t385 + t323 * t316 + t324 * t317 + t335) * t347 + (-t288 * t385 + t323 * t290 + t324 * t292 + t337 * t359 + t357 * t339 + Icges(1,2)) * V_base(5) + (-t289 * t385 + t323 * t291 + t324 * t293 + t338 * t359 + t357 * t340 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
