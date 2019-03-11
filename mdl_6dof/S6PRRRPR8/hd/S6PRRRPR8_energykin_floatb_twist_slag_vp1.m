% Calculate kinetic energy for
% S6PRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPR8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR8_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRPR8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR8_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR8_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPR8_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:48:06
% EndTime: 2019-03-08 23:48:10
% DurationCPUTime: 4.84s
% Computational Cost: add. (4861->401), mult. (12754->585), div. (0->0), fcn. (16200->14), ass. (0->178)
t413 = Icges(5,1) + Icges(6,2);
t412 = Icges(6,1) + Icges(5,3);
t411 = -Icges(5,4) - Icges(6,6);
t410 = -Icges(6,4) + Icges(5,5);
t409 = Icges(6,5) - Icges(5,6);
t408 = Icges(5,2) + Icges(6,3);
t350 = sin(pkin(12));
t352 = cos(pkin(12));
t357 = sin(qJ(2));
t353 = cos(pkin(6));
t359 = cos(qJ(2));
t383 = t353 * t359;
t320 = -t350 * t357 + t352 * t383;
t384 = t353 * t357;
t321 = t350 * t359 + t352 * t384;
t356 = sin(qJ(3));
t351 = sin(pkin(6));
t389 = sin(pkin(7));
t375 = t351 * t389;
t390 = cos(pkin(7));
t393 = cos(qJ(3));
t279 = t321 * t393 + (t320 * t390 - t352 * t375) * t356;
t376 = t351 * t390;
t305 = -t320 * t389 - t352 * t376;
t355 = sin(qJ(4));
t392 = cos(qJ(4));
t262 = t279 * t355 - t305 * t392;
t263 = t279 * t392 + t305 * t355;
t373 = t393 * t389;
t372 = t351 * t373;
t374 = t390 * t393;
t278 = -t320 * t374 + t321 * t356 + t352 * t372;
t406 = t408 * t262 + t411 * t263 + t409 * t278;
t322 = -t350 * t383 - t352 * t357;
t323 = -t350 * t384 + t352 * t359;
t281 = t323 * t393 + (t322 * t390 + t350 * t375) * t356;
t306 = -t322 * t389 + t350 * t376;
t264 = t281 * t355 - t306 * t392;
t265 = t281 * t392 + t306 * t355;
t280 = -t322 * t374 + t323 * t356 - t350 * t372;
t405 = t408 * t264 + t411 * t265 + t409 * t280;
t404 = t409 * t262 + t410 * t263 + t412 * t278;
t403 = t409 * t264 + t410 * t265 + t412 * t280;
t402 = t411 * t262 + t413 * t263 + t410 * t278;
t401 = t411 * t264 + t413 * t265 + t410 * t280;
t304 = t353 * t389 * t356 + (t356 * t359 * t390 + t357 * t393) * t351;
t319 = t353 * t390 - t359 * t375;
t282 = t304 * t355 - t319 * t392;
t283 = t304 * t392 + t319 * t355;
t385 = t351 * t357;
t303 = -t351 * t359 * t374 - t353 * t373 + t356 * t385;
t400 = t408 * t282 + t411 * t283 + t409 * t303;
t399 = t409 * t282 + t410 * t283 + t412 * t303;
t398 = t411 * t282 + t413 * t283 + t410 * t303;
t391 = pkin(8) * t353;
t388 = Icges(2,4) * t350;
t387 = t350 * t351;
t386 = t351 * t352;
t382 = qJD(2) * t351;
t381 = V_base(5) * qJ(1) + V_base(1);
t377 = qJD(1) + V_base(3);
t334 = t350 * t382 + V_base(4);
t344 = qJD(2) * t353 + V_base(6);
t295 = qJD(3) * t306 + t334;
t307 = qJD(3) * t319 + t344;
t257 = qJD(4) * t280 + t295;
t273 = qJD(4) * t303 + t307;
t333 = -t352 * t382 + V_base(5);
t294 = qJD(3) * t305 + t333;
t326 = pkin(1) * t350 - pkin(8) * t386;
t371 = -t326 * V_base(6) + V_base(5) * t391 + t381;
t327 = pkin(1) * t352 + pkin(8) * t387;
t370 = V_base(4) * t326 - t327 * V_base(5) + t377;
t256 = qJD(4) * t278 + t294;
t369 = V_base(6) * t327 + V_base(2) + (-qJ(1) - t391) * V_base(4);
t284 = t321 * pkin(2) + pkin(9) * t305;
t308 = pkin(2) * t385 + pkin(9) * t319;
t368 = -t284 * t344 + t333 * t308 + t371;
t285 = t323 * pkin(2) + pkin(9) * t306;
t367 = t334 * t284 - t285 * t333 + t370;
t366 = t344 * t285 - t308 * t334 + t369;
t252 = pkin(3) * t279 + pkin(10) * t278;
t271 = pkin(3) * t304 + pkin(10) * t303;
t365 = -t252 * t307 + t294 * t271 + t368;
t253 = pkin(3) * t281 + pkin(10) * t280;
t364 = t295 * t252 - t253 * t294 + t367;
t363 = t307 * t253 - t271 * t295 + t366;
t254 = pkin(4) * t283 + qJ(5) * t282;
t362 = qJD(5) * t264 + t256 * t254 + t365;
t222 = pkin(4) * t263 + qJ(5) * t262;
t361 = qJD(5) * t282 + t257 * t222 + t364;
t223 = pkin(4) * t265 + qJ(5) * t264;
t360 = qJD(5) * t262 + t273 * t223 + t363;
t358 = cos(qJ(6));
t354 = sin(qJ(6));
t348 = Icges(2,4) * t352;
t342 = rSges(2,1) * t352 - rSges(2,2) * t350;
t341 = rSges(2,1) * t350 + rSges(2,2) * t352;
t340 = Icges(2,1) * t352 - t388;
t339 = Icges(2,1) * t350 + t348;
t338 = -Icges(2,2) * t350 + t348;
t337 = Icges(2,2) * t352 + t388;
t332 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t331 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t330 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t316 = t353 * rSges(3,3) + (rSges(3,1) * t357 + rSges(3,2) * t359) * t351;
t315 = Icges(3,5) * t353 + (Icges(3,1) * t357 + Icges(3,4) * t359) * t351;
t314 = Icges(3,6) * t353 + (Icges(3,4) * t357 + Icges(3,2) * t359) * t351;
t313 = Icges(3,3) * t353 + (Icges(3,5) * t357 + Icges(3,6) * t359) * t351;
t310 = V_base(5) * rSges(2,3) - t341 * V_base(6) + t381;
t309 = t342 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t302 = t341 * V_base(4) - t342 * V_base(5) + t377;
t293 = rSges(3,1) * t323 + rSges(3,2) * t322 + rSges(3,3) * t387;
t292 = rSges(3,1) * t321 + rSges(3,2) * t320 - rSges(3,3) * t386;
t291 = Icges(3,1) * t323 + Icges(3,4) * t322 + Icges(3,5) * t387;
t290 = Icges(3,1) * t321 + Icges(3,4) * t320 - Icges(3,5) * t386;
t289 = Icges(3,4) * t323 + Icges(3,2) * t322 + Icges(3,6) * t387;
t288 = Icges(3,4) * t321 + Icges(3,2) * t320 - Icges(3,6) * t386;
t287 = Icges(3,5) * t323 + Icges(3,6) * t322 + Icges(3,3) * t387;
t286 = Icges(3,5) * t321 + Icges(3,6) * t320 - Icges(3,3) * t386;
t270 = rSges(4,1) * t304 - rSges(4,2) * t303 + rSges(4,3) * t319;
t269 = Icges(4,1) * t304 - Icges(4,4) * t303 + Icges(4,5) * t319;
t268 = Icges(4,4) * t304 - Icges(4,2) * t303 + Icges(4,6) * t319;
t267 = Icges(4,5) * t304 - Icges(4,6) * t303 + Icges(4,3) * t319;
t266 = pkin(5) * t303 + pkin(11) * t283;
t261 = t282 * t354 + t303 * t358;
t260 = t282 * t358 - t303 * t354;
t251 = -t292 * t344 + t316 * t333 + t371;
t250 = t293 * t344 - t316 * t334 + t369;
t249 = qJD(6) * t283 + t273;
t247 = rSges(5,1) * t283 - rSges(5,2) * t282 + rSges(5,3) * t303;
t246 = rSges(6,1) * t303 - rSges(6,2) * t283 + rSges(6,3) * t282;
t245 = rSges(4,1) * t281 - rSges(4,2) * t280 + rSges(4,3) * t306;
t244 = rSges(4,1) * t279 - rSges(4,2) * t278 + rSges(4,3) * t305;
t237 = Icges(4,1) * t281 - Icges(4,4) * t280 + Icges(4,5) * t306;
t236 = Icges(4,1) * t279 - Icges(4,4) * t278 + Icges(4,5) * t305;
t235 = Icges(4,4) * t281 - Icges(4,2) * t280 + Icges(4,6) * t306;
t234 = Icges(4,4) * t279 - Icges(4,2) * t278 + Icges(4,6) * t305;
t233 = Icges(4,5) * t281 - Icges(4,6) * t280 + Icges(4,3) * t306;
t232 = Icges(4,5) * t279 - Icges(4,6) * t278 + Icges(4,3) * t305;
t231 = t292 * t334 - t293 * t333 + t370;
t230 = pkin(5) * t280 + pkin(11) * t265;
t229 = pkin(5) * t278 + pkin(11) * t263;
t228 = t264 * t354 + t280 * t358;
t227 = t264 * t358 - t280 * t354;
t226 = t262 * t354 + t278 * t358;
t225 = t262 * t358 - t278 * t354;
t221 = qJD(6) * t265 + t257;
t220 = qJD(6) * t263 + t256;
t218 = rSges(5,1) * t265 - rSges(5,2) * t264 + rSges(5,3) * t280;
t217 = rSges(5,1) * t263 - rSges(5,2) * t262 + rSges(5,3) * t278;
t216 = rSges(7,1) * t261 + rSges(7,2) * t260 + rSges(7,3) * t283;
t215 = rSges(6,1) * t280 - rSges(6,2) * t265 + rSges(6,3) * t264;
t214 = rSges(6,1) * t278 - rSges(6,2) * t263 + rSges(6,3) * t262;
t209 = Icges(7,1) * t261 + Icges(7,4) * t260 + Icges(7,5) * t283;
t204 = Icges(7,4) * t261 + Icges(7,2) * t260 + Icges(7,6) * t283;
t199 = Icges(7,5) * t261 + Icges(7,6) * t260 + Icges(7,3) * t283;
t196 = rSges(7,1) * t228 + rSges(7,2) * t227 + rSges(7,3) * t265;
t195 = rSges(7,1) * t226 + rSges(7,2) * t225 + rSges(7,3) * t263;
t194 = Icges(7,1) * t228 + Icges(7,4) * t227 + Icges(7,5) * t265;
t193 = Icges(7,1) * t226 + Icges(7,4) * t225 + Icges(7,5) * t263;
t192 = Icges(7,4) * t228 + Icges(7,2) * t227 + Icges(7,6) * t265;
t191 = Icges(7,4) * t226 + Icges(7,2) * t225 + Icges(7,6) * t263;
t190 = Icges(7,5) * t228 + Icges(7,6) * t227 + Icges(7,3) * t265;
t189 = Icges(7,5) * t226 + Icges(7,6) * t225 + Icges(7,3) * t263;
t188 = -t244 * t307 + t270 * t294 + t368;
t187 = t245 * t307 - t270 * t295 + t366;
t186 = t244 * t295 - t245 * t294 + t367;
t185 = -t217 * t273 + t247 * t256 + t365;
t184 = t218 * t273 - t247 * t257 + t363;
t183 = t217 * t257 - t218 * t256 + t364;
t182 = t246 * t256 + (-t214 - t222) * t273 + t362;
t181 = t215 * t273 + (-t246 - t254) * t257 + t360;
t180 = t214 * t257 + (-t215 - t223) * t256 + t361;
t179 = -t195 * t249 + t216 * t220 + t256 * t266 + (-t222 - t229) * t273 + t362;
t178 = t196 * t249 - t216 * t221 + t230 * t273 + (-t254 - t266) * t257 + t360;
t177 = -t196 * t220 + t229 * t257 + (-t223 - t230) * t256 + t195 * t221 + t361;
t1 = t344 * ((t286 * t333 + t287 * t334 + t313 * t344) * t353 + ((t289 * t359 + t291 * t357) * t334 + (t288 * t359 + t290 * t357) * t333 + (t314 * t359 + t315 * t357) * t344) * t351) / 0.2e1 + m(1) * (t330 ^ 2 + t331 ^ 2 + t332 ^ 2) / 0.2e1 + t307 * ((t233 * t319 - t235 * t303 + t237 * t304) * t295 + (t232 * t319 - t234 * t303 + t236 * t304) * t294 + (t267 * t319 - t268 * t303 + t269 * t304) * t307) / 0.2e1 + t294 * ((t233 * t305 - t235 * t278 + t237 * t279) * t295 + (t232 * t305 - t234 * t278 + t236 * t279) * t294 + (t267 * t305 - t268 * t278 + t269 * t279) * t307) / 0.2e1 + t295 * ((t233 * t306 - t235 * t280 + t237 * t281) * t295 + (t232 * t306 - t234 * t280 + t236 * t281) * t294 + (t267 * t306 - t268 * t280 + t269 * t281) * t307) / 0.2e1 + m(2) * (t302 ^ 2 + t309 ^ 2 + t310 ^ 2) / 0.2e1 + t249 * ((t190 * t283 + t192 * t260 + t194 * t261) * t221 + (t189 * t283 + t191 * t260 + t193 * t261) * t220 + (t283 * t199 + t260 * t204 + t261 * t209) * t249) / 0.2e1 + t221 * ((t265 * t190 + t227 * t192 + t228 * t194) * t221 + (t189 * t265 + t191 * t227 + t193 * t228) * t220 + (t199 * t265 + t204 * t227 + t209 * t228) * t249) / 0.2e1 + t220 * ((t190 * t263 + t192 * t225 + t194 * t226) * t221 + (t263 * t189 + t225 * t191 + t226 * t193) * t220 + (t199 * t263 + t204 * t225 + t209 * t226) * t249) / 0.2e1 + m(3) * (t231 ^ 2 + t250 ^ 2 + t251 ^ 2) / 0.2e1 + m(4) * (t186 ^ 2 + t187 ^ 2 + t188 ^ 2) / 0.2e1 + m(5) * (t183 ^ 2 + t184 ^ 2 + t185 ^ 2) / 0.2e1 + m(6) * (t180 ^ 2 + t181 ^ 2 + t182 ^ 2) / 0.2e1 + m(7) * (t177 ^ 2 + t178 ^ 2 + t179 ^ 2) / 0.2e1 + t333 * ((-t287 * t386 + t289 * t320 + t291 * t321) * t334 + (-t286 * t386 + t288 * t320 + t290 * t321) * t333 + (-t313 * t386 + t314 * t320 + t315 * t321) * t344) / 0.2e1 + t334 * ((t287 * t387 + t289 * t322 + t291 * t323) * t334 + (t286 * t387 + t288 * t322 + t290 * t323) * t333 + (t313 * t387 + t314 * t322 + t315 * t323) * t344) / 0.2e1 + ((t262 * t400 + t263 * t398 + t278 * t399) * t273 + (t262 * t405 + t263 * t401 + t278 * t403) * t257 + (t262 * t406 + t263 * t402 + t278 * t404) * t256) * t256 / 0.2e1 + ((t264 * t400 + t265 * t398 + t280 * t399) * t273 + (t264 * t405 + t265 * t401 + t280 * t403) * t257 + (t264 * t406 + t265 * t402 + t280 * t404) * t256) * t257 / 0.2e1 + ((t282 * t400 + t283 * t398 + t303 * t399) * t273 + (t282 * t405 + t283 * t401 + t303 * t403) * t257 + (t282 * t406 + t283 * t402 + t303 * t404) * t256) * t273 / 0.2e1 + ((-t337 * t350 + t339 * t352 + Icges(1,4)) * V_base(5) + (-t338 * t350 + t340 * t352 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t337 * t352 + t339 * t350 + Icges(1,2)) * V_base(5) + (t338 * t352 + t340 * t350 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t350 + Icges(2,6) * t352 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t352 - Icges(2,6) * t350 + Icges(1,5)) * V_base(4) + (Icges(2,3) / 0.2e1 + Icges(1,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
