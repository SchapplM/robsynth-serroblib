% Calculate kinetic energy for
% S6PRRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRP5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRRP5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRP5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:20:05
% EndTime: 2019-03-09 00:20:10
% DurationCPUTime: 5.03s
% Computational Cost: add. (5622->396), mult. (14811->583), div. (0->0), fcn. (18967->14), ass. (0->180)
t421 = Icges(6,1) + Icges(7,1);
t420 = Icges(6,4) + Icges(7,4);
t419 = Icges(6,5) + Icges(7,5);
t418 = Icges(6,2) + Icges(7,2);
t417 = Icges(6,6) + Icges(7,6);
t416 = Icges(6,3) + Icges(7,3);
t415 = rSges(7,3) + qJ(6);
t345 = sin(pkin(12));
t347 = cos(pkin(12));
t353 = sin(qJ(2));
t348 = cos(pkin(6));
t355 = cos(qJ(2));
t385 = t348 * t355;
t314 = -t345 * t353 + t347 * t385;
t386 = t348 * t353;
t315 = t345 * t355 + t347 * t386;
t352 = sin(qJ(3));
t346 = sin(pkin(6));
t394 = sin(pkin(7));
t374 = t346 * t394;
t395 = cos(pkin(7));
t400 = cos(qJ(3));
t276 = t315 * t400 + (t314 * t395 - t347 * t374) * t352;
t351 = sin(qJ(4));
t375 = t346 * t395;
t368 = -t314 * t394 - t347 * t375;
t399 = cos(qJ(4));
t262 = t276 * t399 + t351 * t368;
t372 = t400 * t394;
t371 = t346 * t372;
t373 = t395 * t400;
t275 = -t314 * t373 + t315 * t352 + t347 * t371;
t350 = sin(qJ(5));
t354 = cos(qJ(5));
t232 = -t262 * t350 + t275 * t354;
t392 = t275 * t350;
t233 = t262 * t354 + t392;
t261 = t276 * t351 - t368 * t399;
t414 = t232 * t417 + t233 * t419 + t261 * t416;
t316 = -t345 * t385 - t347 * t353;
t317 = -t345 * t386 + t347 * t355;
t278 = t317 * t400 + (t316 * t395 + t345 * t374) * t352;
t367 = -t316 * t394 + t345 * t375;
t264 = t278 * t399 + t351 * t367;
t277 = -t316 * t373 + t317 * t352 - t345 * t371;
t234 = -t264 * t350 + t277 * t354;
t391 = t277 * t350;
t235 = t264 * t354 + t391;
t263 = t278 * t351 - t367 * t399;
t413 = t234 * t417 + t235 * t419 + t263 * t416;
t412 = t232 * t418 + t233 * t420 + t261 * t417;
t411 = t234 * t418 + t235 * t420 + t263 * t417;
t410 = t420 * t232 + t233 * t421 + t419 * t261;
t409 = t420 * t234 + t235 * t421 + t419 * t263;
t301 = t348 * t394 * t352 + (t352 * t355 * t395 + t353 * t400) * t346;
t366 = t348 * t395 - t355 * t374;
t280 = t301 * t399 + t351 * t366;
t387 = t346 * t353;
t300 = -t346 * t355 * t373 - t348 * t372 + t352 * t387;
t259 = -t280 * t350 + t300 * t354;
t390 = t300 * t350;
t260 = t280 * t354 + t390;
t279 = t301 * t351 - t366 * t399;
t408 = t259 * t417 + t260 * t419 + t279 * t416;
t407 = t259 * t418 + t260 * t420 + t279 * t417;
t406 = t420 * t259 + t260 * t421 + t419 * t279;
t398 = pkin(8) * t348;
t397 = pkin(5) * t354;
t393 = Icges(2,4) * t345;
t389 = t345 * t346;
t388 = t346 * t347;
t384 = rSges(7,1) * t233 + rSges(7,2) * t232 + pkin(5) * t392 + t261 * t415 + t262 * t397;
t383 = rSges(7,1) * t235 + rSges(7,2) * t234 + pkin(5) * t391 + t263 * t415 + t264 * t397;
t382 = rSges(7,1) * t260 + rSges(7,2) * t259 + pkin(5) * t390 + t279 * t415 + t280 * t397;
t381 = qJD(2) * t346;
t380 = V_base(5) * qJ(1) + V_base(1);
t376 = qJD(1) + V_base(3);
t328 = t345 * t381 + V_base(4);
t338 = qJD(2) * t348 + V_base(6);
t292 = qJD(3) * t367 + t328;
t302 = qJD(3) * t366 + t338;
t258 = qJD(4) * t277 + t292;
t271 = qJD(4) * t300 + t302;
t327 = -t347 * t381 + V_base(5);
t291 = qJD(3) * t368 + t327;
t320 = pkin(1) * t345 - pkin(8) * t388;
t370 = -t320 * V_base(6) + V_base(5) * t398 + t380;
t321 = pkin(1) * t347 + pkin(8) * t389;
t369 = V_base(4) * t320 - t321 * V_base(5) + t376;
t257 = qJD(4) * t275 + t291;
t365 = V_base(6) * t321 + V_base(2) + (-qJ(1) - t398) * V_base(4);
t281 = t315 * pkin(2) + pkin(9) * t368;
t303 = pkin(2) * t387 + pkin(9) * t366;
t364 = -t281 * t338 + t327 * t303 + t370;
t282 = t317 * pkin(2) + pkin(9) * t367;
t363 = t328 * t281 - t282 * t327 + t369;
t362 = t338 * t282 - t303 * t328 + t365;
t253 = pkin(3) * t276 + pkin(10) * t275;
t269 = pkin(3) * t301 + pkin(10) * t300;
t361 = -t253 * t302 + t291 * t269 + t364;
t254 = pkin(3) * t278 + pkin(10) * t277;
t360 = t292 * t253 - t254 * t291 + t363;
t359 = t302 * t254 - t269 * t292 + t362;
t229 = pkin(4) * t262 + pkin(11) * t261;
t255 = pkin(4) * t280 + pkin(11) * t279;
t358 = -t229 * t271 + t257 * t255 + t361;
t230 = pkin(4) * t264 + pkin(11) * t263;
t357 = t258 * t229 - t230 * t257 + t360;
t356 = t271 * t230 - t255 * t258 + t359;
t343 = Icges(2,4) * t347;
t336 = rSges(2,1) * t347 - rSges(2,2) * t345;
t335 = rSges(2,1) * t345 + rSges(2,2) * t347;
t334 = Icges(2,1) * t347 - t393;
t333 = Icges(2,1) * t345 + t343;
t332 = -Icges(2,2) * t345 + t343;
t331 = Icges(2,2) * t347 + t393;
t326 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t325 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t324 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t311 = t348 * rSges(3,3) + (rSges(3,1) * t353 + rSges(3,2) * t355) * t346;
t310 = Icges(3,5) * t348 + (Icges(3,1) * t353 + Icges(3,4) * t355) * t346;
t309 = Icges(3,6) * t348 + (Icges(3,4) * t353 + Icges(3,2) * t355) * t346;
t308 = Icges(3,3) * t348 + (Icges(3,5) * t353 + Icges(3,6) * t355) * t346;
t305 = V_base(5) * rSges(2,3) - t335 * V_base(6) + t380;
t304 = t336 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t299 = t335 * V_base(4) - t336 * V_base(5) + t376;
t290 = rSges(3,1) * t317 + rSges(3,2) * t316 + rSges(3,3) * t389;
t289 = rSges(3,1) * t315 + rSges(3,2) * t314 - rSges(3,3) * t388;
t288 = Icges(3,1) * t317 + Icges(3,4) * t316 + Icges(3,5) * t389;
t287 = Icges(3,1) * t315 + Icges(3,4) * t314 - Icges(3,5) * t388;
t286 = Icges(3,4) * t317 + Icges(3,2) * t316 + Icges(3,6) * t389;
t285 = Icges(3,4) * t315 + Icges(3,2) * t314 - Icges(3,6) * t388;
t284 = Icges(3,5) * t317 + Icges(3,6) * t316 + Icges(3,3) * t389;
t283 = Icges(3,5) * t315 + Icges(3,6) * t314 - Icges(3,3) * t388;
t268 = t301 * rSges(4,1) - t300 * rSges(4,2) + rSges(4,3) * t366;
t267 = Icges(4,1) * t301 - Icges(4,4) * t300 + Icges(4,5) * t366;
t266 = Icges(4,4) * t301 - Icges(4,2) * t300 + Icges(4,6) * t366;
t265 = Icges(4,5) * t301 - Icges(4,6) * t300 + Icges(4,3) * t366;
t252 = -t289 * t338 + t311 * t327 + t370;
t251 = t290 * t338 - t311 * t328 + t365;
t250 = qJD(5) * t279 + t271;
t248 = rSges(5,1) * t280 - rSges(5,2) * t279 + rSges(5,3) * t300;
t247 = t278 * rSges(4,1) - t277 * rSges(4,2) + rSges(4,3) * t367;
t246 = t276 * rSges(4,1) - t275 * rSges(4,2) + rSges(4,3) * t368;
t245 = Icges(5,1) * t280 - Icges(5,4) * t279 + Icges(5,5) * t300;
t244 = Icges(5,4) * t280 - Icges(5,2) * t279 + Icges(5,6) * t300;
t243 = Icges(5,5) * t280 - Icges(5,6) * t279 + Icges(5,3) * t300;
t242 = Icges(4,1) * t278 - Icges(4,4) * t277 + Icges(4,5) * t367;
t241 = Icges(4,1) * t276 - Icges(4,4) * t275 + Icges(4,5) * t368;
t240 = Icges(4,4) * t278 - Icges(4,2) * t277 + Icges(4,6) * t367;
t239 = Icges(4,4) * t276 - Icges(4,2) * t275 + Icges(4,6) * t368;
t238 = Icges(4,5) * t278 - Icges(4,6) * t277 + Icges(4,3) * t367;
t237 = Icges(4,5) * t276 - Icges(4,6) * t275 + Icges(4,3) * t368;
t236 = t289 * t328 - t290 * t327 + t369;
t228 = qJD(5) * t263 + t258;
t227 = qJD(5) * t261 + t257;
t225 = rSges(6,1) * t260 + rSges(6,2) * t259 + rSges(6,3) * t279;
t223 = rSges(5,1) * t264 - rSges(5,2) * t263 + rSges(5,3) * t277;
t222 = rSges(5,1) * t262 - rSges(5,2) * t261 + rSges(5,3) * t275;
t221 = Icges(5,1) * t264 - Icges(5,4) * t263 + Icges(5,5) * t277;
t220 = Icges(5,1) * t262 - Icges(5,4) * t261 + Icges(5,5) * t275;
t217 = Icges(5,4) * t264 - Icges(5,2) * t263 + Icges(5,6) * t277;
t216 = Icges(5,4) * t262 - Icges(5,2) * t261 + Icges(5,6) * t275;
t213 = Icges(5,5) * t264 - Icges(5,6) * t263 + Icges(5,3) * t277;
t212 = Icges(5,5) * t262 - Icges(5,6) * t261 + Icges(5,3) * t275;
t206 = rSges(6,1) * t235 + rSges(6,2) * t234 + rSges(6,3) * t263;
t204 = rSges(6,1) * t233 + rSges(6,2) * t232 + rSges(6,3) * t261;
t190 = -t246 * t302 + t268 * t291 + t364;
t189 = t247 * t302 - t268 * t292 + t362;
t186 = t246 * t292 - t247 * t291 + t363;
t185 = -t222 * t271 + t248 * t257 + t361;
t184 = t223 * t271 - t248 * t258 + t359;
t183 = t222 * t258 - t223 * t257 + t360;
t182 = -t204 * t250 + t225 * t227 + t358;
t181 = t206 * t250 - t225 * t228 + t356;
t180 = t204 * t228 - t206 * t227 + t357;
t179 = qJD(6) * t263 + t227 * t382 - t250 * t384 + t358;
t178 = qJD(6) * t261 - t228 * t382 + t250 * t383 + t356;
t177 = qJD(6) * t279 - t227 * t383 + t228 * t384 + t357;
t1 = m(5) * (t183 ^ 2 + t184 ^ 2 + t185 ^ 2) / 0.2e1 + t291 * ((t238 * t368 - t275 * t240 + t276 * t242) * t292 + (t237 * t368 - t275 * t239 + t276 * t241) * t291 + (t265 * t368 - t275 * t266 + t276 * t267) * t302) / 0.2e1 + t328 * ((t284 * t389 + t286 * t316 + t288 * t317) * t328 + (t283 * t389 + t285 * t316 + t287 * t317) * t327 + (t308 * t389 + t309 * t316 + t310 * t317) * t338) / 0.2e1 + m(6) * (t180 ^ 2 + t181 ^ 2 + t182 ^ 2) / 0.2e1 + t271 * ((t213 * t300 - t217 * t279 + t221 * t280) * t258 + (t212 * t300 - t216 * t279 + t220 * t280) * t257 + (t243 * t300 - t244 * t279 + t245 * t280) * t271) / 0.2e1 + m(2) * (t299 ^ 2 + t304 ^ 2 + t305 ^ 2) / 0.2e1 + t327 * ((-t284 * t388 + t286 * t314 + t288 * t315) * t328 + (-t283 * t388 + t285 * t314 + t287 * t315) * t327 + (-t308 * t388 + t309 * t314 + t310 * t315) * t338) / 0.2e1 + m(7) * (t177 ^ 2 + t178 ^ 2 + t179 ^ 2) / 0.2e1 + m(1) * (t324 ^ 2 + t325 ^ 2 + t326 ^ 2) / 0.2e1 + t258 * ((t213 * t277 - t217 * t263 + t221 * t264) * t258 + (t212 * t277 - t216 * t263 + t220 * t264) * t257 + (t243 * t277 - t244 * t263 + t245 * t264) * t271) / 0.2e1 + t257 * ((t213 * t275 - t217 * t261 + t221 * t262) * t258 + (t212 * t275 - t216 * t261 + t220 * t262) * t257 + (t243 * t275 - t244 * t261 + t245 * t262) * t271) / 0.2e1 + t338 * ((t283 * t327 + t284 * t328 + t308 * t338) * t348 + ((t286 * t355 + t288 * t353) * t328 + (t285 * t355 + t287 * t353) * t327 + (t309 * t355 + t310 * t353) * t338) * t346) / 0.2e1 + m(3) * (t236 ^ 2 + t251 ^ 2 + t252 ^ 2) / 0.2e1 + t302 * ((t238 * t366 - t300 * t240 + t301 * t242) * t292 + (t237 * t366 - t300 * t239 + t301 * t241) * t291 + (t265 * t366 - t300 * t266 + t301 * t267) * t302) / 0.2e1 + t292 * ((t238 * t367 - t277 * t240 + t278 * t242) * t292 + (t237 * t367 - t277 * t239 + t278 * t241) * t291 + (t265 * t367 - t277 * t266 + t278 * t267) * t302) / 0.2e1 + m(4) * (t186 ^ 2 + t189 ^ 2 + t190 ^ 2) / 0.2e1 + ((t232 * t407 + t233 * t406 + t261 * t408) * t250 + (t232 * t411 + t233 * t409 + t261 * t413) * t228 + (t412 * t232 + t410 * t233 + t414 * t261) * t227) * t227 / 0.2e1 + ((t234 * t407 + t235 * t406 + t263 * t408) * t250 + (t411 * t234 + t409 * t235 + t413 * t263) * t228 + (t412 * t234 + t410 * t235 + t263 * t414) * t227) * t228 / 0.2e1 + ((t407 * t259 + t406 * t260 + t408 * t279) * t250 + (t259 * t411 + t260 * t409 + t279 * t413) * t228 + (t412 * t259 + t410 * t260 + t279 * t414) * t227) * t250 / 0.2e1 + ((-t331 * t345 + t333 * t347 + Icges(1,4)) * V_base(5) + (-t332 * t345 + t334 * t347 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t331 * t347 + t333 * t345 + Icges(1,2)) * V_base(5) + (t332 * t347 + t334 * t345 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t347 - Icges(2,6) * t345 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t345 + Icges(2,6) * t347 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
