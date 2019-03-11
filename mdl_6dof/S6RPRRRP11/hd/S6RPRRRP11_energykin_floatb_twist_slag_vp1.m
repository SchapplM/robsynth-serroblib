% Calculate kinetic energy for
% S6RPRRRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP11_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP11_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRP11_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP11_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP11_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP11_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:34:22
% EndTime: 2019-03-09 06:34:27
% DurationCPUTime: 5.00s
% Computational Cost: add. (5547->401), mult. (14586->573), div. (0->0), fcn. (18742->14), ass. (0->181)
t420 = Icges(6,1) + Icges(7,1);
t419 = Icges(6,4) + Icges(7,4);
t418 = Icges(6,5) + Icges(7,5);
t417 = Icges(6,2) + Icges(7,2);
t416 = Icges(6,6) + Icges(7,6);
t415 = Icges(6,3) + Icges(7,3);
t414 = rSges(7,3) + qJ(6);
t347 = cos(pkin(12));
t345 = sin(pkin(12));
t353 = sin(qJ(1));
t385 = t353 * t345;
t348 = cos(pkin(6));
t355 = cos(qJ(1));
t386 = t348 * t355;
t317 = t347 * t386 - t385;
t384 = t353 * t347;
t318 = t345 * t386 + t384;
t352 = sin(qJ(3));
t346 = sin(pkin(6));
t395 = sin(pkin(7));
t374 = t346 * t395;
t396 = cos(pkin(7));
t400 = cos(qJ(3));
t279 = t318 * t400 + (t317 * t396 - t355 * t374) * t352;
t351 = sin(qJ(4));
t375 = t346 * t396;
t367 = -t317 * t395 - t355 * t375;
t399 = cos(qJ(4));
t264 = t279 * t399 + t351 * t367;
t371 = t400 * t395;
t369 = t346 * t371;
t372 = t396 * t400;
t278 = -t317 * t372 + t318 * t352 + t355 * t369;
t350 = sin(qJ(5));
t354 = cos(qJ(5));
t233 = -t264 * t350 + t278 * t354;
t392 = t278 * t350;
t234 = t264 * t354 + t392;
t263 = t279 * t351 - t367 * t399;
t413 = t233 * t416 + t234 * t418 + t263 * t415;
t319 = -t345 * t355 - t348 * t384;
t320 = t347 * t355 - t348 * t385;
t281 = t320 * t400 + (t319 * t396 + t353 * t374) * t352;
t366 = -t319 * t395 + t353 * t375;
t266 = t281 * t399 + t351 * t366;
t280 = -t319 * t372 + t320 * t352 - t353 * t369;
t235 = -t266 * t350 + t280 * t354;
t391 = t280 * t350;
t236 = t266 * t354 + t391;
t265 = t281 * t351 - t366 * t399;
t412 = t235 * t416 + t236 * t418 + t265 * t415;
t411 = t233 * t417 + t234 * t419 + t263 * t416;
t410 = t235 * t417 + t236 * t419 + t265 * t416;
t409 = t233 * t419 + t234 * t420 + t263 * t418;
t408 = t235 * t419 + t236 * t420 + t265 * t418;
t301 = t348 * t395 * t352 + (t347 * t352 * t396 + t345 * t400) * t346;
t365 = -t347 * t374 + t348 * t396;
t276 = t301 * t399 + t351 * t365;
t389 = t345 * t346;
t300 = -t346 * t347 * t372 - t348 * t371 + t352 * t389;
t259 = -t276 * t350 + t300 * t354;
t390 = t300 * t350;
t260 = t276 * t354 + t390;
t275 = t301 * t351 - t365 * t399;
t407 = t259 * t416 + t260 * t418 + t275 * t415;
t406 = t259 * t417 + t260 * t419 + t275 * t416;
t405 = t259 * t419 + t260 * t420 + t275 * t418;
t398 = pkin(5) * t354;
t394 = Icges(2,4) * t353;
t393 = qJ(2) * t348;
t388 = t346 * t353;
t387 = t346 * t355;
t383 = rSges(7,1) * t234 + rSges(7,2) * t233 + pkin(5) * t392 + t263 * t414 + t264 * t398;
t382 = rSges(7,1) * t236 + rSges(7,2) * t235 + pkin(5) * t391 + t265 * t414 + t266 * t398;
t381 = rSges(7,1) * t260 + rSges(7,2) * t259 + pkin(5) * t390 + t275 * t414 + t276 * t398;
t380 = qJD(2) * t346;
t379 = V_base(5) * pkin(8) + V_base(1);
t294 = qJD(3) * t366 + V_base(4);
t293 = qJD(3) * t367 + V_base(5);
t342 = V_base(6) + qJD(1);
t376 = -pkin(8) - t393;
t322 = pkin(1) * t353 - qJ(2) * t387;
t373 = qJD(2) * t348 + t322 * V_base(4) + V_base(3);
t262 = qJD(4) * t280 + t294;
t261 = qJD(4) * t278 + t293;
t307 = qJD(3) * t365 + t342;
t370 = t353 * t380 + t393 * V_base(5) + t379;
t272 = qJD(4) * t300 + t307;
t323 = pkin(1) * t355 + qJ(2) * t388;
t368 = t323 * t342 - t355 * t380 + V_base(2);
t283 = pkin(2) * t318 + pkin(9) * t367;
t304 = pkin(2) * t389 + pkin(9) * t365;
t364 = V_base(5) * t304 + (-t283 - t322) * t342 + t370;
t284 = pkin(2) * t320 + pkin(9) * t366;
t363 = V_base(4) * t283 + (-t284 - t323) * V_base(5) + t373;
t256 = pkin(3) * t279 + pkin(10) * t278;
t271 = pkin(3) * t301 + pkin(10) * t300;
t362 = -t256 * t307 + t271 * t293 + t364;
t257 = pkin(3) * t281 + pkin(10) * t280;
t361 = t256 * t294 - t257 * t293 + t363;
t360 = t342 * t284 + (-t304 + t376) * V_base(4) + t368;
t229 = pkin(4) * t264 + pkin(11) * t263;
t252 = pkin(4) * t276 + pkin(11) * t275;
t359 = -t229 * t272 + t252 * t261 + t362;
t230 = pkin(4) * t266 + pkin(11) * t265;
t358 = t229 * t262 - t230 * t261 + t361;
t357 = t257 * t307 - t271 * t294 + t360;
t356 = t230 * t272 - t252 * t262 + t357;
t343 = Icges(2,4) * t355;
t336 = rSges(2,1) * t355 - rSges(2,2) * t353;
t335 = rSges(2,1) * t353 + rSges(2,2) * t355;
t334 = Icges(2,1) * t355 - t394;
t333 = Icges(2,1) * t353 + t343;
t332 = -Icges(2,2) * t353 + t343;
t331 = Icges(2,2) * t355 + t394;
t330 = Icges(2,5) * t355 - Icges(2,6) * t353;
t329 = Icges(2,5) * t353 + Icges(2,6) * t355;
t328 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t327 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t326 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t313 = rSges(3,3) * t348 + (rSges(3,1) * t345 + rSges(3,2) * t347) * t346;
t312 = Icges(3,5) * t348 + (Icges(3,1) * t345 + Icges(3,4) * t347) * t346;
t311 = Icges(3,6) * t348 + (Icges(3,4) * t345 + Icges(3,2) * t347) * t346;
t310 = Icges(3,3) * t348 + (Icges(3,5) * t345 + Icges(3,6) * t347) * t346;
t306 = V_base(5) * rSges(2,3) - t335 * t342 + t379;
t305 = t336 * t342 + V_base(2) + (-rSges(2,3) - pkin(8)) * V_base(4);
t303 = t335 * V_base(4) - t336 * V_base(5) + V_base(3);
t292 = rSges(3,1) * t320 + rSges(3,2) * t319 + rSges(3,3) * t388;
t291 = rSges(3,1) * t318 + rSges(3,2) * t317 - rSges(3,3) * t387;
t290 = Icges(3,1) * t320 + Icges(3,4) * t319 + Icges(3,5) * t388;
t289 = Icges(3,1) * t318 + Icges(3,4) * t317 - Icges(3,5) * t387;
t288 = Icges(3,4) * t320 + Icges(3,2) * t319 + Icges(3,6) * t388;
t287 = Icges(3,4) * t318 + Icges(3,2) * t317 - Icges(3,6) * t387;
t286 = Icges(3,5) * t320 + Icges(3,6) * t319 + Icges(3,3) * t388;
t285 = Icges(3,5) * t318 + Icges(3,6) * t317 - Icges(3,3) * t387;
t270 = rSges(4,1) * t301 - rSges(4,2) * t300 + rSges(4,3) * t365;
t269 = Icges(4,1) * t301 - Icges(4,4) * t300 + Icges(4,5) * t365;
t268 = Icges(4,4) * t301 - Icges(4,2) * t300 + Icges(4,6) * t365;
t267 = Icges(4,5) * t301 - Icges(4,6) * t300 + Icges(4,3) * t365;
t255 = t313 * V_base(5) + (-t291 - t322) * t342 + t370;
t254 = t342 * t292 + (-t313 + t376) * V_base(4) + t368;
t253 = qJD(5) * t275 + t272;
t251 = t291 * V_base(4) + (-t292 - t323) * V_base(5) + t373;
t249 = rSges(4,1) * t281 - rSges(4,2) * t280 + rSges(4,3) * t366;
t248 = rSges(4,1) * t279 - rSges(4,2) * t278 + rSges(4,3) * t367;
t247 = Icges(4,1) * t281 - Icges(4,4) * t280 + Icges(4,5) * t366;
t246 = Icges(4,1) * t279 - Icges(4,4) * t278 + Icges(4,5) * t367;
t245 = Icges(4,4) * t281 - Icges(4,2) * t280 + Icges(4,6) * t366;
t244 = Icges(4,4) * t279 - Icges(4,2) * t278 + Icges(4,6) * t367;
t243 = Icges(4,5) * t281 - Icges(4,6) * t280 + Icges(4,3) * t366;
t242 = Icges(4,5) * t279 - Icges(4,6) * t278 + Icges(4,3) * t367;
t241 = rSges(5,1) * t276 - rSges(5,2) * t275 + rSges(5,3) * t300;
t239 = Icges(5,1) * t276 - Icges(5,4) * t275 + Icges(5,5) * t300;
t238 = Icges(5,4) * t276 - Icges(5,2) * t275 + Icges(5,6) * t300;
t237 = Icges(5,5) * t276 - Icges(5,6) * t275 + Icges(5,3) * t300;
t232 = qJD(5) * t265 + t262;
t231 = qJD(5) * t263 + t261;
t227 = rSges(5,1) * t266 - rSges(5,2) * t265 + rSges(5,3) * t280;
t226 = rSges(5,1) * t264 - rSges(5,2) * t263 + rSges(5,3) * t278;
t225 = Icges(5,1) * t266 - Icges(5,4) * t265 + Icges(5,5) * t280;
t224 = Icges(5,1) * t264 - Icges(5,4) * t263 + Icges(5,5) * t278;
t223 = Icges(5,4) * t266 - Icges(5,2) * t265 + Icges(5,6) * t280;
t222 = Icges(5,4) * t264 - Icges(5,2) * t263 + Icges(5,6) * t278;
t221 = Icges(5,5) * t266 - Icges(5,6) * t265 + Icges(5,3) * t280;
t220 = Icges(5,5) * t264 - Icges(5,6) * t263 + Icges(5,3) * t278;
t218 = rSges(6,1) * t260 + rSges(6,2) * t259 + rSges(6,3) * t275;
t208 = rSges(6,1) * t236 + rSges(6,2) * t235 + rSges(6,3) * t265;
t206 = rSges(6,1) * t234 + rSges(6,2) * t233 + rSges(6,3) * t263;
t192 = -t248 * t307 + t270 * t293 + t364;
t191 = t307 * t249 - t294 * t270 + t360;
t188 = t248 * t294 - t249 * t293 + t363;
t187 = -t226 * t272 + t241 * t261 + t362;
t186 = t227 * t272 - t241 * t262 + t357;
t185 = t226 * t262 - t227 * t261 + t361;
t184 = -t206 * t253 + t218 * t231 + t359;
t183 = t208 * t253 - t218 * t232 + t356;
t182 = t206 * t232 - t208 * t231 + t358;
t181 = qJD(6) * t265 + t231 * t381 - t253 * t383 + t359;
t180 = qJD(6) * t263 - t232 * t381 + t253 * t382 + t356;
t179 = qJD(6) * t275 - t231 * t382 + t232 * t383 + t358;
t1 = m(1) * (t326 ^ 2 + t327 ^ 2 + t328 ^ 2) / 0.2e1 + m(2) * (t303 ^ 2 + t305 ^ 2 + t306 ^ 2) / 0.2e1 + t272 * ((t221 * t300 - t223 * t275 + t225 * t276) * t262 + (t220 * t300 - t222 * t275 + t224 * t276) * t261 + (t237 * t300 - t238 * t275 + t239 * t276) * t272) / 0.2e1 + t261 * ((t221 * t278 - t223 * t263 + t225 * t264) * t262 + (t220 * t278 - t222 * t263 + t224 * t264) * t261 + (t237 * t278 - t238 * t263 + t239 * t264) * t272) / 0.2e1 + t262 * ((t221 * t280 - t223 * t265 + t225 * t266) * t262 + (t220 * t280 - t222 * t265 + t224 * t266) * t261 + (t237 * t280 - t238 * t265 + t239 * t266) * t272) / 0.2e1 + m(3) * (t251 ^ 2 + t254 ^ 2 + t255 ^ 2) / 0.2e1 + t307 * ((t243 * t365 - t300 * t245 + t301 * t247) * t294 + (t242 * t365 - t300 * t244 + t301 * t246) * t293 + (t267 * t365 - t300 * t268 + t301 * t269) * t307) / 0.2e1 + t294 * ((t243 * t366 - t280 * t245 + t281 * t247) * t294 + (t242 * t366 - t280 * t244 + t281 * t246) * t293 + (t267 * t366 - t280 * t268 + t281 * t269) * t307) / 0.2e1 + t293 * ((t243 * t367 - t278 * t245 + t279 * t247) * t294 + (t242 * t367 - t278 * t244 + t279 * t246) * t293 + (t267 * t367 - t278 * t268 + t279 * t269) * t307) / 0.2e1 + m(7) * (t179 ^ 2 + t180 ^ 2 + t181 ^ 2) / 0.2e1 + m(6) * (t182 ^ 2 + t183 ^ 2 + t184 ^ 2) / 0.2e1 + m(5) * (t185 ^ 2 + t186 ^ 2 + t187 ^ 2) / 0.2e1 + m(4) * (t188 ^ 2 + t191 ^ 2 + t192 ^ 2) / 0.2e1 + ((t233 * t406 + t234 * t405 + t263 * t407) * t253 + (t233 * t410 + t234 * t408 + t263 * t412) * t232 + (t411 * t233 + t409 * t234 + t413 * t263) * t231) * t231 / 0.2e1 + ((t235 * t406 + t236 * t405 + t265 * t407) * t253 + (t410 * t235 + t408 * t236 + t412 * t265) * t232 + (t235 * t411 + t236 * t409 + t265 * t413) * t231) * t232 / 0.2e1 + ((t406 * t259 + t405 * t260 + t407 * t275) * t253 + (t259 * t410 + t260 * t408 + t275 * t412) * t232 + (t259 * t411 + t260 * t409 + t275 * t413) * t231) * t253 / 0.2e1 + ((t285 * V_base(5) + t286 * V_base(4) + t310 * t342) * t348 + ((t288 * t347 + t290 * t345) * V_base(4) + (t287 * t347 + t289 * t345) * V_base(5) + (t311 * t347 + t312 * t345) * t342) * t346 + Icges(2,3) * t342 + t329 * V_base(5) + t330 * V_base(4)) * t342 / 0.2e1 + ((t310 * t388 + t311 * t319 + t312 * t320 + t330) * t342 + (t285 * t388 + t287 * t319 + t289 * t320 - t331 * t353 + t333 * t355 + Icges(1,4)) * V_base(5) + (t286 * t388 + t288 * t319 + t290 * t320 - t353 * t332 + t334 * t355 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t310 * t387 + t311 * t317 + t312 * t318 + t329) * t342 + (-t285 * t387 + t317 * t287 + t318 * t289 + t331 * t355 + t353 * t333 + Icges(1,2)) * V_base(5) + (-t286 * t387 + t288 * t317 + t290 * t318 + t332 * t355 + t334 * t353 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
