% Calculate kinetic energy for
% S6PPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPRRRP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PPRRRP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRP1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRRRP1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:52:14
% EndTime: 2019-03-08 18:52:18
% DurationCPUTime: 4.11s
% Computational Cost: add. (5487->401), mult. (14586->567), div. (0->0), fcn. (18742->14), ass. (0->180)
t419 = Icges(6,1) + Icges(7,1);
t418 = Icges(6,4) + Icges(7,4);
t417 = Icges(6,5) + Icges(7,5);
t416 = Icges(6,2) + Icges(7,2);
t415 = Icges(6,6) + Icges(7,6);
t414 = Icges(6,3) + Icges(7,3);
t413 = rSges(7,3) + qJ(6);
t343 = sin(pkin(11));
t346 = cos(pkin(11));
t412 = Icges(2,5) * t346 - Icges(2,6) * t343 + Icges(1,5);
t411 = Icges(2,5) * t343 + Icges(2,6) * t346 + Icges(1,6);
t342 = sin(pkin(12));
t345 = cos(pkin(12));
t347 = cos(pkin(6));
t383 = t346 * t347;
t314 = -t342 * t343 + t345 * t383;
t315 = t342 * t383 + t343 * t345;
t351 = sin(qJ(3));
t344 = sin(pkin(6));
t393 = sin(pkin(7));
t372 = t344 * t393;
t394 = cos(pkin(7));
t398 = cos(qJ(3));
t273 = t315 * t398 + (t314 * t394 - t346 * t372) * t351;
t350 = sin(qJ(4));
t373 = t344 * t394;
t364 = -t314 * t393 - t346 * t373;
t397 = cos(qJ(4));
t260 = t273 * t397 + t350 * t364;
t369 = t398 * t393;
t366 = t344 * t369;
t370 = t394 * t398;
t272 = -t314 * t370 + t315 * t351 + t346 * t366;
t349 = sin(qJ(5));
t352 = cos(qJ(5));
t231 = -t260 * t349 + t272 * t352;
t390 = t272 * t349;
t232 = t260 * t352 + t390;
t259 = t273 * t350 - t364 * t397;
t410 = t231 * t415 + t232 * t417 + t259 * t414;
t385 = t343 * t347;
t316 = -t342 * t346 - t345 * t385;
t317 = -t342 * t385 + t345 * t346;
t275 = t317 * t398 + (t316 * t394 + t343 * t372) * t351;
t363 = -t316 * t393 + t343 * t373;
t262 = t275 * t397 + t350 * t363;
t274 = -t316 * t370 + t317 * t351 - t343 * t366;
t233 = -t262 * t349 + t274 * t352;
t389 = t274 * t349;
t234 = t262 * t352 + t389;
t261 = t275 * t350 - t363 * t397;
t409 = t233 * t415 + t234 * t417 + t261 * t414;
t408 = t231 * t416 + t232 * t418 + t259 * t415;
t407 = t233 * t416 + t234 * t418 + t261 * t415;
t406 = t418 * t231 + t232 * t419 + t417 * t259;
t405 = t418 * t233 + t234 * t419 + t417 * t261;
t299 = t347 * t393 * t351 + (t345 * t351 * t394 + t342 * t398) * t344;
t362 = -t345 * t372 + t347 * t394;
t278 = t299 * t397 + t350 * t362;
t387 = t342 * t344;
t298 = -t344 * t345 * t370 - t347 * t369 + t351 * t387;
t263 = -t278 * t349 + t298 * t352;
t388 = t298 * t349;
t264 = t278 * t352 + t388;
t277 = t299 * t350 - t362 * t397;
t404 = t263 * t415 + t264 * t417 + t277 * t414;
t403 = t263 * t416 + t264 * t418 + t277 * t415;
t402 = t418 * t263 + t264 * t419 + t417 * t277;
t396 = pkin(5) * t352;
t392 = Icges(2,4) * t343;
t391 = qJ(2) * t347;
t386 = t343 * t344;
t384 = t344 * t346;
t382 = rSges(7,1) * t232 + rSges(7,2) * t231 + pkin(5) * t390 + t413 * t259 + t260 * t396;
t381 = rSges(7,1) * t234 + rSges(7,2) * t233 + pkin(5) * t389 + t413 * t261 + t262 * t396;
t380 = rSges(7,1) * t264 + rSges(7,2) * t263 + pkin(5) * t388 + t413 * t277 + t278 * t396;
t379 = qJD(2) * t344;
t378 = V_base(5) * qJ(1) + V_base(1);
t374 = qJD(1) + V_base(3);
t292 = qJD(3) * t363 + V_base(4);
t291 = qJD(3) * t364 + V_base(5);
t307 = qJD(3) * t362 + V_base(6);
t371 = -qJ(1) - t391;
t258 = qJD(4) * t274 + t292;
t257 = qJD(4) * t272 + t291;
t276 = qJD(4) * t298 + t307;
t368 = t343 * t379 + V_base(5) * t391 + t378;
t320 = pkin(1) * t343 - qJ(2) * t384;
t367 = qJD(2) * t347 + V_base(4) * t320 + t374;
t321 = pkin(1) * t346 + qJ(2) * t386;
t365 = V_base(6) * t321 - t346 * t379 + V_base(2);
t281 = t315 * pkin(2) + pkin(8) * t364;
t302 = pkin(2) * t387 + pkin(8) * t362;
t361 = V_base(5) * t302 + (-t281 - t320) * V_base(6) + t368;
t282 = t317 * pkin(2) + pkin(8) * t363;
t360 = V_base(4) * t281 + (-t282 - t321) * V_base(5) + t367;
t250 = pkin(3) * t273 + pkin(9) * t272;
t269 = pkin(3) * t299 + pkin(9) * t298;
t359 = -t250 * t307 + t291 * t269 + t361;
t251 = pkin(3) * t275 + pkin(9) * t274;
t358 = t292 * t250 - t251 * t291 + t360;
t357 = V_base(6) * t282 + (-t302 + t371) * V_base(4) + t365;
t227 = pkin(4) * t260 + pkin(10) * t259;
t252 = pkin(4) * t278 + pkin(10) * t277;
t356 = -t227 * t276 + t257 * t252 + t359;
t228 = pkin(4) * t262 + pkin(10) * t261;
t355 = t258 * t227 - t228 * t257 + t358;
t354 = t307 * t251 - t269 * t292 + t357;
t353 = t276 * t228 - t252 * t258 + t354;
t340 = Icges(2,4) * t346;
t334 = rSges(2,1) * t346 - rSges(2,2) * t343;
t333 = rSges(2,1) * t343 + rSges(2,2) * t346;
t332 = Icges(2,1) * t346 - t392;
t331 = Icges(2,1) * t343 + t340;
t330 = -Icges(2,2) * t343 + t340;
t329 = Icges(2,2) * t346 + t392;
t326 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t325 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t324 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t311 = rSges(3,3) * t347 + (rSges(3,1) * t342 + rSges(3,2) * t345) * t344;
t310 = Icges(3,5) * t347 + (Icges(3,1) * t342 + Icges(3,4) * t345) * t344;
t309 = Icges(3,6) * t347 + (Icges(3,4) * t342 + Icges(3,2) * t345) * t344;
t308 = Icges(3,3) * t347 + (Icges(3,5) * t342 + Icges(3,6) * t345) * t344;
t304 = V_base(5) * rSges(2,3) - t333 * V_base(6) + t378;
t303 = t334 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t300 = t333 * V_base(4) - t334 * V_base(5) + t374;
t290 = rSges(3,1) * t317 + rSges(3,2) * t316 + rSges(3,3) * t386;
t289 = rSges(3,1) * t315 + rSges(3,2) * t314 - rSges(3,3) * t384;
t288 = Icges(3,1) * t317 + Icges(3,4) * t316 + Icges(3,5) * t386;
t287 = Icges(3,1) * t315 + Icges(3,4) * t314 - Icges(3,5) * t384;
t286 = Icges(3,4) * t317 + Icges(3,2) * t316 + Icges(3,6) * t386;
t285 = Icges(3,4) * t315 + Icges(3,2) * t314 - Icges(3,6) * t384;
t284 = Icges(3,5) * t317 + Icges(3,6) * t316 + Icges(3,3) * t386;
t283 = Icges(3,5) * t315 + Icges(3,6) * t314 - Icges(3,3) * t384;
t268 = t299 * rSges(4,1) - t298 * rSges(4,2) + rSges(4,3) * t362;
t267 = Icges(4,1) * t299 - Icges(4,4) * t298 + Icges(4,5) * t362;
t266 = Icges(4,4) * t299 - Icges(4,2) * t298 + Icges(4,6) * t362;
t265 = Icges(4,5) * t299 - Icges(4,6) * t298 + Icges(4,3) * t362;
t255 = t311 * V_base(5) + (-t289 - t320) * V_base(6) + t368;
t254 = t290 * V_base(6) + (-t311 + t371) * V_base(4) + t365;
t253 = qJD(5) * t277 + t276;
t248 = t289 * V_base(4) + (-t290 - t321) * V_base(5) + t367;
t247 = rSges(5,1) * t278 - rSges(5,2) * t277 + rSges(5,3) * t298;
t246 = Icges(5,1) * t278 - Icges(5,4) * t277 + Icges(5,5) * t298;
t245 = Icges(5,4) * t278 - Icges(5,2) * t277 + Icges(5,6) * t298;
t244 = Icges(5,5) * t278 - Icges(5,6) * t277 + Icges(5,3) * t298;
t243 = t275 * rSges(4,1) - t274 * rSges(4,2) + rSges(4,3) * t363;
t242 = t273 * rSges(4,1) - t272 * rSges(4,2) + rSges(4,3) * t364;
t241 = Icges(4,1) * t275 - Icges(4,4) * t274 + Icges(4,5) * t363;
t240 = Icges(4,1) * t273 - Icges(4,4) * t272 + Icges(4,5) * t364;
t239 = Icges(4,4) * t275 - Icges(4,2) * t274 + Icges(4,6) * t363;
t238 = Icges(4,4) * t273 - Icges(4,2) * t272 + Icges(4,6) * t364;
t237 = Icges(4,5) * t275 - Icges(4,6) * t274 + Icges(4,3) * t363;
t236 = Icges(4,5) * t273 - Icges(4,6) * t272 + Icges(4,3) * t364;
t230 = qJD(5) * t261 + t258;
t229 = qJD(5) * t259 + t257;
t224 = rSges(6,1) * t264 + rSges(6,2) * t263 + rSges(6,3) * t277;
t216 = rSges(5,1) * t262 - rSges(5,2) * t261 + rSges(5,3) * t274;
t215 = rSges(5,1) * t260 - rSges(5,2) * t259 + rSges(5,3) * t272;
t214 = Icges(5,1) * t262 - Icges(5,4) * t261 + Icges(5,5) * t274;
t213 = Icges(5,1) * t260 - Icges(5,4) * t259 + Icges(5,5) * t272;
t212 = Icges(5,4) * t262 - Icges(5,2) * t261 + Icges(5,6) * t274;
t211 = Icges(5,4) * t260 - Icges(5,2) * t259 + Icges(5,6) * t272;
t210 = Icges(5,5) * t262 - Icges(5,6) * t261 + Icges(5,3) * t274;
t209 = Icges(5,5) * t260 - Icges(5,6) * t259 + Icges(5,3) * t272;
t206 = rSges(6,1) * t234 + rSges(6,2) * t233 + rSges(6,3) * t261;
t204 = rSges(6,1) * t232 + rSges(6,2) * t231 + rSges(6,3) * t259;
t190 = -t242 * t307 + t268 * t291 + t361;
t189 = t243 * t307 - t268 * t292 + t357;
t186 = t242 * t292 - t243 * t291 + t360;
t185 = -t215 * t276 + t247 * t257 + t359;
t184 = t216 * t276 - t247 * t258 + t354;
t183 = t215 * t258 - t216 * t257 + t358;
t182 = -t204 * t253 + t224 * t229 + t356;
t181 = t206 * t253 - t224 * t230 + t353;
t180 = t204 * t230 - t206 * t229 + t355;
t179 = qJD(6) * t261 + t229 * t380 - t253 * t382 + t356;
t178 = qJD(6) * t259 - t230 * t380 + t253 * t381 + t353;
t177 = qJD(6) * t277 - t229 * t381 + t230 * t382 + t355;
t1 = m(1) * (t324 ^ 2 + t325 ^ 2 + t326 ^ 2) / 0.2e1 + m(4) * (t186 ^ 2 + t189 ^ 2 + t190 ^ 2) / 0.2e1 + m(5) * (t183 ^ 2 + t184 ^ 2 + t185 ^ 2) / 0.2e1 + m(6) * (t180 ^ 2 + t181 ^ 2 + t182 ^ 2) / 0.2e1 + m(7) * (t177 ^ 2 + t178 ^ 2 + t179 ^ 2) / 0.2e1 + t276 * ((t210 * t298 - t212 * t277 + t214 * t278) * t258 + (t209 * t298 - t211 * t277 + t213 * t278) * t257 + (t244 * t298 - t245 * t277 + t246 * t278) * t276) / 0.2e1 + m(2) * (t300 ^ 2 + t303 ^ 2 + t304 ^ 2) / 0.2e1 + m(3) * (t248 ^ 2 + t254 ^ 2 + t255 ^ 2) / 0.2e1 + t307 * ((t237 * t362 - t298 * t239 + t299 * t241) * t292 + (t236 * t362 - t298 * t238 + t299 * t240) * t291 + (t362 * t265 - t298 * t266 + t299 * t267) * t307) / 0.2e1 + t292 * ((t363 * t237 - t274 * t239 + t275 * t241) * t292 + (t236 * t363 - t274 * t238 + t275 * t240) * t291 + (t265 * t363 - t274 * t266 + t275 * t267) * t307) / 0.2e1 + t291 * ((t237 * t364 - t272 * t239 + t273 * t241) * t292 + (t364 * t236 - t272 * t238 + t273 * t240) * t291 + (t265 * t364 - t272 * t266 + t273 * t267) * t307) / 0.2e1 + t257 * ((t210 * t272 - t212 * t259 + t214 * t260) * t258 + (t209 * t272 - t211 * t259 + t213 * t260) * t257 + (t244 * t272 - t245 * t259 + t246 * t260) * t276) / 0.2e1 + t258 * ((t210 * t274 - t212 * t261 + t214 * t262) * t258 + (t209 * t274 - t211 * t261 + t213 * t262) * t257 + (t244 * t274 - t245 * t261 + t246 * t262) * t276) / 0.2e1 + ((t403 * t231 + t402 * t232 + t404 * t259) * t253 + (t407 * t231 + t405 * t232 + t409 * t259) * t230 + (t408 * t231 + t406 * t232 + t410 * t259) * t229) * t229 / 0.2e1 + ((t403 * t233 + t402 * t234 + t404 * t261) * t253 + (t407 * t233 + t405 * t234 + t409 * t261) * t230 + (t408 * t233 + t406 * t234 + t410 * t261) * t229) * t230 / 0.2e1 + ((t403 * t263 + t402 * t264 + t404 * t277) * t253 + (t407 * t263 + t405 * t264 + t409 * t277) * t230 + (t408 * t263 + t406 * t264 + t410 * t277) * t229) * t253 / 0.2e1 + ((t308 * t386 + t309 * t316 + t310 * t317 + t412) * V_base(6) + (t283 * t386 + t285 * t316 + t287 * t317 - t329 * t343 + t331 * t346 + Icges(1,4)) * V_base(5) + (t284 * t386 + t286 * t316 + t288 * t317 - t330 * t343 + t332 * t346 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t308 * t384 + t309 * t314 + t310 * t315 + t411) * V_base(6) + (-t283 * t384 + t285 * t314 + t287 * t315 + t329 * t346 + t331 * t343 + Icges(1,2)) * V_base(5) + (-t284 * t384 + t286 * t314 + t288 * t315 + t330 * t346 + t332 * t343 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t308 * t347 + (t309 * t345 + t310 * t342) * t344 + Icges(2,3) + Icges(1,3)) * V_base(6) + (t283 * t347 + (t285 * t345 + t287 * t342) * t344 + t411) * V_base(5) + (t284 * t347 + (t286 * t345 + t288 * t342) * t344 + t412) * V_base(4)) * V_base(6) / 0.2e1;
T  = t1;
