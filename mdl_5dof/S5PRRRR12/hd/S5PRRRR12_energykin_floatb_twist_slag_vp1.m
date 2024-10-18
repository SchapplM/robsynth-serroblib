% Calculate kinetic energy for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% m [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR12_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR12_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRRR12_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR12_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR12_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR12_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:07:10
% EndTime: 2024-09-28 18:07:13
% DurationCPUTime: 0.75s
% Computational Cost: add. (3472->440), mult. (3964->638), div. (0->0), fcn. (4118->26), ass. (0->221)
t374 = pkin(7) + pkin(8);
t361 = sin(pkin(6));
t422 = pkin(10) * t361;
t365 = cos(pkin(5));
t421 = t365 * pkin(7);
t372 = cos(qJ(3));
t343 = pkin(3) * t372 + pkin(2);
t420 = -pkin(2) + t343;
t373 = cos(qJ(2));
t344 = t373 * pkin(2) + pkin(1);
t360 = sin(pkin(11));
t419 = Icges(2,4) * t360;
t359 = qJ(2) + qJ(3);
t356 = qJ(4) + t359;
t341 = sin(t356);
t418 = t341 * t360;
t362 = sin(pkin(5));
t417 = t360 * t362;
t416 = t360 * t365;
t415 = t361 * t362;
t363 = cos(pkin(11));
t414 = t361 * t363;
t413 = t362 * t363;
t364 = cos(pkin(6));
t412 = t362 * t364;
t411 = t363 * t365;
t366 = sin(qJ(5));
t410 = t364 * t366;
t370 = cos(qJ(5));
t409 = t364 * t370;
t408 = t365 * t366;
t369 = sin(qJ(2));
t407 = t365 * t369;
t406 = t365 * t370;
t405 = t365 * t373;
t358 = pkin(9) + t374;
t368 = sin(qJ(3));
t396 = t373 * t368 * pkin(3);
t268 = -t358 * t362 + (t369 * t343 + t396) * t365;
t397 = pkin(2) * t407;
t304 = -t362 * t374 + t397;
t404 = t268 - t304;
t352 = cos(t359);
t310 = pkin(3) * t352 + t344;
t403 = t310 - t344;
t402 = qJD(2) * t362;
t401 = -qJD(2) - qJD(3);
t400 = V_base(5) * qJ(1) + V_base(1);
t399 = t360 * t422;
t398 = pkin(10) * t414;
t392 = qJD(1) + V_base(3);
t391 = t366 * t415;
t390 = t370 * t415;
t389 = t364 * t408;
t388 = t364 * t406;
t315 = t360 * t402 + V_base(4);
t326 = qJD(2) * t365 + V_base(6);
t350 = pkin(5) - t359;
t349 = pkin(5) + t359;
t387 = pkin(7) * t362 + t304;
t283 = qJD(3) * t417 + t315;
t305 = qJD(3) * t365 + t326;
t270 = qJD(4) * t417 + t283;
t284 = qJD(4) * t365 + t305;
t330 = pkin(10) * t364 + t358;
t386 = -t330 * t362 - t268 + t397;
t353 = t360 * pkin(1);
t307 = -pkin(7) * t413 + t353;
t385 = -t307 * V_base(6) + V_base(5) * t421 + t400;
t296 = pkin(4) * t363 + t365 * t399;
t298 = pkin(4) * t416 - t398;
t367 = sin(qJ(4));
t371 = cos(qJ(4));
t384 = pkin(3) * t363 + t296 * t371 - t298 * t367;
t297 = -pkin(4) * t360 + t365 * t398;
t299 = pkin(4) * t411 + t399;
t383 = pkin(3) * t360 - t297 * t371 + t299 * t367;
t354 = t363 * pkin(1);
t308 = pkin(7) * t417 + t354;
t382 = V_base(4) * t307 - t308 * V_base(5) + t392;
t269 = V_base(5) + (-qJD(4) + t401) * t413;
t381 = V_base(6) * t308 + V_base(2) + (-qJ(1) - t421) * V_base(4);
t252 = t344 * t360 + t363 * t387 - t353;
t285 = t362 * t369 * pkin(2) + (-pkin(7) + t374) * t365;
t314 = -t363 * t402 + V_base(5);
t380 = -t252 * t326 + t314 * t285 + t385;
t253 = t344 * t363 - t360 * t387 - t354;
t379 = t315 * t252 - t253 * t314 + t382;
t378 = t326 * t253 - t285 * t315 + t381;
t209 = t360 * t403 + t363 * t404;
t249 = (t358 - t374) * t365 + (t369 * t420 + t396) * t362;
t282 = t401 * t413 + V_base(5);
t377 = -t209 * t305 + t282 * t249 + t380;
t210 = -t360 * t404 + t363 * t403;
t376 = t283 * t209 - t210 * t282 + t379;
t375 = t305 * t210 - t249 * t283 + t378;
t351 = sin(t359);
t348 = Icges(2,4) * t363;
t342 = cos(t356);
t340 = -qJ(4) + t350;
t339 = qJ(4) + t349;
t338 = cos(t349);
t337 = sin(t350);
t335 = cos(t339);
t334 = sin(t340);
t332 = cos(t350) / 0.2e1;
t331 = sin(t349) / 0.2e1;
t325 = cos(t340) / 0.2e1;
t324 = sin(t339) / 0.2e1;
t323 = rSges(2,1) * t363 - rSges(2,2) * t360;
t322 = rSges(2,1) * t360 + rSges(2,2) * t363;
t321 = Icges(2,1) * t363 - t419;
t320 = Icges(2,1) * t360 + t348;
t319 = -Icges(2,2) * t360 + t348;
t318 = Icges(2,2) * t363 + t419;
t313 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t312 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t311 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t309 = -pkin(4) * t367 + t371 * t422;
t306 = pkin(4) * t371 + t367 * t422 + pkin(3);
t303 = t332 - t338 / 0.2e1;
t302 = t332 + t338 / 0.2e1;
t301 = t331 - t337 / 0.2e1;
t300 = t331 + t337 / 0.2e1;
t293 = -t360 * t407 + t363 * t373;
t292 = -t360 * t405 - t363 * t369;
t291 = t360 * t373 + t363 * t407;
t290 = -t360 * t369 + t363 * t405;
t289 = t325 - t335 / 0.2e1;
t288 = t325 + t335 / 0.2e1;
t287 = t324 - t334 / 0.2e1;
t286 = t324 + t334 / 0.2e1;
t281 = -t342 * t415 + t364 * t365;
t280 = rSges(3,3) * t365 + (rSges(3,1) * t369 + rSges(3,2) * t373) * t362;
t279 = Icges(3,5) * t365 + (Icges(3,1) * t369 + Icges(3,4) * t373) * t362;
t278 = Icges(3,6) * t365 + (Icges(3,4) * t369 + Icges(3,2) * t373) * t362;
t277 = Icges(3,3) * t365 + (Icges(3,5) * t369 + Icges(3,6) * t373) * t362;
t276 = V_base(5) * rSges(2,3) - t322 * V_base(6) + t400;
t275 = t323 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t274 = -t301 * t360 + t352 * t363;
t273 = -t302 * t360 - t351 * t363;
t272 = t301 * t363 + t352 * t360;
t271 = t302 * t363 - t351 * t360;
t266 = -t287 * t360 + t342 * t363;
t265 = -t288 * t360 - t341 * t363;
t264 = t287 * t363 + t342 * t360;
t263 = t288 * t363 - t418;
t262 = t322 * V_base(4) - t323 * V_base(5) + t392;
t261 = t361 * t408 + (t341 * t370 + t342 * t410) * t362;
t260 = t361 * t406 + (-t341 * t366 + t342 * t409) * t362;
t259 = rSges(4,1) * t303 + rSges(4,2) * t300 + rSges(4,3) * t365;
t258 = Icges(4,1) * t303 + Icges(4,4) * t300 + Icges(4,5) * t365;
t257 = Icges(4,4) * t303 + Icges(4,2) * t300 + Icges(4,6) * t365;
t256 = Icges(4,5) * t303 + Icges(4,6) * t300 + Icges(4,3) * t365;
t255 = -t363 * t412 + (-t342 * t411 + t418) * t361;
t254 = t341 * t414 + (t342 * t361 * t365 + t412) * t360;
t251 = qJD(5) * t281 + t284;
t250 = rSges(5,1) * t289 + rSges(5,2) * t286 + rSges(5,3) * t365;
t248 = Icges(5,1) * t289 + Icges(5,4) * t286 + Icges(5,5) * t365;
t247 = Icges(5,4) * t289 + Icges(5,2) * t286 + Icges(5,6) * t365;
t246 = Icges(5,5) * t289 + Icges(5,6) * t286 + Icges(5,3) * t365;
t245 = rSges(3,1) * t293 + rSges(3,2) * t292 + rSges(3,3) * t417;
t244 = rSges(3,1) * t291 + rSges(3,2) * t290 - rSges(3,3) * t413;
t243 = Icges(3,1) * t293 + Icges(3,4) * t292 + Icges(3,5) * t417;
t242 = Icges(3,1) * t291 + Icges(3,4) * t290 - Icges(3,5) * t413;
t241 = Icges(3,4) * t293 + Icges(3,2) * t292 + Icges(3,6) * t417;
t240 = Icges(3,4) * t291 + Icges(3,2) * t290 - Icges(3,6) * t413;
t239 = Icges(3,5) * t293 + Icges(3,6) * t292 + Icges(3,3) * t417;
t238 = Icges(3,5) * t291 + Icges(3,6) * t290 - Icges(3,3) * t413;
t237 = pkin(3) * t411 + t297 * t367 + t299 * t371;
t236 = -pkin(3) * t416 - t296 * t367 - t298 * t371;
t232 = (-t360 * t366 + t363 * t388) * t342 + (-t360 * t409 - t363 * t408) * t341 - t363 * t390;
t231 = (t360 * t370 + t363 * t389) * t342 + (-t360 * t410 + t363 * t406) * t341 - t363 * t391;
t230 = (-t360 * t388 - t363 * t366) * t342 + (t360 * t408 - t363 * t409) * t341 + t360 * t390;
t229 = (-t360 * t389 + t363 * t370) * t342 + (-t360 * t406 - t363 * t410) * t341 + t360 * t391;
t228 = rSges(4,1) * t274 + rSges(4,2) * t273 + rSges(4,3) * t417;
t227 = rSges(4,1) * t272 + rSges(4,2) * t271 - rSges(4,3) * t413;
t226 = Icges(4,1) * t274 + Icges(4,4) * t273 + Icges(4,5) * t417;
t225 = Icges(4,1) * t272 + Icges(4,4) * t271 - Icges(4,5) * t413;
t224 = Icges(4,4) * t274 + Icges(4,2) * t273 + Icges(4,6) * t417;
t223 = Icges(4,4) * t272 + Icges(4,2) * t271 - Icges(4,6) * t413;
t222 = Icges(4,5) * t274 + Icges(4,6) * t273 + Icges(4,3) * t417;
t221 = Icges(4,5) * t272 + Icges(4,6) * t271 - Icges(4,3) * t413;
t220 = qJD(5) * t254 + t270;
t219 = qJD(5) * t255 + t269;
t218 = rSges(5,1) * t266 + rSges(5,2) * t265 + rSges(5,3) * t417;
t217 = rSges(5,1) * t264 + rSges(5,2) * t263 - rSges(5,3) * t413;
t216 = Icges(5,1) * t266 + Icges(5,4) * t265 + Icges(5,5) * t417;
t215 = Icges(5,1) * t264 + Icges(5,4) * t263 - Icges(5,5) * t413;
t214 = Icges(5,4) * t266 + Icges(5,2) * t265 + Icges(5,6) * t417;
t213 = Icges(5,4) * t264 + Icges(5,2) * t263 - Icges(5,6) * t413;
t212 = Icges(5,5) * t266 + Icges(5,6) * t265 + Icges(5,3) * t417;
t211 = Icges(5,5) * t264 + Icges(5,6) * t263 - Icges(5,3) * t413;
t206 = rSges(6,1) * t261 + rSges(6,2) * t260 + rSges(6,3) * t281;
t205 = Icges(6,1) * t261 + Icges(6,4) * t260 + Icges(6,5) * t281;
t204 = Icges(6,4) * t261 + Icges(6,2) * t260 + Icges(6,6) * t281;
t203 = Icges(6,5) * t261 + Icges(6,6) * t260 + Icges(6,3) * t281;
t202 = -t244 * t326 + t280 * t314 + t385;
t201 = t245 * t326 - t280 * t315 + t381;
t200 = (-t358 + t330) * t365 + ((-t309 * t372 + (-pkin(3) + t306) * t368) * t373 + (t306 * t372 + t309 * t368 - t420) * t369) * t362;
t199 = t244 * t315 - t245 * t314 + t382;
t198 = rSges(6,1) * t231 + rSges(6,2) * t232 + rSges(6,3) * t255;
t197 = rSges(6,1) * t229 + rSges(6,2) * t230 + rSges(6,3) * t254;
t196 = Icges(6,1) * t231 + Icges(6,4) * t232 + Icges(6,5) * t255;
t195 = Icges(6,1) * t229 + Icges(6,4) * t230 + Icges(6,5) * t254;
t194 = Icges(6,4) * t231 + Icges(6,2) * t232 + Icges(6,6) * t255;
t193 = Icges(6,4) * t229 + Icges(6,2) * t230 + Icges(6,6) * t254;
t192 = Icges(6,5) * t231 + Icges(6,6) * t232 + Icges(6,3) * t255;
t191 = Icges(6,5) * t229 + Icges(6,6) * t230 + Icges(6,3) * t254;
t190 = -t227 * t305 + t259 * t282 + t380;
t189 = t228 * t305 - t259 * t283 + t378;
t188 = -t363 * t310 + (t363 * pkin(2) + t236 * t368 + t372 * t384) * t373 + (t236 * t372 - t368 * t384) * t369 + t354 - t386 * t360;
t187 = -t360 * t310 + (t360 * pkin(2) + t237 * t368 + t372 * t383) * t373 + (t237 * t372 - t368 * t383) * t369 + t353 + t386 * t363;
t186 = t227 * t283 - t228 * t282 + t379;
t185 = -t217 * t284 + t250 * t269 + t377;
t184 = t218 * t284 - t250 * t270 + t375;
t183 = t217 * t270 - t218 * t269 + t376;
t182 = -t187 * t284 - t198 * t251 + t200 * t269 + t206 * t219 + t377;
t181 = t188 * t284 + t197 * t251 - t200 * t270 - t206 * t220 + t375;
t180 = t187 * t270 - t188 * t269 - t197 * t219 + t198 * t220 + t376;
t1 = m(1) * (t311 ^ 2 + t312 ^ 2 + t313 ^ 2) / 0.2e1 + m(2) * (t262 ^ 2 + t275 ^ 2 + t276 ^ 2) / 0.2e1 + m(3) * (t199 ^ 2 + t201 ^ 2 + t202 ^ 2) / 0.2e1 + t315 * ((t239 * t417 + t241 * t292 + t243 * t293) * t315 + (t238 * t417 + t240 * t292 + t242 * t293) * t314 + (t277 * t417 + t278 * t292 + t279 * t293) * t326) / 0.2e1 + t314 * ((-t239 * t413 + t241 * t290 + t243 * t291) * t315 + (-t238 * t413 + t240 * t290 + t242 * t291) * t314 + (-t277 * t413 + t278 * t290 + t279 * t291) * t326) / 0.2e1 + t326 * ((t238 * t314 + t239 * t315 + t277 * t326) * t365 + ((t241 * t373 + t243 * t369) * t315 + (t240 * t373 + t242 * t369) * t314 + (t278 * t373 + t279 * t369) * t326) * t362) / 0.2e1 + m(4) * (t186 ^ 2 + t189 ^ 2 + t190 ^ 2) / 0.2e1 + t283 * ((t222 * t417 + t273 * t224 + t274 * t226) * t283 + (t221 * t417 + t223 * t273 + t225 * t274) * t282 + (t256 * t417 + t257 * t273 + t258 * t274) * t305) / 0.2e1 + t282 * ((-t222 * t413 + t224 * t271 + t226 * t272) * t283 + (-t221 * t413 + t271 * t223 + t272 * t225) * t282 + (-t256 * t413 + t257 * t271 + t258 * t272) * t305) / 0.2e1 + t305 * ((t222 * t365 + t224 * t300 + t226 * t303) * t283 + (t221 * t365 + t223 * t300 + t225 * t303) * t282 + (t256 * t365 + t257 * t300 + t258 * t303) * t305) / 0.2e1 + m(5) * (t183 ^ 2 + t184 ^ 2 + t185 ^ 2) / 0.2e1 + t270 * ((t212 * t417 + t265 * t214 + t266 * t216) * t270 + (t211 * t417 + t213 * t265 + t215 * t266) * t269 + (t246 * t417 + t247 * t265 + t248 * t266) * t284) / 0.2e1 + t269 * ((-t212 * t413 + t214 * t263 + t216 * t264) * t270 + (-t211 * t413 + t263 * t213 + t264 * t215) * t269 + (-t246 * t413 + t247 * t263 + t248 * t264) * t284) / 0.2e1 + t284 * ((t212 * t365 + t214 * t286 + t216 * t289) * t270 + (t211 * t365 + t213 * t286 + t215 * t289) * t269 + (t365 * t246 + t286 * t247 + t289 * t248) * t284) / 0.2e1 + m(6) * (t180 ^ 2 + t181 ^ 2 + t182 ^ 2) / 0.2e1 + t220 * ((t254 * t191 + t230 * t193 + t229 * t195) * t220 + (t192 * t254 + t194 * t230 + t196 * t229) * t219 + (t203 * t254 + t204 * t230 + t205 * t229) * t251) / 0.2e1 + t219 * ((t191 * t255 + t193 * t232 + t195 * t231) * t220 + (t255 * t192 + t232 * t194 + t231 * t196) * t219 + (t203 * t255 + t204 * t232 + t205 * t231) * t251) / 0.2e1 + t251 * ((t191 * t281 + t193 * t260 + t195 * t261) * t220 + (t192 * t281 + t194 * t260 + t196 * t261) * t219 + (t281 * t203 + t260 * t204 + t261 * t205) * t251) / 0.2e1 + ((-t318 * t360 + t320 * t363 + Icges(1,4)) * V_base(5) + (-t319 * t360 + t321 * t363 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t318 * t363 + t320 * t360 + Icges(1,2)) * V_base(5) + (t319 * t363 + t321 * t360 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t360 + Icges(2,6) * t363 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t363 - Icges(2,6) * t360 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t1;
