% Calculate kinetic energy for
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRR8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRPRR8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR8_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR8_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRR8_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:35:55
% EndTime: 2019-03-08 22:35:59
% DurationCPUTime: 4.83s
% Computational Cost: add. (4366->401), mult. (11404->584), div. (0->0), fcn. (14382->14), ass. (0->180)
t416 = Icges(4,1) + Icges(5,2);
t415 = Icges(5,1) + Icges(4,3);
t414 = -Icges(4,4) - Icges(5,6);
t413 = -Icges(5,4) + Icges(4,5);
t412 = Icges(5,5) - Icges(4,6);
t411 = Icges(4,2) + Icges(5,3);
t351 = sin(pkin(12));
t353 = cos(pkin(12));
t358 = sin(qJ(2));
t354 = cos(pkin(6));
t360 = cos(qJ(2));
t386 = t354 * t360;
t320 = -t351 * t358 + t353 * t386;
t387 = t354 * t358;
t321 = t351 * t360 + t353 * t387;
t357 = sin(qJ(3));
t352 = sin(pkin(6));
t392 = sin(pkin(7));
t396 = cos(qJ(3));
t374 = t396 * t392;
t373 = t352 * t374;
t393 = cos(pkin(7));
t375 = t393 * t396;
t278 = -t320 * t375 + t321 * t357 + t353 * t373;
t376 = t357 * t392;
t377 = t357 * t393;
t389 = t352 * t353;
t279 = t320 * t377 + t321 * t396 - t376 * t389;
t379 = t352 * t393;
t306 = -t320 * t392 - t353 * t379;
t408 = t411 * t278 + t414 * t279 + t412 * t306;
t322 = -t351 * t386 - t353 * t358;
t323 = -t351 * t387 + t353 * t360;
t280 = -t322 * t375 + t323 * t357 - t351 * t373;
t378 = t352 * t392;
t281 = t323 * t396 + (t322 * t393 + t351 * t378) * t357;
t307 = -t322 * t392 + t351 * t379;
t407 = t411 * t280 + t414 * t281 + t412 * t307;
t406 = t412 * t278 + t413 * t279 + t415 * t306;
t405 = t412 * t280 + t413 * t281 + t415 * t307;
t404 = t414 * t278 + t416 * t279 + t413 * t306;
t403 = t414 * t280 + t416 * t281 + t413 * t307;
t388 = t352 * t358;
t304 = -t352 * t360 * t375 - t354 * t374 + t357 * t388;
t305 = t354 * t376 + (t358 * t396 + t360 * t377) * t352;
t319 = t354 * t393 - t360 * t378;
t402 = t411 * t304 + t414 * t305 + t412 * t319;
t401 = t412 * t304 + t413 * t305 + t415 * t319;
t400 = t414 * t304 + t416 * t305 + t413 * t319;
t395 = cos(qJ(5));
t394 = pkin(8) * t354;
t391 = Icges(2,4) * t351;
t390 = t351 * t352;
t385 = qJD(2) * t352;
t384 = V_base(5) * qJ(1) + V_base(1);
t380 = qJD(1) + V_base(3);
t335 = t351 * t385 + V_base(4);
t345 = qJD(2) * t354 + V_base(6);
t296 = qJD(3) * t307 + t335;
t308 = qJD(3) * t319 + t345;
t250 = qJD(5) * t281 + t296;
t270 = qJD(5) * t305 + t308;
t334 = -t353 * t385 + V_base(5);
t295 = qJD(3) * t306 + t334;
t326 = pkin(1) * t351 - pkin(8) * t389;
t372 = -t326 * V_base(6) + V_base(5) * t394 + t384;
t327 = pkin(1) * t353 + pkin(8) * t390;
t371 = V_base(4) * t326 - t327 * V_base(5) + t380;
t249 = qJD(5) * t279 + t295;
t370 = V_base(6) * t327 + V_base(2) + (-qJ(1) - t394) * V_base(4);
t285 = t321 * pkin(2) + pkin(9) * t306;
t309 = pkin(2) * t388 + pkin(9) * t319;
t369 = -t285 * t345 + t334 * t309 + t372;
t286 = t323 * pkin(2) + pkin(9) * t307;
t368 = t335 * t285 - t286 * t334 + t371;
t367 = t345 * t286 - t309 * t335 + t370;
t268 = pkin(3) * t305 + qJ(4) * t304;
t366 = qJD(4) * t280 + t295 * t268 + t369;
t245 = pkin(3) * t279 + qJ(4) * t278;
t365 = qJD(4) * t304 + t296 * t245 + t368;
t246 = pkin(3) * t281 + qJ(4) * t280;
t364 = qJD(4) * t278 + t308 * t246 + t367;
t258 = pkin(4) * t306 + pkin(10) * t279;
t284 = pkin(4) * t319 + pkin(10) * t305;
t363 = t295 * t284 + (-t245 - t258) * t308 + t366;
t259 = pkin(4) * t307 + pkin(10) * t281;
t362 = t296 * t258 + (-t246 - t259) * t295 + t365;
t361 = t308 * t259 + (-t268 - t284) * t296 + t364;
t359 = cos(qJ(6));
t356 = sin(qJ(5));
t355 = sin(qJ(6));
t349 = Icges(2,4) * t353;
t343 = rSges(2,1) * t353 - rSges(2,2) * t351;
t342 = rSges(2,1) * t351 + rSges(2,2) * t353;
t341 = Icges(2,1) * t353 - t391;
t340 = Icges(2,1) * t351 + t349;
t339 = -Icges(2,2) * t351 + t349;
t338 = Icges(2,2) * t353 + t391;
t333 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t332 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t331 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t317 = t354 * rSges(3,3) + (rSges(3,1) * t358 + rSges(3,2) * t360) * t352;
t316 = Icges(3,5) * t354 + (Icges(3,1) * t358 + Icges(3,4) * t360) * t352;
t315 = Icges(3,6) * t354 + (Icges(3,4) * t358 + Icges(3,2) * t360) * t352;
t314 = Icges(3,3) * t354 + (Icges(3,5) * t358 + Icges(3,6) * t360) * t352;
t311 = V_base(5) * rSges(2,3) - t342 * V_base(6) + t384;
t310 = t343 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t303 = t342 * V_base(4) - t343 * V_base(5) + t380;
t294 = rSges(3,1) * t323 + rSges(3,2) * t322 + rSges(3,3) * t390;
t293 = rSges(3,1) * t321 + rSges(3,2) * t320 - rSges(3,3) * t389;
t292 = Icges(3,1) * t323 + Icges(3,4) * t322 + Icges(3,5) * t390;
t291 = Icges(3,1) * t321 + Icges(3,4) * t320 - Icges(3,5) * t389;
t290 = Icges(3,4) * t323 + Icges(3,2) * t322 + Icges(3,6) * t390;
t289 = Icges(3,4) * t321 + Icges(3,2) * t320 - Icges(3,6) * t389;
t288 = Icges(3,5) * t323 + Icges(3,6) * t322 + Icges(3,3) * t390;
t287 = Icges(3,5) * t321 + Icges(3,6) * t320 - Icges(3,3) * t389;
t283 = t304 * t356 + t319 * t395;
t282 = -t304 * t395 + t319 * t356;
t267 = rSges(5,1) * t319 - rSges(5,2) * t305 + rSges(5,3) * t304;
t266 = rSges(4,1) * t305 - rSges(4,2) * t304 + rSges(4,3) * t319;
t257 = t280 * t356 + t307 * t395;
t256 = -t280 * t395 + t307 * t356;
t255 = t278 * t356 + t306 * t395;
t254 = -t278 * t395 + t306 * t356;
t253 = t283 * t359 + t305 * t355;
t252 = -t283 * t355 + t305 * t359;
t247 = pkin(5) * t283 + pkin(11) * t282;
t244 = -t293 * t345 + t317 * t334 + t372;
t243 = t294 * t345 - t317 * t335 + t370;
t242 = qJD(6) * t282 + t270;
t238 = rSges(6,1) * t283 - rSges(6,2) * t282 + rSges(6,3) * t305;
t237 = rSges(5,1) * t307 - rSges(5,2) * t281 + rSges(5,3) * t280;
t236 = rSges(5,1) * t306 - rSges(5,2) * t279 + rSges(5,3) * t278;
t235 = rSges(4,1) * t281 - rSges(4,2) * t280 + rSges(4,3) * t307;
t234 = rSges(4,1) * t279 - rSges(4,2) * t278 + rSges(4,3) * t306;
t233 = Icges(6,1) * t283 - Icges(6,4) * t282 + Icges(6,5) * t305;
t232 = Icges(6,4) * t283 - Icges(6,2) * t282 + Icges(6,6) * t305;
t231 = Icges(6,5) * t283 - Icges(6,6) * t282 + Icges(6,3) * t305;
t218 = t293 * t335 - t294 * t334 + t371;
t217 = t257 * t359 + t281 * t355;
t216 = -t257 * t355 + t281 * t359;
t215 = t255 * t359 + t279 * t355;
t214 = -t255 * t355 + t279 * t359;
t212 = pkin(5) * t257 + pkin(11) * t256;
t211 = pkin(5) * t255 + pkin(11) * t254;
t210 = qJD(6) * t256 + t250;
t209 = qJD(6) * t254 + t249;
t208 = rSges(6,1) * t257 - rSges(6,2) * t256 + rSges(6,3) * t281;
t207 = rSges(6,1) * t255 - rSges(6,2) * t254 + rSges(6,3) * t279;
t206 = rSges(7,1) * t253 + rSges(7,2) * t252 + rSges(7,3) * t282;
t205 = Icges(6,1) * t257 - Icges(6,4) * t256 + Icges(6,5) * t281;
t204 = Icges(6,1) * t255 - Icges(6,4) * t254 + Icges(6,5) * t279;
t203 = Icges(7,1) * t253 + Icges(7,4) * t252 + Icges(7,5) * t282;
t202 = Icges(6,4) * t257 - Icges(6,2) * t256 + Icges(6,6) * t281;
t201 = Icges(6,4) * t255 - Icges(6,2) * t254 + Icges(6,6) * t279;
t200 = Icges(7,4) * t253 + Icges(7,2) * t252 + Icges(7,6) * t282;
t199 = Icges(6,5) * t257 - Icges(6,6) * t256 + Icges(6,3) * t281;
t198 = Icges(6,5) * t255 - Icges(6,6) * t254 + Icges(6,3) * t279;
t197 = Icges(7,5) * t253 + Icges(7,6) * t252 + Icges(7,3) * t282;
t196 = rSges(7,1) * t217 + rSges(7,2) * t216 + rSges(7,3) * t256;
t195 = rSges(7,1) * t215 + rSges(7,2) * t214 + rSges(7,3) * t254;
t194 = Icges(7,1) * t217 + Icges(7,4) * t216 + Icges(7,5) * t256;
t193 = Icges(7,1) * t215 + Icges(7,4) * t214 + Icges(7,5) * t254;
t192 = Icges(7,4) * t217 + Icges(7,2) * t216 + Icges(7,6) * t256;
t191 = Icges(7,4) * t215 + Icges(7,2) * t214 + Icges(7,6) * t254;
t190 = Icges(7,5) * t217 + Icges(7,6) * t216 + Icges(7,3) * t256;
t189 = Icges(7,5) * t215 + Icges(7,6) * t214 + Icges(7,3) * t254;
t188 = -t234 * t308 + t266 * t295 + t369;
t187 = t235 * t308 - t266 * t296 + t367;
t186 = t234 * t296 - t235 * t295 + t368;
t185 = t267 * t295 + (-t236 - t245) * t308 + t366;
t184 = t237 * t308 + (-t267 - t268) * t296 + t364;
t183 = t236 * t296 + (-t237 - t246) * t295 + t365;
t182 = -t207 * t270 + t238 * t249 + t363;
t181 = t208 * t270 - t238 * t250 + t361;
t180 = t207 * t250 - t208 * t249 + t362;
t179 = -t195 * t242 + t206 * t209 - t211 * t270 + t247 * t249 + t363;
t178 = t196 * t242 - t206 * t210 + t212 * t270 - t247 * t250 + t361;
t177 = t195 * t210 - t196 * t209 + t211 * t250 - t212 * t249 + t362;
t1 = t345 * ((t287 * t334 + t288 * t335 + t314 * t345) * t354 + ((t290 * t360 + t292 * t358) * t335 + (t289 * t360 + t291 * t358) * t334 + (t315 * t360 + t316 * t358) * t345) * t352) / 0.2e1 + m(1) * (t331 ^ 2 + t332 ^ 2 + t333 ^ 2) / 0.2e1 + m(2) * (t303 ^ 2 + t310 ^ 2 + t311 ^ 2) / 0.2e1 + t270 * ((t199 * t305 - t202 * t282 + t205 * t283) * t250 + (t198 * t305 - t201 * t282 + t204 * t283) * t249 + (t231 * t305 - t232 * t282 + t233 * t283) * t270) / 0.2e1 + t242 * ((t190 * t282 + t192 * t252 + t194 * t253) * t210 + (t189 * t282 + t191 * t252 + t193 * t253) * t209 + (t197 * t282 + t200 * t252 + t203 * t253) * t242) / 0.2e1 + t249 * ((t199 * t279 - t202 * t254 + t205 * t255) * t250 + (t198 * t279 - t201 * t254 + t204 * t255) * t249 + (t231 * t279 - t232 * t254 + t233 * t255) * t270) / 0.2e1 + t250 * ((t199 * t281 - t202 * t256 + t205 * t257) * t250 + (t198 * t281 - t201 * t256 + t204 * t257) * t249 + (t231 * t281 - t232 * t256 + t233 * t257) * t270) / 0.2e1 + t210 * ((t190 * t256 + t192 * t216 + t194 * t217) * t210 + (t189 * t256 + t191 * t216 + t193 * t217) * t209 + (t197 * t256 + t200 * t216 + t203 * t217) * t242) / 0.2e1 + t209 * ((t190 * t254 + t192 * t214 + t194 * t215) * t210 + (t189 * t254 + t191 * t214 + t193 * t215) * t209 + (t197 * t254 + t200 * t214 + t203 * t215) * t242) / 0.2e1 + m(3) * (t218 ^ 2 + t243 ^ 2 + t244 ^ 2) / 0.2e1 + m(4) * (t186 ^ 2 + t187 ^ 2 + t188 ^ 2) / 0.2e1 + m(5) * (t183 ^ 2 + t184 ^ 2 + t185 ^ 2) / 0.2e1 + m(6) * (t180 ^ 2 + t181 ^ 2 + t182 ^ 2) / 0.2e1 + m(7) * (t177 ^ 2 + t178 ^ 2 + t179 ^ 2) / 0.2e1 + t334 * ((-t288 * t389 + t290 * t320 + t292 * t321) * t335 + (-t287 * t389 + t289 * t320 + t291 * t321) * t334 + (-t314 * t389 + t315 * t320 + t316 * t321) * t345) / 0.2e1 + t335 * ((t288 * t390 + t290 * t322 + t292 * t323) * t335 + (t287 * t390 + t289 * t322 + t291 * t323) * t334 + (t314 * t390 + t315 * t322 + t316 * t323) * t345) / 0.2e1 + ((t278 * t402 + t279 * t400 + t306 * t401) * t308 + (t278 * t407 + t279 * t403 + t306 * t405) * t296 + (t278 * t408 + t279 * t404 + t306 * t406) * t295) * t295 / 0.2e1 + ((t280 * t402 + t281 * t400 + t307 * t401) * t308 + (t280 * t407 + t281 * t403 + t307 * t405) * t296 + (t280 * t408 + t281 * t404 + t307 * t406) * t295) * t296 / 0.2e1 + ((t304 * t402 + t305 * t400 + t319 * t401) * t308 + (t304 * t407 + t305 * t403 + t319 * t405) * t296 + (t304 * t408 + t305 * t404 + t319 * t406) * t295) * t308 / 0.2e1 + ((-t338 * t351 + t340 * t353 + Icges(1,4)) * V_base(5) + (-t339 * t351 + t341 * t353 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t338 * t353 + t340 * t351 + Icges(1,2)) * V_base(5) + (t339 * t353 + t341 * t351 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t351 + Icges(2,6) * t353 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t353 - Icges(2,6) * t351 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
