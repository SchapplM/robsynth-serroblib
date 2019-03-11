% Calculate kinetic energy for
% S6PPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
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
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPRRPR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PPRRPR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRPR2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRRPR2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:48:56
% EndTime: 2019-03-08 18:49:00
% DurationCPUTime: 4.03s
% Computational Cost: add. (4726->406), mult. (12529->569), div. (0->0), fcn. (15975->14), ass. (0->178)
t411 = Icges(5,1) + Icges(6,2);
t410 = Icges(6,1) + Icges(5,3);
t409 = -Icges(5,4) - Icges(6,6);
t408 = -Icges(6,4) + Icges(5,5);
t407 = Icges(6,5) - Icges(5,6);
t406 = Icges(5,2) + Icges(6,3);
t348 = sin(pkin(11));
t351 = cos(pkin(11));
t405 = Icges(2,5) * t351 - Icges(2,6) * t348 + Icges(1,5);
t404 = Icges(2,5) * t348 + Icges(2,6) * t351 + Icges(1,6);
t347 = sin(pkin(12));
t350 = cos(pkin(12));
t352 = cos(pkin(6));
t381 = t351 * t352;
t320 = -t347 * t348 + t350 * t381;
t321 = t347 * t381 + t348 * t350;
t355 = sin(qJ(3));
t349 = sin(pkin(6));
t388 = sin(pkin(7));
t373 = t349 * t388;
t389 = cos(pkin(7));
t391 = cos(qJ(3));
t276 = t321 * t391 + (t320 * t389 - t351 * t373) * t355;
t374 = t349 * t389;
t304 = -t320 * t388 - t351 * t374;
t354 = sin(qJ(4));
t390 = cos(qJ(4));
t260 = t276 * t354 - t304 * t390;
t261 = t276 * t390 + t304 * t354;
t370 = t391 * t388;
t367 = t349 * t370;
t371 = t389 * t391;
t275 = -t320 * t371 + t321 * t355 + t351 * t367;
t403 = t406 * t260 + t409 * t261 + t407 * t275;
t383 = t348 * t352;
t322 = -t347 * t351 - t350 * t383;
t323 = -t347 * t383 + t350 * t351;
t278 = t323 * t391 + (t322 * t389 + t348 * t373) * t355;
t305 = -t322 * t388 + t348 * t374;
t262 = t278 * t354 - t305 * t390;
t263 = t278 * t390 + t305 * t354;
t277 = -t322 * t371 + t323 * t355 - t348 * t367;
t402 = t406 * t262 + t409 * t263 + t407 * t277;
t401 = t407 * t260 + t408 * t261 + t410 * t275;
t400 = t407 * t262 + t408 * t263 + t410 * t277;
t399 = t409 * t260 + t411 * t261 + t408 * t275;
t398 = t409 * t262 + t411 * t263 + t408 * t277;
t302 = t352 * t388 * t355 + (t350 * t355 * t389 + t347 * t391) * t349;
t319 = -t350 * t373 + t352 * t389;
t280 = t302 * t354 - t319 * t390;
t281 = t302 * t390 + t319 * t354;
t385 = t347 * t349;
t301 = -t349 * t350 * t371 - t352 * t370 + t355 * t385;
t397 = t406 * t280 + t409 * t281 + t407 * t301;
t396 = t407 * t280 + t408 * t281 + t410 * t301;
t395 = t409 * t280 + t411 * t281 + t408 * t301;
t387 = Icges(2,4) * t348;
t386 = qJ(2) * t352;
t384 = t348 * t349;
t382 = t349 * t351;
t380 = qJD(2) * t349;
t379 = V_base(5) * qJ(1) + V_base(1);
t375 = qJD(1) + V_base(3);
t295 = qJD(3) * t305 + V_base(4);
t294 = qJD(3) * t304 + V_base(5);
t312 = qJD(3) * t319 + V_base(6);
t372 = -qJ(1) - t386;
t259 = qJD(4) * t277 + t295;
t258 = qJD(4) * t275 + t294;
t279 = qJD(4) * t301 + t312;
t369 = t348 * t380 + V_base(5) * t386 + t379;
t326 = pkin(1) * t348 - qJ(2) * t382;
t368 = qJD(2) * t352 + V_base(4) * t326 + t375;
t327 = pkin(1) * t351 + qJ(2) * t384;
t366 = V_base(6) * t327 - t351 * t380 + V_base(2);
t284 = t321 * pkin(2) + pkin(8) * t304;
t307 = pkin(2) * t385 + pkin(8) * t319;
t365 = V_base(5) * t307 + (-t284 - t326) * V_base(6) + t369;
t285 = t323 * pkin(2) + pkin(8) * t305;
t364 = V_base(4) * t284 + (-t285 - t327) * V_base(5) + t368;
t249 = pkin(3) * t276 + pkin(9) * t275;
t271 = pkin(3) * t302 + pkin(9) * t301;
t363 = -t249 * t312 + t294 * t271 + t365;
t250 = pkin(3) * t278 + pkin(9) * t277;
t362 = t295 * t249 - t250 * t294 + t364;
t361 = V_base(6) * t285 + (-t307 + t372) * V_base(4) + t366;
t251 = pkin(4) * t281 + qJ(5) * t280;
t360 = qJD(5) * t262 + t258 * t251 + t363;
t220 = pkin(4) * t261 + qJ(5) * t260;
t359 = qJD(5) * t280 + t259 * t220 + t362;
t358 = t312 * t250 - t271 * t295 + t361;
t221 = pkin(4) * t263 + qJ(5) * t262;
t357 = qJD(5) * t260 + t279 * t221 + t358;
t356 = cos(qJ(6));
t353 = sin(qJ(6));
t345 = Icges(2,4) * t351;
t340 = rSges(2,1) * t351 - rSges(2,2) * t348;
t339 = rSges(2,1) * t348 + rSges(2,2) * t351;
t338 = Icges(2,1) * t351 - t387;
t337 = Icges(2,1) * t348 + t345;
t336 = -Icges(2,2) * t348 + t345;
t335 = Icges(2,2) * t351 + t387;
t332 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t331 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t330 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t316 = rSges(3,3) * t352 + (rSges(3,1) * t347 + rSges(3,2) * t350) * t349;
t315 = Icges(3,5) * t352 + (Icges(3,1) * t347 + Icges(3,4) * t350) * t349;
t314 = Icges(3,6) * t352 + (Icges(3,4) * t347 + Icges(3,2) * t350) * t349;
t313 = Icges(3,3) * t352 + (Icges(3,5) * t347 + Icges(3,6) * t350) * t349;
t309 = V_base(5) * rSges(2,3) - t339 * V_base(6) + t379;
t308 = t340 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t303 = t339 * V_base(4) - t340 * V_base(5) + t375;
t293 = rSges(3,1) * t323 + rSges(3,2) * t322 + rSges(3,3) * t384;
t292 = rSges(3,1) * t321 + rSges(3,2) * t320 - rSges(3,3) * t382;
t291 = Icges(3,1) * t323 + Icges(3,4) * t322 + Icges(3,5) * t384;
t290 = Icges(3,1) * t321 + Icges(3,4) * t320 - Icges(3,5) * t382;
t289 = Icges(3,4) * t323 + Icges(3,2) * t322 + Icges(3,6) * t384;
t288 = Icges(3,4) * t321 + Icges(3,2) * t320 - Icges(3,6) * t382;
t287 = Icges(3,5) * t323 + Icges(3,6) * t322 + Icges(3,3) * t384;
t286 = Icges(3,5) * t321 + Icges(3,6) * t320 - Icges(3,3) * t382;
t270 = rSges(4,1) * t302 - rSges(4,2) * t301 + rSges(4,3) * t319;
t269 = Icges(4,1) * t302 - Icges(4,4) * t301 + Icges(4,5) * t319;
t268 = Icges(4,4) * t302 - Icges(4,2) * t301 + Icges(4,6) * t319;
t267 = Icges(4,5) * t302 - Icges(4,6) * t301 + Icges(4,3) * t319;
t266 = pkin(5) * t301 + pkin(10) * t281;
t265 = t280 * t353 + t301 * t356;
t264 = t280 * t356 - t301 * t353;
t254 = t316 * V_base(5) + (-t292 - t326) * V_base(6) + t369;
t253 = t293 * V_base(6) + (-t316 + t372) * V_base(4) + t366;
t252 = qJD(6) * t281 + t279;
t247 = t292 * V_base(4) + (-t293 - t327) * V_base(5) + t368;
t246 = rSges(5,1) * t281 - rSges(5,2) * t280 + rSges(5,3) * t301;
t245 = rSges(6,1) * t301 - rSges(6,2) * t281 + rSges(6,3) * t280;
t238 = rSges(4,1) * t278 - rSges(4,2) * t277 + rSges(4,3) * t305;
t237 = rSges(4,1) * t276 - rSges(4,2) * t275 + rSges(4,3) * t304;
t236 = Icges(4,1) * t278 - Icges(4,4) * t277 + Icges(4,5) * t305;
t235 = Icges(4,1) * t276 - Icges(4,4) * t275 + Icges(4,5) * t304;
t234 = Icges(4,4) * t278 - Icges(4,2) * t277 + Icges(4,6) * t305;
t233 = Icges(4,4) * t276 - Icges(4,2) * t275 + Icges(4,6) * t304;
t232 = Icges(4,5) * t278 - Icges(4,6) * t277 + Icges(4,3) * t305;
t231 = Icges(4,5) * t276 - Icges(4,6) * t275 + Icges(4,3) * t304;
t229 = pkin(5) * t277 + pkin(10) * t263;
t228 = pkin(5) * t275 + pkin(10) * t261;
t227 = t262 * t353 + t277 * t356;
t226 = t262 * t356 - t277 * t353;
t225 = t260 * t353 + t275 * t356;
t224 = t260 * t356 - t275 * t353;
t223 = qJD(6) * t263 + t259;
t222 = qJD(6) * t261 + t258;
t217 = rSges(7,1) * t265 + rSges(7,2) * t264 + rSges(7,3) * t281;
t216 = Icges(7,1) * t265 + Icges(7,4) * t264 + Icges(7,5) * t281;
t215 = Icges(7,4) * t265 + Icges(7,2) * t264 + Icges(7,6) * t281;
t214 = Icges(7,5) * t265 + Icges(7,6) * t264 + Icges(7,3) * t281;
t213 = rSges(5,1) * t263 - rSges(5,2) * t262 + rSges(5,3) * t277;
t212 = rSges(5,1) * t261 - rSges(5,2) * t260 + rSges(5,3) * t275;
t211 = rSges(6,1) * t277 - rSges(6,2) * t263 + rSges(6,3) * t262;
t210 = rSges(6,1) * t275 - rSges(6,2) * t261 + rSges(6,3) * t260;
t196 = rSges(7,1) * t227 + rSges(7,2) * t226 + rSges(7,3) * t263;
t195 = rSges(7,1) * t225 + rSges(7,2) * t224 + rSges(7,3) * t261;
t194 = Icges(7,1) * t227 + Icges(7,4) * t226 + Icges(7,5) * t263;
t193 = Icges(7,1) * t225 + Icges(7,4) * t224 + Icges(7,5) * t261;
t192 = Icges(7,4) * t227 + Icges(7,2) * t226 + Icges(7,6) * t263;
t191 = Icges(7,4) * t225 + Icges(7,2) * t224 + Icges(7,6) * t261;
t190 = Icges(7,5) * t227 + Icges(7,6) * t226 + Icges(7,3) * t263;
t189 = Icges(7,5) * t225 + Icges(7,6) * t224 + Icges(7,3) * t261;
t188 = -t237 * t312 + t270 * t294 + t365;
t187 = t238 * t312 - t270 * t295 + t361;
t186 = t237 * t295 - t238 * t294 + t364;
t185 = -t212 * t279 + t246 * t258 + t363;
t184 = t213 * t279 - t246 * t259 + t358;
t183 = t212 * t259 - t213 * t258 + t362;
t182 = t245 * t258 + (-t210 - t220) * t279 + t360;
t181 = t211 * t279 + (-t245 - t251) * t259 + t357;
t180 = t210 * t259 + (-t211 - t221) * t258 + t359;
t179 = -t195 * t252 + t217 * t222 + t258 * t266 + (-t220 - t228) * t279 + t360;
t178 = t196 * t252 - t217 * t223 + t229 * t279 + (-t251 - t266) * t259 + t357;
t177 = (-t221 - t229) * t258 + t195 * t223 - t196 * t222 + t228 * t259 + t359;
t1 = m(1) * (t330 ^ 2 + t331 ^ 2 + t332 ^ 2) / 0.2e1 + t312 * ((t232 * t319 - t234 * t301 + t236 * t302) * t295 + (t231 * t319 - t233 * t301 + t235 * t302) * t294 + (t267 * t319 - t268 * t301 + t269 * t302) * t312) / 0.2e1 + m(2) * (t303 ^ 2 + t308 ^ 2 + t309 ^ 2) / 0.2e1 + t294 * ((t232 * t304 - t234 * t275 + t236 * t276) * t295 + (t231 * t304 - t233 * t275 + t235 * t276) * t294 + (t267 * t304 - t268 * t275 + t269 * t276) * t312) / 0.2e1 + t295 * ((t232 * t305 - t234 * t277 + t236 * t278) * t295 + (t231 * t305 - t233 * t277 + t235 * t278) * t294 + (t267 * t305 - t268 * t277 + t269 * t278) * t312) / 0.2e1 + t252 * ((t190 * t281 + t192 * t264 + t194 * t265) * t223 + (t189 * t281 + t191 * t264 + t193 * t265) * t222 + (t281 * t214 + t264 * t215 + t265 * t216) * t252) / 0.2e1 + t223 * ((t263 * t190 + t226 * t192 + t227 * t194) * t223 + (t189 * t263 + t191 * t226 + t193 * t227) * t222 + (t214 * t263 + t215 * t226 + t216 * t227) * t252) / 0.2e1 + t222 * ((t190 * t261 + t192 * t224 + t194 * t225) * t223 + (t261 * t189 + t224 * t191 + t225 * t193) * t222 + (t214 * t261 + t215 * t224 + t216 * t225) * t252) / 0.2e1 + m(3) * (t247 ^ 2 + t253 ^ 2 + t254 ^ 2) / 0.2e1 + m(4) * (t186 ^ 2 + t187 ^ 2 + t188 ^ 2) / 0.2e1 + m(5) * (t183 ^ 2 + t184 ^ 2 + t185 ^ 2) / 0.2e1 + m(6) * (t180 ^ 2 + t181 ^ 2 + t182 ^ 2) / 0.2e1 + m(7) * (t177 ^ 2 + t178 ^ 2 + t179 ^ 2) / 0.2e1 + ((t260 * t397 + t261 * t395 + t275 * t396) * t279 + (t260 * t402 + t261 * t398 + t275 * t400) * t259 + (t403 * t260 + t399 * t261 + t401 * t275) * t258) * t258 / 0.2e1 + ((t262 * t397 + t263 * t395 + t277 * t396) * t279 + (t402 * t262 + t398 * t263 + t400 * t277) * t259 + (t262 * t403 + t263 * t399 + t277 * t401) * t258) * t259 / 0.2e1 + ((t280 * t397 + t281 * t395 + t301 * t396) * t279 + (t280 * t402 + t281 * t398 + t301 * t400) * t259 + (t280 * t403 + t281 * t399 + t301 * t401) * t258) * t279 / 0.2e1 + ((t313 * t384 + t314 * t322 + t315 * t323 + t405) * V_base(6) + (t286 * t384 + t288 * t322 + t290 * t323 - t335 * t348 + t337 * t351 + Icges(1,4)) * V_base(5) + (t287 * t384 + t289 * t322 + t291 * t323 - t336 * t348 + t338 * t351 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t313 * t382 + t314 * t320 + t315 * t321 + t404) * V_base(6) + (-t286 * t382 + t288 * t320 + t290 * t321 + t335 * t351 + t337 * t348 + Icges(1,2)) * V_base(5) + (-t287 * t382 + t289 * t320 + t291 * t321 + t336 * t351 + t338 * t348 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t313 * t352 + (t314 * t350 + t315 * t347) * t349 + Icges(2,3) + Icges(1,3)) * V_base(6) + (t286 * t352 + (t288 * t350 + t290 * t347) * t349 + t404) * V_base(5) + (t287 * t352 + (t289 * t350 + t291 * t347) * t349 + t405) * V_base(4)) * V_base(6) / 0.2e1;
T  = t1;
