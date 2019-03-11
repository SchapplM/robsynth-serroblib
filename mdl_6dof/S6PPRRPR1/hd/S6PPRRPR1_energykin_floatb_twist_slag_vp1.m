% Calculate kinetic energy for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
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
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPRRPR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PPRRPR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_energykin_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRPR1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRRPR1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:45:12
% EndTime: 2019-03-08 18:45:17
% DurationCPUTime: 4.36s
% Computational Cost: add. (5508->447), mult. (14201->627), div. (0->0), fcn. (18231->16), ass. (0->192)
t412 = Icges(5,2) + Icges(6,3);
t352 = sin(pkin(11));
t356 = cos(pkin(11));
t411 = Icges(2,5) * t356 - Icges(2,6) * t352 + Icges(1,5);
t410 = Icges(2,5) * t352 + Icges(2,6) * t356 + Icges(1,6);
t351 = sin(pkin(12));
t355 = cos(pkin(12));
t357 = cos(pkin(6));
t389 = t356 * t357;
t319 = -t351 * t352 + t355 * t389;
t320 = t351 * t389 + t352 * t355;
t360 = sin(qJ(3));
t353 = sin(pkin(6));
t399 = sin(pkin(7));
t380 = t353 * t399;
t400 = cos(pkin(7));
t403 = cos(qJ(3));
t278 = t320 * t403 + (t319 * t400 - t356 * t380) * t360;
t359 = sin(qJ(4));
t381 = t353 * t400;
t372 = -t319 * t399 - t356 * t381;
t402 = cos(qJ(4));
t266 = t278 * t402 + t359 * t372;
t377 = t403 * t399;
t374 = t353 * t377;
t378 = t400 * t403;
t277 = -t319 * t378 + t320 * t360 + t356 * t374;
t350 = sin(pkin(13));
t354 = cos(pkin(13));
t231 = -t266 * t350 + t277 * t354;
t396 = t277 * t350;
t232 = t266 * t354 + t396;
t265 = t278 * t359 - t372 * t402;
t409 = -Icges(5,4) * t266 + Icges(6,5) * t232 - Icges(5,6) * t277 + Icges(6,6) * t231 + t412 * t265;
t391 = t352 * t357;
t321 = -t351 * t356 - t355 * t391;
t322 = -t351 * t391 + t355 * t356;
t280 = t322 * t403 + (t321 * t400 + t352 * t380) * t360;
t371 = -t321 * t399 + t352 * t381;
t268 = t280 * t402 + t359 * t371;
t279 = -t321 * t378 + t322 * t360 - t352 * t374;
t233 = -t268 * t350 + t279 * t354;
t395 = t279 * t350;
t234 = t268 * t354 + t395;
t267 = t280 * t359 - t371 * t402;
t408 = -Icges(5,4) * t268 + Icges(6,5) * t234 - Icges(5,6) * t279 + Icges(6,6) * t233 + t412 * t267;
t304 = t357 * t399 * t360 + (t355 * t360 * t400 + t351 * t403) * t353;
t370 = -t355 * t380 + t357 * t400;
t283 = t304 * t402 + t359 * t370;
t393 = t351 * t353;
t303 = -t353 * t355 * t378 - t357 * t377 + t360 * t393;
t261 = -t283 * t350 + t303 * t354;
t394 = t303 * t350;
t262 = t283 * t354 + t394;
t282 = t304 * t359 - t370 * t402;
t407 = -Icges(5,4) * t283 + Icges(6,5) * t262 - Icges(5,6) * t303 + Icges(6,6) * t261 + t412 * t282;
t401 = pkin(5) * t354;
t398 = Icges(2,4) * t352;
t397 = qJ(2) * t357;
t392 = t352 * t353;
t390 = t353 * t356;
t387 = qJD(2) * t353;
t386 = V_base(5) * qJ(1) + V_base(1);
t382 = qJD(1) + V_base(3);
t297 = qJD(3) * t371 + V_base(4);
t296 = qJD(3) * t372 + V_base(5);
t312 = qJD(3) * t370 + V_base(6);
t379 = -qJ(1) - t397;
t264 = qJD(4) * t279 + t297;
t263 = qJD(4) * t277 + t296;
t281 = qJD(4) * t303 + t312;
t376 = t352 * t387 + V_base(5) * t397 + t386;
t325 = pkin(1) * t352 - qJ(2) * t390;
t375 = qJD(2) * t357 + V_base(4) * t325 + t382;
t326 = pkin(1) * t356 + qJ(2) * t392;
t373 = V_base(6) * t326 - t356 * t387 + V_base(2);
t286 = t320 * pkin(2) + pkin(8) * t372;
t307 = pkin(2) * t393 + pkin(8) * t370;
t369 = V_base(5) * t307 + (-t286 - t325) * V_base(6) + t376;
t287 = t322 * pkin(2) + pkin(8) * t371;
t368 = V_base(4) * t286 + (-t287 - t326) * V_base(5) + t375;
t250 = pkin(3) * t278 + pkin(9) * t277;
t273 = pkin(3) * t304 + pkin(9) * t303;
t367 = -t250 * t312 + t296 * t273 + t369;
t251 = pkin(3) * t280 + pkin(9) * t279;
t366 = t297 * t250 - t251 * t296 + t368;
t365 = V_base(6) * t287 + (-t307 + t379) * V_base(4) + t373;
t252 = pkin(4) * t283 + qJ(5) * t282;
t364 = qJD(5) * t267 + t263 * t252 + t367;
t223 = pkin(4) * t266 + qJ(5) * t265;
t363 = qJD(5) * t282 + t264 * t223 + t366;
t362 = t312 * t251 - t273 * t297 + t365;
t224 = pkin(4) * t268 + qJ(5) * t267;
t361 = qJD(5) * t265 + t281 * t224 + t362;
t349 = pkin(13) + qJ(6);
t347 = Icges(2,4) * t356;
t346 = cos(t349);
t345 = sin(t349);
t339 = rSges(2,1) * t356 - rSges(2,2) * t352;
t338 = rSges(2,1) * t352 + rSges(2,2) * t356;
t337 = Icges(2,1) * t356 - t398;
t336 = Icges(2,1) * t352 + t347;
t335 = -Icges(2,2) * t352 + t347;
t334 = Icges(2,2) * t356 + t398;
t331 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t330 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t329 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t316 = rSges(3,3) * t357 + (rSges(3,1) * t351 + rSges(3,2) * t355) * t353;
t315 = Icges(3,5) * t357 + (Icges(3,1) * t351 + Icges(3,4) * t355) * t353;
t314 = Icges(3,6) * t357 + (Icges(3,4) * t351 + Icges(3,2) * t355) * t353;
t313 = Icges(3,3) * t357 + (Icges(3,5) * t351 + Icges(3,6) * t355) * t353;
t309 = V_base(5) * rSges(2,3) - t338 * V_base(6) + t386;
t308 = t339 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t305 = t338 * V_base(4) - t339 * V_base(5) + t382;
t295 = rSges(3,1) * t322 + rSges(3,2) * t321 + rSges(3,3) * t392;
t294 = rSges(3,1) * t320 + rSges(3,2) * t319 - rSges(3,3) * t390;
t293 = Icges(3,1) * t322 + Icges(3,4) * t321 + Icges(3,5) * t392;
t292 = Icges(3,1) * t320 + Icges(3,4) * t319 - Icges(3,5) * t390;
t291 = Icges(3,4) * t322 + Icges(3,2) * t321 + Icges(3,6) * t392;
t290 = Icges(3,4) * t320 + Icges(3,2) * t319 - Icges(3,6) * t390;
t289 = Icges(3,5) * t322 + Icges(3,6) * t321 + Icges(3,3) * t392;
t288 = Icges(3,5) * t320 + Icges(3,6) * t319 - Icges(3,3) * t390;
t272 = t304 * rSges(4,1) - t303 * rSges(4,2) + rSges(4,3) * t370;
t271 = Icges(4,1) * t304 - Icges(4,4) * t303 + Icges(4,5) * t370;
t270 = Icges(4,4) * t304 - Icges(4,2) * t303 + Icges(4,6) * t370;
t269 = Icges(4,5) * t304 - Icges(4,6) * t303 + Icges(4,3) * t370;
t258 = t283 * t346 + t303 * t345;
t257 = -t283 * t345 + t303 * t346;
t255 = t316 * V_base(5) + (-t294 - t325) * V_base(6) + t376;
t254 = t295 * V_base(6) + (-t316 + t379) * V_base(4) + t373;
t253 = qJD(6) * t282 + t281;
t248 = t294 * V_base(4) + (-t295 - t326) * V_base(5) + t375;
t247 = rSges(5,1) * t283 - rSges(5,2) * t282 + rSges(5,3) * t303;
t246 = Icges(5,1) * t283 - Icges(5,4) * t282 + Icges(5,5) * t303;
t244 = Icges(5,5) * t283 - Icges(5,6) * t282 + Icges(5,3) * t303;
t243 = t280 * rSges(4,1) - t279 * rSges(4,2) + rSges(4,3) * t371;
t242 = t278 * rSges(4,1) - t277 * rSges(4,2) + rSges(4,3) * t372;
t241 = Icges(4,1) * t280 - Icges(4,4) * t279 + Icges(4,5) * t371;
t240 = Icges(4,1) * t278 - Icges(4,4) * t277 + Icges(4,5) * t372;
t239 = Icges(4,4) * t280 - Icges(4,2) * t279 + Icges(4,6) * t371;
t238 = Icges(4,4) * t278 - Icges(4,2) * t277 + Icges(4,6) * t372;
t237 = Icges(4,5) * t280 - Icges(4,6) * t279 + Icges(4,3) * t371;
t236 = Icges(4,5) * t278 - Icges(4,6) * t277 + Icges(4,3) * t372;
t230 = t268 * t346 + t279 * t345;
t229 = -t268 * t345 + t279 * t346;
t228 = t266 * t346 + t277 * t345;
t227 = -t266 * t345 + t277 * t346;
t226 = qJD(6) * t267 + t264;
t225 = qJD(6) * t265 + t263;
t220 = rSges(5,1) * t268 - rSges(5,2) * t267 + rSges(5,3) * t279;
t219 = rSges(5,1) * t266 - rSges(5,2) * t265 + rSges(5,3) * t277;
t218 = rSges(6,1) * t262 + rSges(6,2) * t261 + rSges(6,3) * t282;
t217 = Icges(5,1) * t268 - Icges(5,4) * t267 + Icges(5,5) * t279;
t216 = Icges(5,1) * t266 - Icges(5,4) * t265 + Icges(5,5) * t277;
t213 = Icges(5,5) * t268 - Icges(5,6) * t267 + Icges(5,3) * t279;
t212 = Icges(5,5) * t266 - Icges(5,6) * t265 + Icges(5,3) * t277;
t211 = Icges(6,1) * t262 + Icges(6,4) * t261 + Icges(6,5) * t282;
t210 = Icges(6,4) * t262 + Icges(6,2) * t261 + Icges(6,6) * t282;
t208 = rSges(7,1) * t258 + rSges(7,2) * t257 + rSges(7,3) * t282;
t207 = Icges(7,1) * t258 + Icges(7,4) * t257 + Icges(7,5) * t282;
t206 = Icges(7,4) * t258 + Icges(7,2) * t257 + Icges(7,6) * t282;
t205 = Icges(7,5) * t258 + Icges(7,6) * t257 + Icges(7,3) * t282;
t204 = pkin(5) * t394 + pkin(10) * t282 + t283 * t401;
t202 = -t242 * t312 + t272 * t296 + t369;
t201 = t243 * t312 - t272 * t297 + t365;
t200 = rSges(6,1) * t234 + rSges(6,2) * t233 + rSges(6,3) * t267;
t199 = rSges(6,1) * t232 + rSges(6,2) * t231 + rSges(6,3) * t265;
t198 = Icges(6,1) * t234 + Icges(6,4) * t233 + Icges(6,5) * t267;
t197 = Icges(6,1) * t232 + Icges(6,4) * t231 + Icges(6,5) * t265;
t196 = Icges(6,4) * t234 + Icges(6,2) * t233 + Icges(6,6) * t267;
t195 = Icges(6,4) * t232 + Icges(6,2) * t231 + Icges(6,6) * t265;
t192 = rSges(7,1) * t230 + rSges(7,2) * t229 + rSges(7,3) * t267;
t191 = rSges(7,1) * t228 + rSges(7,2) * t227 + rSges(7,3) * t265;
t190 = Icges(7,1) * t230 + Icges(7,4) * t229 + Icges(7,5) * t267;
t189 = Icges(7,1) * t228 + Icges(7,4) * t227 + Icges(7,5) * t265;
t188 = Icges(7,4) * t230 + Icges(7,2) * t229 + Icges(7,6) * t267;
t187 = Icges(7,4) * t228 + Icges(7,2) * t227 + Icges(7,6) * t265;
t186 = Icges(7,5) * t230 + Icges(7,6) * t229 + Icges(7,3) * t267;
t185 = Icges(7,5) * t228 + Icges(7,6) * t227 + Icges(7,3) * t265;
t184 = pkin(5) * t395 + pkin(10) * t267 + t268 * t401;
t183 = pkin(5) * t396 + pkin(10) * t265 + t266 * t401;
t182 = t242 * t297 - t243 * t296 + t368;
t181 = -t219 * t281 + t247 * t263 + t367;
t180 = t220 * t281 - t247 * t264 + t362;
t179 = t219 * t264 - t220 * t263 + t366;
t178 = t218 * t263 + (-t199 - t223) * t281 + t364;
t177 = t200 * t281 + (-t218 - t252) * t264 + t361;
t176 = t199 * t264 + (-t200 - t224) * t263 + t363;
t175 = (-t183 - t223) * t281 - t191 * t253 + t204 * t263 + t208 * t225 + t364;
t174 = t184 * t281 + t192 * t253 - t208 * t226 + (-t204 - t252) * t264 + t361;
t173 = t363 + (-t184 - t224) * t263 + t183 * t264 + t191 * t226 - t192 * t225;
t1 = t312 * ((t237 * t370 - t303 * t239 + t304 * t241) * t297 + (t236 * t370 - t303 * t238 + t304 * t240) * t296 + (t269 * t370 - t303 * t270 + t304 * t271) * t312) / 0.2e1 + t297 * ((t237 * t371 - t279 * t239 + t280 * t241) * t297 + (t236 * t371 - t279 * t238 + t280 * t240) * t296 + (t269 * t371 - t279 * t270 + t280 * t271) * t312) / 0.2e1 + t296 * ((t237 * t372 - t277 * t239 + t278 * t241) * t297 + (t236 * t372 - t277 * t238 + t278 * t240) * t296 + (t269 * t372 - t277 * t270 + t278 * t271) * t312) / 0.2e1 + m(1) * (t329 ^ 2 + t330 ^ 2 + t331 ^ 2) / 0.2e1 + m(2) * (t305 ^ 2 + t308 ^ 2 + t309 ^ 2) / 0.2e1 + t253 * ((t186 * t282 + t188 * t257 + t190 * t258) * t226 + (t185 * t282 + t187 * t257 + t189 * t258) * t225 + (t205 * t282 + t206 * t257 + t207 * t258) * t253) / 0.2e1 + t226 * ((t186 * t267 + t229 * t188 + t230 * t190) * t226 + (t185 * t267 + t187 * t229 + t189 * t230) * t225 + (t205 * t267 + t206 * t229 + t207 * t230) * t253) / 0.2e1 + t225 * ((t186 * t265 + t188 * t227 + t190 * t228) * t226 + (t185 * t265 + t227 * t187 + t228 * t189) * t225 + (t205 * t265 + t206 * t227 + t207 * t228) * t253) / 0.2e1 + m(3) * (t248 ^ 2 + t254 ^ 2 + t255 ^ 2) / 0.2e1 + m(4) * (t182 ^ 2 + t201 ^ 2 + t202 ^ 2) / 0.2e1 + m(7) * (t173 ^ 2 + t174 ^ 2 + t175 ^ 2) / 0.2e1 + m(6) * (t176 ^ 2 + t177 ^ 2 + t178 ^ 2) / 0.2e1 + m(5) * (t179 ^ 2 + t180 ^ 2 + t181 ^ 2) / 0.2e1 + ((t210 * t231 + t211 * t232 + t244 * t277 + t246 * t266 + t265 * t407) * t281 + (t196 * t231 + t198 * t232 + t213 * t277 + t217 * t266 + t265 * t408) * t264 + (t231 * t195 + t197 * t232 + t212 * t277 + t216 * t266 + t409 * t265) * t263) * t263 / 0.2e1 + ((t210 * t233 + t211 * t234 + t244 * t279 + t246 * t268 + t267 * t407) * t281 + (t196 * t233 + t198 * t234 + t213 * t279 + t217 * t268 + t408 * t267) * t264 + (t195 * t233 + t197 * t234 + t212 * t279 + t216 * t268 + t267 * t409) * t263) * t264 / 0.2e1 + ((t210 * t261 + t211 * t262 + t244 * t303 + t246 * t283 + t407 * t282) * t281 + (t196 * t261 + t198 * t262 + t213 * t303 + t217 * t283 + t282 * t408) * t264 + (t195 * t261 + t197 * t262 + t212 * t303 + t216 * t283 + t282 * t409) * t263) * t281 / 0.2e1 + ((t313 * t392 + t314 * t321 + t315 * t322 + t411) * V_base(6) + (t288 * t392 + t290 * t321 + t292 * t322 - t334 * t352 + t336 * t356 + Icges(1,4)) * V_base(5) + (t289 * t392 + t291 * t321 + t293 * t322 - t335 * t352 + t337 * t356 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t313 * t390 + t314 * t319 + t315 * t320 + t410) * V_base(6) + (-t288 * t390 + t290 * t319 + t292 * t320 + t334 * t356 + t336 * t352 + Icges(1,2)) * V_base(5) + (-t289 * t390 + t291 * t319 + t293 * t320 + t335 * t356 + t337 * t352 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t313 * t357 + (t314 * t355 + t315 * t351) * t353 + Icges(2,3) + Icges(1,3)) * V_base(6) + (t288 * t357 + (t290 * t355 + t292 * t351) * t353 + t410) * V_base(5) + (t289 * t357 + (t291 * t355 + t293 * t351) * t353 + t411) * V_base(4)) * V_base(6) / 0.2e1;
T  = t1;
