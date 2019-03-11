% Calculate kinetic energy for
% S6PRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR7_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRPR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_energykin_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR7_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR7_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPR7_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:38:53
% EndTime: 2019-03-08 23:38:58
% DurationCPUTime: 5.17s
% Computational Cost: add. (5643->442), mult. (14426->643), div. (0->0), fcn. (18456->16), ass. (0->192)
t414 = Icges(5,2) + Icges(6,3);
t354 = sin(pkin(12));
t357 = cos(pkin(12));
t362 = sin(qJ(2));
t358 = cos(pkin(6));
t363 = cos(qJ(2));
t391 = t358 * t363;
t319 = -t354 * t362 + t357 * t391;
t392 = t358 * t362;
t320 = t354 * t363 + t357 * t392;
t361 = sin(qJ(3));
t355 = sin(pkin(6));
t400 = sin(pkin(7));
t382 = t355 * t400;
t401 = cos(pkin(7));
t405 = cos(qJ(3));
t281 = t320 * t405 + (t319 * t401 - t357 * t382) * t361;
t360 = sin(qJ(4));
t383 = t355 * t401;
t376 = -t319 * t400 - t357 * t383;
t404 = cos(qJ(4));
t266 = t281 * t404 + t360 * t376;
t380 = t405 * t400;
t379 = t355 * t380;
t381 = t401 * t405;
t280 = -t319 * t381 + t320 * t361 + t357 * t379;
t353 = sin(pkin(13));
t356 = cos(pkin(13));
t232 = -t266 * t353 + t280 * t356;
t398 = t280 * t353;
t233 = t266 * t356 + t398;
t265 = t281 * t360 - t376 * t404;
t411 = -Icges(5,4) * t266 + Icges(6,5) * t233 - Icges(5,6) * t280 + Icges(6,6) * t232 + t414 * t265;
t321 = -t354 * t391 - t357 * t362;
t322 = -t354 * t392 + t357 * t363;
t283 = t322 * t405 + (t321 * t401 + t354 * t382) * t361;
t375 = -t321 * t400 + t354 * t383;
t268 = t283 * t404 + t360 * t375;
t282 = -t321 * t381 + t322 * t361 - t354 * t379;
t234 = -t268 * t353 + t282 * t356;
t397 = t282 * t353;
t235 = t268 * t356 + t397;
t267 = t283 * t360 - t375 * t404;
t410 = -Icges(5,4) * t268 + Icges(6,5) * t235 - Icges(5,6) * t282 + Icges(6,6) * t234 + t414 * t267;
t306 = t358 * t400 * t361 + (t361 * t363 * t401 + t362 * t405) * t355;
t374 = t358 * t401 - t363 * t382;
t285 = t306 * t404 + t360 * t374;
t393 = t355 * t362;
t305 = -t355 * t363 * t381 - t358 * t380 + t361 * t393;
t263 = -t285 * t353 + t305 * t356;
t396 = t305 * t353;
t264 = t285 * t356 + t396;
t284 = t306 * t360 - t374 * t404;
t409 = -Icges(5,4) * t285 + Icges(6,5) * t264 - Icges(5,6) * t305 + Icges(6,6) * t263 + t414 * t284;
t403 = pkin(8) * t358;
t402 = pkin(5) * t356;
t399 = Icges(2,4) * t354;
t395 = t354 * t355;
t394 = t355 * t357;
t389 = qJD(2) * t355;
t388 = V_base(5) * qJ(1) + V_base(1);
t384 = qJD(1) + V_base(3);
t333 = t354 * t389 + V_base(4);
t343 = qJD(2) * t358 + V_base(6);
t297 = qJD(3) * t375 + t333;
t307 = qJD(3) * t374 + t343;
t258 = qJD(4) * t282 + t297;
t275 = qJD(4) * t305 + t307;
t332 = -t357 * t389 + V_base(5);
t296 = qJD(3) * t376 + t332;
t325 = pkin(1) * t354 - pkin(8) * t394;
t378 = -t325 * V_base(6) + V_base(5) * t403 + t388;
t326 = pkin(1) * t357 + pkin(8) * t395;
t377 = V_base(4) * t325 - t326 * V_base(5) + t384;
t257 = qJD(4) * t280 + t296;
t373 = V_base(6) * t326 + V_base(2) + (-qJ(1) - t403) * V_base(4);
t286 = t320 * pkin(2) + pkin(9) * t376;
t308 = pkin(2) * t393 + pkin(9) * t374;
t372 = -t286 * t343 + t332 * t308 + t378;
t287 = t322 * pkin(2) + pkin(9) * t375;
t371 = t333 * t286 - t287 * t332 + t377;
t370 = t343 * t287 - t308 * t333 + t373;
t253 = pkin(3) * t281 + pkin(10) * t280;
t273 = pkin(3) * t306 + pkin(10) * t305;
t369 = -t253 * t307 + t296 * t273 + t372;
t254 = pkin(3) * t283 + pkin(10) * t282;
t368 = t297 * t253 - t254 * t296 + t371;
t367 = t307 * t254 - t273 * t297 + t370;
t255 = pkin(4) * t285 + qJ(5) * t284;
t366 = qJD(5) * t267 + t257 * t255 + t369;
t225 = pkin(4) * t266 + qJ(5) * t265;
t365 = qJD(5) * t284 + t258 * t225 + t368;
t226 = pkin(4) * t268 + qJ(5) * t267;
t364 = qJD(5) * t265 + t275 * t226 + t367;
t352 = pkin(13) + qJ(6);
t350 = Icges(2,4) * t357;
t349 = cos(t352);
t348 = sin(t352);
t341 = rSges(2,1) * t357 - rSges(2,2) * t354;
t340 = rSges(2,1) * t354 + rSges(2,2) * t357;
t339 = Icges(2,1) * t357 - t399;
t338 = Icges(2,1) * t354 + t350;
t337 = -Icges(2,2) * t354 + t350;
t336 = Icges(2,2) * t357 + t399;
t331 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t330 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t329 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t316 = t358 * rSges(3,3) + (rSges(3,1) * t362 + rSges(3,2) * t363) * t355;
t315 = Icges(3,5) * t358 + (Icges(3,1) * t362 + Icges(3,4) * t363) * t355;
t314 = Icges(3,6) * t358 + (Icges(3,4) * t362 + Icges(3,2) * t363) * t355;
t313 = Icges(3,3) * t358 + (Icges(3,5) * t362 + Icges(3,6) * t363) * t355;
t310 = V_base(5) * rSges(2,3) - t340 * V_base(6) + t388;
t309 = t341 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t304 = t340 * V_base(4) - t341 * V_base(5) + t384;
t295 = rSges(3,1) * t322 + rSges(3,2) * t321 + rSges(3,3) * t395;
t294 = rSges(3,1) * t320 + rSges(3,2) * t319 - rSges(3,3) * t394;
t293 = Icges(3,1) * t322 + Icges(3,4) * t321 + Icges(3,5) * t395;
t292 = Icges(3,1) * t320 + Icges(3,4) * t319 - Icges(3,5) * t394;
t291 = Icges(3,4) * t322 + Icges(3,2) * t321 + Icges(3,6) * t395;
t290 = Icges(3,4) * t320 + Icges(3,2) * t319 - Icges(3,6) * t394;
t289 = Icges(3,5) * t322 + Icges(3,6) * t321 + Icges(3,3) * t395;
t288 = Icges(3,5) * t320 + Icges(3,6) * t319 - Icges(3,3) * t394;
t272 = t306 * rSges(4,1) - t305 * rSges(4,2) + rSges(4,3) * t374;
t271 = Icges(4,1) * t306 - Icges(4,4) * t305 + Icges(4,5) * t374;
t270 = Icges(4,4) * t306 - Icges(4,2) * t305 + Icges(4,6) * t374;
t269 = Icges(4,5) * t306 - Icges(4,6) * t305 + Icges(4,3) * t374;
t260 = t285 * t349 + t305 * t348;
t259 = -t285 * t348 + t305 * t349;
t252 = -t294 * t343 + t316 * t332 + t378;
t251 = t295 * t343 - t316 * t333 + t373;
t250 = qJD(6) * t284 + t275;
t248 = rSges(5,1) * t285 - rSges(5,2) * t284 + rSges(5,3) * t305;
t247 = t283 * rSges(4,1) - t282 * rSges(4,2) + rSges(4,3) * t375;
t246 = t281 * rSges(4,1) - t280 * rSges(4,2) + rSges(4,3) * t376;
t245 = Icges(5,1) * t285 - Icges(5,4) * t284 + Icges(5,5) * t305;
t243 = Icges(5,5) * t285 - Icges(5,6) * t284 + Icges(5,3) * t305;
t242 = Icges(4,1) * t283 - Icges(4,4) * t282 + Icges(4,5) * t375;
t241 = Icges(4,1) * t281 - Icges(4,4) * t280 + Icges(4,5) * t376;
t240 = Icges(4,4) * t283 - Icges(4,2) * t282 + Icges(4,6) * t375;
t239 = Icges(4,4) * t281 - Icges(4,2) * t280 + Icges(4,6) * t376;
t238 = Icges(4,5) * t283 - Icges(4,6) * t282 + Icges(4,3) * t375;
t237 = Icges(4,5) * t281 - Icges(4,6) * t280 + Icges(4,3) * t376;
t236 = t294 * t333 - t295 * t332 + t377;
t230 = t268 * t349 + t282 * t348;
t229 = -t268 * t348 + t282 * t349;
t228 = t266 * t349 + t280 * t348;
t227 = -t266 * t348 + t280 * t349;
t224 = qJD(6) * t267 + t258;
t223 = qJD(6) * t265 + t257;
t221 = rSges(5,1) * t268 - rSges(5,2) * t267 + rSges(5,3) * t282;
t220 = rSges(5,1) * t266 - rSges(5,2) * t265 + rSges(5,3) * t280;
t219 = Icges(5,1) * t268 - Icges(5,4) * t267 + Icges(5,5) * t282;
t218 = Icges(5,1) * t266 - Icges(5,4) * t265 + Icges(5,5) * t280;
t215 = Icges(5,5) * t268 - Icges(5,6) * t267 + Icges(5,3) * t282;
t214 = Icges(5,5) * t266 - Icges(5,6) * t265 + Icges(5,3) * t280;
t213 = rSges(6,1) * t264 + rSges(6,2) * t263 + rSges(6,3) * t284;
t212 = Icges(6,1) * t264 + Icges(6,4) * t263 + Icges(6,5) * t284;
t211 = Icges(6,4) * t264 + Icges(6,2) * t263 + Icges(6,6) * t284;
t208 = rSges(7,1) * t260 + rSges(7,2) * t259 + rSges(7,3) * t284;
t207 = Icges(7,1) * t260 + Icges(7,4) * t259 + Icges(7,5) * t284;
t206 = Icges(7,4) * t260 + Icges(7,2) * t259 + Icges(7,6) * t284;
t205 = Icges(7,5) * t260 + Icges(7,6) * t259 + Icges(7,3) * t284;
t204 = pkin(5) * t396 + pkin(11) * t284 + t285 * t402;
t202 = rSges(6,1) * t235 + rSges(6,2) * t234 + rSges(6,3) * t267;
t201 = rSges(6,1) * t233 + rSges(6,2) * t232 + rSges(6,3) * t265;
t200 = Icges(6,1) * t235 + Icges(6,4) * t234 + Icges(6,5) * t267;
t199 = Icges(6,1) * t233 + Icges(6,4) * t232 + Icges(6,5) * t265;
t198 = Icges(6,4) * t235 + Icges(6,2) * t234 + Icges(6,6) * t267;
t197 = Icges(6,4) * t233 + Icges(6,2) * t232 + Icges(6,6) * t265;
t194 = rSges(7,1) * t230 + rSges(7,2) * t229 + rSges(7,3) * t267;
t193 = rSges(7,1) * t228 + rSges(7,2) * t227 + rSges(7,3) * t265;
t192 = Icges(7,1) * t230 + Icges(7,4) * t229 + Icges(7,5) * t267;
t191 = Icges(7,1) * t228 + Icges(7,4) * t227 + Icges(7,5) * t265;
t190 = Icges(7,4) * t230 + Icges(7,2) * t229 + Icges(7,6) * t267;
t189 = Icges(7,4) * t228 + Icges(7,2) * t227 + Icges(7,6) * t265;
t188 = Icges(7,5) * t230 + Icges(7,6) * t229 + Icges(7,3) * t267;
t187 = Icges(7,5) * t228 + Icges(7,6) * t227 + Icges(7,3) * t265;
t186 = -t246 * t307 + t272 * t296 + t372;
t185 = t247 * t307 - t272 * t297 + t370;
t184 = pkin(5) * t397 + pkin(11) * t267 + t268 * t402;
t183 = pkin(5) * t398 + pkin(11) * t265 + t266 * t402;
t182 = t246 * t297 - t247 * t296 + t371;
t181 = -t220 * t275 + t248 * t257 + t369;
t180 = t221 * t275 - t248 * t258 + t367;
t179 = t220 * t258 - t221 * t257 + t368;
t178 = t213 * t257 + (-t201 - t225) * t275 + t366;
t177 = t202 * t275 + (-t213 - t255) * t258 + t364;
t176 = t201 * t258 + (-t202 - t226) * t257 + t365;
t175 = -t193 * t250 + t204 * t257 + t208 * t223 + (-t183 - t225) * t275 + t366;
t174 = t184 * t275 + t194 * t250 - t208 * t224 + (-t204 - t255) * t258 + t364;
t173 = (-t184 - t226) * t257 + t183 * t258 + t193 * t224 - t194 * t223 + t365;
t1 = t343 * ((t288 * t332 + t289 * t333 + t313 * t343) * t358 + ((t291 * t363 + t293 * t362) * t333 + (t290 * t363 + t292 * t362) * t332 + (t314 * t363 + t315 * t362) * t343) * t355) / 0.2e1 + t296 * ((t238 * t376 - t280 * t240 + t281 * t242) * t297 + (t237 * t376 - t280 * t239 + t281 * t241) * t296 + (t269 * t376 - t280 * t270 + t281 * t271) * t307) / 0.2e1 + t307 * ((t238 * t374 - t305 * t240 + t306 * t242) * t297 + (t237 * t374 - t305 * t239 + t306 * t241) * t296 + (t269 * t374 - t305 * t270 + t306 * t271) * t307) / 0.2e1 + t297 * ((t238 * t375 - t282 * t240 + t283 * t242) * t297 + (t237 * t375 - t282 * t239 + t283 * t241) * t296 + (t269 * t375 - t282 * t270 + t283 * t271) * t307) / 0.2e1 + t332 * ((-t289 * t394 + t291 * t319 + t293 * t320) * t333 + (-t288 * t394 + t290 * t319 + t292 * t320) * t332 + (-t313 * t394 + t314 * t319 + t315 * t320) * t343) / 0.2e1 + t333 * ((t289 * t395 + t291 * t321 + t293 * t322) * t333 + (t288 * t395 + t290 * t321 + t292 * t322) * t332 + (t313 * t395 + t314 * t321 + t315 * t322) * t343) / 0.2e1 + m(1) * (t329 ^ 2 + t330 ^ 2 + t331 ^ 2) / 0.2e1 + m(2) * (t304 ^ 2 + t309 ^ 2 + t310 ^ 2) / 0.2e1 + t250 * ((t188 * t284 + t190 * t259 + t192 * t260) * t224 + (t187 * t284 + t189 * t259 + t191 * t260) * t223 + (t205 * t284 + t206 * t259 + t207 * t260) * t250) / 0.2e1 + t224 * ((t267 * t188 + t229 * t190 + t230 * t192) * t224 + (t187 * t267 + t189 * t229 + t191 * t230) * t223 + (t205 * t267 + t206 * t229 + t207 * t230) * t250) / 0.2e1 + t223 * ((t188 * t265 + t190 * t227 + t192 * t228) * t224 + (t265 * t187 + t227 * t189 + t228 * t191) * t223 + (t205 * t265 + t206 * t227 + t207 * t228) * t250) / 0.2e1 + m(3) * (t236 ^ 2 + t251 ^ 2 + t252 ^ 2) / 0.2e1 + m(7) * (t173 ^ 2 + t174 ^ 2 + t175 ^ 2) / 0.2e1 + m(6) * (t176 ^ 2 + t177 ^ 2 + t178 ^ 2) / 0.2e1 + m(5) * (t179 ^ 2 + t180 ^ 2 + t181 ^ 2) / 0.2e1 + m(4) * (t182 ^ 2 + t185 ^ 2 + t186 ^ 2) / 0.2e1 + ((t211 * t232 + t212 * t233 + t243 * t280 + t245 * t266 + t265 * t409) * t275 + (t198 * t232 + t200 * t233 + t215 * t280 + t219 * t266 + t265 * t410) * t258 + (t197 * t232 + t199 * t233 + t214 * t280 + t218 * t266 + t411 * t265) * t257) * t257 / 0.2e1 + ((t211 * t234 + t212 * t235 + t243 * t282 + t245 * t268 + t267 * t409) * t275 + (t198 * t234 + t200 * t235 + t215 * t282 + t219 * t268 + t410 * t267) * t258 + (t197 * t234 + t199 * t235 + t214 * t282 + t218 * t268 + t267 * t411) * t257) * t258 / 0.2e1 + ((t211 * t263 + t212 * t264 + t243 * t305 + t245 * t285 + t409 * t284) * t275 + (t198 * t263 + t200 * t264 + t215 * t305 + t219 * t285 + t284 * t410) * t258 + (t197 * t263 + t199 * t264 + t214 * t305 + t218 * t285 + t284 * t411) * t257) * t275 / 0.2e1 + ((-t336 * t354 + t338 * t357 + Icges(1,4)) * V_base(5) + (-t337 * t354 + t339 * t357 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t336 * t357 + t338 * t354 + Icges(1,2)) * V_base(5) + (t337 * t357 + t339 * t354 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t354 + Icges(2,6) * t357 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t357 - Icges(2,6) * t354 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
