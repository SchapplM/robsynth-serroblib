% Calculate kinetic energy for
% S6RRRRPR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-10 00:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR15_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR15_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPR15_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR15_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR15_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR15_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:34:06
% EndTime: 2019-03-10 00:34:11
% DurationCPUTime: 4.81s
% Computational Cost: add. (4921->401), mult. (12754->587), div. (0->0), fcn. (16200->14), ass. (0->180)
t414 = Icges(5,1) + Icges(6,2);
t413 = Icges(6,1) + Icges(5,3);
t412 = -Icges(5,4) - Icges(6,6);
t411 = -Icges(6,4) + Icges(5,5);
t410 = Icges(6,5) - Icges(5,6);
t409 = Icges(5,2) + Icges(6,3);
t354 = cos(pkin(6));
t361 = cos(qJ(2));
t362 = cos(qJ(1));
t384 = t361 * t362;
t358 = sin(qJ(2));
t359 = sin(qJ(1));
t386 = t359 * t358;
t323 = t354 * t384 - t386;
t385 = t359 * t361;
t387 = t358 * t362;
t324 = t354 * t387 + t385;
t357 = sin(qJ(3));
t353 = sin(pkin(6));
t392 = sin(pkin(7));
t378 = t353 * t392;
t393 = cos(pkin(7));
t396 = cos(qJ(3));
t283 = t324 * t396 + (t323 * t393 - t362 * t378) * t357;
t379 = t353 * t393;
t306 = -t323 * t392 - t362 * t379;
t356 = sin(qJ(4));
t395 = cos(qJ(4));
t264 = t283 * t356 - t306 * t395;
t265 = t283 * t395 + t306 * t356;
t376 = t396 * t392;
t375 = t353 * t376;
t377 = t393 * t396;
t282 = -t323 * t377 + t324 * t357 + t362 * t375;
t408 = t409 * t264 + t412 * t265 + t410 * t282;
t325 = -t354 * t385 - t387;
t326 = -t354 * t386 + t384;
t285 = t326 * t396 + (t325 * t393 + t359 * t378) * t357;
t307 = -t325 * t392 + t359 * t379;
t266 = t285 * t356 - t307 * t395;
t267 = t285 * t395 + t307 * t356;
t284 = -t325 * t377 + t326 * t357 - t359 * t375;
t407 = t409 * t266 + t412 * t267 + t410 * t284;
t406 = t410 * t264 + t411 * t265 + t413 * t282;
t405 = t410 * t266 + t411 * t267 + t413 * t284;
t404 = t412 * t264 + t414 * t265 + t411 * t282;
t403 = t412 * t266 + t414 * t267 + t411 * t284;
t305 = t354 * t392 * t357 + (t357 * t361 * t393 + t358 * t396) * t353;
t322 = t354 * t393 - t361 * t378;
t280 = t305 * t356 - t322 * t395;
t281 = t305 * t395 + t322 * t356;
t390 = t353 * t358;
t304 = -t353 * t361 * t377 - t354 * t376 + t357 * t390;
t402 = t409 * t280 + t412 * t281 + t410 * t304;
t401 = t410 * t280 + t411 * t281 + t413 * t304;
t400 = t412 * t280 + t414 * t281 + t411 * t304;
t394 = pkin(9) * t354;
t391 = Icges(2,4) * t359;
t389 = t353 * t359;
t388 = t353 * t362;
t383 = qJD(2) * t353;
t382 = V_base(5) * pkin(8) + V_base(1);
t336 = t359 * t383 + V_base(4);
t350 = V_base(6) + qJD(1);
t297 = qJD(3) * t307 + t336;
t337 = qJD(2) * t354 + t350;
t259 = qJD(4) * t284 + t297;
t308 = qJD(3) * t322 + t337;
t335 = -t362 * t383 + V_base(5);
t328 = t359 * pkin(1) - pkin(9) * t388;
t374 = -t328 * t350 + V_base(5) * t394 + t382;
t274 = qJD(4) * t304 + t308;
t329 = pkin(1) * t362 + pkin(9) * t389;
t373 = V_base(4) * t328 - t329 * V_base(5) + V_base(3);
t296 = qJD(3) * t306 + t335;
t258 = qJD(4) * t282 + t296;
t372 = t350 * t329 + V_base(2) + (-pkin(8) - t394) * V_base(4);
t286 = t324 * pkin(2) + pkin(10) * t306;
t312 = pkin(2) * t390 + pkin(10) * t322;
t371 = -t286 * t337 + t335 * t312 + t374;
t287 = t326 * pkin(2) + pkin(10) * t307;
t370 = t336 * t286 - t287 * t335 + t373;
t369 = t337 * t287 - t312 * t336 + t372;
t255 = pkin(3) * t283 + pkin(11) * t282;
t273 = pkin(3) * t305 + pkin(11) * t304;
t368 = -t255 * t308 + t296 * t273 + t371;
t256 = pkin(3) * t285 + pkin(11) * t284;
t367 = t297 * t255 - t256 * t296 + t370;
t254 = pkin(4) * t281 + qJ(5) * t280;
t366 = qJD(5) * t266 + t258 * t254 + t368;
t224 = pkin(4) * t265 + qJ(5) * t264;
t365 = qJD(5) * t280 + t259 * t224 + t367;
t364 = t308 * t256 - t273 * t297 + t369;
t225 = pkin(4) * t267 + qJ(5) * t266;
t363 = qJD(5) * t264 + t274 * t225 + t364;
t360 = cos(qJ(6));
t355 = sin(qJ(6));
t351 = Icges(2,4) * t362;
t345 = rSges(2,1) * t362 - t359 * rSges(2,2);
t344 = t359 * rSges(2,1) + rSges(2,2) * t362;
t343 = Icges(2,1) * t362 - t391;
t342 = Icges(2,1) * t359 + t351;
t341 = -Icges(2,2) * t359 + t351;
t340 = Icges(2,2) * t362 + t391;
t334 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t333 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t332 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t318 = rSges(3,3) * t354 + (rSges(3,1) * t358 + rSges(3,2) * t361) * t353;
t317 = Icges(3,5) * t354 + (Icges(3,1) * t358 + Icges(3,4) * t361) * t353;
t316 = Icges(3,6) * t354 + (Icges(3,4) * t358 + Icges(3,2) * t361) * t353;
t315 = Icges(3,3) * t354 + (Icges(3,5) * t358 + Icges(3,6) * t361) * t353;
t311 = V_base(5) * rSges(2,3) - t344 * t350 + t382;
t310 = t345 * t350 + V_base(2) + (-rSges(2,3) - pkin(8)) * V_base(4);
t309 = t344 * V_base(4) - t345 * V_base(5) + V_base(3);
t295 = rSges(3,1) * t326 + rSges(3,2) * t325 + rSges(3,3) * t389;
t294 = t324 * rSges(3,1) + t323 * rSges(3,2) - rSges(3,3) * t388;
t293 = Icges(3,1) * t326 + Icges(3,4) * t325 + Icges(3,5) * t389;
t292 = Icges(3,1) * t324 + Icges(3,4) * t323 - Icges(3,5) * t388;
t291 = Icges(3,4) * t326 + Icges(3,2) * t325 + Icges(3,6) * t389;
t290 = Icges(3,4) * t324 + Icges(3,2) * t323 - Icges(3,6) * t388;
t289 = Icges(3,5) * t326 + Icges(3,6) * t325 + Icges(3,3) * t389;
t288 = Icges(3,5) * t324 + Icges(3,6) * t323 - Icges(3,3) * t388;
t272 = rSges(4,1) * t305 - rSges(4,2) * t304 + rSges(4,3) * t322;
t271 = Icges(4,1) * t305 - Icges(4,4) * t304 + Icges(4,5) * t322;
t270 = Icges(4,4) * t305 - Icges(4,2) * t304 + Icges(4,6) * t322;
t269 = Icges(4,5) * t305 - Icges(4,6) * t304 + Icges(4,3) * t322;
t268 = pkin(5) * t304 + pkin(12) * t281;
t261 = t280 * t355 + t304 * t360;
t260 = t280 * t360 - t304 * t355;
t253 = -t294 * t337 + t318 * t335 + t374;
t252 = t295 * t337 - t318 * t336 + t372;
t251 = qJD(6) * t281 + t274;
t249 = rSges(4,1) * t285 - rSges(4,2) * t284 + rSges(4,3) * t307;
t248 = rSges(4,1) * t283 - rSges(4,2) * t282 + rSges(4,3) * t306;
t247 = t294 * t336 - t295 * t335 + t373;
t246 = Icges(4,1) * t285 - Icges(4,4) * t284 + Icges(4,5) * t307;
t245 = Icges(4,1) * t283 - Icges(4,4) * t282 + Icges(4,5) * t306;
t244 = Icges(4,4) * t285 - Icges(4,2) * t284 + Icges(4,6) * t307;
t243 = Icges(4,4) * t283 - Icges(4,2) * t282 + Icges(4,6) * t306;
t242 = Icges(4,5) * t285 - Icges(4,6) * t284 + Icges(4,3) * t307;
t241 = Icges(4,5) * t283 - Icges(4,6) * t282 + Icges(4,3) * t306;
t240 = rSges(5,1) * t281 - rSges(5,2) * t280 + rSges(5,3) * t304;
t239 = rSges(6,1) * t304 - rSges(6,2) * t281 + rSges(6,3) * t280;
t232 = pkin(5) * t284 + pkin(12) * t267;
t231 = pkin(5) * t282 + pkin(12) * t265;
t230 = t266 * t355 + t284 * t360;
t229 = t266 * t360 - t284 * t355;
t228 = t264 * t355 + t282 * t360;
t227 = t264 * t360 - t282 * t355;
t223 = qJD(6) * t267 + t259;
t222 = qJD(6) * t265 + t258;
t220 = rSges(5,1) * t267 - rSges(5,2) * t266 + rSges(5,3) * t284;
t219 = rSges(5,1) * t265 - rSges(5,2) * t264 + rSges(5,3) * t282;
t218 = rSges(6,1) * t284 - rSges(6,2) * t267 + rSges(6,3) * t266;
t217 = rSges(6,1) * t282 - rSges(6,2) * t265 + rSges(6,3) * t264;
t203 = rSges(7,1) * t261 + rSges(7,2) * t260 + rSges(7,3) * t281;
t202 = Icges(7,1) * t261 + Icges(7,4) * t260 + Icges(7,5) * t281;
t201 = Icges(7,4) * t261 + Icges(7,2) * t260 + Icges(7,6) * t281;
t200 = Icges(7,5) * t261 + Icges(7,6) * t260 + Icges(7,3) * t281;
t198 = rSges(7,1) * t230 + rSges(7,2) * t229 + rSges(7,3) * t267;
t197 = rSges(7,1) * t228 + rSges(7,2) * t227 + rSges(7,3) * t265;
t196 = Icges(7,1) * t230 + Icges(7,4) * t229 + Icges(7,5) * t267;
t195 = Icges(7,1) * t228 + Icges(7,4) * t227 + Icges(7,5) * t265;
t194 = Icges(7,4) * t230 + Icges(7,2) * t229 + Icges(7,6) * t267;
t193 = Icges(7,4) * t228 + Icges(7,2) * t227 + Icges(7,6) * t265;
t192 = Icges(7,5) * t230 + Icges(7,6) * t229 + Icges(7,3) * t267;
t191 = Icges(7,5) * t228 + Icges(7,6) * t227 + Icges(7,3) * t265;
t190 = -t248 * t308 + t272 * t296 + t371;
t189 = t249 * t308 - t272 * t297 + t369;
t188 = t248 * t297 - t249 * t296 + t370;
t187 = -t219 * t274 + t240 * t258 + t368;
t186 = t220 * t274 - t240 * t259 + t364;
t185 = t219 * t259 - t220 * t258 + t367;
t184 = t239 * t258 + (-t217 - t224) * t274 + t366;
t183 = t218 * t274 + (-t239 - t254) * t259 + t363;
t182 = t217 * t259 + (-t218 - t225) * t258 + t365;
t181 = t366 - t197 * t251 + t203 * t222 + t258 * t268 + (-t224 - t231) * t274;
t180 = t198 * t251 - t203 * t223 + t232 * t274 + (-t254 - t268) * t259 + t363;
t179 = t197 * t223 - t198 * t222 + t231 * t259 + (-t225 - t232) * t258 + t365;
t1 = t337 * ((t288 * t335 + t289 * t336 + t315 * t337) * t354 + ((t291 * t361 + t293 * t358) * t336 + (t290 * t361 + t292 * t358) * t335 + (t316 * t361 + t317 * t358) * t337) * t353) / 0.2e1 + m(1) * (t332 ^ 2 + t333 ^ 2 + t334 ^ 2) / 0.2e1 + t308 * ((t242 * t322 - t244 * t304 + t246 * t305) * t297 + (t241 * t322 - t243 * t304 + t245 * t305) * t296 + (t269 * t322 - t270 * t304 + t271 * t305) * t308) / 0.2e1 + t296 * ((t242 * t306 - t244 * t282 + t246 * t283) * t297 + (t241 * t306 - t243 * t282 + t245 * t283) * t296 + (t269 * t306 - t270 * t282 + t271 * t283) * t308) / 0.2e1 + t297 * ((t242 * t307 - t244 * t284 + t246 * t285) * t297 + (t241 * t307 - t243 * t284 + t245 * t285) * t296 + (t269 * t307 - t270 * t284 + t271 * t285) * t308) / 0.2e1 + m(2) * (t309 ^ 2 + t310 ^ 2 + t311 ^ 2) / 0.2e1 + t251 * ((t192 * t281 + t194 * t260 + t196 * t261) * t223 + (t191 * t281 + t193 * t260 + t195 * t261) * t222 + (t281 * t200 + t260 * t201 + t261 * t202) * t251) / 0.2e1 + t223 * ((t267 * t192 + t229 * t194 + t230 * t196) * t223 + (t191 * t267 + t193 * t229 + t195 * t230) * t222 + (t200 * t267 + t201 * t229 + t202 * t230) * t251) / 0.2e1 + t222 * ((t192 * t265 + t194 * t227 + t196 * t228) * t223 + (t265 * t191 + t227 * t193 + t228 * t195) * t222 + (t200 * t265 + t201 * t227 + t202 * t228) * t251) / 0.2e1 + m(3) * (t247 ^ 2 + t252 ^ 2 + t253 ^ 2) / 0.2e1 + m(4) * (t188 ^ 2 + t189 ^ 2 + t190 ^ 2) / 0.2e1 + m(5) * (t185 ^ 2 + t186 ^ 2 + t187 ^ 2) / 0.2e1 + m(6) * (t182 ^ 2 + t183 ^ 2 + t184 ^ 2) / 0.2e1 + m(7) * (t179 ^ 2 + t180 ^ 2 + t181 ^ 2) / 0.2e1 + t335 * ((-t289 * t388 + t323 * t291 + t324 * t293) * t336 + (-t288 * t388 + t323 * t290 + t324 * t292) * t335 + (-t315 * t388 + t323 * t316 + t324 * t317) * t337) / 0.2e1 + t336 * ((t289 * t389 + t291 * t325 + t293 * t326) * t336 + (t288 * t389 + t290 * t325 + t292 * t326) * t335 + (t315 * t389 + t316 * t325 + t317 * t326) * t337) / 0.2e1 + ((t264 * t402 + t265 * t400 + t282 * t401) * t274 + (t264 * t407 + t265 * t403 + t282 * t405) * t259 + (t264 * t408 + t404 * t265 + t406 * t282) * t258) * t258 / 0.2e1 + ((t266 * t402 + t267 * t400 + t284 * t401) * t274 + (t266 * t407 + t267 * t403 + t284 * t405) * t259 + (t266 * t408 + t404 * t267 + t406 * t284) * t258) * t259 / 0.2e1 + ((t280 * t402 + t281 * t400 + t304 * t401) * t274 + (t280 * t407 + t281 * t403 + t304 * t405) * t259 + (t280 * t408 + t404 * t281 + t406 * t304) * t258) * t274 / 0.2e1 + ((-t359 * t340 + t342 * t362 + Icges(1,4)) * V_base(5) + (-t359 * t341 + t343 * t362 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t340 * t362 + t359 * t342 + Icges(1,2)) * V_base(5) + (t341 * t362 + t359 * t343 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t359 + Icges(2,6) * t362) * V_base(5) + (Icges(2,5) * t362 - Icges(2,6) * t359) * V_base(4) + Icges(2,3) * t350 / 0.2e1) * t350;
T  = t1;
