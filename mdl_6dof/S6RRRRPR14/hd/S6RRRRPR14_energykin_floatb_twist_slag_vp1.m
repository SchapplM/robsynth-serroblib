% Calculate kinetic energy for
% S6RRRRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-10 00:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR14_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR14_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPR14_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_energykin_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR14_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR14_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR14_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:08:35
% EndTime: 2019-03-10 00:08:41
% DurationCPUTime: 5.35s
% Computational Cost: add. (5703->442), mult. (14426->645), div. (0->0), fcn. (18456->16), ass. (0->194)
t415 = Icges(5,2) + Icges(6,3);
t359 = cos(pkin(6));
t365 = cos(qJ(2));
t366 = cos(qJ(1));
t392 = t365 * t366;
t363 = sin(qJ(2));
t364 = sin(qJ(1));
t394 = t364 * t363;
t322 = t359 * t392 - t394;
t393 = t364 * t365;
t395 = t363 * t366;
t323 = t359 * t395 + t393;
t362 = sin(qJ(3));
t357 = sin(pkin(6));
t403 = sin(pkin(7));
t385 = t357 * t403;
t404 = cos(pkin(7));
t408 = cos(qJ(3));
t285 = t323 * t408 + (t322 * t404 - t366 * t385) * t362;
t361 = sin(qJ(4));
t386 = t357 * t404;
t378 = -t322 * t403 - t366 * t386;
t407 = cos(qJ(4));
t268 = t285 * t407 + t361 * t378;
t383 = t408 * t403;
t382 = t357 * t383;
t384 = t404 * t408;
t284 = -t322 * t384 + t323 * t362 + t366 * t382;
t356 = sin(pkin(13));
t358 = cos(pkin(13));
t234 = -t268 * t356 + t284 * t358;
t401 = t284 * t356;
t235 = t268 * t358 + t401;
t267 = t285 * t361 - t378 * t407;
t414 = -Icges(5,4) * t268 + Icges(6,5) * t235 - Icges(5,6) * t284 + Icges(6,6) * t234 + t267 * t415;
t324 = -t359 * t393 - t395;
t325 = -t359 * t394 + t392;
t287 = t325 * t408 + (t324 * t404 + t364 * t385) * t362;
t377 = -t324 * t403 + t364 * t386;
t270 = t287 * t407 + t361 * t377;
t286 = -t324 * t384 + t325 * t362 - t364 * t382;
t236 = -t270 * t356 + t286 * t358;
t400 = t286 * t356;
t237 = t270 * t358 + t400;
t269 = t287 * t361 - t377 * t407;
t413 = -Icges(5,4) * t270 + Icges(6,5) * t237 - Icges(5,6) * t286 + Icges(6,6) * t236 + t269 * t415;
t307 = t359 * t403 * t362 + (t362 * t365 * t404 + t363 * t408) * t357;
t376 = t359 * t404 - t365 * t385;
t283 = t307 * t407 + t361 * t376;
t398 = t357 * t363;
t306 = -t357 * t365 * t384 - t359 * t383 + t362 * t398;
t263 = -t283 * t356 + t306 * t358;
t399 = t306 * t356;
t264 = t283 * t358 + t399;
t282 = t307 * t361 - t376 * t407;
t412 = -Icges(5,4) * t283 + Icges(6,5) * t264 - Icges(5,6) * t306 + Icges(6,6) * t263 + t282 * t415;
t406 = pkin(9) * t359;
t405 = pkin(5) * t358;
t402 = Icges(2,4) * t364;
t397 = t357 * t364;
t396 = t357 * t366;
t390 = qJD(2) * t357;
t389 = V_base(5) * pkin(8) + V_base(1);
t335 = t364 * t390 + V_base(4);
t352 = V_base(6) + qJD(1);
t299 = qJD(3) * t377 + t335;
t336 = qJD(2) * t359 + t352;
t260 = qJD(4) * t286 + t299;
t308 = qJD(3) * t376 + t336;
t334 = -t366 * t390 + V_base(5);
t327 = pkin(1) * t364 - pkin(9) * t396;
t381 = -t327 * t352 + t406 * V_base(5) + t389;
t276 = qJD(4) * t306 + t308;
t328 = pkin(1) * t366 + pkin(9) * t397;
t380 = t327 * V_base(4) - t328 * V_base(5) + V_base(3);
t298 = qJD(3) * t378 + t334;
t259 = qJD(4) * t284 + t298;
t379 = t352 * t328 + V_base(2) + (-pkin(8) - t406) * V_base(4);
t288 = pkin(2) * t323 + pkin(10) * t378;
t312 = pkin(2) * t398 + pkin(10) * t376;
t375 = -t288 * t336 + t312 * t334 + t381;
t289 = pkin(2) * t325 + pkin(10) * t377;
t374 = t288 * t335 - t289 * t334 + t380;
t373 = t289 * t336 - t312 * t335 + t379;
t256 = pkin(3) * t285 + pkin(11) * t284;
t275 = pkin(3) * t307 + pkin(11) * t306;
t372 = -t256 * t308 + t275 * t298 + t375;
t257 = pkin(3) * t287 + pkin(11) * t286;
t371 = t256 * t299 - t257 * t298 + t374;
t255 = pkin(4) * t283 + qJ(5) * t282;
t370 = qJD(5) * t269 + t255 * t259 + t372;
t227 = pkin(4) * t268 + qJ(5) * t267;
t369 = qJD(5) * t282 + t227 * t260 + t371;
t368 = t257 * t308 - t275 * t299 + t373;
t228 = pkin(4) * t270 + qJ(5) * t269;
t367 = qJD(5) * t267 + t228 * t276 + t368;
t355 = pkin(13) + qJ(6);
t353 = Icges(2,4) * t366;
t351 = cos(t355);
t350 = sin(t355);
t344 = rSges(2,1) * t366 - rSges(2,2) * t364;
t343 = rSges(2,1) * t364 + rSges(2,2) * t366;
t342 = Icges(2,1) * t366 - t402;
t341 = Icges(2,1) * t364 + t353;
t340 = -Icges(2,2) * t364 + t353;
t339 = Icges(2,2) * t366 + t402;
t333 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t332 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t331 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t318 = rSges(3,3) * t359 + (rSges(3,1) * t363 + rSges(3,2) * t365) * t357;
t317 = Icges(3,5) * t359 + (Icges(3,1) * t363 + Icges(3,4) * t365) * t357;
t316 = Icges(3,6) * t359 + (Icges(3,4) * t363 + Icges(3,2) * t365) * t357;
t315 = Icges(3,3) * t359 + (Icges(3,5) * t363 + Icges(3,6) * t365) * t357;
t311 = V_base(5) * rSges(2,3) - t343 * t352 + t389;
t310 = t344 * t352 + V_base(2) + (-rSges(2,3) - pkin(8)) * V_base(4);
t309 = t343 * V_base(4) - t344 * V_base(5) + V_base(3);
t297 = rSges(3,1) * t325 + rSges(3,2) * t324 + rSges(3,3) * t397;
t296 = rSges(3,1) * t323 + rSges(3,2) * t322 - rSges(3,3) * t396;
t295 = Icges(3,1) * t325 + Icges(3,4) * t324 + Icges(3,5) * t397;
t294 = Icges(3,1) * t323 + Icges(3,4) * t322 - Icges(3,5) * t396;
t293 = Icges(3,4) * t325 + Icges(3,2) * t324 + Icges(3,6) * t397;
t292 = Icges(3,4) * t323 + Icges(3,2) * t322 - Icges(3,6) * t396;
t291 = Icges(3,5) * t325 + Icges(3,6) * t324 + Icges(3,3) * t397;
t290 = Icges(3,5) * t323 + Icges(3,6) * t322 - Icges(3,3) * t396;
t274 = rSges(4,1) * t307 - rSges(4,2) * t306 + rSges(4,3) * t376;
t273 = Icges(4,1) * t307 - Icges(4,4) * t306 + Icges(4,5) * t376;
t272 = Icges(4,4) * t307 - Icges(4,2) * t306 + Icges(4,6) * t376;
t271 = Icges(4,5) * t307 - Icges(4,6) * t306 + Icges(4,3) * t376;
t262 = t283 * t351 + t306 * t350;
t261 = -t283 * t350 + t306 * t351;
t254 = -t296 * t336 + t318 * t334 + t381;
t253 = t297 * t336 - t318 * t335 + t379;
t252 = qJD(6) * t282 + t276;
t250 = rSges(4,1) * t287 - rSges(4,2) * t286 + rSges(4,3) * t377;
t249 = rSges(4,1) * t285 - rSges(4,2) * t284 + rSges(4,3) * t378;
t248 = t296 * t335 - t297 * t334 + t380;
t247 = Icges(4,1) * t287 - Icges(4,4) * t286 + Icges(4,5) * t377;
t246 = Icges(4,1) * t285 - Icges(4,4) * t284 + Icges(4,5) * t378;
t245 = Icges(4,4) * t287 - Icges(4,2) * t286 + Icges(4,6) * t377;
t244 = Icges(4,4) * t285 - Icges(4,2) * t284 + Icges(4,6) * t378;
t243 = Icges(4,5) * t287 - Icges(4,6) * t286 + Icges(4,3) * t377;
t242 = Icges(4,5) * t285 - Icges(4,6) * t284 + Icges(4,3) * t378;
t241 = rSges(5,1) * t283 - rSges(5,2) * t282 + rSges(5,3) * t306;
t240 = Icges(5,1) * t283 - Icges(5,4) * t282 + Icges(5,5) * t306;
t238 = Icges(5,5) * t283 - Icges(5,6) * t282 + Icges(5,3) * t306;
t232 = t270 * t351 + t286 * t350;
t231 = -t270 * t350 + t286 * t351;
t230 = t268 * t351 + t284 * t350;
t229 = -t268 * t350 + t284 * t351;
t226 = qJD(6) * t269 + t260;
t225 = qJD(6) * t267 + t259;
t223 = rSges(5,1) * t270 - rSges(5,2) * t269 + rSges(5,3) * t286;
t222 = rSges(5,1) * t268 - rSges(5,2) * t267 + rSges(5,3) * t284;
t221 = Icges(5,1) * t270 - Icges(5,4) * t269 + Icges(5,5) * t286;
t220 = Icges(5,1) * t268 - Icges(5,4) * t267 + Icges(5,5) * t284;
t217 = Icges(5,5) * t270 - Icges(5,6) * t269 + Icges(5,3) * t286;
t216 = Icges(5,5) * t268 - Icges(5,6) * t267 + Icges(5,3) * t284;
t214 = rSges(6,1) * t264 + rSges(6,2) * t263 + rSges(6,3) * t282;
t213 = Icges(6,1) * t264 + Icges(6,4) * t263 + Icges(6,5) * t282;
t212 = Icges(6,4) * t264 + Icges(6,2) * t263 + Icges(6,6) * t282;
t210 = rSges(7,1) * t262 + rSges(7,2) * t261 + rSges(7,3) * t282;
t209 = Icges(7,1) * t262 + Icges(7,4) * t261 + Icges(7,5) * t282;
t208 = Icges(7,4) * t262 + Icges(7,2) * t261 + Icges(7,6) * t282;
t207 = Icges(7,5) * t262 + Icges(7,6) * t261 + Icges(7,3) * t282;
t206 = pkin(5) * t399 + pkin(12) * t282 + t283 * t405;
t204 = rSges(6,1) * t237 + rSges(6,2) * t236 + rSges(6,3) * t269;
t203 = rSges(6,1) * t235 + rSges(6,2) * t234 + rSges(6,3) * t267;
t202 = Icges(6,1) * t237 + Icges(6,4) * t236 + Icges(6,5) * t269;
t201 = Icges(6,1) * t235 + Icges(6,4) * t234 + Icges(6,5) * t267;
t200 = Icges(6,4) * t237 + Icges(6,2) * t236 + Icges(6,6) * t269;
t199 = Icges(6,4) * t235 + Icges(6,2) * t234 + Icges(6,6) * t267;
t196 = rSges(7,1) * t232 + rSges(7,2) * t231 + rSges(7,3) * t269;
t195 = rSges(7,1) * t230 + rSges(7,2) * t229 + rSges(7,3) * t267;
t194 = Icges(7,1) * t232 + Icges(7,4) * t231 + Icges(7,5) * t269;
t193 = Icges(7,1) * t230 + Icges(7,4) * t229 + Icges(7,5) * t267;
t192 = Icges(7,4) * t232 + Icges(7,2) * t231 + Icges(7,6) * t269;
t191 = Icges(7,4) * t230 + Icges(7,2) * t229 + Icges(7,6) * t267;
t190 = Icges(7,5) * t232 + Icges(7,6) * t231 + Icges(7,3) * t269;
t189 = Icges(7,5) * t230 + Icges(7,6) * t229 + Icges(7,3) * t267;
t188 = -t249 * t308 + t274 * t298 + t375;
t187 = t250 * t308 - t274 * t299 + t373;
t186 = pkin(5) * t400 + pkin(12) * t269 + t270 * t405;
t185 = pkin(5) * t401 + pkin(12) * t267 + t268 * t405;
t184 = t249 * t299 - t250 * t298 + t374;
t183 = -t222 * t276 + t241 * t259 + t372;
t182 = t223 * t276 - t241 * t260 + t368;
t181 = t222 * t260 - t223 * t259 + t371;
t180 = t214 * t259 + (-t203 - t227) * t276 + t370;
t179 = t204 * t276 + (-t214 - t255) * t260 + t367;
t178 = t203 * t260 + (-t204 - t228) * t259 + t369;
t177 = t370 - t195 * t252 + t206 * t259 + t210 * t225 + (-t185 - t227) * t276;
t176 = t186 * t276 + t196 * t252 - t210 * t226 + (-t206 - t255) * t260 + t367;
t175 = t185 * t260 + t195 * t226 - t196 * t225 + (-t186 - t228) * t259 + t369;
t1 = t334 * ((-t291 * t396 + t293 * t322 + t295 * t323) * t335 + (-t290 * t396 + t322 * t292 + t323 * t294) * t334 + (-t315 * t396 + t316 * t322 + t317 * t323) * t336) / 0.2e1 + t335 * ((t291 * t397 + t293 * t324 + t295 * t325) * t335 + (t290 * t397 + t292 * t324 + t294 * t325) * t334 + (t315 * t397 + t316 * t324 + t317 * t325) * t336) / 0.2e1 + t308 * ((t243 * t376 - t306 * t245 + t307 * t247) * t299 + (t242 * t376 - t306 * t244 + t307 * t246) * t298 + (t271 * t376 - t306 * t272 + t307 * t273) * t308) / 0.2e1 + t299 * ((t243 * t377 - t286 * t245 + t287 * t247) * t299 + (t242 * t377 - t286 * t244 + t287 * t246) * t298 + (t271 * t377 - t286 * t272 + t287 * t273) * t308) / 0.2e1 + t298 * ((t243 * t378 - t284 * t245 + t285 * t247) * t299 + (t242 * t378 - t284 * t244 + t285 * t246) * t298 + (t271 * t378 - t284 * t272 + t285 * t273) * t308) / 0.2e1 + m(1) * (t331 ^ 2 + t332 ^ 2 + t333 ^ 2) / 0.2e1 + m(2) * (t309 ^ 2 + t310 ^ 2 + t311 ^ 2) / 0.2e1 + t252 * ((t190 * t282 + t192 * t261 + t194 * t262) * t226 + (t189 * t282 + t191 * t261 + t193 * t262) * t225 + (t207 * t282 + t208 * t261 + t209 * t262) * t252) / 0.2e1 + t226 * ((t190 * t269 + t192 * t231 + t194 * t232) * t226 + (t189 * t269 + t191 * t231 + t193 * t232) * t225 + (t207 * t269 + t208 * t231 + t209 * t232) * t252) / 0.2e1 + t225 * ((t190 * t267 + t192 * t229 + t194 * t230) * t226 + (t189 * t267 + t191 * t229 + t193 * t230) * t225 + (t207 * t267 + t208 * t229 + t209 * t230) * t252) / 0.2e1 + m(3) * (t248 ^ 2 + t253 ^ 2 + t254 ^ 2) / 0.2e1 + m(7) * (t175 ^ 2 + t176 ^ 2 + t177 ^ 2) / 0.2e1 + m(6) * (t178 ^ 2 + t179 ^ 2 + t180 ^ 2) / 0.2e1 + m(5) * (t181 ^ 2 + t182 ^ 2 + t183 ^ 2) / 0.2e1 + m(4) * (t184 ^ 2 + t187 ^ 2 + t188 ^ 2) / 0.2e1 + t336 * ((t290 * t334 + t291 * t335 + t315 * t336) * t359 + ((t293 * t365 + t295 * t363) * t335 + (t292 * t365 + t294 * t363) * t334 + (t316 * t365 + t317 * t363) * t336) * t357) / 0.2e1 + ((t212 * t234 + t213 * t235 + t238 * t284 + t240 * t268 + t267 * t412) * t276 + (t200 * t234 + t202 * t235 + t217 * t284 + t221 * t268 + t267 * t413) * t260 + (t199 * t234 + t201 * t235 + t216 * t284 + t220 * t268 + t414 * t267) * t259) * t259 / 0.2e1 + ((t212 * t236 + t213 * t237 + t238 * t286 + t240 * t270 + t269 * t412) * t276 + (t200 * t236 + t202 * t237 + t217 * t286 + t221 * t270 + t413 * t269) * t260 + (t199 * t236 + t201 * t237 + t216 * t286 + t220 * t270 + t269 * t414) * t259) * t260 / 0.2e1 + ((t212 * t263 + t213 * t264 + t238 * t306 + t240 * t283 + t412 * t282) * t276 + (t200 * t263 + t202 * t264 + t217 * t306 + t221 * t283 + t282 * t413) * t260 + (t199 * t263 + t201 * t264 + t216 * t306 + t220 * t283 + t282 * t414) * t259) * t276 / 0.2e1 + ((-t339 * t364 + t341 * t366 + Icges(1,4)) * V_base(5) + (-t364 * t340 + t342 * t366 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t339 * t366 + t364 * t341 + Icges(1,2)) * V_base(5) + (t340 * t366 + t342 * t364 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t364 + Icges(2,6) * t366) * V_base(5) + (Icges(2,5) * t366 - Icges(2,6) * t364) * V_base(4) + Icges(2,3) * t352 / 0.2e1) * t352;
T  = t1;
