% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP10
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRP10_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP10_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:14:24
% EndTime: 2019-02-26 22:14:27
% DurationCPUTime: 2.78s
% Computational Cost: add. (14858->179), mult. (28394->347), div. (951->12), fcn. (36086->13), ass. (0->146)
t310 = cos(pkin(6));
t313 = sin(qJ(1));
t312 = sin(qJ(2));
t371 = qJD(2) * t312;
t349 = t313 * t371;
t372 = qJD(1) * t313;
t352 = t312 * t372;
t315 = cos(qJ(2));
t316 = cos(qJ(1));
t373 = t316 * t315;
t278 = -t310 * t352 - t349 + (qJD(2) * t310 + qJD(1)) * t373;
t374 = t316 * t312;
t375 = t313 * t315;
t299 = t310 * t374 + t375;
t311 = sin(qJ(3));
t314 = cos(qJ(3));
t309 = sin(pkin(6));
t377 = t309 * t316;
t289 = t299 * t311 + t314 * t377;
t354 = t309 * t372;
t257 = t289 * qJD(3) - t278 * t314 - t311 * t354;
t290 = -t299 * t314 + t311 * t377;
t308 = pkin(11) + qJ(5);
t306 = sin(t308);
t307 = cos(t308);
t376 = t313 * t312;
t345 = -t310 * t373 + t376;
t267 = t290 * t306 + t345 * t307;
t333 = t310 * t375 + t374;
t323 = t333 * qJD(1) + t299 * qJD(2);
t240 = t267 * qJD(5) - t257 * t307 + t323 * t306;
t339 = t345 * t306;
t268 = t290 * t307 - t339;
t414 = t268 * qJD(5) + t257 * t306 + t323 * t307;
t259 = t267 ^ 2;
t379 = t309 * t314;
t298 = t310 * t311 + t312 * t379;
t378 = t309 * t315;
t282 = t298 * t306 + t307 * t378;
t280 = 0.1e1 / t282 ^ 2;
t246 = t259 * t280 + 0.1e1;
t244 = 0.1e1 / t246;
t283 = t298 * t307 - t306 * t378;
t380 = t309 * t311;
t297 = t310 * t314 - t312 * t380;
t350 = qJD(2) * t378;
t286 = t297 * qJD(3) + t314 * t350;
t351 = t309 * t371;
t251 = t283 * qJD(5) + t286 * t306 - t307 * t351;
t279 = 0.1e1 / t282;
t383 = t267 * t280;
t337 = -t251 * t383 + t279 * t414;
t224 = t337 * t244;
t247 = atan2(t267, t282);
t242 = sin(t247);
t243 = cos(t247);
t338 = -t242 * t282 + t243 * t267;
t219 = t338 * t224 + t242 * t414 + t243 * t251;
t236 = t242 * t267 + t243 * t282;
t234 = 0.1e1 / t236 ^ 2;
t413 = t219 * t234;
t277 = -t299 * qJD(1) - t333 * qJD(2);
t300 = -t310 * t376 + t373;
t291 = -t300 * t311 + t313 * t379;
t353 = qJD(1) * t377;
t254 = t291 * qJD(3) + t277 * t314 + t311 * t353;
t292 = t300 * t314 + t313 * t380;
t270 = t292 * t307 + t333 * t306;
t276 = t310 * t349 + t352 + (-qJD(1) * t310 - qJD(2)) * t373;
t237 = t270 * qJD(5) + t254 * t306 + t276 * t307;
t269 = t292 * t306 - t333 * t307;
t398 = 0.2e1 * t269;
t233 = 0.1e1 / t236;
t408 = t233 * t413;
t347 = t398 * t408;
t412 = -t234 * t237 + t347;
t262 = 0.1e1 / t270 ^ 2;
t284 = t291 ^ 2;
t382 = t284 * t262;
t250 = 0.1e1 + t382;
t253 = -t292 * qJD(3) - t277 * t311 + t314 * t353;
t238 = -t269 * qJD(5) + t254 * t307 - t276 * t306;
t261 = 0.1e1 / t270;
t391 = t238 * t261 * t262;
t355 = t284 * t391;
t385 = t262 * t291;
t395 = (t253 * t385 - t355) / t250 ^ 2;
t409 = 0.2e1 * t395;
t407 = -0.2e1 * t269;
t406 = t251 * t280;
t334 = t279 * t289 - t297 * t383;
t405 = t306 * t334;
t248 = 0.1e1 / t250;
t387 = t248 * t262;
t403 = t238 * t387 + t261 * t409;
t260 = t269 ^ 2;
t232 = t260 * t234 + 0.1e1;
t230 = 0.1e1 / t232;
t392 = t234 * t269;
t397 = (t237 * t392 - t260 * t408) / t232 ^ 2;
t402 = -t230 * t413 - 0.2e1 * t233 * t397;
t344 = t385 * t395;
t358 = t291 * t391;
t401 = 0.2e1 * t248 * t358 - t253 * t387 + 0.2e1 * t344;
t367 = 0.2e1 * t397;
t400 = t412 * t230 + t367 * t392;
t255 = t290 * qJD(3) - t278 * t311 + t314 * t354;
t399 = 0.2e1 * t267;
t386 = t279 * t406;
t396 = (-t259 * t386 + t383 * t414) / t246 ^ 2;
t394 = t230 * t233;
t390 = t242 * t269;
t389 = t243 * t269;
t388 = t248 * t261;
t384 = t267 * t279;
t381 = t291 * t306;
t370 = qJD(3) * t311;
t369 = qJD(5) * t307;
t368 = t314 * qJD(5);
t366 = -0.2e1 * t396;
t363 = t279 * t396;
t360 = t230 * t392;
t356 = t248 * t385;
t346 = t386 * t399;
t336 = t268 * t279 - t283 * t383;
t271 = -t299 * t307 - t314 * t339;
t293 = (t306 * t314 * t315 - t307 * t312) * t309;
t335 = -t271 * t279 - t293 * t383;
t329 = t314 * t333;
t328 = -t242 + (-t243 * t384 + t242) * t244;
t327 = qJD(3) * t333;
t326 = -qJD(5) * t329 - t277;
t325 = t300 * qJD(5) + t276 * t314 + t311 * t327;
t285 = -t298 * qJD(3) - t311 * t350;
t258 = ((-qJD(2) + t368) * t315 * t307 + (-t315 * t370 + (-qJD(2) * t314 + qJD(5)) * t312) * t306) * t309;
t252 = -t282 * qJD(5) + t286 * t307 + t306 * t351;
t241 = (-t345 * t368 - t278) * t307 + (t299 * qJD(5) - t314 * t323 + t345 * t370) * t306;
t229 = t244 * t405;
t228 = t335 * t244;
t227 = t336 * t244;
t222 = (t242 * t289 + t243 * t297) * t306 + t338 * t229;
t220 = t338 * t227 + t242 * t268 + t243 * t283;
t218 = t335 * t366 + (t293 * t346 - t241 * t279 + (t251 * t271 - t258 * t267 - t293 * t414) * t280) * t244;
t216 = t336 * t366 + (t283 * t346 - t240 * t279 + (-t251 * t268 - t252 * t267 - t283 * t414) * t280) * t244;
t215 = t366 * t405 + (t334 * t369 + (t297 * t346 - t255 * t279 + (-t251 * t289 - t267 * t285 - t297 * t414) * t280) * t306) * t244;
t1 = [t363 * t398 + (-t237 * t279 + t269 * t406) * t244, t218, t215, 0, t216, 0; t414 * t394 - (t328 * t237 + ((t224 * t244 * t384 + t366) * t242 + (t363 * t399 - t224 + (t224 - t337) * t244) * t243) * t269) * t360 + t402 * t267 + t400 * t328 * t269 (t325 * t306 + t326 * t307) * t394 - ((t218 * t267 + t228 * t414 + t258 + (-t228 * t282 - t271) * t224) * t243 + (-t218 * t282 - t228 * t251 - t241 + (-t228 * t267 - t293) * t224) * t242) * t360 + t402 * (-t300 * t307 - t306 * t329) + t400 * (t338 * t228 - t242 * t271 + t243 * t293) (t222 * t392 - t233 * t381) * t367 + ((t253 * t306 + t291 * t369) * t233 + t412 * t222 + (-t381 * t219 - (t297 * t369 + t215 * t267 + t229 * t414 + t285 * t306 + (-t229 * t282 + t289 * t306) * t224) * t389 - (t289 * t369 - t215 * t282 - t229 * t251 - t255 * t306 + (-t229 * t267 - t297 * t306) * t224) * t390) * t234) * t230, 0 (t220 * t392 - t233 * t270) * t367 + (t220 * t347 + t238 * t233 + (-t270 * t219 - t220 * t237 - (t216 * t267 + t227 * t414 + t252 + (-t227 * t282 + t268) * t224) * t389 - (-t216 * t282 - t227 * t251 - t240 + (-t227 * t267 - t283) * t224) * t390) * t234) * t230, 0; t240 * t356 - t255 * t388 + t401 * t268 - t403 * t289 -(-t326 * t306 + t325 * t307) * t356 + (-t276 * t311 + t314 * t327) * t388 - t403 * t311 * t333 + t401 * (t300 * t306 - t307 * t329) (t261 * t292 + t307 * t382) * t409 + (0.2e1 * t307 * t355 - t254 * t261 + (qJD(5) * t284 * t306 - 0.2e1 * t253 * t291 * t307 + t238 * t292) * t262) * t248, 0, t344 * t407 + (t358 * t407 + (t237 * t291 + t253 * t269) * t262) * t248, 0;];
JaD_rot  = t1;
