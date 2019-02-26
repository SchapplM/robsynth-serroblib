% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPP8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPP8_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP8_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:29:35
% EndTime: 2019-02-26 22:29:38
% DurationCPUTime: 2.84s
% Computational Cost: add. (9833->178), mult. (28394->346), div. (951->12), fcn. (36086->13), ass. (0->147)
t300 = cos(pkin(6));
t304 = sin(qJ(1));
t303 = sin(qJ(2));
t363 = qJD(2) * t303;
t341 = t304 * t363;
t364 = qJD(1) * t304;
t344 = t303 * t364;
t307 = cos(qJ(2));
t308 = cos(qJ(1));
t365 = t307 * t308;
t271 = -t300 * t344 - t341 + (qJD(2) * t300 + qJD(1)) * t365;
t367 = t304 * t307;
t368 = t303 * t308;
t292 = t300 * t368 + t367;
t302 = sin(qJ(3));
t306 = cos(qJ(3));
t299 = sin(pkin(6));
t371 = t299 * t308;
t282 = t292 * t302 + t306 * t371;
t346 = t299 * t364;
t248 = t282 * qJD(3) - t271 * t306 - t302 * t346;
t283 = -t292 * t306 + t302 * t371;
t301 = sin(qJ(4));
t305 = cos(qJ(4));
t369 = t303 * t304;
t337 = -t300 * t365 + t369;
t260 = t283 * t301 + t337 * t305;
t325 = t300 * t367 + t368;
t315 = t325 * qJD(1) + t292 * qJD(2);
t233 = t260 * qJD(4) - t248 * t305 + t315 * t301;
t331 = t337 * t301;
t261 = t283 * t305 - t331;
t407 = t261 * qJD(4) + t248 * t301 + t315 * t305;
t252 = t260 ^ 2;
t372 = t299 * t306;
t291 = t300 * t302 + t303 * t372;
t366 = t305 * t307;
t278 = t291 * t301 + t299 * t366;
t273 = 0.1e1 / t278 ^ 2;
t241 = t252 * t273 + 0.1e1;
t237 = 0.1e1 / t241;
t373 = t299 * t302;
t290 = t300 * t306 - t303 * t373;
t342 = qJD(2) * t299 * t307;
t277 = t290 * qJD(3) + t306 * t342;
t370 = t301 * t307;
t279 = t291 * t305 - t299 * t370;
t343 = t299 * t363;
t249 = t279 * qJD(4) + t277 * t301 - t305 * t343;
t272 = 0.1e1 / t278;
t375 = t260 * t273;
t329 = -t249 * t375 + t272 * t407;
t217 = t329 * t237;
t242 = atan2(t260, t278);
t235 = sin(t242);
t236 = cos(t242);
t330 = -t235 * t278 + t236 * t260;
t212 = t330 * t217 + t235 * t407 + t236 * t249;
t229 = t235 * t260 + t236 * t278;
t227 = 0.1e1 / t229 ^ 2;
t406 = t212 * t227;
t293 = -t300 * t369 + t365;
t285 = t293 * t306 + t304 * t373;
t262 = t285 * t301 - t325 * t305;
t391 = 0.2e1 * t262;
t226 = 0.1e1 / t229;
t401 = t226 * t406;
t339 = t391 * t401;
t270 = -t292 * qJD(1) - t325 * qJD(2);
t284 = -t293 * t302 + t304 * t372;
t345 = qJD(1) * t371;
t245 = t284 * qJD(3) + t270 * t306 + t302 * t345;
t263 = t285 * t305 + t325 * t301;
t269 = t300 * t341 + t344 + (-qJD(1) * t300 - qJD(2)) * t365;
t230 = t263 * qJD(4) + t245 * t301 + t269 * t305;
t385 = t230 * t227;
t405 = -t385 + t339;
t255 = 0.1e1 / t263 ^ 2;
t275 = t284 ^ 2;
t378 = t255 * t275;
t243 = 0.1e1 + t378;
t244 = -t285 * qJD(3) - t270 * t302 + t306 * t345;
t231 = -t262 * qJD(4) + t245 * t305 - t269 * t301;
t254 = 0.1e1 / t263;
t384 = t231 * t254 * t255;
t350 = t275 * t384;
t377 = t255 * t284;
t389 = (t244 * t377 - t350) / t243 ^ 2;
t402 = 0.2e1 * t389;
t400 = -0.2e1 * t262;
t399 = t249 * t273;
t326 = t272 * t282 - t290 * t375;
t398 = t301 * t326;
t239 = 0.1e1 / t243;
t380 = t239 * t255;
t396 = t231 * t380 + t254 * t402;
t253 = t262 ^ 2;
t225 = t227 * t253 + 0.1e1;
t223 = 0.1e1 / t225;
t390 = (-t253 * t401 + t262 * t385) / t225 ^ 2;
t395 = -t223 * t406 - 0.2e1 * t226 * t390;
t336 = t377 * t389;
t349 = t284 * t384;
t394 = 0.2e1 * t239 * t349 - t244 * t380 + 0.2e1 * t336;
t359 = 0.2e1 * t390;
t386 = t227 * t262;
t393 = t405 * t223 + t359 * t386;
t246 = t283 * qJD(3) - t271 * t302 + t306 * t346;
t392 = 0.2e1 * t260;
t379 = t272 * t399;
t388 = (-t252 * t379 + t375 * t407) / t241 ^ 2;
t387 = t223 * t226;
t383 = t235 * t262;
t382 = t236 * t262;
t381 = t239 * t254;
t376 = t260 * t272;
t374 = t284 * t301;
t362 = qJD(3) * t302;
t361 = qJD(4) * t305;
t360 = qJD(4) * t306;
t358 = -0.2e1 * t388;
t354 = t272 * t388;
t352 = t223 * t386;
t347 = t239 * t377;
t338 = t379 * t392;
t328 = t261 * t272 - t279 * t375;
t264 = -t292 * t305 - t306 * t331;
t286 = (-t303 * t305 + t306 * t370) * t299;
t327 = -t264 * t272 - t286 * t375;
t321 = t306 * t325;
t320 = -t235 + (-t236 * t376 + t235) * t237;
t319 = qJD(3) * t325;
t318 = -qJD(4) * t321 - t270;
t317 = t293 * qJD(4) + t269 * t306 + t302 * t319;
t276 = -t291 * qJD(3) - t302 * t342;
t251 = ((-qJD(2) + t360) * t366 + (-t307 * t362 + (-qJD(2) * t306 + qJD(4)) * t303) * t301) * t299;
t250 = -t278 * qJD(4) + t277 * t305 + t301 * t343;
t234 = (-t337 * t360 - t271) * t305 + (t292 * qJD(4) - t306 * t315 + t337 * t362) * t301;
t222 = t237 * t398;
t221 = t327 * t237;
t220 = t328 * t237;
t215 = (t235 * t282 + t236 * t290) * t301 + t330 * t222;
t213 = t330 * t220 + t235 * t261 + t236 * t279;
t211 = t327 * t358 + (t286 * t338 - t234 * t272 + (t249 * t264 - t251 * t260 - t286 * t407) * t273) * t237;
t209 = t328 * t358 + (t279 * t338 - t233 * t272 + (-t249 * t261 - t250 * t260 - t279 * t407) * t273) * t237;
t208 = t358 * t398 + (t326 * t361 + (t290 * t338 - t246 * t272 + (-t249 * t282 - t260 * t276 - t290 * t407) * t273) * t301) * t237;
t1 = [t354 * t391 + (-t230 * t272 + t262 * t399) * t237, t211, t208, t209, 0, 0; t407 * t387 - (t320 * t230 + ((t217 * t237 * t376 + t358) * t235 + (t354 * t392 - t217 + (t217 - t329) * t237) * t236) * t262) * t352 + t395 * t260 + t393 * t320 * t262 (t317 * t301 + t318 * t305) * t387 - ((t211 * t260 + t221 * t407 + t251 + (-t221 * t278 - t264) * t217) * t236 + (-t211 * t278 - t221 * t249 - t234 + (-t221 * t260 - t286) * t217) * t235) * t352 + t395 * (-t293 * t305 - t301 * t321) + t393 * (t330 * t221 - t235 * t264 + t236 * t286) (t215 * t386 - t226 * t374) * t359 + ((t244 * t301 + t284 * t361) * t226 + t405 * t215 + (-t374 * t212 - (t290 * t361 + t208 * t260 + t222 * t407 + t276 * t301 + (-t222 * t278 + t282 * t301) * t217) * t382 - (t282 * t361 - t208 * t278 - t222 * t249 - t246 * t301 + (-t222 * t260 - t290 * t301) * t217) * t383) * t227) * t223 (t213 * t386 - t226 * t263) * t359 + (t213 * t339 + t231 * t226 + (-t263 * t212 - t213 * t230 - (t209 * t260 + t220 * t407 + t250 + (-t220 * t278 + t261) * t217) * t382 - (-t209 * t278 - t220 * t249 - t233 + (-t220 * t260 - t279) * t217) * t383) * t227) * t223, 0, 0; t233 * t347 - t246 * t381 + t394 * t261 - t396 * t282 -(-t318 * t301 + t317 * t305) * t347 + (-t269 * t302 + t306 * t319) * t381 - t396 * t302 * t325 + t394 * (t293 * t301 - t305 * t321) (t254 * t285 + t305 * t378) * t402 + (0.2e1 * t305 * t350 - t245 * t254 + (qJD(4) * t275 * t301 - 0.2e1 * t244 * t284 * t305 + t231 * t285) * t255) * t239, t336 * t400 + (t349 * t400 + (t230 * t284 + t244 * t262) * t255) * t239, 0, 0;];
JaD_rot  = t1;
