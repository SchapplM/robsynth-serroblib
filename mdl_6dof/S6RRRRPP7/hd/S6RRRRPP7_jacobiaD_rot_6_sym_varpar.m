% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPP7_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP7_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:28:55
% EndTime: 2019-02-26 22:28:58
% DurationCPUTime: 2.91s
% Computational Cost: add. (14858->179), mult. (28394->347), div. (951->12), fcn. (36086->13), ass. (0->147)
t311 = cos(pkin(6));
t314 = sin(qJ(1));
t313 = sin(qJ(2));
t372 = qJD(2) * t313;
t350 = t314 * t372;
t373 = qJD(1) * t314;
t353 = t313 * t373;
t316 = cos(qJ(2));
t317 = cos(qJ(1));
t374 = t316 * t317;
t279 = -t311 * t353 - t350 + (qJD(2) * t311 + qJD(1)) * t374;
t375 = t314 * t316;
t376 = t313 * t317;
t300 = t311 * t376 + t375;
t312 = sin(qJ(3));
t315 = cos(qJ(3));
t310 = sin(pkin(6));
t378 = t310 * t317;
t290 = t300 * t312 + t315 * t378;
t355 = t310 * t373;
t258 = qJD(3) * t290 - t279 * t315 - t312 * t355;
t291 = -t300 * t315 + t312 * t378;
t309 = qJ(4) + pkin(11);
t307 = sin(t309);
t308 = cos(t309);
t377 = t313 * t314;
t346 = -t311 * t374 + t377;
t268 = t291 * t307 + t308 * t346;
t334 = t311 * t375 + t376;
t324 = qJD(1) * t334 + qJD(2) * t300;
t241 = qJD(4) * t268 - t258 * t308 + t307 * t324;
t340 = t346 * t307;
t269 = t291 * t308 - t340;
t415 = qJD(4) * t269 + t258 * t307 + t308 * t324;
t260 = t268 ^ 2;
t380 = t310 * t315;
t299 = t311 * t312 + t313 * t380;
t379 = t310 * t316;
t283 = t299 * t307 + t308 * t379;
t281 = 0.1e1 / t283 ^ 2;
t247 = t260 * t281 + 0.1e1;
t245 = 0.1e1 / t247;
t284 = t299 * t308 - t307 * t379;
t381 = t310 * t312;
t298 = t311 * t315 - t313 * t381;
t351 = qJD(2) * t379;
t287 = qJD(3) * t298 + t315 * t351;
t352 = t310 * t372;
t252 = qJD(4) * t284 + t287 * t307 - t308 * t352;
t280 = 0.1e1 / t283;
t383 = t268 * t281;
t338 = -t252 * t383 + t280 * t415;
t225 = t338 * t245;
t248 = atan2(t268, t283);
t243 = sin(t248);
t244 = cos(t248);
t339 = -t243 * t283 + t244 * t268;
t220 = t225 * t339 + t243 * t415 + t244 * t252;
t237 = t243 * t268 + t244 * t283;
t235 = 0.1e1 / t237 ^ 2;
t414 = t220 * t235;
t301 = -t311 * t377 + t374;
t293 = t301 * t315 + t314 * t381;
t270 = t293 * t307 - t308 * t334;
t399 = 0.2e1 * t270;
t234 = 0.1e1 / t237;
t409 = t234 * t414;
t347 = t399 * t409;
t278 = -qJD(1) * t300 - qJD(2) * t334;
t292 = -t301 * t312 + t314 * t380;
t354 = qJD(1) * t378;
t255 = qJD(3) * t292 + t278 * t315 + t312 * t354;
t271 = t293 * t308 + t307 * t334;
t277 = t311 * t350 + t353 + (-qJD(1) * t311 - qJD(2)) * t374;
t238 = qJD(4) * t271 + t255 * t307 + t277 * t308;
t393 = t238 * t235;
t413 = -t393 + t347;
t263 = 0.1e1 / t271 ^ 2;
t285 = t292 ^ 2;
t386 = t263 * t285;
t251 = 0.1e1 + t386;
t254 = -qJD(3) * t293 - t278 * t312 + t315 * t354;
t239 = -qJD(4) * t270 + t255 * t308 - t277 * t307;
t262 = 0.1e1 / t271;
t392 = t239 * t262 * t263;
t360 = t285 * t392;
t385 = t263 * t292;
t396 = (t254 * t385 - t360) / t251 ^ 2;
t410 = 0.2e1 * t396;
t408 = -0.2e1 * t270;
t407 = t252 * t281;
t335 = t280 * t290 - t298 * t383;
t406 = t307 * t335;
t249 = 0.1e1 / t251;
t388 = t249 * t263;
t404 = t239 * t388 + t262 * t410;
t261 = t270 ^ 2;
t233 = t235 * t261 + 0.1e1;
t231 = 0.1e1 / t233;
t398 = (-t261 * t409 + t270 * t393) / t233 ^ 2;
t403 = -t231 * t414 - 0.2e1 * t234 * t398;
t345 = t385 * t396;
t356 = t292 * t392;
t402 = 0.2e1 * t249 * t356 - t254 * t388 + 0.2e1 * t345;
t368 = 0.2e1 * t398;
t394 = t235 * t270;
t401 = t413 * t231 + t368 * t394;
t256 = qJD(3) * t291 - t279 * t312 + t315 * t355;
t400 = 0.2e1 * t268;
t387 = t280 * t407;
t397 = (-t260 * t387 + t383 * t415) / t247 ^ 2;
t395 = t231 * t234;
t391 = t243 * t270;
t390 = t244 * t270;
t389 = t249 * t262;
t384 = t268 * t280;
t382 = t292 * t307;
t371 = qJD(3) * t312;
t370 = qJD(4) * t308;
t369 = qJD(4) * t315;
t367 = -0.2e1 * t397;
t364 = t280 * t397;
t362 = t231 * t394;
t359 = t249 * t385;
t348 = t387 * t400;
t337 = t269 * t280 - t284 * t383;
t272 = -t300 * t308 - t315 * t340;
t294 = (t307 * t315 * t316 - t308 * t313) * t310;
t336 = -t272 * t280 - t294 * t383;
t330 = t315 * t334;
t329 = -t243 + (-t244 * t384 + t243) * t245;
t328 = qJD(3) * t334;
t327 = -qJD(4) * t330 - t278;
t326 = t301 * qJD(4) + t277 * t315 + t312 * t328;
t286 = -qJD(3) * t299 - t312 * t351;
t259 = ((-qJD(2) + t369) * t316 * t308 + (-t316 * t371 + (-qJD(2) * t315 + qJD(4)) * t313) * t307) * t310;
t253 = -qJD(4) * t283 + t287 * t308 + t307 * t352;
t242 = (-t346 * t369 - t279) * t308 + (t300 * qJD(4) - t315 * t324 + t346 * t371) * t307;
t230 = t245 * t406;
t229 = t336 * t245;
t228 = t337 * t245;
t223 = (t243 * t290 + t244 * t298) * t307 + t339 * t230;
t221 = t228 * t339 + t243 * t269 + t244 * t284;
t219 = t336 * t367 + (t294 * t348 - t242 * t280 + (t252 * t272 - t259 * t268 - t294 * t415) * t281) * t245;
t217 = t337 * t367 + (t284 * t348 - t241 * t280 + (-t252 * t269 - t253 * t268 - t284 * t415) * t281) * t245;
t216 = t367 * t406 + (t335 * t370 + (t298 * t348 - t256 * t280 + (-t252 * t290 - t268 * t286 - t298 * t415) * t281) * t307) * t245;
t1 = [t364 * t399 + (-t238 * t280 + t270 * t407) * t245, t219, t216, t217, 0, 0; t415 * t395 - (t329 * t238 + ((t225 * t245 * t384 + t367) * t243 + (t364 * t400 - t225 + (t225 - t338) * t245) * t244) * t270) * t362 + t403 * t268 + t401 * t329 * t270 (t307 * t326 + t308 * t327) * t395 - ((t219 * t268 + t229 * t415 + t259 + (-t229 * t283 - t272) * t225) * t244 + (-t219 * t283 - t229 * t252 - t242 + (-t229 * t268 - t294) * t225) * t243) * t362 + t403 * (-t301 * t308 - t307 * t330) + t401 * (t229 * t339 - t243 * t272 + t244 * t294) (t223 * t394 - t234 * t382) * t368 + ((t254 * t307 + t292 * t370) * t234 + t413 * t223 + (-t382 * t220 - (t298 * t370 + t216 * t268 + t230 * t415 + t286 * t307 + (-t230 * t283 + t290 * t307) * t225) * t390 - (t290 * t370 - t216 * t283 - t230 * t252 - t256 * t307 + (-t230 * t268 - t298 * t307) * t225) * t391) * t235) * t231 (t221 * t394 - t234 * t271) * t368 + (t221 * t347 + t239 * t234 + (-t271 * t220 - t221 * t238 - (t217 * t268 + t228 * t415 + t253 + (-t228 * t283 + t269) * t225) * t390 - (-t217 * t283 - t228 * t252 - t241 + (-t228 * t268 - t284) * t225) * t391) * t235) * t231, 0, 0; t241 * t359 - t256 * t389 + t402 * t269 - t404 * t290 -(-t307 * t327 + t308 * t326) * t359 + (-t277 * t312 + t315 * t328) * t389 - t404 * t312 * t334 + t402 * (t301 * t307 - t308 * t330) (t262 * t293 + t308 * t386) * t410 + (0.2e1 * t308 * t360 - t255 * t262 + (qJD(4) * t285 * t307 - 0.2e1 * t254 * t292 * t308 + t239 * t293) * t263) * t249, t345 * t408 + (t356 * t408 + (t238 * t292 + t254 * t270) * t263) * t249, 0, 0;];
JaD_rot  = t1;
