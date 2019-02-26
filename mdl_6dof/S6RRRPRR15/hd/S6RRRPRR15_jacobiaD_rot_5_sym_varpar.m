% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR15
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR15_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR15_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:24:26
% EndTime: 2019-02-26 22:24:27
% DurationCPUTime: 1.83s
% Computational Cost: add. (8022->178), mult. (25165->336), div. (705->12), fcn. (31370->15), ass. (0->146)
t342 = sin(qJ(3));
t346 = cos(qJ(3));
t340 = cos(pkin(6));
t343 = sin(qJ(2));
t344 = sin(qJ(1));
t396 = t344 * t343;
t381 = t340 * t396;
t347 = cos(qJ(2));
t348 = cos(qJ(1));
t392 = t348 * t347;
t358 = t381 - t392;
t393 = t348 * t343;
t395 = t344 * t347;
t359 = t340 * t395 + t393;
t337 = sin(pkin(7));
t338 = sin(pkin(6));
t403 = t338 * t344;
t383 = t337 * t403;
t339 = cos(pkin(7));
t400 = t339 * t346;
t300 = -t342 * t358 - t346 * t383 + t359 * t400;
t318 = t337 * t359 + t339 * t403;
t341 = sin(qJ(5));
t345 = cos(qJ(5));
t370 = t300 * t345 - t318 * t341;
t423 = t370 * qJD(5);
t394 = t346 * t347;
t399 = t342 * t343;
t360 = -t339 * t399 + t394;
t363 = t339 * t394 - t399;
t406 = t337 * t340;
t378 = qJD(3) * t406;
t291 = t346 * t378 + (t360 * qJD(2) + t363 * qJD(3)) * t338;
t397 = t343 * t346;
t398 = t342 * t347;
t361 = t339 * t398 + t397;
t316 = t361 * t338 + t342 * t406;
t313 = 0.1e1 / t316 ^ 2;
t422 = t291 * t313;
t325 = t340 * t393 + t395;
t310 = t359 * qJD(1) + t325 * qJD(2);
t311 = -qJD(1) * t381 - qJD(2) * t396 + (qJD(2) * t340 + qJD(1)) * t392;
t324 = -t340 * t392 + t396;
t401 = t339 * t342;
t365 = t324 * t401 - t325 * t346;
t402 = t338 * t348;
t382 = t337 * t402;
t373 = qJD(3) * t382;
t391 = qJD(1) * t338;
t380 = t344 * t391;
t374 = t337 * t380;
t271 = t365 * qJD(3) - t310 * t400 + t346 * t374 + (-t311 + t373) * t342;
t299 = t342 * t382 + t365;
t287 = atan2(t299, t316);
t282 = sin(t287);
t283 = cos(t287);
t263 = t282 * t299 + t283 * t316;
t260 = 0.1e1 / t263;
t281 = t300 * t341 + t318 * t345;
t275 = 0.1e1 / t281;
t312 = 0.1e1 / t316;
t261 = 0.1e1 / t263 ^ 2;
t276 = 0.1e1 / t281 ^ 2;
t421 = 0.2e1 * t299;
t364 = -t339 * t359 + t383;
t301 = t364 * t342 - t346 * t358;
t294 = t301 ^ 2;
t257 = t294 * t261 + 0.1e1;
t309 = t325 * qJD(1) + t359 * qJD(2);
t308 = t324 * qJD(1) + t358 * qJD(2);
t379 = t348 * t391;
t357 = t308 * t339 + t337 * t379;
t389 = qJD(3) * t342;
t268 = t358 * t389 + t357 * t342 + (t364 * qJD(3) - t309) * t346;
t414 = t268 * t261;
t293 = t299 ^ 2;
t286 = t293 * t313 + 0.1e1;
t284 = 0.1e1 / t286;
t329 = t346 * t373;
t390 = qJD(3) * t324;
t270 = t342 * t374 + t311 * t346 - t325 * t389 - t329 + (-t310 * t342 - t346 * t390) * t339;
t408 = t299 * t313;
t369 = -t270 * t312 - t291 * t408;
t251 = t369 * t284;
t372 = -t282 * t316 + t283 * t299;
t246 = t372 * t251 - t282 * t270 + t283 * t291;
t419 = t246 * t260 * t261;
t420 = (-t294 * t419 + t301 * t414) / t257 ^ 2;
t267 = t301 * qJD(3) - t309 * t342 - t357 * t346;
t302 = -t308 * t337 + t339 * t379;
t258 = t281 * qJD(5) - t267 * t345 + t302 * t341;
t274 = t370 ^ 2;
t266 = t274 * t276 + 0.1e1;
t413 = t276 * t370;
t259 = t267 * t341 + t302 * t345 + t423;
t416 = t259 * t275 * t276;
t418 = (-t258 * t413 - t274 * t416) / t266 ^ 2;
t410 = t312 * t422;
t417 = (-t270 * t408 - t293 * t410) / t286 ^ 2;
t415 = t261 * t301;
t412 = t282 * t301;
t411 = t283 * t301;
t409 = t299 * t312;
t407 = t325 * t342;
t405 = t337 * t341;
t404 = t337 * t345;
t388 = 0.2e1 * t420;
t387 = 0.2e1 * t419;
t386 = 0.2e1 * t418;
t385 = -0.2e1 * t417;
t384 = t312 * t417;
t377 = t301 * t387;
t376 = -0.2e1 * t370 * t416;
t375 = t410 * t421;
t298 = -t407 + (-t324 * t339 - t382) * t346;
t317 = -t324 * t337 + t339 * t402;
t371 = t298 * t345 - t317 * t341;
t279 = t298 * t341 + t317 * t345;
t368 = t345 * t275 - t341 * t413;
t295 = t324 * t400 + t346 * t382 + t407;
t315 = t363 * t338 + t346 * t406;
t367 = t295 * t312 - t315 * t408;
t305 = -t324 * t346 - t325 * t401;
t321 = t360 * t338;
t366 = -t305 * t312 - t321 * t408;
t306 = -t342 * t359 - t358 * t400;
t289 = t306 * t341 - t358 * t404;
t288 = -t306 * t345 - t358 * t405;
t307 = -t346 * t359 + t358 * t401;
t362 = -t339 * t397 - t398;
t356 = -t282 + (-t283 * t409 + t282) * t284;
t304 = (-t361 * qJD(2) + t362 * qJD(3)) * t338;
t303 = -t310 * t337 - t339 * t380;
t290 = -t342 * t378 + (t362 * qJD(2) - t361 * qJD(3)) * t338;
t273 = -t311 * t401 - t310 * t346 + (t324 * t342 - t325 * t400) * qJD(3);
t272 = t307 * qJD(3) + t308 * t342 - t309 * t400;
t264 = 0.1e1 / t266;
t255 = 0.1e1 / t257;
t254 = t366 * t284;
t253 = t367 * t284;
t250 = t356 * t301;
t248 = t372 * t254 - t282 * t305 + t283 * t321;
t247 = t372 * t253 + t282 * t295 + t283 * t315;
t245 = t366 * t385 + (t321 * t375 - t273 * t312 + (t270 * t321 + t291 * t305 - t299 * t304) * t313) * t284;
t244 = t367 * t385 + (t315 * t375 - t271 * t312 + (t270 * t315 - t290 * t299 - t291 * t295) * t313) * t284;
t1 = [0.2e1 * t301 * t384 + (-t268 * t312 + t301 * t422) * t284, t245, t244, 0, 0, 0; -0.2e1 * t299 * t260 * t420 + ((t329 + (t339 * t390 - t311) * t346 + (qJD(3) * t325 + t310 * t339 - t374) * t342) * t260 + (-t299 * t246 - t250 * t268) * t261) * t255 + ((t250 * t387 - t356 * t414) * t255 + (t250 * t388 + (-(t251 * t284 * t409 + t385) * t412 - (t384 * t421 - t251 + (t251 - t369) * t284) * t411) * t255) * t261) * t301 (t248 * t415 - t260 * t307) * t388 + ((-t306 * qJD(3) + t308 * t346 + t309 * t401) * t260 + t248 * t377 + (-t307 * t246 - t248 * t268 - (t245 * t299 - t254 * t270 + t304 + (-t254 * t316 - t305) * t251) * t411 - (-t245 * t316 - t254 * t291 - t273 + (-t254 * t299 - t321) * t251) * t412) * t261) * t255 (t247 * t415 + t260 * t300) * t388 + (t247 * t377 - t267 * t260 + (t300 * t246 - t247 * t268 - (t244 * t299 - t253 * t270 + t290 + (-t253 * t316 + t295) * t251) * t411 - (-t244 * t316 - t253 * t291 - t271 + (-t253 * t299 - t315) * t251) * t412) * t261) * t255, 0, 0, 0; (t275 * t371 - t279 * t413) * t386 + ((t279 * qJD(5) - t271 * t345 + t303 * t341) * t275 + t279 * t376 + (t371 * t259 + (t371 * qJD(5) + t271 * t341 + t303 * t345) * t370 - t279 * t258) * t276) * t264 (-t275 * t288 - t289 * t413) * t386 + ((t289 * qJD(5) - t272 * t345 - t309 * t405) * t275 + t289 * t376 + (-t288 * t259 + (-t288 * qJD(5) + t272 * t341 - t309 * t404) * t370 - t289 * t258) * t276) * t264, t368 * t301 * t386 + (-t368 * t268 + ((qJD(5) * t275 + t376) * t341 + (-t258 * t341 + (t259 + t423) * t345) * t276) * t301) * t264, 0, -0.2e1 * t418 - 0.2e1 * (t258 * t276 * t264 - (-t264 * t416 - t276 * t418) * t370) * t370, 0;];
JaD_rot  = t1;
