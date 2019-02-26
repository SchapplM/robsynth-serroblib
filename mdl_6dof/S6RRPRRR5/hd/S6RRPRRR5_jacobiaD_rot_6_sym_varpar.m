% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:56:20
% EndTime: 2019-02-26 21:56:22
% DurationCPUTime: 1.97s
% Computational Cost: add. (9451->155), mult. (25311->308), div. (744->12), fcn. (32679->15), ass. (0->136)
t341 = sin(pkin(12));
t342 = cos(pkin(12));
t345 = sin(qJ(2));
t347 = cos(qJ(2));
t327 = t341 * t345 - t347 * t342;
t343 = cos(pkin(6));
t324 = t327 * t343;
t319 = qJD(2) * t324;
t369 = t341 * t347 + t342 * t345;
t326 = t369 * qJD(2);
t348 = cos(qJ(1));
t416 = sin(pkin(6));
t382 = t348 * t416;
t325 = t369 * t343;
t417 = sin(qJ(1));
t385 = t417 * t325;
t422 = -t417 * t326 - qJD(1) * t385 + (-qJD(1) * t327 - t319) * t348 - qJD(4) * t382;
t344 = sin(qJ(4));
t346 = cos(qJ(4));
t362 = -t348 * t325 + t417 * t327;
t297 = -t344 * t362 + t346 * t382;
t383 = t347 * t416;
t384 = t345 * t416;
t356 = -t341 * t383 - t342 * t384;
t313 = -t343 * t346 - t344 * t356;
t292 = atan2(-t297, t313);
t287 = sin(t292);
t288 = cos(t292);
t263 = -t287 * t297 + t288 * t313;
t261 = 0.1e1 / t263 ^ 2;
t363 = -t348 * t327 - t385;
t374 = t417 * t416;
t359 = -t344 * t363 + t346 * t374;
t296 = t359 ^ 2;
t259 = t261 * t296 + 0.1e1;
t303 = t344 * t374 + t346 * t363;
t355 = t362 * qJD(1) + t417 * t319 - t348 * t326;
t373 = qJD(1) * t382;
t267 = t303 * qJD(4) + t344 * t355 - t346 * t373;
t409 = t261 * t359;
t295 = t297 ^ 2;
t311 = 0.1e1 / t313 ^ 2;
t291 = t295 * t311 + 0.1e1;
t289 = 0.1e1 / t291;
t368 = qJD(1) * t374;
t395 = qJD(4) * t346;
t269 = t422 * t344 - t346 * t368 - t362 * t395;
t314 = t343 * t344 - t346 * t356;
t322 = -t341 * t384 + t342 * t383;
t318 = t322 * qJD(2);
t293 = t314 * qJD(4) + t318 * t344;
t310 = 0.1e1 / t313;
t401 = t297 * t311;
t367 = -t269 * t310 + t293 * t401;
t251 = t367 * t289;
t370 = -t287 * t313 - t288 * t297;
t246 = t370 * t251 - t287 * t269 + t288 * t293;
t260 = 0.1e1 / t263;
t262 = t260 * t261;
t414 = t246 * t262;
t394 = 0.2e1 * (-t267 * t409 - t296 * t414) / t259 ^ 2;
t421 = t293 * t311;
t305 = -t348 * t324 - t369 * t417;
t364 = -t305 * t310 + t322 * t401;
t420 = t344 * t364;
t270 = (qJD(4) * t362 + t368) * t344 + t422 * t346;
t308 = t417 * t324 - t348 * t369;
t340 = qJ(5) + qJ(6);
t337 = sin(t340);
t338 = cos(t340);
t279 = t303 * t338 - t308 * t337;
t273 = 0.1e1 / t279;
t274 = 0.1e1 / t279 ^ 2;
t419 = -0.2e1 * t297;
t418 = -0.2e1 * t359;
t320 = t343 * t326;
t361 = t327 * qJD(2);
t283 = t305 * qJD(1) - t417 * t320 - t348 * t361;
t339 = qJD(5) + qJD(6);
t377 = t303 * t339 - t283;
t268 = t359 * qJD(4) + t344 * t373 + t346 * t355;
t379 = -t308 * t339 + t268;
t253 = t379 * t337 + t377 * t338;
t278 = t303 * t337 + t308 * t338;
t272 = t278 ^ 2;
t266 = t272 * t274 + 0.1e1;
t407 = t274 * t278;
t254 = -t377 * t337 + t379 * t338;
t411 = t254 * t273 * t274;
t413 = (t253 * t407 - t272 * t411) / t266 ^ 2;
t403 = t310 * t421;
t412 = (t269 * t401 - t295 * t403) / t291 ^ 2;
t410 = t261 * t267;
t408 = t273 * t337;
t406 = t278 * t338;
t405 = t287 * t359;
t404 = t288 * t359;
t402 = t297 * t310;
t400 = t308 * t344;
t399 = t308 * t346;
t393 = -0.2e1 * t413;
t392 = 0.2e1 * t413;
t391 = -0.2e1 * t412;
t390 = t262 * t418;
t389 = t310 * t412;
t388 = t278 * t411;
t387 = t261 * t405;
t386 = t261 * t404;
t381 = 0.2e1 * t388;
t380 = t403 * t419;
t378 = t305 * t339 - t270;
t285 = t308 * qJD(1) - t348 * t320 + t417 * t361;
t299 = -t344 * t382 - t346 * t362;
t376 = -t299 * t339 - t285;
t371 = t339 * t399 - t355;
t366 = t274 * t406 - t408;
t365 = -t299 * t310 + t314 * t401;
t360 = -t287 + (t288 * t402 + t287) * t289;
t358 = -qJD(4) * t400 - t283 * t346 + t339 * t363;
t317 = t356 * qJD(2);
t294 = -t313 * qJD(4) + t318 * t346;
t281 = t337 * t363 + t338 * t399;
t280 = t337 * t399 - t338 * t363;
t277 = -t299 * t338 + t305 * t337;
t276 = -t299 * t337 - t305 * t338;
t264 = 0.1e1 / t266;
t257 = 0.1e1 / t259;
t256 = t289 * t420;
t255 = t365 * t289;
t250 = t360 * t359;
t248 = (-t287 * t305 + t288 * t322) * t344 + t370 * t256;
t247 = t370 * t255 - t287 * t299 + t288 * t314;
t245 = t365 * t391 + (t314 * t380 - t270 * t310 + (t269 * t314 + t293 * t299 + t294 * t297) * t311) * t289;
t243 = t391 * t420 + (t364 * t395 + (t322 * t380 - t285 * t310 + (t269 * t322 + t293 * t305 + t297 * t317) * t311) * t344) * t289;
t242 = t393 + 0.2e1 * (t253 * t274 * t264 + (-t264 * t411 - t274 * t413) * t278) * t278;
t1 = [t389 * t418 + (-t267 * t310 - t359 * t421) * t289, t243, 0, t245, 0, 0; t297 * t260 * t394 + (-t269 * t260 + (t246 * t297 + t250 * t267) * t261) * t257 - (-t250 * t261 * t394 + (-0.2e1 * t250 * t414 + (-t251 * t289 * t402 + t391) * t387 + (t389 * t419 - t251 + (t251 - t367) * t289) * t386 - t360 * t410) * t257) * t359 (-t248 * t409 - t260 * t400) * t394 + (-t248 * t410 + (-t283 * t344 + t308 * t395) * t260 + (t248 * t390 - t261 * t400) * t246 + (t322 * t395 - t243 * t297 - t256 * t269 + t317 * t344 + (-t256 * t313 - t305 * t344) * t251) * t386 + (-t305 * t395 - t243 * t313 - t256 * t293 - t285 * t344 + (t256 * t297 - t322 * t344) * t251) * t387) * t257, 0 (-t247 * t409 - t260 * t303) * t394 + (t247 * t246 * t390 + t268 * t260 + (-t303 * t246 - t247 * t267 + (-t245 * t297 - t255 * t269 + t294 + (-t255 * t313 - t299) * t251) * t404 + (-t245 * t313 - t255 * t293 - t270 + (t255 * t297 - t314) * t251) * t405) * t261) * t257, 0, 0; (-t273 * t276 + t277 * t407) * t392 + ((t378 * t337 + t376 * t338) * t273 + t277 * t381 + (-t276 * t254 - (-t376 * t337 + t378 * t338) * t278 - t277 * t253) * t274) * t264 (-t273 * t280 + t281 * t407) * t392 + (t281 * t381 + t371 * t273 * t338 + t358 * t408 + (t371 * t278 * t337 - t281 * t253 - t280 * t254 - t358 * t406) * t274) * t264, 0, -t366 * t359 * t393 + (t366 * t267 - ((-t273 * t339 - 0.2e1 * t388) * t338 + (t253 * t338 + (-t278 * t339 + t254) * t337) * t274) * t359) * t264, t242, t242;];
JaD_rot  = t1;
