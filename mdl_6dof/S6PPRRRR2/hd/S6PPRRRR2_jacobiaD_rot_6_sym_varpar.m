% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRR2
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPRRRR2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_jacobiaD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:43:16
% EndTime: 2019-02-26 19:43:18
% DurationCPUTime: 1.77s
% Computational Cost: add. (10307->121), mult. (29136->245), div. (577->12), fcn. (38233->17), ass. (0->123)
t401 = sin(pkin(13));
t402 = sin(pkin(12));
t369 = t402 * t401;
t405 = cos(pkin(13));
t406 = cos(pkin(12));
t375 = t406 * t405;
t408 = cos(pkin(6));
t338 = -t408 * t369 + t375;
t345 = sin(qJ(3));
t347 = cos(qJ(3));
t371 = t402 * t405;
t374 = t406 * t401;
t361 = t408 * t371 + t374;
t404 = sin(pkin(6));
t370 = t402 * t404;
t403 = sin(pkin(7));
t407 = cos(pkin(7));
t410 = t361 * t407 - t403 * t370;
t322 = t338 * t345 + t410 * t347;
t337 = t408 * t374 + t371;
t360 = -t408 * t375 + t369;
t373 = t404 * t403;
t355 = -t360 * t407 - t406 * t373;
t321 = t337 * t347 + t355 * t345;
t344 = sin(qJ(4));
t346 = cos(qJ(4));
t376 = t407 * t404;
t354 = t360 * t403 - t406 * t376;
t312 = t321 * t346 + t354 * t344;
t320 = -t337 * t345 + t355 * t347;
t316 = t320 * qJD(3);
t294 = t312 * qJD(4) + t316 * t344;
t310 = t321 * t344 - t354 * t346;
t308 = t310 ^ 2;
t359 = t405 * t376 + t408 * t403;
t372 = t404 * t401;
t334 = t359 * t345 + t347 * t372;
t336 = -t405 * t373 + t408 * t407;
t327 = t334 * t344 - t336 * t346;
t325 = 0.1e1 / t327 ^ 2;
t302 = t308 * t325 + 0.1e1;
t300 = 0.1e1 / t302;
t328 = t334 * t346 + t336 * t344;
t333 = -t345 * t372 + t359 * t347;
t329 = t333 * qJD(3);
t306 = t328 * qJD(4) + t329 * t344;
t324 = 0.1e1 / t327;
t390 = t310 * t325;
t272 = (-t294 * t324 + t306 * t390) * t300;
t303 = atan2(-t310, t327);
t298 = sin(t303);
t299 = cos(t303);
t368 = -t298 * t327 - t299 * t310;
t268 = t368 * t272 - t294 * t298 + t299 * t306;
t282 = -t298 * t310 + t299 * t327;
t279 = 0.1e1 / t282;
t280 = 0.1e1 / t282 ^ 2;
t413 = t268 * t279 * t280;
t323 = t338 * t347 - t345 * t410;
t356 = t361 * t403 + t407 * t370;
t313 = t323 * t344 - t356 * t346;
t412 = 0.2e1 * t313 * t413;
t364 = -t320 * t324 + t333 * t390;
t411 = t344 * t364;
t391 = t306 * t324 * t325;
t409 = -0.2e1 * (t294 * t390 - t308 * t391) / t302 ^ 2;
t314 = t323 * t346 + t356 * t344;
t343 = qJ(5) + qJ(6);
t340 = sin(t343);
t341 = cos(t343);
t293 = t314 * t341 + t322 * t340;
t289 = 0.1e1 / t293;
t290 = 0.1e1 / t293 ^ 2;
t319 = t323 * qJD(3);
t342 = qJD(5) + qJD(6);
t378 = t314 * t342 - t319;
t318 = t322 * qJD(3);
t297 = -t313 * qJD(4) - t318 * t346;
t379 = t322 * t342 + t297;
t283 = t379 * t340 + t378 * t341;
t292 = t314 * t340 - t322 * t341;
t288 = t292 ^ 2;
t287 = t288 * t290 + 0.1e1;
t396 = t290 * t292;
t284 = -t378 * t340 + t379 * t341;
t398 = t284 * t289 * t290;
t400 = (t283 * t396 - t288 * t398) / t287 ^ 2;
t399 = t280 * t313;
t397 = t289 * t340;
t395 = t292 * t341;
t296 = t314 * qJD(4) - t318 * t344;
t394 = t296 * t280;
t393 = t298 * t313;
t392 = t299 * t313;
t389 = t322 * t344;
t388 = t322 * t346;
t385 = qJD(4) * t346;
t309 = t313 ^ 2;
t278 = t280 * t309 + 0.1e1;
t384 = 0.2e1 * (-t309 * t413 + t313 * t394) / t278 ^ 2;
t383 = -0.2e1 * t400;
t381 = t292 * t398;
t380 = -0.2e1 * t310 * t391;
t377 = t342 * t388 - t318;
t366 = t290 * t395 - t397;
t365 = -t312 * t324 + t328 * t390;
t362 = qJD(4) * t389 - t319 * t346 + t323 * t342;
t330 = t334 * qJD(3);
t317 = t321 * qJD(3);
t307 = -t327 * qJD(4) + t329 * t346;
t305 = t323 * t340 - t341 * t388;
t304 = -t323 * t341 - t340 * t388;
t295 = -t310 * qJD(4) + t316 * t346;
t285 = 0.1e1 / t287;
t276 = 0.1e1 / t278;
t274 = t300 * t411;
t273 = t365 * t300;
t270 = (-t298 * t320 + t299 * t333) * t344 + t368 * t274;
t269 = t368 * t273 - t298 * t312 + t299 * t328;
t266 = t365 * t409 + (t328 * t380 - t295 * t324 + (t294 * t328 + t306 * t312 + t307 * t310) * t325) * t300;
t265 = t409 * t411 + (t364 * t385 + (t333 * t380 + t317 * t324 + (t294 * t333 + t306 * t320 - t310 * t330) * t325) * t344) * t300;
t264 = t383 + 0.2e1 * (t283 * t285 * t290 + (-t285 * t398 - t290 * t400) * t292) * t292;
t1 = [0, 0, t265, t266, 0, 0; 0, 0 (t270 * t399 + t279 * t389) * t384 + ((-t319 * t344 - t322 * t385) * t279 + (-t394 + t412) * t270 + (t389 * t268 - (t333 * t385 - t265 * t310 - t274 * t294 - t330 * t344 + (-t274 * t327 - t320 * t344) * t272) * t392 - (-t320 * t385 - t265 * t327 - t274 * t306 + t317 * t344 + (t274 * t310 - t333 * t344) * t272) * t393) * t280) * t276 (t269 * t399 - t279 * t314) * t384 + (t269 * t412 + t297 * t279 + (-t314 * t268 - t269 * t296 - (-t266 * t310 - t273 * t294 + t307 + (-t273 * t327 - t312) * t272) * t392 - (-t266 * t327 - t273 * t306 - t295 + (t273 * t310 - t328) * t272) * t393) * t280) * t276, 0, 0; 0, 0, 0.2e1 * (-t289 * t304 + t305 * t396) * t400 + (0.2e1 * t305 * t381 - t377 * t289 * t341 + t362 * t397 + (-t377 * t292 * t340 - t305 * t283 - t304 * t284 - t362 * t395) * t290) * t285, t366 * t313 * t383 + (t366 * t296 + ((-t289 * t342 - 0.2e1 * t381) * t341 + (t283 * t341 + (-t292 * t342 + t284) * t340) * t290) * t313) * t285, t264, t264;];
JaD_rot  = t1;
