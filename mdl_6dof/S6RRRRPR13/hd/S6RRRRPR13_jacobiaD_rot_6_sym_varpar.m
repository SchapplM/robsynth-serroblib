% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR13_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR13_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:37:36
% EndTime: 2019-02-26 22:37:39
% DurationCPUTime: 2.11s
% Computational Cost: add. (6795->177), mult. (19689->350), div. (745->12), fcn. (24766->15), ass. (0->145)
t427 = qJD(4) - qJD(6);
t333 = cos(pkin(6));
t337 = sin(qJ(2));
t420 = sin(qJ(1));
t386 = t420 * t337;
t379 = t333 * t386;
t340 = cos(qJ(2));
t341 = cos(qJ(1));
t400 = t341 * t340;
t303 = -qJD(1) * t379 - qJD(2) * t386 + (qJD(2) * t333 + qJD(1)) * t400;
t336 = sin(qJ(3));
t419 = sin(pkin(6));
t331 = t341 * t419 * t336;
t385 = t420 * t340;
t401 = t341 * t337;
t358 = -t333 * t401 - t385;
t375 = t420 * t419;
t421 = cos(qJ(3));
t359 = t421 * t375;
t384 = qJD(3) * t421;
t281 = -qJD(1) * t359 - qJD(3) * t331 + t303 * t336 - t358 * t384;
t376 = t421 * t419;
t364 = t341 * t376;
t352 = -t336 * t358 + t364;
t304 = t352 ^ 2;
t383 = t337 * t419;
t351 = -t333 * t421 + t336 * t383;
t317 = 0.1e1 / t351 ^ 2;
t295 = t304 * t317 + 0.1e1;
t405 = t352 * t317;
t316 = 0.1e1 / t351;
t320 = -t333 * t336 - t337 * t376;
t382 = t340 * t419;
t377 = t336 * t382;
t306 = -qJD(2) * t377 + t320 * qJD(3);
t425 = t306 * t317;
t407 = t316 * t425;
t414 = (t281 * t405 + t304 * t407) / t295 ^ 2;
t426 = -0.2e1 * t414;
t378 = t352 * t382;
t357 = -t333 * t400 + t386;
t403 = t316 * t357;
t424 = t336 * (t317 * t378 + t403);
t339 = cos(qJ(4));
t423 = t427 * t339;
t335 = sin(qJ(4));
t422 = t427 * t335;
t363 = qJD(3) * t376;
t282 = (qJD(1) * t375 + qJD(3) * t358) * t336 + t303 * t421 - t341 * t363;
t355 = -t333 * t385 - t401;
t301 = t358 * qJD(1) + t355 * qJD(2);
t356 = t379 - t400;
t349 = t336 * t356 + t359;
t280 = qJD(1) * t331 + t349 * qJD(3) + t301 * t421;
t314 = t336 * t375 - t356 * t421;
t290 = t314 * t339 - t335 * t355;
t300 = t357 * qJD(1) + t356 * qJD(2);
t254 = t290 * qJD(4) + t280 * t335 + t300 * t339;
t289 = t314 * t335 + t339 * t355;
t255 = -t289 * qJD(4) + t280 * t339 - t300 * t335;
t334 = sin(qJ(6));
t338 = cos(qJ(6));
t370 = t289 * t338 - t290 * t334;
t245 = t370 * qJD(6) + t254 * t334 + t255 * t338;
t296 = atan2(t352, -t351);
t291 = sin(t296);
t292 = cos(t296);
t266 = t291 * t352 - t292 * t351;
t263 = 0.1e1 / t266;
t276 = t289 * t334 + t290 * t338;
t268 = 0.1e1 / t276;
t264 = 0.1e1 / t266 ^ 2;
t269 = 0.1e1 / t276 ^ 2;
t244 = t276 * qJD(6) - t254 * t338 + t255 * t334;
t267 = t370 ^ 2;
t250 = t267 * t269 + 0.1e1;
t412 = t269 * t370;
t270 = t268 * t269;
t415 = t245 * t270;
t418 = (-t244 * t412 - t267 * t415) / t250 ^ 2;
t305 = t349 ^ 2;
t260 = t264 * t305 + 0.1e1;
t279 = -qJD(1) * t364 + t314 * qJD(3) + t301 * t336;
t410 = t279 * t264;
t293 = 0.1e1 / t295;
t361 = -t281 * t316 - t306 * t405;
t247 = t361 * t293;
t369 = t291 * t351 + t292 * t352;
t241 = t369 * t247 + t281 * t291 + t292 * t306;
t265 = t263 * t264;
t416 = t241 * t265;
t417 = (-t305 * t416 - t349 * t410) / t260 ^ 2;
t413 = t264 * t349;
t366 = -t334 * t335 - t338 * t339;
t285 = t366 * t349;
t411 = t269 * t285;
t409 = t291 * t349;
t408 = t292 * t349;
t406 = t352 * t316;
t404 = t352 * t320;
t402 = t355 * t336;
t395 = 0.2e1 * t418;
t394 = -0.2e1 * t417;
t393 = 0.2e1 * t417;
t392 = 0.2e1 * t414;
t391 = 0.2e1 * t265 * t349;
t390 = -0.2e1 * t270 * t370;
t389 = t264 * t409;
t388 = t264 * t408;
t387 = t355 * t421;
t381 = t245 * t390;
t380 = t316 * t426;
t374 = qJD(2) * t383;
t309 = -t358 * t421 - t331;
t287 = -t309 * t335 + t339 * t357;
t288 = -t309 * t339 - t335 * t357;
t371 = t287 * t338 - t288 * t334;
t272 = t287 * t334 + t288 * t338;
t297 = t335 * t387 + t339 * t356;
t298 = -t335 * t356 + t339 * t387;
t368 = t297 * t338 - t298 * t334;
t278 = t297 * t334 + t298 * t338;
t367 = -t334 * t339 + t335 * t338;
t365 = t340 * t376;
t362 = -qJD(4) * t387 + t301;
t360 = t309 * t316 + t317 * t404;
t354 = t291 + (-t292 * t406 - t291) * t293;
t348 = -qJD(3) * t402 - qJD(4) * t356 + t421 * t300;
t307 = -qJD(2) * t365 + t351 * qJD(3);
t302 = t355 * qJD(1) + t358 * qJD(2);
t284 = t367 * t349;
t262 = t362 * t335 + t348 * t339;
t261 = t348 * t335 - t362 * t339;
t258 = 0.1e1 / t260;
t257 = -t287 * qJD(4) - t282 * t339 + t302 * t335;
t256 = t288 * qJD(4) - t282 * t335 - t302 * t339;
t253 = t293 * t424;
t252 = t360 * t293;
t248 = 0.1e1 / t250;
t246 = t354 * t349;
t243 = (-t291 * t357 - t292 * t382) * t336 + t369 * t253;
t242 = -t369 * t252 + t291 * t309 + t292 * t320;
t240 = t360 * t392 + (-0.2e1 * t404 * t407 - t282 * t316 + (-t281 * t320 - t306 * t309 - t307 * t352) * t317) * t293;
t238 = t424 * t426 + ((t365 * t405 + t421 * t403) * qJD(3) + (0.2e1 * t378 * t407 - t302 * t316 + (t281 * t382 + t306 * t357 - t352 * t374) * t317) * t336) * t293;
t1 = [t349 * t380 + (-t279 * t316 + t349 * t425) * t293, t238, t240, 0, 0, 0; t352 * t263 * t394 + (t281 * t263 + (-t241 * t352 - t246 * t279) * t264) * t258 - (-t246 * t264 * t394 + (0.2e1 * t246 * t416 - (t247 * t293 * t406 + t392) * t389 - (-t352 * t380 + t247 + (-t247 + t361) * t293) * t388 + t354 * t410) * t258) * t349 (t243 * t413 + t263 * t402) * t393 + (t243 * t410 + (-t300 * t336 - t355 * t384) * t263 + (t243 * t391 + t264 * t402) * t241 - (-t340 * t363 + t336 * t374 + t238 * t352 + t253 * t281 + (t253 * t351 - t336 * t357) * t247) * t388 - (-t357 * t384 + t238 * t351 - t253 * t306 + t302 * t336 + (-t253 * t352 + t377) * t247) * t389) * t258 (t242 * t413 + t263 * t314) * t393 + (t242 * t241 * t391 - t280 * t263 + (t314 * t241 + t242 * t279 - (t240 * t352 - t252 * t281 + t307 + (-t252 * t351 + t309) * t247) * t408 - (t240 * t351 + t252 * t306 + t282 + (t252 * t352 - t320) * t247) * t409) * t264) * t258, 0, 0, 0; (t268 * t371 - t272 * t412) * t395 + ((t272 * qJD(6) - t256 * t338 + t257 * t334) * t268 + t272 * t381 + (t371 * t245 + (t371 * qJD(6) + t256 * t334 + t257 * t338) * t370 - t272 * t244) * t269) * t248 (t268 * t368 - t278 * t412) * t395 + ((t278 * qJD(6) - t261 * t338 + t262 * t334) * t268 + t278 * t381 + (t368 * t245 + (t368 * qJD(6) + t261 * t334 + t262 * t338) * t370 - t278 * t244) * t269) * t248 (t268 * t284 + t370 * t411) * t395 + (t244 * t411 + (t284 * t269 - t285 * t390) * t245 + (t367 * t268 + t366 * t412) * t279 - ((t423 * t268 + t422 * t412) * t338 + (t422 * t268 - t423 * t412) * t334) * t349) * t248 (t268 * t276 + t370 * t412) * t395 + (-t245 * t268 - t370 * t381 + (0.2e1 * t370 * t244 + t276 * t245) * t269) * t248, 0, -0.2e1 * t418 - 0.2e1 * (t244 * t269 * t248 - (-t248 * t415 - t269 * t418) * t370) * t370;];
JaD_rot  = t1;
