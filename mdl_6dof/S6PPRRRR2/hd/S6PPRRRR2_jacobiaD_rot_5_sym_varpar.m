% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function JaD_rot = S6PPRRRR2_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR2_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_jacobiaD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:43:16
% EndTime: 2019-02-26 19:43:18
% DurationCPUTime: 1.68s
% Computational Cost: add. (9293->119), mult. (27507->244), div. (559->12), fcn. (36133->17), ass. (0->119)
t368 = sin(pkin(13));
t369 = sin(pkin(12));
t337 = t369 * t368;
t372 = cos(pkin(13));
t373 = cos(pkin(12));
t343 = t373 * t372;
t375 = cos(pkin(6));
t308 = -t337 * t375 + t343;
t312 = sin(qJ(3));
t315 = cos(qJ(3));
t339 = t369 * t372;
t342 = t373 * t368;
t329 = t339 * t375 + t342;
t371 = sin(pkin(6));
t338 = t369 * t371;
t370 = sin(pkin(7));
t374 = cos(pkin(7));
t377 = t329 * t374 - t370 * t338;
t292 = t308 * t312 + t377 * t315;
t307 = t342 * t375 + t339;
t328 = -t343 * t375 + t337;
t341 = t371 * t370;
t323 = -t328 * t374 - t341 * t373;
t291 = t307 * t315 + t312 * t323;
t311 = sin(qJ(4));
t314 = cos(qJ(4));
t344 = t374 * t371;
t322 = t328 * t370 - t344 * t373;
t282 = t291 * t314 + t311 * t322;
t290 = -t307 * t312 + t315 * t323;
t285 = t290 * qJD(3);
t258 = qJD(4) * t282 + t285 * t311;
t280 = t291 * t311 - t314 * t322;
t278 = t280 ^ 2;
t327 = t344 * t372 + t370 * t375;
t340 = t371 * t368;
t304 = t312 * t327 + t315 * t340;
t306 = -t341 * t372 + t374 * t375;
t297 = t304 * t311 - t306 * t314;
t295 = 0.1e1 / t297 ^ 2;
t272 = t278 * t295 + 0.1e1;
t270 = 0.1e1 / t272;
t298 = t304 * t314 + t306 * t311;
t303 = -t312 * t340 + t315 * t327;
t299 = t303 * qJD(3);
t276 = qJD(4) * t298 + t299 * t311;
t294 = 0.1e1 / t297;
t357 = t280 * t295;
t242 = (-t258 * t294 + t276 * t357) * t270;
t273 = atan2(-t280, t297);
t268 = sin(t273);
t269 = cos(t273);
t336 = -t268 * t297 - t269 * t280;
t238 = t242 * t336 - t258 * t268 + t269 * t276;
t252 = -t268 * t280 + t269 * t297;
t249 = 0.1e1 / t252;
t250 = 0.1e1 / t252 ^ 2;
t380 = t238 * t249 * t250;
t293 = t308 * t315 - t312 * t377;
t324 = t329 * t370 + t338 * t374;
t283 = t293 * t311 - t314 * t324;
t379 = 0.2e1 * t283 * t380;
t332 = -t290 * t294 + t303 * t357;
t378 = t311 * t332;
t358 = t276 * t294 * t295;
t376 = -0.2e1 * (t258 * t357 - t278 * t358) / t272 ^ 2;
t284 = t293 * t314 + t311 * t324;
t310 = sin(qJ(5));
t313 = cos(qJ(5));
t267 = t284 * t313 + t292 * t310;
t263 = 0.1e1 / t267;
t264 = 0.1e1 / t267 ^ 2;
t287 = t292 * qJD(3);
t261 = -qJD(4) * t283 - t287 * t314;
t288 = t293 * qJD(3);
t253 = qJD(5) * t267 + t261 * t310 - t288 * t313;
t266 = t284 * t310 - t292 * t313;
t262 = t266 ^ 2;
t257 = t262 * t264 + 0.1e1;
t362 = t264 * t266;
t351 = qJD(5) * t266;
t254 = t261 * t313 + t288 * t310 - t351;
t365 = t254 * t263 * t264;
t367 = (t253 * t362 - t262 * t365) / t257 ^ 2;
t366 = t250 * t283;
t260 = qJD(4) * t284 - t287 * t311;
t364 = t260 * t250;
t363 = t263 * t310;
t361 = t266 * t313;
t360 = t268 * t283;
t359 = t269 * t283;
t356 = t292 * t311;
t355 = t292 * t314;
t352 = qJD(4) * t314;
t279 = t283 ^ 2;
t248 = t250 * t279 + 0.1e1;
t350 = 0.2e1 * (-t279 * t380 + t283 * t364) / t248 ^ 2;
t349 = -0.2e1 * t367;
t347 = t266 * t365;
t346 = -0.2e1 * t280 * t358;
t345 = qJD(5) * t355 - t287;
t334 = t264 * t361 - t363;
t333 = -t282 * t294 + t298 * t357;
t330 = qJD(4) * t356 + qJD(5) * t293 - t288 * t314;
t300 = t304 * qJD(3);
t286 = t291 * qJD(3);
t277 = -qJD(4) * t297 + t299 * t314;
t275 = t293 * t310 - t313 * t355;
t274 = -t293 * t313 - t310 * t355;
t259 = -qJD(4) * t280 + t285 * t314;
t255 = 0.1e1 / t257;
t246 = 0.1e1 / t248;
t244 = t270 * t378;
t243 = t333 * t270;
t240 = (-t268 * t290 + t269 * t303) * t311 + t336 * t244;
t239 = t243 * t336 - t268 * t282 + t269 * t298;
t236 = t333 * t376 + (t298 * t346 - t259 * t294 + (t258 * t298 + t276 * t282 + t277 * t280) * t295) * t270;
t235 = t376 * t378 + (t332 * t352 + (t303 * t346 + t286 * t294 + (t258 * t303 + t276 * t290 - t280 * t300) * t295) * t311) * t270;
t1 = [0, 0, t235, t236, 0, 0; 0, 0 (t240 * t366 + t249 * t356) * t350 + ((-t288 * t311 - t292 * t352) * t249 + (-t364 + t379) * t240 + (t356 * t238 - (t303 * t352 - t235 * t280 - t244 * t258 - t300 * t311 + (-t244 * t297 - t290 * t311) * t242) * t359 - (-t290 * t352 - t235 * t297 - t244 * t276 + t286 * t311 + (t244 * t280 - t303 * t311) * t242) * t360) * t250) * t246 (t239 * t366 - t249 * t284) * t350 + (t239 * t379 + t261 * t249 + (-t284 * t238 - t239 * t260 - (-t236 * t280 - t243 * t258 + t277 + (-t243 * t297 - t282) * t242) * t359 - (-t236 * t297 - t243 * t276 - t259 + (t243 * t280 - t298) * t242) * t360) * t250) * t246, 0, 0; 0, 0, 0.2e1 * (-t263 * t274 + t275 * t362) * t367 + (0.2e1 * t275 * t347 - t345 * t263 * t313 + t330 * t363 + (-t266 * t310 * t345 - t275 * t253 - t274 * t254 - t330 * t361) * t264) * t255, t334 * t283 * t349 + (t334 * t260 + ((-qJD(5) * t263 - 0.2e1 * t347) * t313 + (t253 * t313 + (t254 - t351) * t310) * t264) * t283) * t255, t349 + 0.2e1 * (t253 * t264 * t255 + (-t255 * t365 - t264 * t367) * t266) * t266, 0;];
JaD_rot  = t1;
