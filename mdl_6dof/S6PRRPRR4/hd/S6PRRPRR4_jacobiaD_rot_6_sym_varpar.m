% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRR4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:53
% EndTime: 2019-02-26 20:05:56
% DurationCPUTime: 2.88s
% Computational Cost: add. (10809->132), mult. (30690->255), div. (801->12), fcn. (39581->15), ass. (0->122)
t351 = sin(qJ(3));
t354 = cos(qJ(3));
t355 = cos(qJ(2));
t417 = cos(pkin(11));
t418 = cos(pkin(6));
t386 = t418 * t417;
t416 = sin(pkin(11));
t419 = sin(qJ(2));
t370 = -t416 * t355 - t419 * t386;
t348 = sin(pkin(6));
t393 = t348 * t417;
t329 = -t351 * t393 - t370 * t354;
t350 = sin(qJ(5));
t353 = cos(qJ(5));
t364 = -t370 * t351 + t354 * t393;
t308 = t329 * t353 + t364 * t350;
t367 = qJD(3) * t370;
t340 = t355 * t386 - t416 * t419;
t424 = t340 * qJD(2) - qJD(3) * t393;
t315 = t351 * t367 + t354 * t424;
t362 = t351 * t424 - t354 * t367;
t275 = t308 * qJD(5) + t315 * t350 - t362 * t353;
t442 = -t329 * t350 + t364 * t353;
t276 = t442 * qJD(5) + t315 * t353 + t362 * t350;
t303 = t442 ^ 2;
t396 = t351 * t419;
t343 = t348 * t396 - t418 * t354;
t395 = t354 * t419;
t344 = t348 * t395 + t418 * t351;
t380 = t343 * t353 - t344 * t350;
t322 = 0.1e1 / t380 ^ 2;
t289 = t303 * t322 + 0.1e1;
t287 = 0.1e1 / t289;
t326 = t343 * t350 + t344 * t353;
t403 = t348 * t355;
t394 = qJD(2) * t403;
t332 = t344 * qJD(3) + t351 * t394;
t333 = -t343 * qJD(3) + t354 * t394;
t292 = t326 * qJD(5) - t332 * t353 + t333 * t350;
t293 = t380 * qJD(5) + t332 * t350 + t333 * t353;
t321 = 0.1e1 / t380;
t406 = t442 * t322;
t375 = t308 * t321 - t326 * t406;
t408 = t292 * t321 * t322;
t390 = -0.2e1 * t442 * t408;
t422 = -0.2e1 * (-t275 * t406 + t303 * t408) / t289 ^ 2;
t252 = (t276 * t321 + (t275 * t326 + t292 * t308 - t293 * t442) * t322 + t326 * t390) * t287 + t375 * t422;
t260 = (t275 * t321 - t292 * t406) * t287;
t290 = atan2(t442, -t380);
t285 = sin(t290);
t286 = cos(t290);
t384 = t285 * t380 + t286 * t442;
t255 = t384 * t260 - t285 * t275 + t286 * t292;
t444 = t375 * t287;
t257 = -t285 * t308 + t286 * t326 + t384 * t444;
t271 = t285 * t442 - t286 * t380;
t269 = 0.1e1 / t271 ^ 2;
t385 = t418 * t416;
t368 = -t417 * t355 + t419 * t385;
t392 = t348 * t416;
t331 = t351 * t392 - t354 * t368;
t371 = t351 * t368 + t354 * t392;
t382 = -t331 * t350 - t353 * t371;
t304 = t382 ^ 2;
t267 = t304 * t269 + 0.1e1;
t265 = 0.1e1 / t267;
t268 = 0.1e1 / t271;
t312 = t331 * t353 - t350 * t371;
t369 = t355 * t385 + t417 * t419;
t338 = t369 * qJD(2);
t316 = t331 * qJD(3) - t338 * t351;
t317 = t371 * qJD(3) - t338 * t354;
t278 = t312 * qJD(5) - t316 * t353 + t317 * t350;
t279 = t382 * qJD(5) + t316 * t350 + t317 * t353;
t415 = t255 * t268 * t269;
t391 = -0.2e1 * t382 * t415;
t413 = t269 * t382;
t399 = 0.2e1 * (-t278 * t413 - t304 * t415) / t267 ^ 2;
t409 = t286 * t382;
t410 = t285 * t382;
t456 = t265 * ((t312 * t255 + t257 * t278 - (t252 * t442 - t444 * t275 + t293 + (t380 * t444 - t308) * t260) * t409 - (t252 * t380 - (t442 * t444 + t326) * t260 - t444 * t292 - t276) * t410) * t269 - t257 * t391 - t279 * t268) + (t257 * t413 + t268 * t312) * t399;
t349 = sin(qJ(6));
t352 = cos(qJ(6));
t299 = t312 * t352 - t349 * t369;
t339 = t368 * qJD(2);
t272 = t299 * qJD(6) + t279 * t349 - t339 * t352;
t298 = t312 * t349 + t352 * t369;
t400 = qJD(6) * t298;
t273 = t279 * t352 + t339 * t349 - t400;
t294 = t298 ^ 2;
t296 = 0.1e1 / t299 ^ 2;
t282 = t294 * t296 + 0.1e1;
t280 = 0.1e1 / t282;
t295 = 0.1e1 / t299;
t407 = t296 * t298;
t377 = -t349 * t295 + t352 * t407;
t411 = t273 * t295 * t296;
t437 = 0.2e1 * t298;
t389 = t411 * t437;
t443 = (t377 * t278 + (t296 * ((-t273 + t400) * t349 - t272 * t352) + (qJD(6) * t295 + t389) * t352) * t382) * t280;
t430 = t377 * t382;
t379 = t350 * t354 - t351 * t353;
t425 = t379 * t369;
t423 = qJD(5) - qJD(3);
t414 = (t272 * t407 - t294 * t411) / t282 ^ 2;
t398 = -0.2e1 * t414;
t397 = 0.2e1 * t414;
t378 = t350 * t351 + t353 * t354;
t320 = t378 * t369;
t302 = -t320 * t352 + t349 * t368;
t301 = -t320 * t349 - t352 * t368;
t318 = t379 * t340;
t334 = t379 * t403;
t374 = t318 * t321 - t334 * t406;
t365 = t423 * t378;
t300 = ((-t350 * t395 + t353 * t396) * qJD(2) + t365 * t355) * t348;
t284 = t378 * t339 + t423 * t425;
t283 = t379 * t370 * qJD(2) + t365 * t340;
t264 = t374 * t287;
t258 = t384 * t264 - t285 * t318 + t286 * t334;
t254 = t374 * t422 + (t334 * t390 + t283 * t321 + (t275 * t334 + t292 * t318 - t300 * t442) * t322) * t287;
t1 = [0, t254, -t252, 0, t252, 0; 0 (-t258 * t413 + t268 * t425) * t399 + (t258 * t391 + (t425 * t255 - t258 * t278 + (t254 * t442 - t264 * t275 + t300 + (t264 * t380 - t318) * t260) * t409 + (t254 * t380 - t264 * t292 - t283 + (-t264 * t442 - t334) * t260) * t410) * t269 + (t379 * t339 - t365 * t369) * t268) * t265, t456, 0, -t456, 0; 0 (-t295 * t301 + t302 * t407) * t397 + ((t302 * qJD(6) + t284 * t349 - t338 * t352) * t295 + t302 * t389 + (-t301 * t273 - (-t301 * qJD(6) + t284 * t352 + t338 * t349) * t298 - t302 * t272) * t296) * t280, -t397 * t430 - t443, 0, -t398 * t430 + t443, t398 + (t272 * t296 * t280 + (-t280 * t411 - t296 * t414) * t298) * t437;];
JaD_rot  = t1;
