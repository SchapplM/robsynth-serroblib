% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:11
% EndTime: 2019-02-26 21:30:13
% DurationCPUTime: 2.08s
% Computational Cost: add. (8421->151), mult. (24081->304), div. (726->12), fcn. (31139->15), ass. (0->131)
t317 = sin(qJ(1));
t320 = cos(qJ(1));
t397 = sin(pkin(11));
t399 = cos(pkin(6));
t352 = t399 * t397;
t398 = cos(pkin(11));
t353 = t399 * t398;
t400 = sin(qJ(2));
t401 = cos(qJ(2));
t334 = t400 * t352 - t401 * t353;
t338 = t401 * t397 + t400 * t398;
t287 = -t317 * t338 - t320 * t334;
t316 = sin(qJ(5));
t319 = cos(qJ(5));
t314 = sin(pkin(6));
t380 = t314 * t320;
t343 = -t287 * t319 + t316 * t380;
t276 = t343 ^ 2;
t308 = t400 * t397 - t401 * t398;
t302 = t308 * t314;
t294 = -t302 * t319 + t399 * t316;
t292 = 0.1e1 / t294 ^ 2;
t271 = t276 * t292 + 0.1e1;
t269 = 0.1e1 / t271;
t336 = qJD(2) * t308;
t407 = qJD(1) * t334 + t336;
t304 = t401 * t352 + t400 * t353;
t408 = -qJD(1) * t338 - qJD(2) * t304;
t328 = t407 * t317 + t408 * t320;
t377 = qJD(1) * t314;
t363 = t317 * t377;
t364 = t319 * t380;
t246 = -qJD(5) * t364 + (-t287 * qJD(5) + t363) * t316 + t328 * t319;
t295 = t302 * t316 + t399 * t319;
t303 = t338 * t314;
t297 = qJD(2) * t303;
t274 = t295 * qJD(5) - t297 * t319;
t291 = 0.1e1 / t294;
t383 = t343 * t292;
t348 = -t246 * t291 - t274 * t383;
t230 = t348 * t269;
t272 = atan2(t343, t294);
t267 = sin(t272);
t268 = cos(t272);
t351 = -t267 * t294 + t268 * t343;
t225 = t351 * t230 - t267 * t246 + t268 * t274;
t242 = t267 * t343 + t268 * t294;
t240 = 0.1e1 / t242 ^ 2;
t411 = t240 * t225;
t331 = t317 * t334 - t320 * t338;
t381 = t314 * t317;
t278 = -t331 * t316 + t319 * t381;
t315 = sin(qJ(6));
t318 = cos(qJ(6));
t349 = t317 * t304 + t308 * t320;
t255 = t278 * t315 + t318 * t349;
t410 = 0.2e1 * t255;
t239 = 0.1e1 / t242;
t409 = t239 * t411;
t277 = t316 * t381 + t331 * t319;
t359 = 0.2e1 * t277 * t409;
t327 = t408 * t317 - t407 * t320;
t362 = t320 * t377;
t249 = t278 * qJD(5) + t316 * t362 - t327 * t319;
t389 = t249 * t240;
t406 = -t389 + t359;
t405 = t274 * t292;
t286 = t304 * t320 - t317 * t308;
t345 = t286 * t291 + t303 * t383;
t404 = t319 * t345;
t403 = -qJD(6) * t316 * t349 + t327;
t382 = t349 * t319;
t402 = -qJD(5) * t382 + t331 * qJD(6);
t299 = t334 * qJD(2);
t307 = t338 * qJD(2);
t266 = t349 * qJD(1) + t299 * t320 + t317 * t307;
t247 = t343 * qJD(5) - t328 * t316 + t319 * t363;
t256 = t278 * t318 - t315 * t349;
t252 = 0.1e1 / t256;
t253 = 0.1e1 / t256 ^ 2;
t275 = t277 ^ 2;
t238 = t240 * t275 + 0.1e1;
t396 = (-t275 * t409 + t277 * t389) / t238 ^ 2;
t250 = -qJD(5) * t277 + t327 * t316 + t319 * t362;
t335 = -t286 * qJD(1) + t317 * t299 - t307 * t320;
t233 = t256 * qJD(6) + t250 * t315 - t318 * t335;
t251 = t255 ^ 2;
t245 = t251 * t253 + 0.1e1;
t388 = t253 * t255;
t375 = qJD(6) * t255;
t234 = t250 * t318 + t315 * t335 - t375;
t392 = t234 * t252 * t253;
t395 = (t233 * t388 - t251 * t392) / t245 ^ 2;
t385 = t291 * t405;
t393 = (-t246 * t383 - t276 * t385) / t271 ^ 2;
t391 = t240 * t277;
t243 = 0.1e1 / t245;
t390 = t243 * t253;
t387 = t267 * t277;
t386 = t268 * t277;
t384 = t343 * t291;
t379 = t316 * t315;
t378 = t316 * t318;
t376 = qJD(5) * t316;
t374 = 0.2e1 * t396;
t373 = -0.2e1 * t395;
t372 = -0.2e1 * t393;
t371 = 0.2e1 * t393;
t369 = t253 * t395;
t368 = t233 * t390;
t367 = t255 * t392;
t366 = t343 * t385;
t358 = t291 * t371;
t357 = 0.2e1 * t367;
t344 = t287 * t316 + t364;
t258 = -t286 * t315 + t318 * t344;
t257 = t286 * t318 + t315 * t344;
t347 = -t315 * t252 + t318 * t388;
t346 = -t291 * t344 + t295 * t383;
t339 = -t267 + (-t268 * t384 + t267) * t269;
t298 = t314 * t336;
t273 = -t294 * qJD(5) + t297 * t316;
t262 = t331 * t315 - t349 * t378;
t236 = 0.1e1 / t238;
t235 = t269 * t404;
t232 = t346 * t269;
t227 = (t267 * t286 - t268 * t303) * t319 + t351 * t235;
t226 = -t351 * t232 + t267 * t344 + t268 * t295;
t224 = t346 * t371 + (0.2e1 * t295 * t366 - t247 * t291 + (t246 * t295 - t273 * t343 - t274 * t344) * t292) * t269;
t222 = t372 * t404 + (-t345 * t376 + (-0.2e1 * t303 * t366 - t266 * t291 + (-t246 * t303 - t274 * t286 - t298 * t343) * t292) * t319) * t269;
t1 = [t277 * t358 + (-t249 * t291 + t277 * t405) * t269, t222, 0, 0, t224, 0; -0.2e1 * t343 * t239 * t396 + (-t246 * t239 - t343 * t411 - (t339 * t249 + ((t230 * t269 * t384 + t372) * t267 + (t343 * t358 - t230 + (t230 - t348) * t269) * t268) * t277) * t391) * t236 + (t406 * t236 + t391 * t374) * t339 * t277 (t227 * t391 - t239 * t382) * t374 + ((-t319 * t335 - t349 * t376) * t239 + t406 * t227 + (-t382 * t225 - (t303 * t376 + t222 * t343 - t235 * t246 + t298 * t319 + (-t235 * t294 + t286 * t319) * t230) * t386 - (-t286 * t376 - t222 * t294 - t235 * t274 - t266 * t319 + (-t235 * t343 + t303 * t319) * t230) * t387) * t240) * t236, 0, 0 (t226 * t391 - t239 * t278) * t374 + (t226 * t359 + t250 * t239 + (-t278 * t225 - t226 * t249 - (t224 * t343 + t232 * t246 + t273 + (t232 * t294 + t344) * t230) * t386 - (-t224 * t294 + t232 * t274 - t247 + (t232 * t343 - t295) * t230) * t387) * t240) * t236, 0; 0.2e1 * (-t252 * t257 + t258 * t388) * t395 + ((t258 * qJD(6) - t247 * t315 - t266 * t318) * t252 + t258 * t357 + (-t257 * t234 - (-t257 * qJD(6) - t247 * t318 + t266 * t315) * t255 - t258 * t233) * t253) * t243 (t369 * t410 - t368) * t262 + (-t234 * t390 + t252 * t373) * (-t331 * t318 - t349 * t379) + (t262 * t357 + (t379 * t252 - t378 * t388) * t335 + (t252 * t403 - t388 * t402) * t318 + (t252 * t402 + t388 * t403) * t315) * t243, 0, 0, t347 * t277 * t373 + (t347 * t249 + ((-qJD(6) * t252 - 0.2e1 * t367) * t318 + (t233 * t318 + (t234 - t375) * t315) * t253) * t277) * t243, t373 + (t368 + (-t243 * t392 - t369) * t255) * t410;];
JaD_rot  = t1;
