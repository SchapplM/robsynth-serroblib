% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR12_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR12_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:00:25
% EndTime: 2019-02-26 22:00:27
% DurationCPUTime: 2.12s
% Computational Cost: add. (11947->152), mult. (18084->300), div. (989->12), fcn. (22873->13), ass. (0->132)
t295 = qJD(4) + qJD(5);
t299 = sin(qJ(2));
t300 = sin(qJ(1));
t302 = cos(qJ(2));
t303 = cos(qJ(1));
t382 = cos(pkin(6));
t333 = t303 * t382;
t319 = -t299 * t333 - t300 * t302;
t334 = t300 * t382;
t320 = t303 * t299 + t302 * t334;
t297 = sin(pkin(6));
t358 = t297 * t303;
t394 = t320 * qJD(1) - t319 * qJD(2) + t295 * t358;
t284 = t300 * t299 - t302 * t333;
t296 = qJ(4) + qJ(5);
t293 = sin(t296);
t294 = cos(t296);
t322 = t284 * t294 + t293 * t358;
t269 = t322 ^ 2;
t359 = t297 * t302;
t340 = t294 * t359;
t281 = t382 * t293 + t340;
t279 = 0.1e1 / t281 ^ 2;
t258 = t269 * t279 + 0.1e1;
t256 = 0.1e1 / t258;
t355 = qJD(1) * t300;
t337 = t297 * t355;
t392 = t394 * t294;
t241 = (t284 * t295 + t337) * t293 - t392;
t354 = qJD(2) * t299;
t317 = -t382 * t295 + t297 * t354;
t341 = t293 * t359;
t265 = -t317 * t294 - t295 * t341;
t278 = 0.1e1 / t281;
t368 = t322 * t279;
t326 = -t241 * t278 - t265 * t368;
t225 = t326 * t256;
t259 = atan2(t322, t281);
t250 = sin(t259);
t251 = cos(t259);
t328 = -t250 * t281 + t251 * t322;
t220 = t328 * t225 - t250 * t241 + t251 * t265;
t237 = t250 * t322 + t251 * t281;
t235 = 0.1e1 / t237 ^ 2;
t393 = t220 * t235;
t362 = t294 * t295;
t242 = t284 * t362 + t394 * t293 + t294 * t337;
t360 = t297 * t300;
t271 = t320 * t293 + t294 * t360;
t329 = t299 * t334;
t357 = t303 * t302;
t286 = -t329 + t357;
t298 = sin(qJ(6));
t301 = cos(qJ(6));
t252 = t271 * t298 - t286 * t301;
t391 = 0.2e1 * t252;
t234 = 0.1e1 / t237;
t390 = t234 * t393;
t270 = t293 * t360 - t320 * t294;
t330 = 0.2e1 * t270 * t390;
t314 = (t382 * qJD(1) + qJD(2)) * t357 - qJD(2) * t329 - t299 * t355;
t383 = t295 * t360 - t314;
t384 = qJD(1) * t358 + t320 * t295;
t244 = t384 * t293 + t383 * t294;
t374 = t244 * t235;
t389 = -t374 + t330;
t388 = t265 * t279;
t361 = t297 * t299;
t343 = t322 * t361;
t321 = -t278 * t319 + t279 * t343;
t387 = t294 * t321;
t386 = -t286 * t293 * qJD(6) - t314;
t385 = -t320 * qJD(6) + t286 * t362;
t253 = t271 * t301 + t286 * t298;
t247 = 0.1e1 / t253;
t248 = 0.1e1 / t253 ^ 2;
t268 = t270 ^ 2;
t231 = t268 * t235 + 0.1e1;
t381 = (-t268 * t390 + t270 * t374) / t231 ^ 2;
t245 = -t383 * t293 + t384 * t294;
t266 = t319 * qJD(1) - t320 * qJD(2);
t232 = t253 * qJD(6) + t245 * t298 - t266 * t301;
t246 = t252 ^ 2;
t240 = t246 * t248 + 0.1e1;
t373 = t248 * t252;
t352 = qJD(6) * t252;
t233 = t245 * t301 + t266 * t298 - t352;
t377 = t233 * t247 * t248;
t380 = (t232 * t373 - t246 * t377) / t240 ^ 2;
t370 = t278 * t388;
t378 = (-t241 * t368 - t269 * t370) / t258 ^ 2;
t376 = t235 * t270;
t238 = 0.1e1 / t240;
t375 = t238 * t248;
t372 = t250 * t270;
t371 = t251 * t270;
t369 = t322 * t278;
t282 = t382 * t294 - t341;
t367 = t322 * t282;
t366 = t286 * t294;
t365 = t293 * t295;
t364 = t293 * t298;
t363 = t293 * t301;
t353 = qJD(2) * t302;
t351 = 0.2e1 * t381;
t350 = -0.2e1 * t380;
t349 = -0.2e1 * t378;
t348 = 0.2e1 * t378;
t346 = t248 * t380;
t345 = t232 * t375;
t344 = t252 * t377;
t332 = t278 * t348;
t331 = 0.2e1 * t344;
t323 = -t284 * t293 + t294 * t358;
t255 = t298 * t319 + t301 * t323;
t254 = t298 * t323 - t301 * t319;
t325 = -t298 * t247 + t301 * t373;
t324 = -t278 * t323 + t279 * t367;
t318 = -t250 + (-t251 * t369 + t250) * t256;
t267 = -qJD(1) * t329 - t300 * t354 + (qJD(2) * t382 + qJD(1)) * t357;
t264 = t317 * t293 - t295 * t340;
t261 = t286 * t363 - t320 * t298;
t229 = 0.1e1 / t231;
t228 = t256 * t387;
t227 = t324 * t256;
t222 = (-t250 * t319 - t251 * t361) * t294 + t328 * t228;
t221 = -t328 * t227 + t250 * t323 + t251 * t282;
t218 = t324 * t348 + (0.2e1 * t367 * t370 - t242 * t278 + (t241 * t282 - t264 * t322 - t265 * t323) * t279) * t256;
t217 = t349 * t387 + (-t321 * t365 + (-0.2e1 * t343 * t370 + t267 * t278 + (t265 * t319 + (-t241 * t299 + t322 * t353) * t297) * t279) * t294) * t256;
t216 = t325 * t270 * t350 + (t325 * t244 + ((-qJD(6) * t247 - 0.2e1 * t344) * t301 + (t232 * t301 + (t233 - t352) * t298) * t248) * t270) * t238;
t215 = (t221 * t376 - t234 * t271) * t351 + (t221 * t330 + t245 * t234 + (-t271 * t220 - t221 * t244 - (t218 * t322 + t227 * t241 + t264 + (t227 * t281 + t323) * t225) * t371 - (-t218 * t281 + t227 * t265 - t242 + (t227 * t322 - t282) * t225) * t372) * t235) * t229;
t1 = [t270 * t332 + (-t244 * t278 + t270 * t388) * t256, t217, 0, t218, t218, 0; -0.2e1 * t322 * t234 * t381 + ((-t284 * t365 - t293 * t337 + t392) * t234 - t322 * t393 - (t318 * t244 + ((t225 * t256 * t369 + t349) * t250 + (t322 * t332 - t225 + (t225 - t326) * t256) * t251) * t270) * t376) * t229 + (t229 * t389 + t376 * t351) * t318 * t270 (t222 * t376 + t234 * t366) * t351 + ((-t266 * t294 + t286 * t365) * t234 + t389 * t222 + (t366 * t220 - (t217 * t322 - t228 * t241 + (-t294 * t353 + t299 * t365) * t297 + (-t228 * t281 - t294 * t319) * t225) * t371 - (t319 * t365 - t217 * t281 - t228 * t265 + t267 * t294 + (-t228 * t322 + t294 * t361) * t225) * t372) * t235) * t229, 0, t215, t215, 0; 0.2e1 * (-t247 * t254 + t255 * t373) * t380 + ((t255 * qJD(6) - t242 * t298 + t267 * t301) * t247 + t255 * t331 + (-t254 * t233 - (-t254 * qJD(6) - t242 * t301 - t267 * t298) * t252 - t255 * t232) * t248) * t238 (t346 * t391 - t345) * t261 + (-t233 * t375 + t247 * t350) * (t286 * t364 + t320 * t301) + (t261 * t331 + (t364 * t247 - t363 * t373) * t266 + (-t386 * t247 - t385 * t373) * t301 + (t385 * t247 - t386 * t373) * t298) * t238, 0, t216, t216, t350 + (t345 + (-t238 * t377 - t346) * t252) * t391;];
JaD_rot  = t1;
