% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRP8_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP8_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:43:49
% EndTime: 2019-02-26 22:43:50
% DurationCPUTime: 1.62s
% Computational Cost: add. (11947->153), mult. (18084->300), div. (989->12), fcn. (22873->13), ass. (0->136)
t301 = qJ(3) + qJ(4);
t298 = sin(t301);
t303 = cos(pkin(6));
t307 = cos(qJ(2));
t379 = sin(qJ(1));
t336 = t379 * t307;
t305 = sin(qJ(2));
t308 = cos(qJ(1));
t356 = t308 * t305;
t319 = -t303 * t356 - t336;
t299 = cos(t301);
t302 = sin(pkin(6));
t358 = t302 * t308;
t340 = t299 * t358;
t271 = -t298 * t319 + t340;
t360 = t302 * t305;
t342 = t298 * t360;
t281 = -t303 * t299 + t342;
t260 = atan2(-t271, t281);
t251 = sin(t260);
t252 = cos(t260);
t238 = -t251 * t271 + t252 * t281;
t236 = 0.1e1 / t238 ^ 2;
t337 = t379 * t305;
t329 = t303 * t337;
t355 = t308 * t307;
t287 = -t329 + t355;
t338 = t302 * t379;
t276 = t287 * t298 - t299 * t338;
t270 = t276 ^ 2;
t232 = t270 * t236 + 0.1e1;
t318 = -t303 * t336 - t356;
t266 = t319 * qJD(1) + t318 * qJD(2);
t300 = qJD(3) + qJD(4);
t325 = t300 * t338 + t266;
t335 = qJD(1) * t358;
t361 = t299 * t300;
t242 = t287 * t361 + t325 * t298 - t299 * t335;
t372 = t242 * t236;
t269 = t271 ^ 2;
t279 = 0.1e1 / t281 ^ 2;
t259 = t269 * t279 + 0.1e1;
t257 = 0.1e1 / t259;
t334 = qJD(2) * t379;
t268 = -qJD(1) * t329 - t305 * t334 + (qJD(2) * t303 + qJD(1)) * t355;
t293 = t298 * t358;
t333 = t379 * qJD(1);
t328 = t302 * t333;
t244 = t268 * t298 - t300 * t293 - t299 * t328 - t319 * t361;
t353 = qJD(2) * t307;
t321 = t300 * t303 + t302 * t353;
t341 = t299 * t360;
t263 = t321 * t298 + t300 * t341;
t278 = 0.1e1 / t281;
t365 = t271 * t279;
t324 = -t244 * t278 + t263 * t365;
t226 = t324 * t257;
t326 = -t251 * t281 - t252 * t271;
t221 = t326 * t226 - t251 * t244 + t252 * t263;
t235 = 0.1e1 / t238;
t237 = t235 * t236;
t377 = t221 * t237;
t351 = 0.2e1 * (-t270 * t377 + t276 * t372) / t232 ^ 2;
t383 = t263 * t279;
t339 = t303 * t355;
t284 = -t337 + t339;
t359 = t302 * t307;
t320 = -t278 * t284 + t359 * t365;
t382 = t298 * t320;
t245 = (t300 * t319 + t328) * t298 + t268 * t299 - t300 * t340;
t277 = t287 * t299 + t298 * t338;
t306 = cos(qJ(5));
t304 = sin(qJ(5));
t363 = t318 * t304;
t256 = t277 * t306 - t363;
t248 = 0.1e1 / t256;
t249 = 0.1e1 / t256 ^ 2;
t381 = -0.2e1 * t271;
t380 = 0.2e1 * t276;
t243 = t325 * t299 + (-t287 * t300 + t335) * t298;
t265 = -qJD(1) * t339 - t308 * t353 + (t303 * t334 + t333) * t305;
t233 = t256 * qJD(5) + t243 * t304 + t265 * t306;
t362 = t318 * t306;
t255 = t277 * t304 + t362;
t247 = t255 ^ 2;
t241 = t247 * t249 + 0.1e1;
t371 = t249 * t255;
t352 = qJD(5) * t255;
t234 = t243 * t306 - t265 * t304 - t352;
t374 = t234 * t248 * t249;
t376 = (t233 * t371 - t247 * t374) / t241 ^ 2;
t367 = t278 * t383;
t375 = (t244 * t365 - t269 * t367) / t259 ^ 2;
t373 = t236 * t276;
t370 = t251 * t276;
t369 = t252 * t276;
t368 = t255 * t306;
t366 = t271 * t278;
t364 = t318 * t298;
t357 = t304 * t248;
t354 = qJD(2) * t305;
t350 = -0.2e1 * t376;
t349 = 0.2e1 * t376;
t348 = -0.2e1 * t375;
t347 = t237 * t380;
t346 = t278 * t375;
t345 = t236 * t370;
t344 = t236 * t369;
t343 = t255 * t374;
t332 = 0.2e1 * t343;
t331 = t367 * t381;
t273 = -t299 * t319 - t293;
t327 = -qJD(5) * t299 * t318 + t266;
t254 = -t273 * t306 + t284 * t304;
t253 = -t273 * t304 - t284 * t306;
t323 = t249 * t368 - t357;
t282 = t303 * t298 + t341;
t322 = -t273 * t278 + t282 * t365;
t316 = -t251 + (t252 * t366 + t251) * t257;
t315 = qJD(5) * t287 + t265 * t299 - t300 * t364;
t267 = t318 * qJD(1) + t319 * qJD(2);
t264 = t321 * t299 - t300 * t342;
t262 = t287 * t304 + t299 * t362;
t261 = -t287 * t306 + t299 * t363;
t239 = 0.1e1 / t241;
t230 = 0.1e1 / t232;
t229 = t257 * t382;
t228 = t322 * t257;
t225 = t316 * t276;
t223 = (-t251 * t284 + t252 * t359) * t298 + t326 * t229;
t222 = t326 * t228 - t251 * t273 + t252 * t282;
t219 = t322 * t348 + (t282 * t331 - t245 * t278 + (t244 * t282 + t263 * t273 + t264 * t271) * t279) * t257;
t218 = t348 * t382 + (t320 * t361 + (t331 * t359 - t267 * t278 + (t263 * t284 + (t244 * t307 - t271 * t354) * t302) * t279) * t298) * t257;
t217 = t323 * t276 * t350 + (t323 * t242 + ((-qJD(5) * t248 - 0.2e1 * t343) * t306 + (t233 * t306 + (t234 - t352) * t304) * t249) * t276) * t239;
t216 = (t222 * t373 - t235 * t277) * t351 + (t222 * t221 * t347 + t243 * t235 + (-t277 * t221 - t222 * t242 - (-t219 * t271 - t228 * t244 + t264 + (-t228 * t281 - t273) * t226) * t369 - (-t219 * t281 - t228 * t263 - t245 + (t228 * t271 - t282) * t226) * t370) * t236) * t230;
t1 = [t346 * t380 + (-t242 * t278 + t276 * t383) * t257, t218, t219, t219, 0, 0; t271 * t235 * t351 + (-t244 * t235 + (t221 * t271 - t225 * t242) * t236) * t230 + (t225 * t236 * t351 + (0.2e1 * t225 * t377 - (-t226 * t257 * t366 + t348) * t345 - (t346 * t381 - t226 + (t226 - t324) * t257) * t344 - t316 * t372) * t230) * t276 (t223 * t373 - t235 * t364) * t351 + (-t223 * t372 + (t265 * t298 + t318 * t361) * t235 + (t223 * t347 - t236 * t364) * t221 - (-t218 * t271 - t229 * t244 + (-t298 * t354 + t307 * t361) * t302 + (-t229 * t281 - t284 * t298) * t226) * t344 - (-t284 * t361 - t218 * t281 - t229 * t263 - t267 * t298 + (t229 * t271 - t298 * t359) * t226) * t345) * t230, t216, t216, 0, 0; (-t248 * t253 + t254 * t371) * t349 + ((t254 * qJD(5) - t245 * t304 - t267 * t306) * t248 + t254 * t332 + (-t253 * t234 - (-t253 * qJD(5) - t245 * t306 + t267 * t304) * t255 - t254 * t233) * t249) * t239 (-t248 * t261 + t262 * t371) * t349 + (t262 * t332 - t327 * t248 * t306 + t315 * t357 + (-t327 * t255 * t304 - t262 * t233 - t261 * t234 - t315 * t368) * t249) * t239, t217, t217, t350 + 0.2e1 * (t233 * t249 * t239 + (-t239 * t374 - t249 * t376) * t255) * t255, 0;];
JaD_rot  = t1;
