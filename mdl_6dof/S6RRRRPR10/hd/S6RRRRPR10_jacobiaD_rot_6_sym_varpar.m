% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR10
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
% Datum: 2019-02-26 22:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR10_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR10_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:35:53
% EndTime: 2019-02-26 22:35:55
% DurationCPUTime: 1.61s
% Computational Cost: add. (11947->154), mult. (18084->300), div. (989->12), fcn. (22873->13), ass. (0->135)
t303 = cos(pkin(6));
t305 = sin(qJ(2));
t306 = sin(qJ(1));
t350 = t306 * t305;
t334 = t303 * t350;
t308 = cos(qJ(2));
t309 = cos(qJ(1));
t348 = t308 * t309;
t288 = -t334 + t348;
t301 = qJ(3) + qJ(4);
t298 = sin(t301);
t299 = cos(t301);
t302 = sin(pkin(6));
t354 = t302 * t306;
t276 = t288 * t298 - t299 * t354;
t307 = cos(qJ(6));
t349 = t306 * t308;
t351 = t305 * t309;
t287 = t303 * t349 + t351;
t304 = sin(qJ(6));
t359 = t287 * t304;
t323 = t276 * t307 - t359;
t380 = t323 * qJD(6);
t286 = t303 * t351 + t349;
t352 = t302 * t309;
t293 = t298 * t352;
t272 = t286 * t299 - t293;
t355 = t302 * t305;
t284 = t298 * t303 + t299 * t355;
t259 = atan2(-t272, t284);
t250 = sin(t259);
t251 = cos(t259);
t237 = -t250 * t272 + t251 * t284;
t235 = 0.1e1 / t237 ^ 2;
t277 = t288 * t299 + t298 * t354;
t270 = t277 ^ 2;
t231 = t235 * t270 + 0.1e1;
t265 = -t286 * qJD(1) - t287 * qJD(2);
t300 = qJD(3) + qJD(4);
t327 = t300 * t354 + t265;
t347 = qJD(1) * t302;
t332 = t309 * t347;
t357 = t298 * t300;
t242 = -t288 * t357 + t298 * t332 + t327 * t299;
t369 = t235 * t277;
t269 = t272 ^ 2;
t281 = 0.1e1 / t284 ^ 2;
t258 = t269 * t281 + 0.1e1;
t256 = 0.1e1 / t258;
t328 = qJD(2) * t303 + qJD(1);
t346 = qJD(2) * t305;
t267 = -qJD(1) * t334 - t306 * t346 + t328 * t348;
t337 = t299 * t352;
t329 = -t267 * t299 + t300 * t337;
t333 = t306 * t347;
t244 = -t286 * t357 + t298 * t333 - t329;
t345 = qJD(2) * t308;
t319 = t300 * t303 + t302 * t345;
t336 = t300 * t355;
t263 = -t298 * t336 + t319 * t299;
t280 = 0.1e1 / t284;
t361 = t272 * t281;
t322 = -t244 * t280 + t263 * t361;
t225 = t322 * t256;
t325 = -t250 * t284 - t251 * t272;
t220 = t325 * t225 - t250 * t244 + t251 * t263;
t234 = 0.1e1 / t237;
t236 = t234 * t235;
t374 = t220 * t236;
t344 = 0.2e1 * (t242 * t369 - t270 * t374) / t231 ^ 2;
t379 = t263 * t281;
t335 = t303 * t348;
t285 = t335 - t350;
t353 = t302 * t308;
t318 = -t280 * t285 + t353 * t361;
t378 = t299 * t318;
t358 = t287 * t307;
t255 = t276 * t304 + t358;
t247 = 0.1e1 / t255;
t248 = 0.1e1 / t255 ^ 2;
t377 = -0.2e1 * t272;
t376 = 0.2e1 * t277;
t241 = (t288 * t300 - t332) * t299 + t327 * t298;
t264 = -qJD(1) * t335 - t309 * t345 + t328 * t350;
t232 = t255 * qJD(6) - t241 * t307 - t264 * t304;
t246 = t323 ^ 2;
t240 = t246 * t248 + 0.1e1;
t367 = t248 * t323;
t233 = t241 * t304 - t264 * t307 + t380;
t371 = t233 * t247 * t248;
t373 = (-t232 * t367 - t246 * t371) / t240 ^ 2;
t363 = t280 * t379;
t372 = (t244 * t361 - t269 * t363) / t258 ^ 2;
t370 = t235 * t242;
t368 = t247 * t307;
t366 = t250 * t277;
t365 = t251 * t277;
t364 = t323 * t304;
t362 = t272 * t280;
t360 = t287 * t299;
t356 = t299 * t300;
t343 = 0.2e1 * t373;
t342 = -0.2e1 * t372;
t341 = t236 * t376;
t340 = t280 * t372;
t339 = t235 * t366;
t338 = t235 * t365;
t331 = -0.2e1 * t323 * t371;
t330 = t363 * t377;
t326 = -qJD(6) * t287 * t298 + t265;
t271 = t286 * t298 + t337;
t324 = -t271 * t307 - t285 * t304;
t253 = -t271 * t304 + t285 * t307;
t321 = -t248 * t364 + t368;
t283 = -t298 * t355 + t299 * t303;
t320 = t271 * t280 + t283 * t361;
t317 = -t250 + (t251 * t362 + t250) * t256;
t243 = t267 * t298 + t286 * t356 - t300 * t293 - t299 * t333;
t316 = qJD(6) * t288 - t264 * t298 + t287 * t356;
t266 = -t287 * qJD(1) - t286 * qJD(2);
t262 = -t319 * t298 - t299 * t336;
t261 = t288 * t307 - t298 * t359;
t260 = t288 * t304 + t298 * t358;
t238 = 0.1e1 / t240;
t229 = 0.1e1 / t231;
t228 = t256 * t378;
t227 = t320 * t256;
t224 = t317 * t277;
t222 = (-t250 * t285 + t251 * t353) * t299 + t325 * t228;
t221 = t325 * t227 + t250 * t271 + t251 * t283;
t218 = t320 * t342 + (t283 * t330 + t243 * t280 + (t244 * t283 + t262 * t272 - t263 * t271) * t281) * t256;
t217 = t342 * t378 + (-t318 * t357 + (t330 * t353 - t266 * t280 + (t263 * t285 + (t244 * t308 - t272 * t346) * t302) * t281) * t299) * t256;
t216 = t321 * t277 * t343 + (-t321 * t242 + ((qJD(6) * t247 + t331) * t304 + (-t232 * t304 + (t233 + t380) * t307) * t248) * t277) * t238;
t215 = (t221 * t369 + t234 * t276) * t344 + (t221 * t220 * t341 - t241 * t234 + (t276 * t220 - t221 * t242 - (-t218 * t272 - t227 * t244 + t262 + (-t227 * t284 + t271) * t225) * t365 - (-t218 * t284 - t227 * t263 + t243 + (t227 * t272 - t283) * t225) * t366) * t235) * t229;
t1 = [t340 * t376 + (-t242 * t280 + t277 * t379) * t256, t217, t218, t218, 0, 0; t272 * t234 * t344 + (((t286 * t300 - t333) * t298 + t329) * t234 + (t220 * t272 - t224 * t242) * t235) * t229 + (t224 * t235 * t344 + (0.2e1 * t224 * t374 - (-t225 * t256 * t362 + t342) * t339 - (t340 * t377 - t225 + (t225 - t322) * t256) * t338 - t317 * t370) * t229) * t277 (t222 * t369 + t234 * t360) * t344 + (-t222 * t370 + (t264 * t299 + t287 * t357) * t234 + (t222 * t341 + t235 * t360) * t220 - (-t217 * t272 - t228 * t244 + (-t299 * t346 - t308 * t357) * t302 + (-t228 * t284 - t285 * t299) * t225) * t338 - (t285 * t357 - t217 * t284 - t228 * t263 - t266 * t299 + (t228 * t272 - t299 * t353) * t225) * t339) * t229, t215, t215, 0, 0; (t247 * t324 - t253 * t367) * t343 + ((t253 * qJD(6) + t243 * t307 + t266 * t304) * t247 + t253 * t331 + (t324 * t233 + (t324 * qJD(6) - t243 * t304 + t266 * t307) * t323 - t253 * t232) * t248) * t238 (-t247 * t260 - t261 * t367) * t343 + (t261 * t331 + t326 * t247 * t304 + t316 * t368 + (t307 * t323 * t326 - t261 * t232 - t260 * t233 - t316 * t364) * t248) * t238, t216, t216, 0, -0.2e1 * t373 - 0.2e1 * (t232 * t238 * t248 - (-t238 * t371 - t248 * t373) * t323) * t323;];
JaD_rot  = t1;
