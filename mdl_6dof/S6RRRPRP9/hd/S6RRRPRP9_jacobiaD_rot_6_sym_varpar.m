% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP9
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRP9_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP9_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:13:43
% EndTime: 2019-02-26 22:13:46
% DurationCPUTime: 2.44s
% Computational Cost: add. (7107->148), mult. (21419->299), div. (940->12), fcn. (26530->11), ass. (0->135)
t288 = cos(qJ(3));
t284 = sin(qJ(3));
t290 = cos(qJ(1));
t356 = t290 * t284;
t286 = sin(qJ(1));
t289 = cos(qJ(2));
t359 = t286 * t289;
t261 = t288 * t359 - t356;
t283 = sin(qJ(5));
t287 = cos(qJ(5));
t338 = t284 * t359;
t310 = t288 * t290 + t338;
t395 = t261 * t283 - t287 * t310;
t421 = -0.2e1 * t395;
t357 = t289 * t290;
t262 = t286 * t284 + t288 * t357;
t285 = sin(qJ(2));
t361 = t285 * t286;
t405 = -qJD(2) * t361 - qJD(3) * t290;
t233 = qJD(1) * t262 - qJD(3) * t338 + t405 * t288;
t243 = t261 * t287 + t283 * t310;
t279 = t286 * t288;
t354 = qJD(1) * t290;
t355 = qJD(1) * t286;
t299 = t289 * (qJD(3) * t279 + t284 * t354) - t288 * t355 + t405 * t284;
t216 = qJD(5) * t243 + t233 * t283 - t287 * t299;
t217 = -qJD(5) * t395 + t233 * t287 + t283 * t299;
t362 = t284 * t287;
t339 = t285 * t362;
t363 = t283 * t288;
t340 = t285 * t363;
t353 = qJD(2) * t289;
t398 = -t283 * t284 - t287 * t288;
t220 = -qJD(5) * t339 - qJD(3) * t340 + (qJD(3) * t362 + qJD(5) * t363) * t285 + t398 * t353;
t397 = qJD(3) - qJD(5);
t297 = t397 * t398;
t318 = -t362 + t363;
t260 = t318 * t289;
t304 = qJD(2) * t260;
t221 = t285 * t297 + t304;
t235 = t395 ^ 2;
t259 = t318 * t285;
t256 = 0.1e1 / t259 ^ 2;
t226 = t235 * t256 + 0.1e1;
t224 = 0.1e1 / t226;
t255 = 0.1e1 / t259;
t258 = t398 * t285;
t400 = t221 * t256;
t382 = t255 * t400;
t417 = t382 * t421;
t420 = (t217 * t255 - (-t216 * t258 - t220 * t395 + t221 * t243) * t256 + t258 * t417) * t224;
t370 = t395 * t256;
t314 = -t216 * t255 + t221 * t370;
t200 = t314 * t224;
t406 = t243 * t255 + t258 * t370;
t410 = t406 * t224;
t419 = (t395 * t410 - t258) * t200 - t410 * t221 + t217;
t418 = -(t259 * t410 - t243) * t200 - t410 * t216 + t220;
t227 = atan2(-t395, t259);
t222 = sin(t227);
t223 = cos(t227);
t319 = -t222 * t259 - t223 * t395;
t411 = t222 * t243 + t223 * t258 + t319 * t410;
t352 = qJD(2) * t290;
t334 = t285 * t352;
t232 = (-qJD(3) * t289 + qJD(1)) * t356 + (-t334 + (-qJD(1) * t289 + qJD(3)) * t286) * t288;
t301 = -qJD(1) * t338 + qJD(3) * t262 - t284 * t334 - t288 * t354;
t324 = t289 * t356 - t279;
t394 = t262 * t283 - t324 * t287;
t214 = -qJD(5) * t394 + t232 * t287 + t283 * t301;
t249 = t262 * t287 + t324 * t283;
t213 = qJD(5) * t249 + t232 * t283 - t287 * t301;
t402 = -0.2e1 * t394;
t195 = t319 * t200 - t216 * t222 + t221 * t223;
t212 = -t222 * t395 + t223 * t259;
t210 = 0.1e1 / t212 ^ 2;
t401 = t195 * t210;
t309 = t398 * t290;
t236 = t394 ^ 2;
t208 = t210 * t236 + 0.1e1;
t206 = 0.1e1 / t208;
t209 = 0.1e1 / t212;
t384 = t213 * t210;
t389 = t209 * t401;
t390 = (-t236 * t389 + t384 * t394) / t208 ^ 2;
t399 = -t206 * t401 - 0.2e1 * t209 * t390;
t392 = 0.2e1 * t394;
t331 = t389 * t392;
t350 = 0.2e1 * t390;
t385 = t210 * t394;
t396 = t350 * t385 + (t331 - t384) * t206;
t238 = 0.1e1 / t249 ^ 2;
t281 = t285 ^ 2;
t282 = t290 ^ 2;
t365 = t281 * t282;
t234 = t238 * t365 + 0.1e1;
t230 = 0.1e1 / t234;
t360 = t285 * t290;
t342 = t238 * t360;
t237 = 0.1e1 / t249;
t383 = t214 * t237 * t238;
t387 = (-t365 * t383 + (-t281 * t286 * t354 + t282 * t285 * t353) * t238) / t234 ^ 2;
t323 = t342 * t387;
t327 = t360 * t383;
t378 = t230 * t238;
t328 = t361 * t378;
t333 = t289 * t352;
t393 = qJD(1) * t328 + 0.2e1 * t230 * t327 - t333 * t378 + 0.2e1 * t323;
t252 = t285 * t309;
t391 = 0.2e1 * t252;
t388 = (t216 * t370 - t235 * t382) / t226 ^ 2;
t386 = t206 * t209;
t381 = t222 * t394;
t380 = t223 * t394;
t372 = t238 * t285;
t371 = t395 * t255;
t349 = -0.2e1 * t388;
t348 = t255 * t388;
t347 = t206 * t385;
t343 = t237 * t387;
t329 = t230 * t342;
t250 = (-t339 + t340) * t286;
t311 = t250 * t255 + t260 * t370;
t308 = qJD(1) * t318;
t307 = qJD(2) * t318;
t305 = -t222 + (t223 * t371 + t222) * t224;
t251 = t318 * t360;
t219 = -t285 * t307 + t289 * t297;
t218 = t286 * t304 + (t286 * t297 + t290 * t308) * t285;
t205 = t311 * t224;
t198 = t319 * t205 + t222 * t250 + t223 * t260;
t194 = t311 * t349 + (t260 * t417 + t218 * t255 + (t216 * t260 + t219 * t395 - t221 * t250) * t256) * t224;
t192 = 0.2e1 * t406 * t388 - t420;
t191 = t406 * t349 + t420;
t1 = [t348 * t392 + (-t213 * t255 + t394 * t400) * t224, t194, t191, 0, t192, 0; -t216 * t386 - (t305 * t213 + ((-t200 * t224 * t371 + t349) * t222 + (t348 * t421 - t200 + (t200 - t314) * t224) * t223) * t394) * t347 - t399 * t395 + t396 * t305 * t394 (t198 * t385 + t209 * t251) * t350 + (t198 * t331 + (t251 * t195 - t198 * t213 - (-t194 * t395 - t205 * t216 + t219 + (-t205 * t259 + t250) * t200) * t380 - (-t194 * t259 - t205 * t221 + t218 + (t205 * t395 - t260) * t200) * t381) * t210 + (-t307 * t357 + (t286 * t308 - t397 * t309) * t285) * t209) * t206, -t214 * t386 - ((-t191 * t395 + t418) * t223 + (-t191 * t259 + t419) * t222) * t347 - t399 * t249 + t396 * t411, 0 (-t209 * t249 - t385 * t411) * t350 + (-t411 * t331 + t214 * t209 + (-t249 * t195 + t411 * t213 - (-t192 * t395 - t418) * t380 - (-t192 * t259 - t419) * t381) * t210) * t206, 0; t214 * t328 + t217 * t329 + 0.2e1 * t343 * t361 + (-t285 * t354 - t286 * t353) * t237 * t230 - t393 * t243 (-t237 * t289 + t252 * t372) * t230 * t355 + ((-0.2e1 * t343 + (-qJD(2) * t252 - t214) * t378) * t289 + (t238 * t387 * t391 + (t383 * t391 + (t397 * t290 * t318 + t355 * t398) * t372 + (-t238 * t289 * t309 - t237) * qJD(2)) * t230) * t285) * t290, -t213 * t329 + t393 * t394, 0, t323 * t402 + (t327 * t402 + (t213 * t360 + (-t285 * t355 + t333) * t394) * t238) * t230, 0;];
JaD_rot  = t1;
