% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR14V3
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR14V3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14V3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobiaD_rot_6_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:12:09
% EndTime: 2019-04-12 15:12:12
% DurationCPUTime: 2.28s
% Computational Cost: add. (6638->203), mult. (19789->388), div. (992->12), fcn. (24542->13), ass. (0->160)
t301 = sin(qJ(5));
t303 = sin(qJ(2));
t302 = sin(qJ(4));
t370 = qJD(4) * t302;
t308 = cos(qJ(2));
t307 = cos(qJ(4));
t337 = qJD(2) * t307 - qJD(5);
t407 = t337 * t308;
t413 = (-t303 * t370 + t407) * t301;
t304 = sin(qJ(1));
t372 = qJD(2) * t308;
t309 = cos(qJ(1));
t379 = t303 * t309;
t412 = qJD(1) * t379 + t304 * t372;
t306 = cos(qJ(5));
t338 = qJD(5) * t307 - qJD(2);
t410 = (t338 * t301 + t306 * t370) * t303 - t306 * t407;
t374 = t309 * t302;
t376 = t307 * t308;
t284 = t304 * t376 - t374;
t381 = t303 * t306;
t267 = t284 * t301 - t304 * t381;
t380 = t303 * t307;
t281 = t301 * t380 + t306 * t308;
t260 = atan2(-t267, t281);
t248 = sin(t260);
t249 = cos(t260);
t229 = -t248 * t267 + t249 * t281;
t227 = 0.1e1 / t229 ^ 2;
t375 = t307 * t309;
t378 = t304 * t302;
t287 = t308 * t375 + t378;
t272 = t287 * t301 - t306 * t379;
t266 = t272 ^ 2;
t225 = t227 * t266 + 0.1e1;
t339 = -qJD(1) * t308 + qJD(4);
t368 = qJD(4) * t308;
t340 = qJD(1) - t368;
t371 = qJD(2) * t309;
t350 = t303 * t371;
t256 = t340 * t374 + (t339 * t304 - t350) * t307;
t273 = t287 * t306 + t301 * t379;
t373 = qJD(1) * t304;
t322 = t303 * t373 - t308 * t371;
t233 = t273 * qJD(5) + t256 * t301 + t322 * t306;
t396 = t227 * t272;
t265 = t267 ^ 2;
t279 = 0.1e1 / t281 ^ 2;
t259 = t265 * t279 + 0.1e1;
t250 = 0.1e1 / t259;
t369 = qJD(4) * t307;
t344 = t309 * t369;
t345 = t304 * t370;
t382 = t303 * t304;
t351 = qJD(2) * t382;
t258 = t287 * qJD(1) - t307 * t351 - t308 * t345 - t344;
t269 = t284 * t306 + t301 * t382;
t235 = t269 * qJD(5) + t258 * t301 - t412 * t306;
t332 = t338 * t306;
t253 = t303 * t332 + t413;
t278 = 0.1e1 / t281;
t386 = t267 * t279;
t328 = -t235 * t278 + t253 * t386;
t216 = t328 * t250;
t333 = -t248 * t281 - t249 * t267;
t210 = t333 * t216 - t235 * t248 + t249 * t253;
t226 = 0.1e1 / t229;
t228 = t226 * t227;
t401 = t210 * t228;
t364 = 0.2e1 * (t233 * t396 - t266 * t401) / t225 ^ 2;
t409 = t253 * t279;
t324 = t308 * t378 + t375;
t383 = t302 * t303;
t355 = t267 * t383;
t323 = -t278 * t324 + t279 * t355;
t408 = t301 * t323;
t282 = -t301 * t308 + t306 * t380;
t276 = t282 * t309;
t405 = qJD(6) * t276 + t322 * t302 - t303 * t344;
t367 = qJD(5) * t303;
t236 = (-qJD(5) * t284 + t412) * t301 + (t304 * t367 + t258) * t306;
t377 = t304 * t307;
t286 = t308 * t374 - t377;
t300 = sin(qJ(6));
t305 = cos(qJ(6));
t247 = t273 * t305 + t286 * t300;
t241 = 0.1e1 / t247;
t242 = 0.1e1 / t247 ^ 2;
t404 = -0.2e1 * t267;
t403 = 0.2e1 * t272;
t234 = (t309 * t367 + t256) * t306 + (-qJD(5) * t287 - t322) * t301;
t255 = t324 * qJD(1) + t302 * t350 - t308 * t344 - t345;
t218 = t247 * qJD(6) + t234 * t300 + t255 * t305;
t246 = t273 * t300 - t286 * t305;
t240 = t246 ^ 2;
t232 = t240 * t242 + 0.1e1;
t393 = t242 * t246;
t365 = qJD(6) * t246;
t219 = t234 * t305 - t255 * t300 - t365;
t398 = t219 * t241 * t242;
t400 = (t218 * t393 - t240 * t398) / t232 ^ 2;
t388 = t278 * t409;
t399 = (t235 * t386 - t265 * t388) / t259 ^ 2;
t397 = t227 * t233;
t395 = t241 * t300;
t394 = t241 * t305;
t392 = t246 * t300;
t391 = t246 * t305;
t390 = t248 * t272;
t389 = t249 * t272;
t387 = t267 * t278;
t385 = t286 * t301;
t384 = t286 * t306;
t366 = qJD(5) * t306;
t363 = -0.2e1 * t400;
t362 = 0.2e1 * t400;
t361 = -0.2e1 * t399;
t360 = t228 * t403;
t359 = t278 * t399;
t358 = t246 * t398;
t357 = t227 * t390;
t356 = t227 * t389;
t354 = t303 * t374;
t343 = t210 * t360;
t342 = 0.2e1 * t358;
t341 = t388 * t404;
t334 = qJD(6) * t384 + t256;
t245 = -t269 * t305 - t300 * t324;
t244 = -t269 * t300 + t305 * t324;
t329 = -qJD(6) * t354 + t282 * t373 + t410 * t309;
t327 = t242 * t391 - t395;
t326 = -t269 * t278 + t282 * t386;
t274 = t281 * t304;
t285 = t301 * t376 - t381;
t325 = t274 * t278 + t285 * t386;
t321 = -t302 * t372 - t303 * t369;
t319 = -t248 + (t249 * t387 + t248) * t250;
t318 = qJD(1) * t281;
t317 = qJD(5) * t385 + qJD(6) * t287 + t255 * t306;
t275 = t281 * t309;
t264 = -t276 * t305 - t300 * t354;
t263 = -t276 * t300 + t305 * t354;
t262 = t287 * t300 - t305 * t384;
t261 = -t287 * t305 - t300 * t384;
t257 = t340 * t377 + (t339 * t309 + t351) * t302;
t252 = t308 * t332 + (-t302 * t368 - t337 * t303) * t301;
t239 = t253 * t304 + t309 * t318;
t230 = 0.1e1 / t232;
t223 = 0.1e1 / t225;
t222 = t250 * t408;
t221 = t325 * t250;
t220 = t326 * t250;
t215 = t319 * t272;
t214 = (t248 * t324 - t249 * t383) * t301 - t333 * t222;
t212 = t333 * t221 + t248 * t274 + t249 * t285;
t211 = t333 * t220 - t248 * t269 + t249 * t282;
t209 = t325 * t361 + (t285 * t341 + t239 * t278 + (t235 * t285 + t252 * t267 - t253 * t274) * t279) * t250;
t207 = t326 * t361 + (t282 * t341 - t236 * t278 + (t235 * t282 + t253 * t269 - t267 * t410) * t279) * t250;
t206 = 0.2e1 * t399 * t408 + (-t323 * t366 + (0.2e1 * t355 * t388 - t257 * t278 + (-t235 * t383 - t253 * t324 + t321 * t267) * t279) * t301) * t250;
t1 = [t359 * t403 + (-t233 * t278 + t272 * t409) * t250, t209, 0, t206, t207, 0; t267 * t226 * t364 + (-t235 * t226 + (t210 * t267 - t215 * t233) * t227) * t223 + (t215 * t227 * t364 + (0.2e1 * t215 * t401 - (-t216 * t250 * t387 + t361) * t357 - (t359 * t404 - t216 + (t216 - t328) * t250) * t356 - t319 * t397) * t223) * t272 (t212 * t396 + t226 * t275) * t364 + (t212 * t343 + (t275 * t210 - t212 * t233 - (-t209 * t267 - t221 * t235 + t252 + (-t221 * t281 + t274) * t216) * t389 - (-t209 * t281 - t221 * t253 + t239 + (t221 * t267 - t285) * t216) * t390) * t227 + (t304 * t318 + (-t338 * t381 - t413) * t309) * t226) * t223, 0 (t214 * t396 + t226 * t385) * t364 + (-t214 * t397 + (t255 * t301 - t286 * t366) * t226 + (t214 * t360 + t227 * t385) * t210 - (t324 * t366 - t206 * t281 + t222 * t253 - t257 * t301 + (-t222 * t267 + t301 * t383) * t216) * t357 - (-t366 * t383 - t206 * t267 - (-t216 * t281 - t235) * t222 + (t216 * t324 + t321) * t301) * t356) * t223 (t211 * t396 - t226 * t273) * t364 + (t211 * t343 + t234 * t226 + (-t273 * t210 - t211 * t233 - (-t207 * t267 - t220 * t235 - t410 + (-t220 * t281 - t269) * t216) * t389 - (-t207 * t281 - t220 * t253 - t236 + (t220 * t267 - t282) * t216) * t390) * t227) * t223, 0; (-t241 * t244 + t245 * t393) * t362 + ((t245 * qJD(6) - t236 * t300 - t257 * t305) * t241 + t245 * t342 + (-t244 * t219 - (-t244 * qJD(6) - t236 * t305 + t257 * t300) * t246 - t245 * t218) * t242) * t230 (-t241 * t263 + t264 * t393) * t362 + (t264 * t342 + t329 * t395 - t405 * t394 + (-t264 * t218 - t263 * t219 - t329 * t391 - t405 * t392) * t242) * t230, 0 (-t241 * t261 + t262 * t393) * t362 + (t262 * t342 - t334 * t394 + t317 * t395 + (-t262 * t218 - t261 * t219 - t317 * t391 - t334 * t392) * t242) * t230, t327 * t272 * t363 + (t327 * t233 + ((-qJD(6) * t241 - 0.2e1 * t358) * t305 + (t218 * t305 + (t219 - t365) * t300) * t242) * t272) * t230, t363 + 0.2e1 * (t218 * t242 * t230 + (-t230 * t398 - t242 * t400) * t246) * t246;];
JaD_rot  = t1;
