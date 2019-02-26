% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP12
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRP12_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP12_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:15:32
% EndTime: 2019-02-26 22:15:34
% DurationCPUTime: 2.67s
% Computational Cost: add. (9860->184), mult. (28472->359), div. (954->12), fcn. (36187->13), ass. (0->148)
t302 = sin(qJ(2));
t305 = cos(qJ(2));
t400 = cos(pkin(6));
t402 = sin(qJ(1));
t338 = t400 * t402;
t403 = cos(qJ(1));
t319 = -t302 * t403 - t305 * t338;
t339 = t400 * t403;
t320 = t302 * t339 + t305 * t402;
t270 = -qJD(1) * t319 + qJD(2) * t320;
t301 = sin(qJ(5));
t303 = cos(qJ(5));
t330 = t302 * t338;
t354 = t402 * t302;
t271 = -qJD(1) * t330 - qJD(2) * t354 + (qJD(1) * t403 + qJD(2) * t339) * t305;
t300 = sin(pkin(6));
t304 = cos(qJ(3));
t401 = sin(qJ(3));
t340 = t403 * t401;
t282 = -t300 * t340 + t304 * t320;
t356 = t304 * t402;
t341 = t300 * t356;
t312 = -qJD(1) * t341 + qJD(3) * t282 + t271 * t401;
t331 = t305 * t339;
t292 = t354 - t331;
t377 = t300 * t304;
t342 = t403 * t377;
t315 = t320 * t401 + t342;
t409 = t292 * t301 - t303 * t315;
t234 = -qJD(5) * t409 + t270 * t303 + t301 * t312;
t260 = t292 * t303 + t301 * t315;
t233 = qJD(5) * t260 + t270 * t301 - t303 * t312;
t358 = t300 * t401;
t332 = t402 * t358;
t248 = qJD(1) * t332 - qJD(3) * t315 + t271 * t304;
t291 = t302 * t377 + t400 * t401;
t352 = t401 * qJD(2);
t376 = t300 * t305;
t278 = qJD(3) * t291 + t352 * t376;
t290 = t302 * t358 - t304 * t400;
t281 = t290 * t301 - t303 * t376;
t374 = qJD(2) * t300;
t353 = t302 * t374;
t250 = qJD(5) * t281 - t278 * t303 + t301 * t353;
t375 = t301 * t305;
t324 = t290 * t303 + t300 * t375;
t275 = 0.1e1 / t324 ^ 2;
t408 = t250 * t275;
t274 = 0.1e1 / t324;
t384 = t409 * t275;
t325 = t274 * t282 + t291 * t384;
t407 = t303 * t325;
t355 = t403 * t305;
t321 = t355 - t330;
t318 = t321 * t401;
t317 = t318 - t341;
t263 = -t301 * t319 - t303 * t317;
t268 = -qJD(1) * t331 - qJD(2) * t355 + (qJD(1) * t402 + qJD(2) * t338) * t302;
t269 = -qJD(1) * t320 + qJD(2) * t319;
t285 = t304 * t321 + t332;
t313 = -qJD(1) * t342 + qJD(3) * t285 + t269 * t401;
t232 = -qJD(5) * t263 - t268 * t303 + t301 * t313;
t264 = t301 * t317 - t303 * t319;
t255 = 0.1e1 / t264;
t256 = 0.1e1 / t264 ^ 2;
t277 = t285 ^ 2;
t246 = t256 * t277 + 0.1e1;
t242 = 0.1e1 / t246;
t388 = t242 * t256;
t247 = t269 * t304 - qJD(3) * t318 + (qJD(1) * t340 + qJD(3) * t356) * t300;
t386 = t256 * t285;
t361 = t247 * t386;
t392 = t232 * t255 * t256;
t396 = (-t277 * t392 + t361) / t246 ^ 2;
t406 = t232 * t388 + 0.2e1 * t255 * t396;
t368 = t256 * t396;
t343 = t285 * t368;
t346 = -0.2e1 * t242 * t392;
t405 = t247 * t388 + t285 * t346 - 0.2e1 * t343;
t245 = atan2(-t409, -t324);
t238 = sin(t245);
t239 = cos(t245);
t230 = -t238 * t409 - t239 * t324;
t227 = 0.1e1 / t230;
t228 = 0.1e1 / t230 ^ 2;
t404 = 0.2e1 * t263;
t254 = t263 ^ 2;
t226 = t228 * t254 + 0.1e1;
t231 = qJD(5) * t264 - t268 * t301 - t303 * t313;
t393 = t228 * t263;
t253 = t409 ^ 2;
t244 = t253 * t275 + 0.1e1;
t240 = 0.1e1 / t244;
t328 = t233 * t274 + t250 * t384;
t218 = t328 * t240;
t335 = t238 * t324 - t239 * t409;
t213 = t218 * t335 - t238 * t233 + t239 * t250;
t229 = t227 * t228;
t398 = t213 * t229;
t399 = (t231 * t393 - t254 * t398) / t226 ^ 2;
t387 = t274 * t408;
t395 = (t233 * t384 + t253 * t387) / t244 ^ 2;
t394 = t228 * t231;
t391 = t238 * t263;
t390 = t239 * t263;
t389 = t242 * t255;
t385 = t409 * t274;
t381 = t285 * t303;
t373 = qJD(3) * t304;
t372 = qJD(5) * t301;
t371 = 0.2e1 * t399;
t370 = -0.2e1 * t395;
t369 = t229 * t404;
t367 = t274 * t395;
t366 = t228 * t391;
t365 = t228 * t390;
t363 = t242 * t386;
t360 = t409 * t387;
t359 = t319 * t401;
t357 = t303 * t401;
t351 = t401 * qJD(5);
t350 = -0.2e1 * t227 * t399;
t349 = t228 * t371;
t348 = t213 * t369;
t345 = 0.2e1 * t360;
t329 = t319 * t351 + t269;
t327 = t260 * t274 + t281 * t384;
t265 = t292 * t357 + t301 * t320;
t288 = (t301 * t302 - t305 * t357) * t300;
t326 = t265 * t274 + t288 * t384;
t323 = -t238 + (-t239 * t385 + t238) * t240;
t314 = qJD(5) * t321 - t268 * t401 - t319 * t373;
t279 = t304 * t305 * t374 - qJD(3) * t290;
t266 = t301 * t321 - t319 * t357;
t252 = ((t351 + qJD(2)) * t375 + (-t305 * t373 + (t352 + qJD(5)) * t302) * t303) * t300;
t251 = qJD(5) * t324 + t278 * t301 + t303 * t353;
t235 = (-t292 * t351 + t271) * t301 + (qJD(5) * t320 + t270 * t401 + t292 * t373) * t303;
t224 = 0.1e1 / t226;
t223 = t240 * t407;
t222 = t326 * t240;
t221 = t327 * t240;
t217 = t323 * t263;
t216 = (t238 * t282 - t239 * t291) * t303 - t335 * t223;
t214 = t221 * t335 - t238 * t260 + t239 * t281;
t212 = t326 * t370 + (t288 * t345 + t235 * t274 + (t233 * t288 + t250 * t265 + t252 * t409) * t275) * t240;
t210 = t327 * t370 + (t281 * t345 + t234 * t274 + (t233 * t281 + t250 * t260 + t251 * t409) * t275) * t240;
t209 = 0.2e1 * t395 * t407 + (t325 * t372 + (-0.2e1 * t291 * t360 - t248 * t274 + (-t233 * t291 - t250 * t282 - t279 * t409) * t275) * t303) * t240;
t1 = [-t367 * t404 + (t231 * t274 + t263 * t408) * t240, t212, t209, 0, t210, 0; -t409 * t350 + (-t233 * t227 + (t213 * t409 - t217 * t231) * t228) * t224 + (t217 * t349 + (0.2e1 * t217 * t398 - (t218 * t240 * t385 + t370) * t366 - (0.2e1 * t409 * t367 - t218 + (t218 - t328) * t240) * t365 - t323 * t394) * t224) * t263, t266 * t350 + ((t301 * t329 + t303 * t314) * t227 - t266 * t228 * t213 - ((-t212 * t409 - t222 * t233 + t252 + (t222 * t324 - t265) * t218) * t239 + (t212 * t324 - t222 * t250 - t235 + (t222 * t409 - t288) * t218) * t238) * t393) * t224 + (t263 * t349 + (-t394 + t348) * t224) * (t222 * t335 - t238 * t265 + t239 * t288) (t216 * t393 + t227 * t381) * t371 + (-t216 * t394 + (-t247 * t303 + t285 * t372) * t227 + (t216 * t369 + t228 * t381) * t213 - (t291 * t372 - t209 * t409 + t223 * t233 - t279 * t303 + (-t223 * t324 + t282 * t303) * t218) * t365 - (-t282 * t372 + t209 * t324 + t223 * t250 + t248 * t303 + (-t223 * t409 + t291 * t303) * t218) * t366) * t224, 0 (t214 * t393 - t227 * t264) * t371 + (t214 * t348 + t232 * t227 + (-t264 * t213 - t214 * t231 - (-t210 * t409 - t221 * t233 + t251 + (t221 * t324 - t260) * t218) * t390 - (t210 * t324 - t221 * t250 - t234 + (t221 * t409 - t281) * t218) * t391) * t228) * t224, 0; -t234 * t363 + t248 * t389 - t260 * t405 - t282 * t406 (-t301 * t314 + t303 * t329) * t363 + (qJD(3) * t359 - t268 * t304) * t389 + t406 * t319 * t304 + t405 * (t301 * t359 + t303 * t321) 0.2e1 * t301 * t242 * t361 + t313 * t389 - t406 * t317 + (qJD(5) * t303 * t388 + (t346 - 0.2e1 * t368) * t301) * t277, 0, t343 * t404 + (t285 * t392 * t404 + (-t231 * t285 - t247 * t263) * t256) * t242, 0;];
JaD_rot  = t1;
