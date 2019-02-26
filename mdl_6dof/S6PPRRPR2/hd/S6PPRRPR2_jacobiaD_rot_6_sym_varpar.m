% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRPR2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPRRPR2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:57
% EndTime: 2019-02-26 19:40:59
% DurationCPUTime: 1.71s
% Computational Cost: add. (9293->121), mult. (27507->247), div. (559->12), fcn. (36133->17), ass. (0->119)
t319 = cos(pkin(12));
t320 = cos(pkin(11));
t322 = cos(pkin(6));
t316 = sin(pkin(12));
t380 = sin(pkin(11));
t353 = t380 * t316;
t314 = t320 * t319 - t322 * t353;
t325 = sin(qJ(3));
t328 = cos(qJ(3));
t317 = sin(pkin(7));
t321 = cos(pkin(7));
t352 = t380 * t319;
t340 = t320 * t316 + t322 * t352;
t318 = sin(pkin(6));
t354 = t318 * t380;
t383 = -t317 * t354 + t340 * t321;
t295 = t314 * t325 + t383 * t328;
t360 = t320 * t322;
t313 = t316 * t360 + t352;
t312 = t319 * t360 - t353;
t362 = t318 * t320;
t342 = -t312 * t321 + t317 * t362;
t294 = -t313 * t325 - t342 * t328;
t290 = t294 * qJD(3);
t308 = -t312 * t317 - t321 * t362;
t327 = cos(qJ(4));
t324 = sin(qJ(4));
t381 = -t313 * t328 + t342 * t325;
t336 = t381 * t324;
t262 = qJD(4) * t336 + (qJD(4) * t308 + t290) * t327;
t285 = t308 * t324 - t327 * t381;
t282 = t285 ^ 2;
t361 = t319 * t321;
t363 = t317 * t322;
t307 = (t316 * t328 + t325 * t361) * t318 + t325 * t363;
t311 = -t317 * t318 * t319 + t321 * t322;
t301 = t307 * t327 + t311 * t324;
t298 = 0.1e1 / t301 ^ 2;
t275 = t282 * t298 + 0.1e1;
t273 = 0.1e1 / t275;
t300 = -t307 * t324 + t311 * t327;
t306 = t328 * t363 + (-t316 * t325 + t328 * t361) * t318;
t302 = t306 * qJD(3);
t281 = t300 * qJD(4) + t302 * t327;
t297 = 0.1e1 / t301;
t369 = t285 * t298;
t245 = (-t262 * t297 + t281 * t369) * t273;
t276 = atan2(-t285, t301);
t271 = sin(t276);
t272 = cos(t276);
t347 = -t271 * t301 - t272 * t285;
t241 = t347 * t245 - t262 * t271 + t272 * t281;
t255 = -t271 * t285 + t272 * t301;
t252 = 0.1e1 / t255;
t253 = 0.1e1 / t255 ^ 2;
t387 = t241 * t252 * t253;
t296 = t314 * t328 - t325 * t383;
t309 = t340 * t317 + t321 * t354;
t287 = t296 * t324 - t309 * t327;
t323 = sin(qJ(6));
t326 = cos(qJ(6));
t346 = t287 * t326 - t295 * t323;
t386 = t346 * qJD(6);
t288 = t296 * t327 + t309 * t324;
t385 = 0.2e1 * t288 * t387;
t343 = -t294 * t297 + t306 * t369;
t384 = t327 * t343;
t370 = t281 * t297 * t298;
t382 = -0.2e1 * (t262 * t369 - t282 * t370) / t275 ^ 2;
t367 = t295 * t326;
t270 = t287 * t323 + t367;
t266 = 0.1e1 / t270;
t267 = 0.1e1 / t270 ^ 2;
t292 = t295 * qJD(3);
t263 = t288 * qJD(4) - t292 * t324;
t293 = t296 * qJD(3);
t256 = t270 * qJD(6) - t263 * t326 + t293 * t323;
t265 = t346 ^ 2;
t260 = t265 * t267 + 0.1e1;
t374 = t267 * t346;
t257 = t263 * t323 + t293 * t326 + t386;
t377 = t257 * t266 * t267;
t379 = (-t256 * t374 - t265 * t377) / t260 ^ 2;
t378 = t253 * t288;
t264 = -t287 * qJD(4) - t292 * t327;
t376 = t264 * t253;
t375 = t266 * t326;
t373 = t346 * t323;
t372 = t271 * t288;
t371 = t272 * t288;
t368 = t295 * t324;
t366 = t295 * t327;
t359 = qJD(4) * t324;
t283 = t288 ^ 2;
t251 = t253 * t283 + 0.1e1;
t358 = 0.2e1 * (-t283 * t387 + t288 * t376) / t251 ^ 2;
t357 = 0.2e1 * t379;
t351 = -0.2e1 * t346 * t377;
t350 = -0.2e1 * t285 * t370;
t348 = -qJD(6) * t368 - t292;
t345 = -t267 * t373 + t375;
t284 = -t308 * t327 - t336;
t344 = t284 * t297 + t300 * t369;
t339 = qJD(4) * t366 + qJD(6) * t296 + t293 * t324;
t303 = t307 * qJD(3);
t291 = t381 * qJD(3);
t280 = -t301 * qJD(4) - t302 * t324;
t278 = t296 * t326 - t323 * t368;
t277 = t296 * t323 + t324 * t367;
t261 = t285 * qJD(4) + t290 * t324;
t258 = 0.1e1 / t260;
t249 = 0.1e1 / t251;
t247 = t273 * t384;
t246 = t344 * t273;
t243 = (-t271 * t294 + t272 * t306) * t327 + t347 * t247;
t242 = t347 * t246 + t271 * t284 + t272 * t300;
t239 = t344 * t382 + (t300 * t350 + t261 * t297 + (t262 * t300 + t280 * t285 - t281 * t284) * t298) * t273;
t238 = t382 * t384 + (-t343 * t359 + (t306 * t350 - t291 * t297 + (t262 * t306 + t281 * t294 - t285 * t303) * t298) * t327) * t273;
t1 = [0, 0, t238, t239, 0, 0; 0, 0 (t243 * t378 + t252 * t366) * t358 + ((-t293 * t327 + t295 * t359) * t252 + (-t376 + t385) * t243 + (t366 * t241 - (-t306 * t359 - t238 * t285 - t247 * t262 - t303 * t327 + (-t247 * t301 - t294 * t327) * t245) * t371 - (t294 * t359 - t238 * t301 - t247 * t281 - t291 * t327 + (t247 * t285 - t306 * t327) * t245) * t372) * t253) * t249 (t242 * t378 + t252 * t287) * t358 + (t242 * t385 - t263 * t252 + (t287 * t241 - t242 * t264 - (-t239 * t285 - t246 * t262 + t280 + (-t246 * t301 + t284) * t245) * t371 - (-t239 * t301 - t246 * t281 + t261 + (t246 * t285 - t300) * t245) * t372) * t253) * t249, 0, 0; 0, 0 (-t266 * t277 - t278 * t374) * t357 + (t278 * t351 + t348 * t266 * t323 + t339 * t375 + (t326 * t346 * t348 - t278 * t256 - t277 * t257 - t339 * t373) * t267) * t258, t345 * t288 * t357 + (-t345 * t264 + ((qJD(6) * t266 + t351) * t323 + (-t256 * t323 + (t257 + t386) * t326) * t267) * t288) * t258, 0, -0.2e1 * t379 - 0.2e1 * (t256 * t258 * t267 - (-t258 * t377 - t267 * t379) * t346) * t346;];
JaD_rot  = t1;
