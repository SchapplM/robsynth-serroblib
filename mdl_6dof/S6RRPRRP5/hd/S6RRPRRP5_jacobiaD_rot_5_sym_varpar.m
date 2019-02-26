% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRP5_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP5_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:48:21
% EndTime: 2019-02-26 21:48:23
% DurationCPUTime: 1.88s
% Computational Cost: add. (8421->153), mult. (24081->307), div. (726->12), fcn. (31139->15), ass. (0->130)
t315 = sin(pkin(11));
t316 = cos(pkin(11));
t320 = sin(qJ(2));
t323 = cos(qJ(2));
t305 = t320 * t315 - t323 * t316;
t317 = cos(pkin(6));
t302 = t305 * t317;
t297 = qJD(2) * t302;
t345 = t323 * t315 + t320 * t316;
t304 = t345 * qJD(2);
t324 = cos(qJ(1));
t389 = sin(pkin(6));
t354 = t324 * t389;
t303 = t345 * t317;
t390 = sin(qJ(1));
t357 = t390 * t303;
t395 = -t390 * t304 - qJD(1) * t357 + (-qJD(1) * t305 - t297) * t324 - qJD(4) * t354;
t319 = sin(qJ(4));
t322 = cos(qJ(4));
t338 = -t324 * t303 + t390 * t305;
t275 = -t319 * t338 + t322 * t354;
t355 = t323 * t389;
t356 = t320 * t389;
t332 = -t315 * t355 - t316 * t356;
t291 = -t317 * t322 - t319 * t332;
t270 = atan2(-t275, t291);
t265 = sin(t270);
t266 = cos(t270);
t241 = -t265 * t275 + t266 * t291;
t239 = 0.1e1 / t241 ^ 2;
t339 = -t324 * t305 - t357;
t350 = t390 * t389;
t335 = -t319 * t339 + t322 * t350;
t274 = t335 ^ 2;
t237 = t274 * t239 + 0.1e1;
t281 = t319 * t350 + t322 * t339;
t331 = t338 * qJD(1) + t390 * t297 - t324 * t304;
t348 = qJD(1) * t354;
t245 = t281 * qJD(4) + t319 * t331 - t322 * t348;
t382 = t245 * t239;
t273 = t275 ^ 2;
t289 = 0.1e1 / t291 ^ 2;
t269 = t273 * t289 + 0.1e1;
t267 = 0.1e1 / t269;
t344 = qJD(1) * t350;
t368 = qJD(4) * t322;
t247 = t319 * t395 - t322 * t344 - t338 * t368;
t292 = t317 * t319 - t322 * t332;
t300 = -t315 * t356 + t316 * t355;
t296 = t300 * qJD(2);
t271 = t292 * qJD(4) + t296 * t319;
t288 = 0.1e1 / t291;
t376 = t275 * t289;
t343 = -t247 * t288 + t271 * t376;
t229 = t343 * t267;
t346 = -t265 * t291 - t266 * t275;
t224 = t346 * t229 - t265 * t247 + t266 * t271;
t238 = 0.1e1 / t241;
t240 = t238 * t239;
t387 = t224 * t240;
t366 = 0.2e1 * (-t274 * t387 - t335 * t382) / t237 ^ 2;
t394 = t271 * t289;
t283 = -t324 * t302 - t345 * t390;
t340 = -t283 * t288 + t300 * t376;
t393 = t319 * t340;
t248 = (qJD(4) * t338 + t344) * t319 + t395 * t322;
t286 = t390 * t302 - t324 * t345;
t318 = sin(qJ(5));
t321 = cos(qJ(5));
t257 = t281 * t321 - t286 * t318;
t251 = 0.1e1 / t257;
t252 = 0.1e1 / t257 ^ 2;
t392 = -0.2e1 * t275;
t391 = -0.2e1 * t335;
t246 = t335 * qJD(4) + t319 * t348 + t322 * t331;
t298 = t317 * t304;
t337 = t305 * qJD(2);
t261 = t283 * qJD(1) - t390 * t298 - t324 * t337;
t232 = t257 * qJD(5) + t246 * t318 - t261 * t321;
t256 = t281 * t318 + t286 * t321;
t250 = t256 ^ 2;
t244 = t250 * t252 + 0.1e1;
t381 = t252 * t256;
t367 = qJD(5) * t256;
t233 = t246 * t321 + t261 * t318 - t367;
t384 = t233 * t251 * t252;
t386 = (t232 * t381 - t250 * t384) / t244 ^ 2;
t378 = t288 * t394;
t385 = (t247 * t376 - t273 * t378) / t269 ^ 2;
t383 = t239 * t335;
t380 = t265 * t335;
t379 = t266 * t335;
t377 = t275 * t288;
t375 = t286 * t319;
t374 = t286 * t322;
t373 = t318 * t251;
t371 = t321 * t256;
t365 = -0.2e1 * t386;
t364 = 0.2e1 * t386;
t363 = -0.2e1 * t385;
t362 = t240 * t391;
t361 = t288 * t385;
t360 = t239 * t380;
t359 = t239 * t379;
t358 = t256 * t384;
t353 = 0.2e1 * t358;
t352 = t378 * t392;
t277 = -t319 * t354 - t322 * t338;
t347 = qJD(5) * t374 - t331;
t255 = -t277 * t321 + t283 * t318;
t254 = -t277 * t318 - t283 * t321;
t342 = t252 * t371 - t373;
t341 = -t277 * t288 + t292 * t376;
t336 = -t265 + (t266 * t377 + t265) * t267;
t333 = -qJD(4) * t375 + qJD(5) * t339 - t261 * t322;
t295 = t332 * qJD(2);
t272 = -t291 * qJD(4) + t296 * t322;
t263 = t286 * qJD(1) - t324 * t298 + t390 * t337;
t259 = t318 * t339 + t321 * t374;
t258 = t318 * t374 - t321 * t339;
t242 = 0.1e1 / t244;
t235 = 0.1e1 / t237;
t234 = t267 * t393;
t231 = t341 * t267;
t228 = t336 * t335;
t226 = (-t265 * t283 + t266 * t300) * t319 + t346 * t234;
t225 = t346 * t231 - t265 * t277 + t266 * t292;
t223 = t341 * t363 + (t292 * t352 - t248 * t288 + (t247 * t292 + t271 * t277 + t272 * t275) * t289) * t267;
t221 = t363 * t393 + (t340 * t368 + (t300 * t352 - t263 * t288 + (t247 * t300 + t271 * t283 + t275 * t295) * t289) * t319) * t267;
t1 = [t361 * t391 + (-t245 * t288 - t335 * t394) * t267, t221, 0, t223, 0, 0; t275 * t238 * t366 + (-t247 * t238 + (t224 * t275 + t228 * t245) * t239) * t235 - (-t228 * t239 * t366 + (-0.2e1 * t228 * t387 + (-t229 * t267 * t377 + t363) * t360 + (t361 * t392 - t229 + (t229 - t343) * t267) * t359 - t336 * t382) * t235) * t335 (-t226 * t383 - t238 * t375) * t366 + (-t226 * t382 + (-t261 * t319 + t286 * t368) * t238 + (t226 * t362 - t239 * t375) * t224 + (t300 * t368 - t221 * t275 - t234 * t247 + t295 * t319 + (-t234 * t291 - t283 * t319) * t229) * t359 + (-t283 * t368 - t221 * t291 - t234 * t271 - t263 * t319 + (t234 * t275 - t300 * t319) * t229) * t360) * t235, 0 (-t225 * t383 - t238 * t281) * t366 + (t225 * t224 * t362 + t246 * t238 + (-t281 * t224 - t225 * t245 + (-t223 * t275 - t231 * t247 + t272 + (-t231 * t291 - t277) * t229) * t379 + (-t223 * t291 - t231 * t271 - t248 + (t231 * t275 - t292) * t229) * t380) * t239) * t235, 0, 0; (-t251 * t254 + t255 * t381) * t364 + ((t255 * qJD(5) - t248 * t318 - t263 * t321) * t251 + t255 * t353 + (-t254 * t233 - (-t254 * qJD(5) - t248 * t321 + t263 * t318) * t256 - t255 * t232) * t252) * t242 (-t251 * t258 + t259 * t381) * t364 + (t259 * t353 + t347 * t251 * t321 + t333 * t373 + (t347 * t256 * t318 - t259 * t232 - t258 * t233 - t333 * t371) * t252) * t242, 0, -t342 * t335 * t365 + (t342 * t245 - ((-qJD(5) * t251 - 0.2e1 * t358) * t321 + (t232 * t321 + (t233 - t367) * t318) * t252) * t335) * t242, t365 + 0.2e1 * (t232 * t252 * t242 + (-t242 * t384 - t252 * t386) * t256) * t256, 0;];
JaD_rot  = t1;
