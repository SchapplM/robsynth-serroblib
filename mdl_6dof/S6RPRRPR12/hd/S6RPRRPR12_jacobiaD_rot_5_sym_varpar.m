% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPR12_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR12_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:07:17
% EndTime: 2019-02-26 21:07:19
% DurationCPUTime: 2.95s
% Computational Cost: add. (11423->140), mult. (34278->290), div. (691->12), fcn. (44248->15), ass. (0->136)
t317 = cos(pkin(6));
t401 = sin(pkin(12));
t404 = cos(qJ(1));
t358 = t404 * t401;
t315 = cos(pkin(12));
t403 = sin(qJ(1));
t366 = t403 * t315;
t338 = t317 * t366 + t358;
t303 = t338 * qJD(1);
t313 = t403 * t401;
t367 = t404 * t315;
t307 = -t317 * t313 + t367;
t304 = t307 * qJD(1);
t316 = cos(pkin(7));
t319 = sin(qJ(3));
t321 = cos(qJ(3));
t314 = sin(pkin(6));
t402 = sin(pkin(7));
t365 = t314 * t402;
t355 = t403 * t365;
t347 = qJD(1) * t355;
t356 = t404 * t365;
t353 = -t317 * t367 + t313;
t411 = t353 * t316;
t331 = t411 + t356;
t339 = -t317 * t358 - t366;
t407 = -t339 * t319 + t331 * t321;
t263 = t407 * qJD(3) - (-t303 * t316 + t347) * t319 - t304 * t321;
t381 = t339 * t321;
t287 = t331 * t319 + t381;
t318 = sin(qJ(4));
t320 = cos(qJ(4));
t380 = t314 * t316;
t332 = t353 * t402 - t404 * t380;
t272 = t287 * t318 + t332 * t320;
t360 = t403 * t380;
t341 = qJD(1) * t360 + t303 * t402;
t249 = t272 * qJD(4) - t263 * t320 + t341 * t318;
t273 = t287 * t320 - t332 * t318;
t419 = t273 * qJD(4) + t263 * t318 + t341 * t320;
t416 = -t338 * t316 + t355;
t378 = qJD(3) * t319;
t379 = t316 * t321;
t261 = (t319 * t411 + t381) * qJD(3) - t303 * t379 - t304 * t319 + t321 * t347 + t356 * t378;
t364 = t402 * t317;
t298 = (t315 * t316 * t319 + t401 * t321) * t314 + t319 * t364;
t305 = -t315 * t365 + t317 * t316;
t280 = t298 * t320 + t305 * t318;
t297 = t321 * t364 + (t315 * t379 - t401 * t319) * t314;
t292 = t297 * qJD(3);
t264 = t280 * qJD(4) + t292 * t318;
t279 = t298 * t318 - t305 * t320;
t277 = 0.1e1 / t279 ^ 2;
t413 = t264 * t277;
t276 = 0.1e1 / t279;
t388 = t272 * t277;
t350 = t276 * t407 - t297 * t388;
t412 = t318 * t350;
t410 = t331 * qJD(1);
t409 = t416 * t321;
t255 = atan2(t272, t279);
t250 = sin(t255);
t251 = cos(t255);
t245 = t250 * t272 + t251 * t279;
t242 = 0.1e1 / t245;
t288 = t307 * t319 - t409;
t281 = 0.1e1 / t288;
t243 = 0.1e1 / t245 ^ 2;
t282 = 0.1e1 / t288 ^ 2;
t406 = 0.2e1 * t272;
t289 = t307 * t321 + t416 * t319;
t330 = t338 * t402 + t360;
t329 = t330 * t320;
t274 = t289 * t318 - t329;
t405 = 0.2e1 * t274;
t267 = t274 ^ 2;
t241 = t267 * t243 + 0.1e1;
t302 = t339 * qJD(1);
t260 = t409 * qJD(3) + t302 * t321 - t307 * t378 + t410 * t319;
t275 = t289 * t320 + t330 * t318;
t328 = t332 * qJD(1);
t246 = t275 * qJD(4) + t260 * t318 + t320 * t328;
t395 = t246 * t243;
t266 = t272 ^ 2;
t254 = t266 * t277 + 0.1e1;
t252 = 0.1e1 / t254;
t352 = -t264 * t388 + t276 * t419;
t234 = t352 * t252;
t357 = -t250 * t279 + t251 * t272;
t230 = t357 * t234 + t250 * t419 + t251 * t264;
t244 = t242 * t243;
t399 = t230 * t244;
t400 = (-t267 * t399 + t274 * t395) / t241 ^ 2;
t377 = qJD(4) * t318;
t247 = qJD(4) * t329 + t260 * t320 - t289 * t377 - t318 * t328;
t268 = t275 ^ 2;
t258 = t268 * t282 + 0.1e1;
t387 = t275 * t282;
t259 = t289 * qJD(3) + t302 * t319 - t410 * t321;
t283 = t281 * t282;
t391 = t259 * t283;
t398 = (t247 * t387 - t268 * t391) / t258 ^ 2;
t390 = t276 * t413;
t397 = (-t266 * t390 + t388 * t419) / t254 ^ 2;
t396 = t243 * t274;
t394 = t250 * t274;
t393 = t251 * t274;
t392 = t259 * t282;
t389 = t272 * t276;
t386 = t275 * t289;
t385 = t281 * t288;
t384 = t288 * t318;
t376 = qJD(4) * t320;
t375 = 0.2e1 * t400;
t374 = -0.2e1 * t397;
t373 = t244 * t405;
t372 = -0.2e1 * t275 * t407;
t371 = t281 * t398;
t370 = t276 * t397;
t369 = t243 * t394;
t368 = t243 * t393;
t363 = t390 * t406;
t351 = t273 * t276 - t280 * t388;
t343 = -t250 + (-t251 * t389 + t250) * t252;
t293 = t298 * qJD(3);
t265 = -t279 * qJD(4) + t292 * t320;
t256 = 0.1e1 / t258;
t239 = 0.1e1 / t241;
t238 = t252 * t412;
t236 = t351 * t252;
t233 = t343 * t274;
t232 = (t250 * t407 + t251 * t297) * t318 + t357 * t238;
t231 = t357 * t236 + t250 * t273 + t251 * t280;
t228 = t351 * t374 + (t280 * t363 - t249 * t276 + (-t264 * t273 - t265 * t272 - t280 * t419) * t277) * t252;
t227 = t374 * t412 + (t350 * t376 + (t297 * t363 - t261 * t276 + (-t264 * t407 + t272 * t293 - t297 * t419) * t277) * t318) * t252;
t1 = [t370 * t405 + (-t246 * t276 + t274 * t413) * t252, 0, t227, t228, 0, 0; -0.2e1 * t272 * t242 * t400 + (t419 * t242 + (-t272 * t230 - t233 * t246) * t243) * t239 + (t233 * t243 * t375 + (0.2e1 * t233 * t399 - (t234 * t252 * t389 + t374) * t369 - (t370 * t406 - t234 + (t234 - t352) * t252) * t368 - t343 * t395) * t239) * t274, 0 (t232 * t396 + t242 * t384) * t375 + (-t232 * t395 + (-t259 * t318 - t288 * t376) * t242 + (t232 * t373 + t243 * t384) * t230 - (t297 * t376 + t227 * t272 + t238 * t419 - t293 * t318 + (-t238 * t279 + t318 * t407) * t234) * t368 - (t407 * t376 - t227 * t279 - t238 * t264 - t261 * t318 + (-t238 * t272 - t297 * t318) * t234) * t369) * t239 (t231 * t396 - t242 * t275) * t375 + (t231 * t230 * t373 + t247 * t242 + (-t275 * t230 - t231 * t246 - (t228 * t272 + t236 * t419 + t265 + (-t236 * t279 + t273) * t234) * t393 - (-t228 * t279 - t236 * t264 - t249 + (-t236 * t272 - t280) * t234) * t394) * t243) * t239, 0, 0; t282 * t372 * t398 - 0.2e1 * t273 * t371 + (t407 * t247 * t282 - t249 * t281 - t261 * t387 - t273 * t392 + t372 * t391) * t256, 0, 0.2e1 * (t282 * t386 + t320 * t385) * t398 + (t377 * t385 + (-t247 * t289 - t260 * t275) * t282 + (0.2e1 * t283 * t386 + (t282 * t288 - t281) * t320) * t259) * t256, t371 * t405 + (-t246 * t281 + t274 * t392) * t256, 0, 0;];
JaD_rot  = t1;
