% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRP2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:51:05
% EndTime: 2019-02-26 19:51:07
% DurationCPUTime: 2.30s
% Computational Cost: add. (13772->162), mult. (39160->316), div. (784->12), fcn. (51130->15), ass. (0->128)
t373 = sin(pkin(11));
t377 = cos(pkin(6));
t337 = t377 * t373;
t375 = cos(pkin(11));
t338 = t377 * t375;
t378 = sin(qJ(2));
t379 = cos(qJ(2));
t302 = t379 * t337 + t378 * t338;
t305 = t378 * t373 - t379 * t375;
t374 = sin(pkin(10));
t376 = cos(pkin(10));
t288 = t376 * t302 - t374 * t305;
t313 = sin(qJ(4));
t315 = cos(qJ(4));
t311 = sin(pkin(6));
t350 = t311 * t376;
t275 = -t288 * t313 - t315 * t350;
t382 = t378 * t337 - t379 * t338;
t299 = t382 * qJD(2);
t325 = t379 * t373 + t378 * t375;
t304 = t325 * qJD(2);
t281 = -t376 * t299 - t374 * t304;
t251 = t275 * qJD(4) + t281 * t315;
t314 = cos(qJ(5));
t329 = -t288 * t315 + t313 * t350;
t312 = sin(qJ(5));
t360 = -t374 * t325 - t376 * t382;
t383 = t360 * t312;
t261 = -t314 * t329 - t383;
t381 = qJD(2) * t288;
t233 = t261 * qJD(5) + t251 * t312 - t314 * t381;
t345 = t360 * t314;
t259 = -t312 * t329 + t345;
t254 = t259 ^ 2;
t300 = t305 * t311;
t301 = t325 * t311;
t332 = -t301 * t315 - t377 * t313;
t272 = -t300 * t314 - t312 * t332;
t270 = 0.1e1 / t272 ^ 2;
t242 = t254 * t270 + 0.1e1;
t240 = 0.1e1 / t242;
t292 = -t301 * t313 + t377 * t315;
t298 = qJD(2) * t300;
t268 = t292 * qJD(4) - t298 * t315;
t273 = t300 * t312 - t314 * t332;
t297 = qJD(2) * t301;
t247 = t273 * qJD(5) + t268 * t312 - t297 * t314;
t269 = 0.1e1 / t272;
t364 = t259 * t270;
t220 = (-t233 * t269 + t247 * t364) * t240;
t243 = atan2(-t259, t272);
t238 = sin(t243);
t239 = cos(t243);
t336 = -t238 * t272 - t239 * t259;
t216 = t336 * t220 - t233 * t238 + t239 * t247;
t232 = -t238 * t259 + t239 * t272;
t229 = 0.1e1 / t232;
t230 = 0.1e1 / t232 ^ 2;
t386 = t216 * t229 * t230;
t290 = -t376 * t325 + t374 * t382;
t326 = -t374 * t302 - t376 * t305;
t349 = t311 * t374;
t328 = -t313 * t349 - t315 * t326;
t262 = t290 * t314 - t312 * t328;
t385 = -0.2e1 * t262;
t344 = 0.2e1 * t262 * t386;
t333 = -t269 * t275 + t292 * t364;
t384 = t312 * t333;
t367 = t247 * t269 * t270;
t380 = -0.2e1 * (t233 * t364 - t254 * t367) / t242 ^ 2;
t263 = -t290 * t312 - t314 * t328;
t256 = 0.1e1 / t263;
t257 = 0.1e1 / t263 ^ 2;
t277 = -t313 * t326 + t315 * t349;
t274 = t277 ^ 2;
t366 = t257 * t274;
t246 = 0.1e1 + t366;
t327 = t374 * t299 - t376 * t304;
t252 = t328 * qJD(4) - t313 * t327;
t253 = t277 * qJD(4) + t315 * t327;
t282 = t326 * qJD(2);
t236 = -t262 * qJD(5) + t253 * t314 + t282 * t312;
t370 = t236 * t256 * t257;
t352 = t274 * t370;
t365 = t257 * t277;
t372 = (t252 * t365 - t352) / t246 ^ 2;
t371 = t230 * t262;
t369 = t238 * t262;
t368 = t239 * t262;
t363 = t277 * t312;
t362 = t312 * t315;
t361 = t314 * t315;
t359 = qJD(4) * t313;
t358 = qJD(5) * t312;
t357 = qJD(5) * t314;
t356 = qJD(5) * t315;
t255 = t262 ^ 2;
t228 = t230 * t255 + 0.1e1;
t235 = t263 * qJD(5) + t253 * t312 - t282 * t314;
t355 = 0.2e1 * (t235 * t371 - t255 * t386) / t228 ^ 2;
t354 = 0.2e1 * t372;
t351 = t277 * t370;
t346 = t315 * t360;
t343 = -0.2e1 * t259 * t367;
t335 = -t261 * t269 + t273 * t364;
t264 = -t288 * t314 + t312 * t346;
t280 = -t300 * t362 - t301 * t314;
t334 = -t264 * t269 + t280 * t364;
t267 = t332 * qJD(4) + t298 * t313;
t266 = t290 * t361 + t312 * t326;
t265 = t290 * t362 - t314 * t326;
t250 = t329 * qJD(4) - t281 * t313;
t249 = (-t300 * t356 + t298) * t314 + (qJD(5) * t301 - t297 * t315 + t300 * t359) * t312;
t248 = -t272 * qJD(5) + t268 * t314 + t297 * t312;
t244 = 0.1e1 / t246;
t237 = -t281 * t314 + t288 * t358 + t346 * t357 - t359 * t383 - t362 * t381;
t234 = -qJD(5) * t345 + t251 * t314 + t312 * t381 + t329 * t358;
t226 = 0.1e1 / t228;
t225 = t240 * t384;
t224 = t334 * t240;
t222 = t335 * t240;
t219 = (-t238 * t275 + t239 * t292) * t312 + t336 * t225;
t218 = t336 * t224 - t238 * t264 + t239 * t280;
t217 = t336 * t222 - t238 * t261 + t239 * t273;
t215 = t334 * t380 + (t280 * t343 - t237 * t269 + (t233 * t280 + t247 * t264 + t249 * t259) * t270) * t240;
t213 = t335 * t380 + (t273 * t343 - t234 * t269 + (t233 * t273 + t247 * t261 + t248 * t259) * t270) * t240;
t212 = t380 * t384 + (t333 * t357 + (t292 * t343 - t250 * t269 + (t233 * t292 + t247 * t275 + t259 * t267) * t270) * t312) * t240;
t1 = [0, t215, 0, t212, t213, 0; 0 (t218 * t371 - t229 * t265) * t355 + (t218 * t344 + (-t265 * t216 - t218 * t235 - (-t215 * t259 - t224 * t233 + t249 + (-t224 * t272 - t264) * t220) * t368 - (-t215 * t272 - t224 * t247 - t237 + (t224 * t259 - t280) * t220) * t369) * t230 + ((t290 * t356 - t327) * t314 + (qJD(5) * t326 - t282 * t315 - t290 * t359) * t312) * t229) * t226, 0 (t219 * t371 - t229 * t363) * t355 + ((t252 * t312 + t277 * t357) * t229 + t219 * t344 + (-t219 * t235 - t363 * t216 - (t292 * t357 - t212 * t259 - t225 * t233 + t267 * t312 + (-t225 * t272 - t275 * t312) * t220) * t368 - (-t275 * t357 - t212 * t272 - t225 * t247 - t250 * t312 + (t225 * t259 - t292 * t312) * t220) * t369) * t230) * t226 (t217 * t371 - t229 * t263) * t355 + (t217 * t344 + t236 * t229 + (-t263 * t216 - t217 * t235 - (-t213 * t259 - t222 * t233 + t248 + (-t222 * t272 - t261) * t220) * t368 - (-t213 * t272 - t222 * t247 - t234 + (t222 * t259 - t273) * t220) * t369) * t230) * t226, 0; 0 (t256 * t290 * t313 + t266 * t365) * t354 + (0.2e1 * t266 * t351 + (-t290 * qJD(4) * t315 + t282 * t313) * t256 + (-(-t282 * t361 + t312 * t327 + t326 * t357) * t277 - t266 * t252 + (t313 * t236 - (-t312 * t356 - t314 * t359) * t277) * t290) * t257) * t244, 0 (-t256 * t328 + t314 * t366) * t354 + (0.2e1 * t314 * t352 - t253 * t256 + (-0.2e1 * t252 * t277 * t314 - t236 * t328 + t274 * t358) * t257) * t244, t365 * t372 * t385 + (t351 * t385 + (t235 * t277 + t252 * t262) * t257) * t244, 0;];
JaD_rot  = t1;
