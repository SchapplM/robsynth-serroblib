% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR10_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR10_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:59:18
% EndTime: 2019-02-26 21:59:19
% DurationCPUTime: 1.54s
% Computational Cost: add. (9346->151), mult. (14155->297), div. (744->12), fcn. (17850->13), ass. (0->135)
t297 = cos(pkin(6));
t298 = sin(qJ(2));
t373 = sin(qJ(1));
t334 = t373 * t298;
t320 = t297 * t334;
t329 = qJD(2) * t373;
t299 = cos(qJ(2));
t300 = cos(qJ(1));
t349 = t300 * t299;
t296 = sin(pkin(6));
t351 = t296 * t300;
t378 = -qJD(1) * t320 - t298 * t329 + (qJD(2) * t297 + qJD(1)) * t349 - qJD(4) * t351;
t293 = pkin(12) + qJ(4);
t289 = sin(t293);
t290 = cos(t293);
t333 = t373 * t299;
t350 = t300 * t298;
t312 = -t297 * t350 - t333;
t262 = -t289 * t312 + t290 * t351;
t353 = t296 * t298;
t273 = t289 * t353 - t290 * t297;
t251 = atan2(-t262, t273);
t246 = sin(t251);
t247 = cos(t251);
t229 = -t246 * t262 + t247 * t273;
t227 = 0.1e1 / t229 ^ 2;
t278 = -t320 + t349;
t335 = t296 * t373;
t311 = -t278 * t289 + t290 * t335;
t261 = t311 ^ 2;
t223 = t227 * t261 + 0.1e1;
t310 = -t297 * t333 - t350;
t255 = t312 * qJD(1) + t310 * qJD(2);
t268 = t278 * t290 + t289 * t335;
t332 = qJD(1) * t351;
t233 = t268 * qJD(4) + t255 * t289 - t290 * t332;
t366 = t233 * t227;
t260 = t262 ^ 2;
t270 = 0.1e1 / t273 ^ 2;
t250 = t260 * t270 + 0.1e1;
t248 = 0.1e1 / t250;
t328 = t373 * qJD(1);
t319 = t296 * t328;
t346 = qJD(4) * t290;
t235 = t378 * t289 - t290 * t319 - t312 * t346;
t274 = t289 * t297 + t290 * t353;
t347 = qJD(2) * t299;
t331 = t296 * t347;
t258 = t274 * qJD(4) + t289 * t331;
t269 = 0.1e1 / t273;
t358 = t262 * t270;
t316 = -t235 * t269 + t258 * t358;
t217 = t316 * t248;
t317 = -t246 * t273 - t247 * t262;
t212 = t317 * t217 - t235 * t246 + t247 * t258;
t226 = 0.1e1 / t229;
t228 = t226 * t227;
t371 = t212 * t228;
t345 = 0.2e1 * (-t261 * t371 - t311 * t366) / t223 ^ 2;
t377 = t258 * t270;
t336 = t297 * t349;
t275 = -t334 + t336;
t352 = t296 * t299;
t313 = -t269 * t275 + t352 * t358;
t376 = t289 * t313;
t236 = (qJD(4) * t312 + t319) * t289 + t378 * t290;
t295 = qJ(5) + qJ(6);
t292 = cos(t295);
t291 = sin(t295);
t356 = t310 * t291;
t245 = t268 * t292 - t356;
t239 = 0.1e1 / t245;
t240 = 0.1e1 / t245 ^ 2;
t375 = -0.2e1 * t262;
t374 = -0.2e1 * t311;
t254 = -qJD(1) * t336 - t300 * t347 + (t297 * t329 + t328) * t298;
t294 = qJD(5) + qJD(6);
t323 = t268 * t294 + t254;
t234 = t311 * qJD(4) + t255 * t290 + t289 * t332;
t354 = t310 * t294;
t325 = t234 - t354;
t224 = t325 * t291 + t323 * t292;
t355 = t310 * t292;
t244 = t268 * t291 + t355;
t238 = t244 ^ 2;
t232 = t238 * t240 + 0.1e1;
t364 = t240 * t244;
t225 = -t323 * t291 + t325 * t292;
t368 = t225 * t239 * t240;
t370 = (t224 * t364 - t238 * t368) / t232 ^ 2;
t360 = t269 * t377;
t369 = (t235 * t358 - t260 * t360) / t250 ^ 2;
t367 = t227 * t311;
t365 = t239 * t291;
t363 = t244 * t292;
t362 = t246 * t311;
t361 = t247 * t311;
t359 = t262 * t269;
t357 = t310 * t289;
t348 = qJD(2) * t298;
t344 = -0.2e1 * t370;
t343 = 0.2e1 * t370;
t342 = -0.2e1 * t369;
t341 = t228 * t374;
t340 = t269 * t369;
t339 = t244 * t368;
t338 = t227 * t362;
t337 = t227 * t361;
t327 = 0.2e1 * t339;
t326 = t360 * t375;
t324 = t275 * t294 - t236;
t256 = t310 * qJD(1) + t312 * qJD(2);
t264 = -t289 * t351 - t290 * t312;
t322 = -t264 * t294 - t256;
t318 = -t290 * t354 + t255;
t315 = t240 * t363 - t365;
t314 = -t264 * t269 + t274 * t358;
t308 = -t246 + (t247 * t359 + t246) * t248;
t307 = -qJD(4) * t357 + t254 * t290 + t278 * t294;
t259 = -t273 * qJD(4) + t290 * t331;
t253 = t278 * t291 + t290 * t355;
t252 = -t278 * t292 + t290 * t356;
t243 = -t264 * t292 + t275 * t291;
t242 = -t264 * t291 - t275 * t292;
t230 = 0.1e1 / t232;
t221 = 0.1e1 / t223;
t220 = t248 * t376;
t218 = t314 * t248;
t216 = t308 * t311;
t214 = (-t246 * t275 + t247 * t352) * t289 + t317 * t220;
t213 = t317 * t218 - t246 * t264 + t247 * t274;
t211 = t314 * t342 + (t274 * t326 - t236 * t269 + (t235 * t274 + t258 * t264 + t259 * t262) * t270) * t248;
t209 = t342 * t376 + (t313 * t346 + (t326 * t352 - t256 * t269 + (t258 * t275 + (t235 * t299 - t262 * t348) * t296) * t270) * t289) * t248;
t208 = t344 + 0.2e1 * (t224 * t230 * t240 + (-t230 * t368 - t240 * t370) * t244) * t244;
t1 = [t340 * t374 + (-t233 * t269 - t311 * t377) * t248, t209, 0, t211, 0, 0; t262 * t226 * t345 + (-t235 * t226 + (t212 * t262 + t216 * t233) * t227) * t221 - (-t216 * t227 * t345 + (-0.2e1 * t216 * t371 + (-t217 * t248 * t359 + t342) * t338 + (t340 * t375 - t217 + (t217 - t316) * t248) * t337 - t308 * t366) * t221) * t311 (-t214 * t367 - t226 * t357) * t345 + (-t214 * t366 + (t254 * t289 + t310 * t346) * t226 + (t214 * t341 - t227 * t357) * t212 + (-t209 * t262 - t220 * t235 + (-t289 * t348 + t299 * t346) * t296 + (-t220 * t273 - t275 * t289) * t217) * t337 + (-t275 * t346 - t209 * t273 - t220 * t258 - t256 * t289 + (t220 * t262 - t289 * t352) * t217) * t338) * t221, 0 (-t213 * t367 - t226 * t268) * t345 + (t213 * t212 * t341 + t234 * t226 + (-t268 * t212 - t213 * t233 + (-t211 * t262 - t218 * t235 + t259 + (-t218 * t273 - t264) * t217) * t361 + (-t211 * t273 - t218 * t258 - t236 + (t218 * t262 - t274) * t217) * t362) * t227) * t221, 0, 0; (-t239 * t242 + t243 * t364) * t343 + ((t324 * t291 + t322 * t292) * t239 + t243 * t327 + (-t242 * t225 - (-t322 * t291 + t324 * t292) * t244 - t243 * t224) * t240) * t230 (-t239 * t252 + t253 * t364) * t343 + (t253 * t327 - t318 * t239 * t292 + t307 * t365 + (-t318 * t244 * t291 - t253 * t224 - t252 * t225 - t307 * t363) * t240) * t230, 0, -t315 * t311 * t344 + (t315 * t233 - ((-t239 * t294 - 0.2e1 * t339) * t292 + (t224 * t292 + (-t244 * t294 + t225) * t291) * t240) * t311) * t230, t208, t208;];
JaD_rot  = t1;
