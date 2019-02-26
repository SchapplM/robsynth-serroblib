% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRRP4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:16:54
% EndTime: 2019-02-26 20:16:56
% DurationCPUTime: 2.03s
% Computational Cost: add. (15558->158), mult. (28006->308), div. (1030->12), fcn. (36105->13), ass. (0->128)
t290 = qJ(4) + qJ(5);
t287 = sin(t290);
t288 = cos(t290);
t293 = sin(qJ(2));
t295 = cos(qJ(2));
t355 = cos(pkin(11));
t356 = cos(pkin(6));
t317 = t356 * t355;
t354 = sin(pkin(11));
t304 = -t293 * t317 - t354 * t295;
t302 = t304 * qJD(2);
t292 = sin(qJ(3));
t294 = cos(qJ(3));
t291 = sin(pkin(6));
t325 = t291 * t355;
t306 = t292 * t325 + t294 * t304;
t289 = qJD(4) + qJD(5);
t339 = t288 * t289;
t266 = t292 * t304 - t294 * t325;
t315 = t354 * t293 - t295 * t317;
t277 = t315 * qJD(2);
t358 = t266 * qJD(3) - t277 * t294 + t315 * t289;
t226 = t358 * t287 + t288 * t302 - t306 * t339;
t248 = -t287 * t306 - t315 * t288;
t243 = t248 ^ 2;
t337 = t291 * t293;
t284 = t356 * t292 + t294 * t337;
t336 = t291 * t295;
t263 = t284 * t287 + t288 * t336;
t261 = 0.1e1 / t263 ^ 2;
t235 = t243 * t261 + 0.1e1;
t233 = 0.1e1 / t235;
t334 = qJD(2) * t291;
t309 = -t284 * t289 + t293 * t334;
t283 = -t292 * t337 + t356 * t294;
t326 = t295 * t334;
t318 = t283 * qJD(3) - t289 * t336 + t294 * t326;
t240 = t318 * t287 - t309 * t288;
t260 = 0.1e1 / t263;
t345 = t248 * t261;
t213 = (-t226 * t260 + t240 * t345) * t233;
t236 = atan2(-t248, t263);
t231 = sin(t236);
t232 = cos(t236);
t314 = -t231 * t263 - t232 * t248;
t208 = t314 * t213 - t231 * t226 + t232 * t240;
t225 = -t231 * t248 + t232 * t263;
t222 = 0.1e1 / t225;
t223 = 0.1e1 / t225 ^ 2;
t361 = t208 * t222 * t223;
t316 = t356 * t354;
t303 = -t355 * t293 - t295 * t316;
t282 = -t293 * t316 + t355 * t295;
t324 = t291 * t354;
t305 = -t282 * t294 - t292 * t324;
t251 = -t287 * t305 + t288 * t303;
t360 = -0.2e1 * t251;
t321 = 0.2e1 * t251 * t361;
t310 = -t260 * t266 + t283 * t345;
t359 = t287 * t310;
t347 = t240 * t260 * t261;
t357 = -0.2e1 * (t226 * t345 - t243 * t347) / t235 ^ 2;
t252 = -t287 * t303 - t288 * t305;
t245 = 0.1e1 / t252;
t246 = 0.1e1 / t252 ^ 2;
t268 = -t282 * t292 + t294 * t324;
t265 = t268 ^ 2;
t343 = t265 * t246;
t239 = 0.1e1 + t343;
t279 = t282 * qJD(2);
t319 = -t289 * t305 - t279;
t278 = t303 * qJD(2);
t256 = t268 * qJD(3) + t278 * t294;
t320 = -t289 * t303 + t256;
t229 = -t319 * t287 + t320 * t288;
t350 = t229 * t245 * t246;
t328 = t265 * t350;
t255 = t305 * qJD(3) - t278 * t292;
t344 = t255 * t268;
t353 = (t246 * t344 - t328) / t239 ^ 2;
t352 = t223 * t251;
t228 = t320 * t287 + t319 * t288;
t351 = t228 * t223;
t349 = t231 * t251;
t348 = t232 * t251;
t346 = t246 * t268;
t342 = t268 * t287;
t341 = t287 * t289;
t340 = t287 * t294;
t338 = t288 * t294;
t335 = t294 * t289;
t333 = qJD(2) * t294;
t332 = qJD(3) * t292;
t244 = t251 ^ 2;
t221 = t244 * t223 + 0.1e1;
t331 = 0.2e1 * (-t244 * t361 + t251 * t351) / t221 ^ 2;
t330 = 0.2e1 * t353;
t327 = t268 * t350;
t322 = -0.2e1 * t248 * t347;
t250 = t315 * t287 - t288 * t306;
t264 = t284 * t288 - t287 * t336;
t312 = -t250 * t260 + t264 * t345;
t308 = t294 * t315;
t257 = -t287 * t308 + t288 * t304;
t272 = (-t288 * t293 + t295 * t340) * t291;
t311 = -t257 * t260 + t272 * t345;
t270 = -t284 * qJD(3) - t292 * t326;
t259 = t282 * t287 + t303 * t338;
t258 = -t282 * t288 + t303 * t340;
t253 = t306 * qJD(3) + t277 * t292;
t242 = ((-qJD(2) + t335) * t295 * t288 + (-t295 * t332 + (t289 - t333) * t293) * t287) * t291;
t241 = t309 * t287 + t318 * t288;
t237 = 0.1e1 / t239;
t230 = t277 * t288 - t304 * t341 - t308 * t339 + (t304 * t333 + t315 * t332) * t287;
t227 = -t287 * t302 + t358 * t288 + t306 * t341;
t219 = 0.1e1 / t221;
t218 = t233 * t359;
t216 = t311 * t233;
t215 = t312 * t233;
t212 = (-t231 * t266 + t232 * t283) * t287 + t314 * t218;
t211 = t314 * t216 - t231 * t257 + t232 * t272;
t210 = t314 * t215 - t231 * t250 + t232 * t264;
t209 = t346 * t353 * t360 + (t327 * t360 + (t228 * t268 + t251 * t255) * t246) * t237;
t207 = t311 * t357 + (t272 * t322 - t230 * t260 + (t226 * t272 + t240 * t257 + t242 * t248) * t261) * t233;
t205 = t312 * t357 + (t264 * t322 - t227 * t260 + (t226 * t264 + t240 * t250 + t241 * t248) * t261) * t233;
t204 = t357 * t359 + (t310 * t339 + (t283 * t322 - t253 * t260 + (t226 * t283 + t240 * t266 + t248 * t270) * t261) * t287) * t233;
t203 = (t210 * t352 - t222 * t252) * t331 + (t210 * t321 + t229 * t222 + (-t252 * t208 - t210 * t228 - (-t205 * t248 - t215 * t226 + t241 + (-t215 * t263 - t250) * t213) * t348 - (-t205 * t263 - t215 * t240 - t227 + (t215 * t248 - t264) * t213) * t349) * t223) * t219;
t1 = [0, t207, t204, t205, t205, 0; 0 (t211 * t352 - t222 * t258) * t331 + (t211 * t321 + (-t258 * t208 - t211 * t228 - (-t207 * t248 - t216 * t226 + t242 + (-t216 * t263 - t257) * t213) * t348 - (-t207 * t263 - t216 * t240 - t230 + (t216 * t248 - t272) * t213) * t349) * t223 + ((t303 * t335 - t278) * t288 + (-t279 * t294 + t282 * t289 - t303 * t332) * t287) * t222) * t219 (t212 * t352 - t222 * t342) * t331 + ((t255 * t287 + t268 * t339) * t222 + (-t351 + t321) * t212 + (-t342 * t208 - (t283 * t339 - t204 * t248 - t218 * t226 + t270 * t287 + (-t218 * t263 - t266 * t287) * t213) * t348 - (-t266 * t339 - t204 * t263 - t218 * t240 - t253 * t287 + (t218 * t248 - t283 * t287) * t213) * t349) * t223) * t219, t203, t203, 0; 0 (t245 * t292 * t303 + t259 * t346) * t330 + (0.2e1 * t259 * t327 + (-qJD(3) * t294 * t303 + t279 * t292) * t245 + (-(t278 * t287 - t279 * t338 + t282 * t339) * t268 - t259 * t255 - (-t292 * t229 - (t287 * t335 + t288 * t332) * t268) * t303) * t246) * t237 (-t245 * t305 + t288 * t343) * t330 + (0.2e1 * t288 * t328 - t245 * t256 + (-t229 * t305 + t265 * t341 - 0.2e1 * t288 * t344) * t246) * t237, t209, t209, 0;];
JaD_rot  = t1;
