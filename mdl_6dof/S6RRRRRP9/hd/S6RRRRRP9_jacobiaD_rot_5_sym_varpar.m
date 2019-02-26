% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRP9_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP9_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:44:21
% EndTime: 2019-02-26 22:44:22
% DurationCPUTime: 1.40s
% Computational Cost: add. (5342->150), mult. (14155->298), div. (744->12), fcn. (17850->13), ass. (0->132)
t283 = cos(pkin(6));
t285 = sin(qJ(2));
t359 = sin(qJ(1));
t322 = t359 * t285;
t308 = t283 * t322;
t317 = qJD(2) * t359;
t287 = cos(qJ(2));
t288 = cos(qJ(1));
t337 = t288 * t287;
t282 = sin(pkin(6));
t339 = t282 * t288;
t364 = -qJD(1) * t308 - t285 * t317 + (qJD(2) * t283 + qJD(1)) * t337 - qJD(3) * t339;
t284 = sin(qJ(3));
t286 = cos(qJ(3));
t321 = t359 * t287;
t338 = t288 * t285;
t300 = -t283 * t338 - t321;
t251 = -t284 * t300 + t286 * t339;
t341 = t282 * t285;
t262 = -t283 * t286 + t284 * t341;
t240 = atan2(-t251, t262);
t235 = sin(t240);
t236 = cos(t240);
t218 = -t235 * t251 + t236 * t262;
t216 = 0.1e1 / t218 ^ 2;
t267 = -t308 + t337;
t323 = t282 * t359;
t299 = -t267 * t284 + t286 * t323;
t248 = t299 ^ 2;
t214 = t216 * t248 + 0.1e1;
t298 = -t283 * t321 - t338;
t244 = t300 * qJD(1) + t298 * qJD(2);
t257 = t267 * t286 + t284 * t323;
t320 = qJD(1) * t339;
t222 = t257 * qJD(3) + t244 * t284 - t286 * t320;
t353 = t216 * t299;
t247 = t251 ^ 2;
t260 = 0.1e1 / t262 ^ 2;
t239 = t247 * t260 + 0.1e1;
t237 = 0.1e1 / t239;
t316 = t359 * qJD(1);
t307 = t282 * t316;
t334 = qJD(3) * t286;
t224 = t364 * t284 - t286 * t307 - t300 * t334;
t263 = t283 * t284 + t286 * t341;
t335 = qJD(2) * t287;
t319 = t282 * t335;
t249 = t263 * qJD(3) + t284 * t319;
t259 = 0.1e1 / t262;
t344 = t251 * t260;
t304 = -t224 * t259 + t249 * t344;
t206 = t304 * t237;
t305 = -t235 * t262 - t236 * t251;
t201 = t305 * t206 - t235 * t224 + t236 * t249;
t215 = 0.1e1 / t218;
t217 = t215 * t216;
t357 = t201 * t217;
t333 = 0.2e1 * (-t222 * t353 - t248 * t357) / t214 ^ 2;
t363 = t249 * t260;
t324 = t283 * t337;
t264 = -t322 + t324;
t340 = t282 * t287;
t301 = -t259 * t264 + t340 * t344;
t362 = t284 * t301;
t225 = (qJD(3) * t300 + t307) * t284 + t364 * t286;
t281 = qJ(4) + qJ(5);
t278 = sin(t281);
t279 = cos(t281);
t234 = t257 * t279 - t278 * t298;
t228 = 0.1e1 / t234;
t229 = 0.1e1 / t234 ^ 2;
t361 = -0.2e1 * t251;
t360 = -0.2e1 * t299;
t346 = t259 * t363;
t356 = (t224 * t344 - t247 * t346) / t239 ^ 2;
t243 = -qJD(1) * t324 - t288 * t335 + (t283 * t317 + t316) * t285;
t280 = qJD(4) + qJD(5);
t311 = t257 * t280 + t243;
t223 = t299 * qJD(3) + t244 * t286 + t284 * t320;
t313 = -t280 * t298 + t223;
t210 = -t311 * t278 + t313 * t279;
t355 = t210 * t228 * t229;
t354 = t216 * t222;
t209 = t313 * t278 + t311 * t279;
t233 = t257 * t278 + t279 * t298;
t227 = t233 ^ 2;
t221 = t227 * t229 + 0.1e1;
t350 = t229 * t233;
t352 = 0.1e1 / t221 ^ 2 * (t209 * t350 - t227 * t355);
t351 = t228 * t278;
t349 = t233 * t279;
t348 = t235 * t299;
t347 = t236 * t299;
t345 = t251 * t259;
t343 = t298 * t284;
t342 = t298 * t286;
t336 = qJD(2) * t285;
t332 = -0.2e1 * t356;
t331 = t217 * t360;
t330 = -0.2e1 * t352;
t329 = 0.2e1 * t352;
t328 = t259 * t356;
t327 = t233 * t355;
t326 = t216 * t348;
t325 = t216 * t347;
t315 = 0.2e1 * t327;
t314 = t346 * t361;
t312 = t264 * t280 - t225;
t245 = t298 * qJD(1) + t300 * qJD(2);
t253 = -t284 * t339 - t286 * t300;
t310 = -t253 * t280 - t245;
t306 = -t280 * t342 + t244;
t303 = t229 * t349 - t351;
t302 = -t253 * t259 + t263 * t344;
t296 = -t235 + (t236 * t345 + t235) * t237;
t295 = -qJD(3) * t343 + t243 * t286 + t267 * t280;
t250 = -t262 * qJD(3) + t286 * t319;
t242 = t267 * t278 + t279 * t342;
t241 = -t267 * t279 + t278 * t342;
t232 = -t253 * t279 + t264 * t278;
t231 = -t253 * t278 - t264 * t279;
t219 = 0.1e1 / t221;
t212 = 0.1e1 / t214;
t211 = t237 * t362;
t208 = t302 * t237;
t205 = t296 * t299;
t203 = (-t235 * t264 + t236 * t340) * t284 + t305 * t211;
t202 = t305 * t208 - t235 * t253 + t236 * t263;
t200 = t302 * t332 + (t263 * t314 - t225 * t259 + (t224 * t263 + t249 * t253 + t250 * t251) * t260) * t237;
t198 = t332 * t362 + (t301 * t334 + (t314 * t340 - t245 * t259 + (t249 * t264 + (t224 * t287 - t251 * t336) * t282) * t260) * t284) * t237;
t197 = t330 + 0.2e1 * (t209 * t229 * t219 + (-t219 * t355 - t229 * t352) * t233) * t233;
t1 = [t328 * t360 + (-t222 * t259 - t299 * t363) * t237, t198, t200, 0, 0, 0; t251 * t215 * t333 + (-t224 * t215 + (t201 * t251 + t205 * t222) * t216) * t212 - (-t205 * t216 * t333 + (-0.2e1 * t205 * t357 + (-t206 * t237 * t345 + t332) * t326 + (t328 * t361 - t206 + (t206 - t304) * t237) * t325 - t296 * t354) * t212) * t299 (-t203 * t353 - t215 * t343) * t333 + (-t203 * t354 + (t243 * t284 + t298 * t334) * t215 + (t203 * t331 - t216 * t343) * t201 + (-t198 * t251 - t211 * t224 + (-t284 * t336 + t287 * t334) * t282 + (-t211 * t262 - t264 * t284) * t206) * t325 + (-t264 * t334 - t198 * t262 - t211 * t249 - t245 * t284 + (t211 * t251 - t284 * t340) * t206) * t326) * t212 (-t202 * t353 - t215 * t257) * t333 + (t202 * t201 * t331 + t223 * t215 + (-t257 * t201 - t202 * t222 + (-t200 * t251 - t208 * t224 + t250 + (-t208 * t262 - t253) * t206) * t347 + (-t200 * t262 - t208 * t249 - t225 + (t208 * t251 - t263) * t206) * t348) * t216) * t212, 0, 0, 0; (-t228 * t231 + t232 * t350) * t329 + ((t312 * t278 + t310 * t279) * t228 + t232 * t315 + (-t231 * t210 - (-t310 * t278 + t312 * t279) * t233 - t232 * t209) * t229) * t219 (-t228 * t241 + t242 * t350) * t329 + (t242 * t315 - t306 * t228 * t279 + t295 * t351 + (-t306 * t233 * t278 - t242 * t209 - t241 * t210 - t295 * t349) * t229) * t219, -t303 * t299 * t330 + (t303 * t222 - ((-t228 * t280 - 0.2e1 * t327) * t279 + (t209 * t279 + (-t233 * t280 + t210) * t278) * t229) * t299) * t219, t197, t197, 0;];
JaD_rot  = t1;
