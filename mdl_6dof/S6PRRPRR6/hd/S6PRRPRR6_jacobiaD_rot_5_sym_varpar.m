% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR6
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRR6_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_jacobiaD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:05
% EndTime: 2019-02-26 20:07:07
% DurationCPUTime: 1.16s
% Computational Cost: add. (5650->131), mult. (16984->252), div. (538->12), fcn. (21480->15), ass. (0->123)
t297 = sin(qJ(2));
t299 = cos(qJ(2));
t363 = cos(pkin(12));
t364 = cos(pkin(6));
t329 = t364 * t363;
t362 = sin(pkin(12));
t315 = -t362 * t297 + t299 * t329;
t279 = t315 * qJD(2);
t293 = sin(pkin(7));
t294 = sin(pkin(6));
t332 = t363 * t294 * t293;
t367 = -qJD(3) * t332 + t279;
t296 = sin(qJ(3));
t314 = -t297 * t329 - t299 * t362;
t309 = t314 * qJD(2);
t298 = cos(qJ(3));
t310 = t315 * t298;
t366 = qJD(3) * t310 + t296 * t309;
t295 = cos(pkin(7));
t311 = t315 * t296;
t306 = t295 * t311 - t298 * t314;
t347 = t295 * t298;
t233 = qJD(3) * t306 + t367 * t296 - t309 * t347;
t352 = t314 * t296;
t257 = -t295 * t310 + t298 * t332 - t352;
t255 = t257 ^ 2;
t343 = t298 * t299;
t346 = t296 * t297;
t320 = t295 * t343 - t346;
t336 = t364 * t293;
t270 = -t294 * t320 - t298 * t336;
t268 = 0.1e1 / t270 ^ 2;
t249 = t255 * t268 + 0.1e1;
t353 = t257 * t268;
t344 = t297 * t298;
t345 = t296 * t299;
t318 = t295 * t345 + t344;
t319 = t295 * t344 + t345;
t330 = qJD(3) * t336;
t253 = t296 * t330 + (qJD(2) * t319 + qJD(3) * t318) * t294;
t267 = 0.1e1 / t270;
t354 = t253 * t267 * t268;
t365 = -0.2e1 * (t233 * t353 - t255 * t354) / t249 ^ 2;
t250 = atan2(-t257, t270);
t245 = sin(t250);
t246 = cos(t250);
t227 = -t245 * t257 + t246 * t270;
t224 = 0.1e1 / t227;
t328 = t364 * t362;
t313 = t297 * t328 - t299 * t363;
t312 = t297 * t363 + t299 * t328;
t337 = t294 * t362;
t331 = t293 * t337;
t316 = -t295 * t312 + t331;
t261 = t296 * t316 - t298 * t313;
t272 = t293 * t312 + t295 * t337;
t292 = pkin(13) + qJ(5);
t290 = sin(t292);
t291 = cos(t292);
t242 = t261 * t291 + t272 * t290;
t238 = 0.1e1 / t242;
t225 = 0.1e1 / t227 ^ 2;
t239 = 0.1e1 / t242 ^ 2;
t247 = 0.1e1 / t249;
t217 = (-t233 * t267 + t253 * t353) * t247;
t327 = -t245 * t270 - t246 * t257;
t213 = t217 * t327 - t245 * t233 + t246 * t253;
t361 = t213 * t224 * t225;
t280 = t312 * qJD(2);
t281 = t313 * qJD(2);
t348 = t295 * t296;
t351 = t313 * t296;
t236 = t281 * t348 - t280 * t298 + (t298 * t316 + t351) * qJD(3);
t349 = t291 * t293;
t228 = qJD(5) * t242 + t236 * t290 + t281 * t349;
t241 = t261 * t290 - t272 * t291;
t237 = t241 ^ 2;
t232 = t237 * t239 + 0.1e1;
t357 = t239 * t241;
t342 = qJD(5) * t241;
t350 = t290 * t293;
t229 = t236 * t291 - t281 * t350 - t342;
t358 = t229 * t238 * t239;
t360 = (t228 * t357 - t237 * t358) / t232 ^ 2;
t260 = -t298 * t331 + t312 * t347 - t351;
t359 = t225 * t260;
t356 = t245 * t260;
t355 = t246 * t260;
t256 = t260 ^ 2;
t223 = t256 * t225 + 0.1e1;
t235 = qJD(3) * t261 - t280 * t296 - t281 * t347;
t341 = 0.2e1 * (t235 * t359 - t256 * t361) / t223 ^ 2;
t340 = -0.2e1 * t360;
t339 = t241 * t358;
t338 = qJD(3) * t352;
t334 = -0.2e1 * t257 * t354;
t333 = 0.2e1 * t260 * t361;
t324 = -t290 * t238 + t291 * t357;
t259 = -t296 * t332 + t306;
t271 = t294 * t318 + t296 * t336;
t323 = -t259 * t267 + t271 * t353;
t263 = -t314 * t347 + t311;
t278 = t319 * t294;
t322 = -t263 * t267 + t278 * t353;
t265 = -t298 * t312 + t313 * t348;
t321 = -t265 * t290 - t313 * t349;
t252 = t265 * t291 - t313 * t350;
t264 = -t296 * t312 - t313 * t347;
t317 = -t295 * t346 + t343;
t262 = (qJD(2) * t320 + qJD(3) * t317) * t294;
t254 = t298 * t330 + (qJD(2) * t317 + qJD(3) * t320) * t294;
t244 = -qJD(3) * t264 + t280 * t348 + t281 * t298;
t243 = t279 * t347 + t295 * t338 + t366;
t234 = t366 * t295 + t367 * t298 + t338;
t230 = 0.1e1 / t232;
t221 = 0.1e1 / t223;
t219 = t322 * t247;
t218 = t323 * t247;
t215 = t219 * t327 - t245 * t263 + t246 * t278;
t214 = t218 * t327 - t245 * t259 + t246 * t271;
t212 = t322 * t365 + (t278 * t334 - t243 * t267 + (t233 * t278 + t253 * t263 + t257 * t262) * t268) * t247;
t211 = t323 * t365 + (t271 * t334 - t234 * t267 + (t233 * t271 + t253 * t259 + t254 * t257) * t268) * t247;
t1 = [0, t212, t211, 0, 0, 0; 0 (t215 * t359 - t224 * t264) * t341 + ((t265 * qJD(3) - t280 * t347 + t281 * t296) * t224 + t215 * t333 + (-t264 * t213 - t215 * t235 - (-t212 * t257 - t219 * t233 + t262 + (-t219 * t270 - t263) * t217) * t355 - (-t212 * t270 - t219 * t253 - t243 + (t219 * t257 - t278) * t217) * t356) * t225) * t221 (t214 * t359 - t224 * t261) * t341 + (t214 * t333 + t236 * t224 + (-t261 * t213 - t214 * t235 - (-t211 * t257 - t218 * t233 + t254 + (-t218 * t270 - t259) * t217) * t355 - (-t211 * t270 - t218 * t253 - t234 + (t218 * t257 - t271) * t217) * t356) * t225) * t221, 0, 0, 0; 0, 0.2e1 * (t238 * t321 + t252 * t357) * t360 + ((qJD(5) * t252 + t244 * t290 + t280 * t349) * t238 + 0.2e1 * t252 * t339 + (t321 * t229 - (qJD(5) * t321 + t244 * t291 - t280 * t350) * t241 - t252 * t228) * t239) * t230, t324 * t260 * t340 + (t324 * t235 + ((-qJD(5) * t238 - 0.2e1 * t339) * t291 + (t228 * t291 + (t229 - t342) * t290) * t239) * t260) * t230, 0, t340 + 0.2e1 * (t228 * t239 * t230 + (-t230 * t358 - t239 * t360) * t241) * t241, 0;];
JaD_rot  = t1;
