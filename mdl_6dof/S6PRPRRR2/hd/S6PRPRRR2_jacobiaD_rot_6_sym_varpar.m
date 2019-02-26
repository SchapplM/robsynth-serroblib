% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRR2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:54:12
% EndTime: 2019-02-26 19:54:13
% DurationCPUTime: 1.30s
% Computational Cost: add. (6481->116), mult. (17394->236), div. (577->12), fcn. (22601->15), ass. (0->113)
t304 = sin(pkin(12));
t307 = cos(pkin(12));
t311 = sin(qJ(2));
t313 = cos(qJ(2));
t293 = t304 * t311 - t313 * t307;
t309 = cos(pkin(6));
t322 = t293 * t309;
t287 = qJD(2) * t322;
t327 = t304 * t313 + t311 * t307;
t292 = t327 * qJD(2);
t305 = sin(pkin(11));
t308 = cos(pkin(11));
t268 = -t287 * t308 - t292 * t305;
t290 = t327 * t309;
t274 = t290 * t308 - t293 * t305;
t310 = sin(qJ(4));
t306 = sin(pkin(6));
t344 = t306 * t310;
t335 = t308 * t344;
t312 = cos(qJ(4));
t340 = qJD(4) * t312;
t244 = -qJD(4) * t335 + t268 * t310 + t274 * t340;
t343 = t306 * t312;
t262 = t274 * t310 + t308 * t343;
t260 = t262 ^ 2;
t289 = t327 * t306;
t281 = t289 * t310 - t309 * t312;
t279 = 0.1e1 / t281 ^ 2;
t256 = t260 * t279 + 0.1e1;
t254 = 0.1e1 / t256;
t282 = t289 * t312 + t309 * t310;
t288 = t293 * t306;
t286 = qJD(2) * t288;
t258 = t282 * qJD(4) - t286 * t310;
t278 = 0.1e1 / t281;
t348 = t262 * t279;
t224 = (-t244 * t278 + t258 * t348) * t254;
t257 = atan2(-t262, t281);
t252 = sin(t257);
t253 = cos(t257);
t330 = -t252 * t281 - t253 * t262;
t220 = t330 * t224 - t244 * t252 + t253 * t258;
t236 = -t252 * t262 + t253 * t281;
t233 = 0.1e1 / t236;
t234 = 0.1e1 / t236 ^ 2;
t361 = t220 * t233 * t234;
t328 = -t290 * t305 - t293 * t308;
t323 = t305 * t343 - t310 * t328;
t360 = -0.2e1 * t323 * t361;
t273 = -t305 * t327 - t308 * t322;
t324 = -t273 * t278 - t288 * t348;
t359 = t310 * t324;
t349 = t258 * t278 * t279;
t358 = -0.2e1 * (t244 * t348 - t260 * t349) / t256 ^ 2;
t266 = t305 * t344 + t312 * t328;
t276 = t305 * t322 - t308 * t327;
t303 = qJ(5) + qJ(6);
t300 = sin(t303);
t301 = cos(t303);
t249 = t266 * t301 - t276 * t300;
t241 = 0.1e1 / t249;
t242 = 0.1e1 / t249 ^ 2;
t291 = t293 * qJD(2);
t321 = t309 * t292;
t269 = t291 * t308 + t305 * t321;
t302 = qJD(5) + qJD(6);
t332 = t266 * t302 + t269;
t329 = t287 * t305 - t292 * t308;
t247 = t323 * qJD(4) + t312 * t329;
t333 = -t276 * t302 + t247;
t231 = t333 * t300 + t332 * t301;
t248 = t266 * t300 + t276 * t301;
t240 = t248 ^ 2;
t239 = t240 * t242 + 0.1e1;
t353 = t242 * t248;
t232 = -t332 * t300 + t333 * t301;
t356 = t232 * t241 * t242;
t357 = (t231 * t353 - t240 * t356) / t239 ^ 2;
t355 = t234 * t323;
t354 = t241 * t300;
t352 = t248 * t301;
t351 = t252 * t323;
t350 = t253 * t323;
t347 = t276 * t310;
t346 = t276 * t312;
t261 = t323 ^ 2;
t230 = t234 * t261 + 0.1e1;
t246 = t266 * qJD(4) + t310 * t329;
t339 = 0.2e1 * (-t246 * t355 - t261 * t361) / t230 ^ 2;
t338 = -0.2e1 * t357;
t336 = t248 * t356;
t334 = -0.2e1 * t262 * t349;
t331 = t302 * t346 - t329;
t326 = t242 * t352 - t354;
t264 = t274 * t312 - t335;
t325 = -t264 * t278 + t282 * t348;
t320 = -qJD(4) * t347 + t269 * t312 + t302 * t328;
t285 = t306 * t292;
t267 = t305 * t291 - t308 * t321;
t259 = -t281 * qJD(4) - t286 * t312;
t251 = t300 * t328 + t301 * t346;
t250 = t300 * t346 - t301 * t328;
t245 = -t262 * qJD(4) + t268 * t312;
t237 = 0.1e1 / t239;
t228 = 0.1e1 / t230;
t226 = t254 * t359;
t225 = t325 * t254;
t222 = (-t252 * t273 - t253 * t288) * t310 + t330 * t226;
t221 = t330 * t225 - t252 * t264 + t253 * t282;
t219 = t325 * t358 + (t282 * t334 - t245 * t278 + (t244 * t282 + t258 * t264 + t259 * t262) * t279) * t254;
t217 = t358 * t359 + (t324 * t340 + (-t288 * t334 - t267 * t278 + (-t244 * t288 + t258 * t273 - t262 * t285) * t279) * t310) * t254;
t216 = t338 + 0.2e1 * (t231 * t237 * t242 + (-t237 * t356 - t242 * t357) * t248) * t248;
t1 = [0, t217, 0, t219, 0, 0; 0 (-t222 * t355 - t233 * t347) * t339 + ((t269 * t310 + t276 * t340) * t233 + t222 * t360 + (-t222 * t246 - t347 * t220 + (-t288 * t340 - t217 * t262 - t226 * t244 - t285 * t310 + (-t226 * t281 - t273 * t310) * t224) * t350 + (-t273 * t340 - t217 * t281 - t226 * t258 - t267 * t310 + (t226 * t262 + t288 * t310) * t224) * t351) * t234) * t228, 0 (-t221 * t355 - t233 * t266) * t339 + (t221 * t360 + t247 * t233 + (-t266 * t220 - t221 * t246 + (-t219 * t262 - t225 * t244 + t259 + (-t225 * t281 - t264) * t224) * t350 + (-t219 * t281 - t225 * t258 - t245 + (t225 * t262 - t282) * t224) * t351) * t234) * t228, 0, 0; 0, 0.2e1 * (-t241 * t250 + t251 * t353) * t357 + (0.2e1 * t251 * t336 + t331 * t241 * t301 + t320 * t354 + (t331 * t248 * t300 - t251 * t231 - t250 * t232 - t320 * t352) * t242) * t237, 0, -t326 * t323 * t338 + (t326 * t246 - ((-t241 * t302 - 0.2e1 * t336) * t301 + (t231 * t301 + (-t248 * t302 + t232) * t300) * t242) * t323) * t237, t216, t216;];
JaD_rot  = t1;
