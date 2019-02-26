% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRR7_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR7_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:50:42
% EndTime: 2019-02-26 22:50:43
% DurationCPUTime: 1.50s
% Computational Cost: add. (6358->150), mult. (14832->298), div. (762->12), fcn. (18655->13), ass. (0->132)
t288 = cos(pkin(6));
t290 = sin(qJ(2));
t364 = sin(qJ(1));
t327 = t364 * t290;
t313 = t288 * t327;
t322 = qJD(2) * t364;
t292 = cos(qJ(2));
t293 = cos(qJ(1));
t342 = t293 * t292;
t287 = sin(pkin(6));
t344 = t287 * t293;
t369 = -qJD(1) * t313 - t290 * t322 + (qJD(2) * t288 + qJD(1)) * t342 - qJD(3) * t344;
t289 = sin(qJ(3));
t291 = cos(qJ(3));
t326 = t364 * t292;
t343 = t293 * t290;
t305 = -t288 * t343 - t326;
t256 = -t289 * t305 + t291 * t344;
t346 = t287 * t290;
t267 = -t288 * t291 + t289 * t346;
t247 = atan2(-t256, t267);
t240 = sin(t247);
t241 = cos(t247);
t223 = -t240 * t256 + t241 * t267;
t221 = 0.1e1 / t223 ^ 2;
t272 = -t313 + t342;
t328 = t287 * t364;
t304 = -t272 * t289 + t291 * t328;
t253 = t304 ^ 2;
t219 = t221 * t253 + 0.1e1;
t303 = -t288 * t326 - t343;
t249 = t305 * qJD(1) + t303 * qJD(2);
t262 = t272 * t291 + t289 * t328;
t325 = qJD(1) * t344;
t227 = t262 * qJD(3) + t249 * t289 - t291 * t325;
t358 = t221 * t304;
t252 = t256 ^ 2;
t265 = 0.1e1 / t267 ^ 2;
t246 = t252 * t265 + 0.1e1;
t244 = 0.1e1 / t246;
t321 = t364 * qJD(1);
t312 = t287 * t321;
t339 = qJD(3) * t291;
t229 = t369 * t289 - t291 * t312 - t305 * t339;
t268 = t288 * t289 + t291 * t346;
t340 = qJD(2) * t292;
t324 = t287 * t340;
t254 = t268 * qJD(3) + t289 * t324;
t264 = 0.1e1 / t267;
t349 = t256 * t265;
t309 = -t229 * t264 + t254 * t349;
t211 = t309 * t244;
t310 = -t240 * t267 - t241 * t256;
t206 = t310 * t211 - t240 * t229 + t241 * t254;
t220 = 0.1e1 / t223;
t222 = t220 * t221;
t362 = t206 * t222;
t338 = 0.2e1 * (-t227 * t358 - t253 * t362) / t219 ^ 2;
t368 = t254 * t265;
t329 = t288 * t342;
t269 = -t327 + t329;
t345 = t287 * t292;
t306 = -t264 * t269 + t345 * t349;
t367 = t289 * t306;
t230 = (qJD(3) * t305 + t312) * t289 + t369 * t291;
t286 = qJ(4) + qJ(5) + qJ(6);
t283 = sin(t286);
t284 = cos(t286);
t239 = t262 * t284 - t283 * t303;
t233 = 0.1e1 / t239;
t234 = 0.1e1 / t239 ^ 2;
t366 = -0.2e1 * t256;
t365 = -0.2e1 * t304;
t351 = t264 * t368;
t361 = (t229 * t349 - t252 * t351) / t246 ^ 2;
t248 = -qJD(1) * t329 - t293 * t340 + (t288 * t322 + t321) * t290;
t285 = qJD(4) + qJD(5) + qJD(6);
t316 = t262 * t285 + t248;
t228 = t304 * qJD(3) + t249 * t291 + t289 * t325;
t318 = -t285 * t303 + t228;
t215 = -t316 * t283 + t318 * t284;
t360 = t215 * t233 * t234;
t359 = t221 * t227;
t214 = t318 * t283 + t316 * t284;
t238 = t262 * t283 + t284 * t303;
t232 = t238 ^ 2;
t226 = t232 * t234 + 0.1e1;
t355 = t234 * t238;
t357 = 0.1e1 / t226 ^ 2 * (t214 * t355 - t232 * t360);
t356 = t233 * t283;
t354 = t238 * t284;
t353 = t240 * t304;
t352 = t241 * t304;
t350 = t256 * t264;
t348 = t303 * t289;
t347 = t303 * t291;
t341 = qJD(2) * t290;
t337 = -0.2e1 * t361;
t336 = t222 * t365;
t335 = -0.2e1 * t357;
t334 = 0.2e1 * t357;
t333 = t264 * t361;
t332 = t238 * t360;
t331 = t221 * t353;
t330 = t221 * t352;
t320 = 0.2e1 * t332;
t319 = t351 * t366;
t317 = t269 * t285 - t230;
t250 = t303 * qJD(1) + t305 * qJD(2);
t258 = -t289 * t344 - t291 * t305;
t315 = -t258 * t285 - t250;
t311 = -t285 * t347 + t249;
t308 = t234 * t354 - t356;
t307 = -t258 * t264 + t268 * t349;
t301 = -t240 + (t241 * t350 + t240) * t244;
t300 = -qJD(3) * t348 + t248 * t291 + t272 * t285;
t255 = -t267 * qJD(3) + t291 * t324;
t243 = t272 * t283 + t284 * t347;
t242 = -t272 * t284 + t283 * t347;
t237 = -t258 * t284 + t269 * t283;
t236 = -t258 * t283 - t269 * t284;
t224 = 0.1e1 / t226;
t217 = 0.1e1 / t219;
t216 = t244 * t367;
t213 = t307 * t244;
t210 = t301 * t304;
t208 = (-t240 * t269 + t241 * t345) * t289 + t310 * t216;
t207 = t310 * t213 - t240 * t258 + t241 * t268;
t205 = t307 * t337 + (t268 * t319 - t230 * t264 + (t229 * t268 + t254 * t258 + t255 * t256) * t265) * t244;
t203 = t337 * t367 + (t306 * t339 + (t319 * t345 - t250 * t264 + (t254 * t269 + (t229 * t292 - t256 * t341) * t287) * t265) * t289) * t244;
t202 = t335 + 0.2e1 * (t214 * t234 * t224 + (-t224 * t360 - t234 * t357) * t238) * t238;
t1 = [t333 * t365 + (-t227 * t264 - t304 * t368) * t244, t203, t205, 0, 0, 0; t256 * t220 * t338 + (-t229 * t220 + (t206 * t256 + t210 * t227) * t221) * t217 - (-t210 * t221 * t338 + (-0.2e1 * t210 * t362 + (-t211 * t244 * t350 + t337) * t331 + (t333 * t366 - t211 + (t211 - t309) * t244) * t330 - t301 * t359) * t217) * t304 (-t208 * t358 - t220 * t348) * t338 + (-t208 * t359 + (t248 * t289 + t303 * t339) * t220 + (t208 * t336 - t221 * t348) * t206 + (-t203 * t256 - t216 * t229 + (-t289 * t341 + t292 * t339) * t287 + (-t216 * t267 - t269 * t289) * t211) * t330 + (-t269 * t339 - t203 * t267 - t216 * t254 - t250 * t289 + (t216 * t256 - t289 * t345) * t211) * t331) * t217 (-t207 * t358 - t220 * t262) * t338 + (t207 * t206 * t336 + t228 * t220 + (-t262 * t206 - t207 * t227 + (-t205 * t256 - t213 * t229 + t255 + (-t213 * t267 - t258) * t211) * t352 + (-t205 * t267 - t213 * t254 - t230 + (t213 * t256 - t268) * t211) * t353) * t221) * t217, 0, 0, 0; (-t233 * t236 + t237 * t355) * t334 + ((t317 * t283 + t315 * t284) * t233 + t237 * t320 + (-t236 * t215 - (-t315 * t283 + t317 * t284) * t238 - t237 * t214) * t234) * t224 (-t233 * t242 + t243 * t355) * t334 + (t243 * t320 - t311 * t233 * t284 + t300 * t356 + (-t311 * t238 * t283 - t243 * t214 - t242 * t215 - t300 * t354) * t234) * t224, -t308 * t304 * t335 + (t308 * t227 - ((-t233 * t285 - 0.2e1 * t332) * t284 + (t214 * t284 + (-t238 * t285 + t215) * t283) * t234) * t304) * t224, t202, t202, t202;];
JaD_rot  = t1;
