% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR9_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR9_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:42:23
% EndTime: 2019-02-26 21:42:25
% DurationCPUTime: 1.44s
% Computational Cost: add. (8849->150), mult. (13478->297), div. (726->12), fcn. (17045->13), ass. (0->128)
t275 = cos(pkin(6));
t276 = sin(qJ(2));
t346 = sin(qJ(1));
t308 = t346 * t276;
t298 = t275 * t308;
t302 = t346 * qJD(2);
t277 = cos(qJ(2));
t278 = cos(qJ(1));
t324 = t278 * t277;
t274 = sin(pkin(6));
t326 = t274 * t278;
t351 = -qJD(1) * t298 - t276 * t302 + (qJD(2) * t275 + qJD(1)) * t324 - qJD(4) * t326;
t273 = pkin(11) + qJ(4);
t269 = sin(t273);
t271 = cos(t273);
t307 = t346 * t277;
t325 = t278 * t276;
t290 = -t275 * t325 - t307;
t241 = -t269 * t290 + t271 * t326;
t328 = t274 * t276;
t252 = t269 * t328 - t275 * t271;
t230 = atan2(-t241, t252);
t225 = sin(t230);
t226 = cos(t230);
t208 = -t225 * t241 + t226 * t252;
t206 = 0.1e1 / t208 ^ 2;
t257 = -t298 + t324;
t309 = t274 * t346;
t289 = -t257 * t269 + t271 * t309;
t240 = t289 ^ 2;
t202 = t240 * t206 + 0.1e1;
t288 = -t275 * t307 - t325;
t234 = t290 * qJD(1) + t288 * qJD(2);
t247 = t257 * t271 + t269 * t309;
t306 = qJD(1) * t326;
t212 = t247 * qJD(4) + t234 * t269 - t271 * t306;
t339 = t212 * t206;
t239 = t241 ^ 2;
t250 = 0.1e1 / t252 ^ 2;
t229 = t239 * t250 + 0.1e1;
t227 = 0.1e1 / t229;
t303 = t346 * qJD(1);
t297 = t274 * t303;
t321 = qJD(4) * t271;
t214 = t269 * t351 - t271 * t297 - t290 * t321;
t253 = t275 * t269 + t271 * t328;
t322 = qJD(2) * t277;
t305 = t274 * t322;
t237 = t253 * qJD(4) + t269 * t305;
t249 = 0.1e1 / t252;
t333 = t241 * t250;
t294 = -t214 * t249 + t237 * t333;
t196 = t294 * t227;
t295 = -t225 * t252 - t226 * t241;
t191 = t295 * t196 - t225 * t214 + t226 * t237;
t205 = 0.1e1 / t208;
t207 = t205 * t206;
t344 = t191 * t207;
t319 = 0.2e1 * (-t240 * t344 - t289 * t339) / t202 ^ 2;
t350 = t237 * t250;
t310 = t275 * t324;
t254 = -t308 + t310;
t327 = t274 * t277;
t291 = -t249 * t254 + t327 * t333;
t349 = t269 * t291;
t215 = (qJD(4) * t290 + t297) * t269 + t351 * t271;
t272 = pkin(12) + qJ(6);
t268 = sin(t272);
t270 = cos(t272);
t224 = t247 * t270 - t268 * t288;
t218 = 0.1e1 / t224;
t219 = 0.1e1 / t224 ^ 2;
t348 = -0.2e1 * t241;
t347 = -0.2e1 * t289;
t213 = t289 * qJD(4) + t234 * t271 + t269 * t306;
t233 = -qJD(1) * t310 - t278 * t322 + (t275 * t302 + t303) * t276;
t203 = t224 * qJD(6) + t213 * t268 + t233 * t270;
t223 = t247 * t268 + t270 * t288;
t217 = t223 ^ 2;
t211 = t217 * t219 + 0.1e1;
t338 = t219 * t223;
t320 = qJD(6) * t223;
t204 = t213 * t270 - t233 * t268 - t320;
t341 = t204 * t218 * t219;
t343 = (t203 * t338 - t217 * t341) / t211 ^ 2;
t335 = t249 * t350;
t342 = (t214 * t333 - t239 * t335) / t229 ^ 2;
t340 = t206 * t289;
t337 = t225 * t289;
t336 = t226 * t289;
t334 = t241 * t249;
t332 = t288 * t269;
t331 = t288 * t271;
t330 = t268 * t218;
t329 = t270 * t223;
t323 = qJD(2) * t276;
t318 = -0.2e1 * t343;
t317 = 0.2e1 * t343;
t316 = -0.2e1 * t342;
t315 = t207 * t347;
t314 = t249 * t342;
t313 = t206 * t337;
t312 = t206 * t336;
t311 = t223 * t341;
t301 = 0.2e1 * t311;
t300 = t335 * t348;
t243 = -t269 * t326 - t271 * t290;
t296 = -qJD(6) * t331 + t234;
t222 = -t243 * t270 + t254 * t268;
t221 = -t243 * t268 - t254 * t270;
t293 = t219 * t329 - t330;
t292 = -t243 * t249 + t253 * t333;
t286 = -t225 + (t226 * t334 + t225) * t227;
t285 = -qJD(4) * t332 + qJD(6) * t257 + t233 * t271;
t238 = -t252 * qJD(4) + t271 * t305;
t235 = t288 * qJD(1) + t290 * qJD(2);
t232 = t257 * t268 + t270 * t331;
t231 = -t257 * t270 + t268 * t331;
t209 = 0.1e1 / t211;
t200 = 0.1e1 / t202;
t199 = t227 * t349;
t197 = t292 * t227;
t195 = t286 * t289;
t193 = (-t225 * t254 + t226 * t327) * t269 + t295 * t199;
t192 = t295 * t197 - t225 * t243 + t226 * t253;
t190 = t292 * t316 + (t253 * t300 - t215 * t249 + (t214 * t253 + t237 * t243 + t238 * t241) * t250) * t227;
t188 = t316 * t349 + (t291 * t321 + (t300 * t327 - t235 * t249 + (t237 * t254 + (t214 * t277 - t241 * t323) * t274) * t250) * t269) * t227;
t1 = [t314 * t347 + (-t212 * t249 - t289 * t350) * t227, t188, 0, t190, 0, 0; t241 * t205 * t319 + (-t214 * t205 + (t191 * t241 + t195 * t212) * t206) * t200 - (-t195 * t206 * t319 + (-0.2e1 * t195 * t344 + (-t196 * t227 * t334 + t316) * t313 + (t314 * t348 - t196 + (t196 - t294) * t227) * t312 - t286 * t339) * t200) * t289 (-t193 * t340 - t205 * t332) * t319 + (-t193 * t339 + (t233 * t269 + t288 * t321) * t205 + (t193 * t315 - t206 * t332) * t191 + (-t188 * t241 - t199 * t214 + (-t269 * t323 + t277 * t321) * t274 + (-t199 * t252 - t254 * t269) * t196) * t312 + (-t254 * t321 - t188 * t252 - t199 * t237 - t235 * t269 + (t199 * t241 - t269 * t327) * t196) * t313) * t200, 0 (-t192 * t340 - t205 * t247) * t319 + (t192 * t191 * t315 + t213 * t205 + (-t247 * t191 - t192 * t212 + (-t190 * t241 - t197 * t214 + t238 + (-t197 * t252 - t243) * t196) * t336 + (-t190 * t252 - t197 * t237 - t215 + (t197 * t241 - t253) * t196) * t337) * t206) * t200, 0, 0; (-t218 * t221 + t222 * t338) * t317 + ((t222 * qJD(6) - t215 * t268 - t235 * t270) * t218 + t222 * t301 + (-t221 * t204 - (-t221 * qJD(6) - t215 * t270 + t235 * t268) * t223 - t222 * t203) * t219) * t209 (-t218 * t231 + t232 * t338) * t317 + (t232 * t301 - t296 * t218 * t270 + t285 * t330 + (-t296 * t223 * t268 - t232 * t203 - t231 * t204 - t285 * t329) * t219) * t209, 0, -t293 * t289 * t318 + (t293 * t212 - ((-qJD(6) * t218 - 0.2e1 * t311) * t270 + (t203 * t270 + (t204 - t320) * t268) * t219) * t289) * t209, 0, t318 + 0.2e1 * (t203 * t219 * t209 + (-t209 * t341 - t219 * t343) * t223) * t223;];
JaD_rot  = t1;
