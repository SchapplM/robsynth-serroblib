% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP9
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
% Datum: 2019-02-26 21:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRP9_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP9_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:50:28
% EndTime: 2019-02-26 21:50:30
% DurationCPUTime: 1.48s
% Computational Cost: add. (8428->149), mult. (13478->297), div. (726->12), fcn. (17045->13), ass. (0->128)
t272 = cos(pkin(6));
t274 = sin(qJ(2));
t346 = sin(qJ(1));
t307 = t346 * t274;
t297 = t272 * t307;
t302 = qJD(2) * t346;
t276 = cos(qJ(2));
t277 = cos(qJ(1));
t323 = t277 * t276;
t271 = sin(pkin(6));
t326 = t271 * t277;
t351 = -qJD(1) * t297 - t274 * t302 + (qJD(2) * t272 + qJD(1)) * t323 - qJD(4) * t326;
t270 = pkin(11) + qJ(4);
t268 = sin(t270);
t269 = cos(t270);
t306 = t346 * t276;
t324 = t277 * t274;
t289 = -t272 * t324 - t306;
t241 = -t268 * t289 + t269 * t326;
t328 = t271 * t274;
t251 = t268 * t328 - t272 * t269;
t230 = atan2(-t241, t251);
t217 = sin(t230);
t218 = cos(t230);
t208 = -t217 * t241 + t218 * t251;
t206 = 0.1e1 / t208 ^ 2;
t257 = -t297 + t323;
t308 = t271 * t346;
t288 = -t257 * t268 + t269 * t308;
t240 = t288 ^ 2;
t202 = t240 * t206 + 0.1e1;
t287 = -t272 * t306 - t324;
t234 = t289 * qJD(1) + t287 * qJD(2);
t247 = t257 * t269 + t268 * t308;
t305 = qJD(1) * t326;
t212 = t247 * qJD(4) + t234 * t268 - t269 * t305;
t339 = t206 * t288;
t239 = t241 ^ 2;
t249 = 0.1e1 / t251 ^ 2;
t229 = t239 * t249 + 0.1e1;
t223 = 0.1e1 / t229;
t301 = t346 * qJD(1);
t296 = t271 * t301;
t320 = qJD(4) * t269;
t214 = t351 * t268 - t269 * t296 - t289 * t320;
t252 = t272 * t268 + t269 * t328;
t321 = qJD(2) * t276;
t304 = t271 * t321;
t237 = t252 * qJD(4) + t268 * t304;
t248 = 0.1e1 / t251;
t332 = t241 * t249;
t293 = -t214 * t248 + t237 * t332;
t196 = t293 * t223;
t294 = -t217 * t251 - t218 * t241;
t191 = t294 * t196 - t217 * t214 + t218 * t237;
t205 = 0.1e1 / t208;
t207 = t205 * t206;
t344 = t191 * t207;
t318 = 0.2e1 * (-t212 * t339 - t240 * t344) / t202 ^ 2;
t350 = t237 * t249;
t309 = t272 * t323;
t254 = -t307 + t309;
t327 = t271 * t276;
t290 = -t248 * t254 + t327 * t332;
t349 = t268 * t290;
t215 = (qJD(4) * t289 + t296) * t268 + t351 * t269;
t275 = cos(qJ(5));
t273 = sin(qJ(5));
t330 = t287 * t273;
t228 = t247 * t275 - t330;
t220 = 0.1e1 / t228;
t221 = 0.1e1 / t228 ^ 2;
t348 = -0.2e1 * t241;
t347 = -0.2e1 * t288;
t213 = t288 * qJD(4) + t234 * t269 + t268 * t305;
t233 = -qJD(1) * t309 - t277 * t321 + (t272 * t302 + t301) * t274;
t203 = t228 * qJD(5) + t213 * t273 + t233 * t275;
t329 = t287 * t275;
t227 = t247 * t273 + t329;
t219 = t227 ^ 2;
t211 = t219 * t221 + 0.1e1;
t336 = t221 * t227;
t319 = qJD(5) * t227;
t204 = t213 * t275 - t233 * t273 - t319;
t341 = t204 * t220 * t221;
t343 = (t203 * t336 - t219 * t341) / t211 ^ 2;
t334 = t248 * t350;
t342 = (t214 * t332 - t239 * t334) / t229 ^ 2;
t340 = t206 * t212;
t338 = t217 * t288;
t337 = t218 * t288;
t335 = t227 * t275;
t333 = t241 * t248;
t331 = t287 * t268;
t325 = t273 * t220;
t322 = qJD(2) * t274;
t317 = -0.2e1 * t343;
t316 = 0.2e1 * t343;
t315 = -0.2e1 * t342;
t314 = t207 * t347;
t313 = t248 * t342;
t312 = t206 * t338;
t311 = t206 * t337;
t310 = t227 * t341;
t300 = 0.2e1 * t310;
t299 = t334 * t348;
t243 = -t268 * t326 - t269 * t289;
t295 = -qJD(5) * t269 * t287 + t234;
t226 = -t243 * t275 + t254 * t273;
t225 = -t243 * t273 - t254 * t275;
t292 = t221 * t335 - t325;
t291 = -t243 * t248 + t252 * t332;
t285 = -t217 + (t218 * t333 + t217) * t223;
t284 = -qJD(4) * t331 + qJD(5) * t257 + t233 * t269;
t238 = -t251 * qJD(4) + t269 * t304;
t235 = t287 * qJD(1) + t289 * qJD(2);
t232 = t257 * t273 + t269 * t329;
t231 = -t257 * t275 + t269 * t330;
t209 = 0.1e1 / t211;
t200 = 0.1e1 / t202;
t199 = t223 * t349;
t197 = t291 * t223;
t195 = t285 * t288;
t193 = (-t217 * t254 + t218 * t327) * t268 + t294 * t199;
t192 = t294 * t197 - t217 * t243 + t218 * t252;
t190 = t291 * t315 + (t252 * t299 - t215 * t248 + (t214 * t252 + t237 * t243 + t238 * t241) * t249) * t223;
t188 = t315 * t349 + (t290 * t320 + (t299 * t327 - t235 * t248 + (t237 * t254 + (t214 * t276 - t241 * t322) * t271) * t249) * t268) * t223;
t1 = [t313 * t347 + (-t212 * t248 - t288 * t350) * t223, t188, 0, t190, 0, 0; t241 * t205 * t318 + (-t214 * t205 + (t191 * t241 + t195 * t212) * t206) * t200 - (-t195 * t206 * t318 + (-0.2e1 * t195 * t344 + (-t196 * t223 * t333 + t315) * t312 + (t313 * t348 - t196 + (t196 - t293) * t223) * t311 - t285 * t340) * t200) * t288 (-t193 * t339 - t205 * t331) * t318 + (-t193 * t340 + (t233 * t268 + t287 * t320) * t205 + (t193 * t314 - t206 * t331) * t191 + (-t188 * t241 - t199 * t214 + (-t268 * t322 + t276 * t320) * t271 + (-t199 * t251 - t254 * t268) * t196) * t311 + (-t254 * t320 - t188 * t251 - t199 * t237 - t235 * t268 + (t199 * t241 - t268 * t327) * t196) * t312) * t200, 0 (-t192 * t339 - t205 * t247) * t318 + (t192 * t191 * t314 + t213 * t205 + (-t247 * t191 - t192 * t212 + (-t190 * t241 - t197 * t214 + t238 + (-t197 * t251 - t243) * t196) * t337 + (-t190 * t251 - t197 * t237 - t215 + (t197 * t241 - t252) * t196) * t338) * t206) * t200, 0, 0; (-t220 * t225 + t226 * t336) * t316 + ((t226 * qJD(5) - t215 * t273 - t235 * t275) * t220 + t226 * t300 + (-t225 * t204 - (-t225 * qJD(5) - t215 * t275 + t235 * t273) * t227 - t226 * t203) * t221) * t209 (-t220 * t231 + t232 * t336) * t316 + (t232 * t300 - t295 * t220 * t275 + t284 * t325 + (-t295 * t227 * t273 - t232 * t203 - t231 * t204 - t284 * t335) * t221) * t209, 0, -t292 * t288 * t317 + (t292 * t212 - ((-qJD(5) * t220 - 0.2e1 * t310) * t275 + (t203 * t275 + (t204 - t319) * t273) * t221) * t288) * t209, t317 + 0.2e1 * (t203 * t221 * t209 + (-t209 * t341 - t221 * t343) * t227) * t227, 0;];
JaD_rot  = t1;
