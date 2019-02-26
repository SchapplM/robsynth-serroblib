% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPP1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPPP1_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPP1_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:02:48
% EndTime: 2019-02-26 22:02:50
% DurationCPUTime: 2.07s
% Computational Cost: add. (6418->168), mult. (20749->341), div. (681->12), fcn. (25776->13), ass. (0->150)
t278 = cos(qJ(2));
t275 = sin(qJ(2));
t271 = sin(pkin(10));
t273 = cos(pkin(6));
t277 = cos(qJ(3));
t274 = sin(qJ(3));
t362 = cos(pkin(10));
t316 = t274 * t362;
t293 = t271 * t277 + t273 * t316;
t290 = t293 * t275;
t272 = sin(pkin(6));
t318 = t272 * t362;
t243 = t278 * t318 + t290;
t307 = t275 * t318;
t244 = t278 * t293 - t307;
t315 = t277 * t362;
t292 = -t271 * t274 + t273 * t315;
t287 = t292 * qJD(3);
t217 = qJD(2) * t244 + t275 * t287;
t241 = 0.1e1 / t243 ^ 2;
t368 = t217 * t241;
t279 = cos(qJ(1));
t270 = t279 * t277;
t276 = sin(qJ(1));
t256 = t270 * t278 + t276 * t274;
t367 = t256 * t272;
t313 = qJD(2) * t362;
t304 = t272 * t313;
t366 = -t304 + t287;
t333 = qJD(3) * t279;
t319 = t277 * t333;
t340 = qJD(1) * t276;
t365 = -t274 * t340 + t319;
t338 = qJD(1) * t279;
t343 = t276 * t277;
t364 = qJD(3) * t343 + t274 * t338;
t342 = t276 * t278;
t346 = t274 * t279;
t255 = t277 * t342 - t346;
t308 = -t255 * t271 + t276 * t307;
t317 = t273 * t362;
t341 = t274 * t342 + t270;
t221 = t317 * t341 - t308;
t211 = atan2(-t221, t243);
t206 = sin(t211);
t207 = cos(t211);
t200 = -t206 * t221 + t207 * t243;
t197 = 0.1e1 / t200;
t240 = 0.1e1 / t243;
t296 = t278 * t346 - t343;
t344 = t275 * t279;
t249 = t272 * t296 + t273 * t344;
t245 = 0.1e1 / t249;
t198 = 0.1e1 / t200 ^ 2;
t246 = 0.1e1 / t249 ^ 2;
t363 = -0.2e1 * t221;
t289 = t296 * t362;
t225 = t256 * t271 + t273 * t289 - t279 * t307;
t219 = t225 ^ 2;
t196 = t198 * t219 + 0.1e1;
t335 = qJD(2) * t279;
t322 = t275 * t335;
t334 = qJD(3) * t278;
t339 = qJD(1) * t278;
t230 = (qJD(1) - t334) * t270 + (t322 + (-qJD(3) + t339) * t276) * t274;
t321 = t274 * t334;
t231 = t279 * t321 + (t276 * t339 + t322) * t277 - t364;
t314 = qJD(1) * t362;
t305 = t276 * t314;
t301 = t272 * t305;
t302 = t278 * t304;
t201 = -t230 * t317 - t231 * t271 + t275 * t301 - t279 * t302;
t357 = t198 * t225;
t218 = t221 ^ 2;
t210 = t218 * t241 + 0.1e1;
t208 = 0.1e1 / t210;
t337 = qJD(2) * t276;
t323 = t275 * t337;
t232 = -t277 * t340 + (-t323 - t333) * t274 + t364 * t278;
t233 = qJD(1) * t256 - t276 * t321 - t277 * t323 - t319;
t300 = t272 * t279 * t314;
t203 = t232 * t317 + t233 * t271 - t275 * t300 - t276 * t302;
t352 = t221 * t241;
t299 = -t203 * t240 + t217 * t352;
t189 = t299 * t208;
t303 = -t206 * t243 - t207 * t221;
t185 = t189 * t303 - t203 * t206 + t207 * t217;
t360 = t185 * t197 * t198;
t361 = (t201 * t357 - t219 * t360) / t196 ^ 2;
t294 = -t275 * t340 + t278 * t335;
t202 = -t231 * t362 + (t230 * t273 + t272 * t294) * t271;
t226 = t256 * t362 + (t272 * t344 - t273 * t296) * t271;
t220 = t226 ^ 2;
t214 = t220 * t246 + 0.1e1;
t215 = -t230 * t272 + t273 * t294;
t247 = t245 * t246;
t351 = t226 * t246;
t359 = (-t215 * t220 * t247 + t202 * t351) / t214 ^ 2;
t354 = t240 * t368;
t358 = (t203 * t352 - t218 * t354) / t210 ^ 2;
t356 = t206 * t225;
t355 = t207 * t225;
t353 = t221 * t240;
t345 = t275 * t276;
t248 = -t272 * t341 - t273 * t345;
t350 = t246 * t248;
t347 = t274 * t275;
t254 = (-t272 * t347 + t273 * t278) * t279;
t349 = t246 * t254;
t348 = t271 * t273;
t336 = qJD(2) * t278;
t332 = 0.2e1 * t361;
t331 = 0.2e1 * t360;
t330 = 0.2e1 * t359;
t329 = -0.2e1 * t358;
t328 = 0.2e1 * t226 * t247;
t327 = t240 * t358;
t324 = t274 * t335;
t312 = t341 * t273;
t311 = t225 * t331;
t310 = t354 * t363;
t309 = t215 * t328;
t234 = t255 * t317 - t271 * t341;
t253 = t292 * t275;
t298 = -t234 * t240 + t253 * t352;
t237 = t243 * t276;
t297 = t237 * t240 + t244 * t352;
t295 = -t275 * t338 - t276 * t336;
t291 = -t206 + (t207 * t353 + t206) * t208;
t288 = qJD(1) * t293;
t239 = (-t275 * t315 + (t272 * t278 + t273 * t347) * t271) * t279;
t238 = t243 * t279;
t236 = -t256 * t348 - t289;
t235 = t256 * t317 - t271 * t296;
t229 = -qJD(3) * t290 + t292 * t336;
t224 = -t255 * t362 + (-t272 * t345 + t312) * t271;
t223 = -t312 * t362 + t308;
t216 = -qJD(2) * t243 + t278 * t287;
t212 = 0.1e1 / t214;
t205 = -t232 * t271 + t233 * t317;
t204 = (t293 * t337 + t300) * t278 + (t366 * t276 + t279 * t288) * t275;
t194 = 0.1e1 / t196;
t193 = t298 * t208;
t191 = t297 * t208;
t188 = t291 * t225;
t187 = t193 * t303 - t206 * t234 + t207 * t253;
t186 = t191 * t303 + t206 * t237 + t207 * t244;
t184 = t298 * t329 + (t253 * t310 - t205 * t240 + (t203 * t253 + t217 * t234 + t221 * t229) * t241) * t208;
t183 = t297 * t329 + (t244 * t310 + t204 * t240 + (t203 * t244 + t216 * t221 - t217 * t237) * t241) * t208;
t1 = [0.2e1 * t225 * t327 + (-t201 * t240 + t225 * t368) * t208, t183, t184, 0, 0, 0; -0.2e1 * t223 * t197 * t361 + (-t203 * t197 + (-t185 * t223 - t188 * t201) * t198) * t194 + (t188 * t331 * t194 + (t188 * t332 + (-(-t189 * t208 * t353 + t329) * t356 - (t327 * t363 - t189 + (t189 - t299) * t208) * t355 - t291 * t201) * t194) * t198) * t225 (t186 * t357 + t197 * t238) * t332 + (t186 * t311 + (t238 * t185 - t186 * t201 - (-t183 * t221 - t191 * t203 + t216 + (-t191 * t243 + t237) * t189) * t355 - (-t183 * t243 - t191 * t217 + t204 + (t191 * t221 - t244) * t189) * t356) * t198 + ((-t293 * t335 + t301) * t278 + (t276 * t288 - t366 * t279) * t275) * t197) * t194 (t187 * t357 - t197 * t235) * t332 + ((t230 * t271 - t231 * t317) * t197 + t187 * t311 + (-t235 * t185 - t187 * t201 - (-t184 * t221 - t193 * t203 + t229 + (-t193 * t243 - t234) * t189) * t355 - (-t184 * t243 - t193 * t217 - t205 + (t193 * t221 - t253) * t189) * t356) * t198) * t194, 0, 0, 0; -0.2e1 * t224 * t245 * t359 + t226 * t330 * t350 + ((-t233 * t362 + (t232 * t273 + t272 * t295) * t271) * t245 - t224 * t246 * t215 - (-t232 * t272 + t273 * t295) * t351 - t202 * t350 + t248 * t309) * t212 (t226 * t349 - t239 * t245) * t330 + (-t202 * t349 + (-t239 * t246 + t254 * t328) * t215 + ((-t271 * t272 * t340 - t270 * t313 + t324 * t348) * t245 - (-t272 * t324 - t273 * t340) * t351) * t278 + (-(-t365 * t272 - t273 * t335) * t351 + (t277 * t305 + t316 * t333 + (-t272 * t335 + t365 * t273) * t271) * t245) * t275) * t212 (-t236 * t245 + t351 * t367) * t330 + ((t230 * t362 + t231 * t348) * t245 + t309 * t367 + (-t236 * t215 + (-t202 * t256 + t226 * t231) * t272) * t246) * t212, 0, 0, 0;];
JaD_rot  = t1;
