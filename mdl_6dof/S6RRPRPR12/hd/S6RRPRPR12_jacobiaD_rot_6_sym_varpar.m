% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR12_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR12_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:44:22
% EndTime: 2019-02-26 21:44:24
% DurationCPUTime: 1.73s
% Computational Cost: add. (8428->153), mult. (13478->302), div. (726->12), fcn. (17045->13), ass. (0->128)
t267 = sin(qJ(2));
t268 = sin(qJ(1));
t270 = cos(qJ(2));
t271 = cos(qJ(1));
t349 = cos(pkin(6));
t301 = t271 * t349;
t253 = t267 * t268 - t270 * t301;
t264 = qJ(4) + pkin(11);
t262 = sin(t264);
t263 = cos(t264);
t265 = sin(pkin(6));
t327 = t265 * t271;
t290 = t253 * t263 + t262 * t327;
t238 = t290 ^ 2;
t328 = t265 * t270;
t287 = -t349 * t262 - t263 * t328;
t248 = 0.1e1 / t287 ^ 2;
t227 = t238 * t248 + 0.1e1;
t221 = 0.1e1 / t227;
t324 = qJD(1) * t268;
t308 = t265 * t324;
t286 = -t267 * t301 - t268 * t270;
t302 = t268 * t349;
t288 = t271 * t267 + t270 * t302;
t279 = t288 * qJD(1) - t286 * qJD(2);
t309 = t263 * t327;
t357 = -qJD(4) * t309 - t279 * t263;
t210 = (qJD(4) * t253 + t308) * t262 + t357;
t251 = -t262 * t328 + t349 * t263;
t323 = qJD(2) * t267;
t306 = t265 * t323;
t236 = t251 * qJD(4) - t263 * t306;
t247 = 0.1e1 / t287;
t335 = t290 * t248;
t294 = t210 * t247 - t236 * t335;
t194 = t294 * t221;
t228 = atan2(t290, -t287);
t215 = sin(t228);
t216 = cos(t228);
t296 = t215 * t287 + t216 * t290;
t189 = t296 * t194 - t210 * t215 + t216 * t236;
t206 = t215 * t290 - t216 * t287;
t204 = 0.1e1 / t206 ^ 2;
t358 = t189 * t204;
t321 = qJD(4) * t262;
t303 = t265 * t321;
t320 = qJD(4) * t263;
t211 = t253 * t320 + t279 * t262 + t263 * t308 + t271 * t303;
t329 = t265 * t268;
t240 = t288 * t262 + t263 * t329;
t297 = t267 * t302;
t326 = t271 * t270;
t255 = -t297 + t326;
t266 = sin(qJ(6));
t269 = cos(qJ(6));
t223 = t240 * t266 - t255 * t269;
t356 = 0.2e1 * t223;
t203 = 0.1e1 / t206;
t355 = t203 * t358;
t282 = (t349 * qJD(1) + qJD(2)) * t326 - qJD(2) * t297 - t267 * t324;
t307 = qJD(1) * t327;
t213 = t240 * qJD(4) + t262 * t307 - t282 * t263;
t284 = t288 * t263;
t239 = t262 * t329 - t284;
t300 = 0.2e1 * t239 * t355;
t354 = -t204 * t213 + t300;
t353 = t236 * t248;
t330 = t265 * t267;
t310 = t290 * t330;
t289 = t247 * t286 + t248 * t310;
t352 = t263 * t289;
t351 = -t255 * t262 * qJD(6) - t282;
t350 = -qJD(6) * t288 + t255 * t320;
t224 = t240 * t269 + t255 * t266;
t218 = 0.1e1 / t224;
t219 = 0.1e1 / t224 ^ 2;
t237 = t239 ^ 2;
t200 = t204 * t237 + 0.1e1;
t342 = t204 * t239;
t348 = (t213 * t342 - t237 * t355) / t200 ^ 2;
t214 = qJD(4) * t284 + t282 * t262 + t263 * t307 - t268 * t303;
t233 = t286 * qJD(1) - t288 * qJD(2);
t201 = t224 * qJD(6) + t214 * t266 - t233 * t269;
t217 = t223 ^ 2;
t209 = t217 * t219 + 0.1e1;
t338 = t219 * t223;
t319 = qJD(6) * t223;
t202 = t214 * t269 + t233 * t266 - t319;
t344 = t202 * t218 * t219;
t347 = (t201 * t338 - t217 * t344) / t209 ^ 2;
t337 = t247 * t353;
t345 = (-t210 * t335 + t238 * t337) / t227 ^ 2;
t207 = 0.1e1 / t209;
t341 = t207 * t219;
t340 = t215 * t239;
t339 = t216 * t239;
t336 = t290 * t247;
t334 = t290 * t251;
t333 = t255 * t263;
t332 = t262 * t266;
t331 = t262 * t269;
t322 = qJD(2) * t270;
t318 = 0.2e1 * t348;
t317 = -0.2e1 * t347;
t316 = -0.2e1 * t345;
t315 = 0.2e1 * t345;
t313 = t219 * t347;
t312 = t201 * t341;
t311 = t223 * t344;
t299 = t247 * t315;
t298 = 0.2e1 * t311;
t291 = -t253 * t262 + t309;
t226 = t266 * t286 + t269 * t291;
t225 = t266 * t291 - t269 * t286;
t293 = -t218 * t266 + t269 * t338;
t292 = t247 * t291 + t248 * t334;
t285 = -t215 + (t216 * t336 + t215) * t221;
t235 = t287 * qJD(4) + t262 * t306;
t234 = -qJD(1) * t297 - t268 * t323 + (qJD(2) * t349 + qJD(1)) * t326;
t230 = t255 * t331 - t288 * t266;
t198 = 0.1e1 / t200;
t197 = t221 * t352;
t195 = t292 * t221;
t191 = (-t215 * t286 - t216 * t330) * t263 + t296 * t197;
t190 = -t296 * t195 + t215 * t291 + t216 * t251;
t188 = t292 * t315 + (-0.2e1 * t334 * t337 + t211 * t247 + (t210 * t251 - t235 * t290 - t236 * t291) * t248) * t221;
t186 = t316 * t352 + (-t289 * t321 + (0.2e1 * t310 * t337 - t234 * t247 + (t236 * t286 + (-t210 * t267 + t290 * t322) * t265) * t248) * t263) * t221;
t1 = [-t239 * t299 + (t213 * t247 + t239 * t353) * t221, t186, 0, t188, 0, 0; -0.2e1 * t290 * t203 * t348 + ((-t253 * t321 - t262 * t308 - t357) * t203 - t290 * t358 - (t285 * t213 + ((-t194 * t221 * t336 + t316) * t215 + (-t290 * t299 - t194 + (t194 - t294) * t221) * t216) * t239) * t342) * t198 + (t354 * t198 + t342 * t318) * t285 * t239 (t191 * t342 + t203 * t333) * t318 + ((-t233 * t263 + t255 * t321) * t203 + t354 * t191 + (t333 * t189 - (t186 * t290 - t197 * t210 + (-t263 * t322 + t267 * t321) * t265 + (t197 * t287 - t263 * t286) * t194) * t339 - (t286 * t321 + t186 * t287 - t197 * t236 + t234 * t263 + (-t197 * t290 + t263 * t330) * t194) * t340) * t204) * t198, 0 (t190 * t342 - t203 * t240) * t318 + (t190 * t300 + t214 * t203 + (-t240 * t189 - t190 * t213 - (t188 * t290 + t195 * t210 + t235 + (-t195 * t287 + t291) * t194) * t339 - (t188 * t287 + t195 * t236 - t211 + (t195 * t290 - t251) * t194) * t340) * t204) * t198, 0, 0; 0.2e1 * (-t218 * t225 + t226 * t338) * t347 + ((t226 * qJD(6) - t211 * t266 + t234 * t269) * t218 + t226 * t298 + (-t225 * t202 - (-t225 * qJD(6) - t211 * t269 - t234 * t266) * t223 - t226 * t201) * t219) * t207 (t313 * t356 - t312) * t230 + (-t202 * t341 + t218 * t317) * (t255 * t332 + t288 * t269) + (t230 * t298 + (t332 * t218 - t331 * t338) * t233 + (-t351 * t218 - t350 * t338) * t269 + (t350 * t218 - t351 * t338) * t266) * t207, 0, t293 * t239 * t317 + (t293 * t213 + ((-qJD(6) * t218 - 0.2e1 * t311) * t269 + (t201 * t269 + (t202 - t319) * t266) * t219) * t239) * t207, 0, t317 + (t312 + (-t207 * t344 - t313) * t223) * t356;];
JaD_rot  = t1;
