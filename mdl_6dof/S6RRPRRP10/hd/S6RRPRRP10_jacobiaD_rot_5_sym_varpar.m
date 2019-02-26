% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRP10
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
% Datum: 2019-02-26 21:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRP10_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP10_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:51:04
% EndTime: 2019-02-26 21:51:06
% DurationCPUTime: 1.48s
% Computational Cost: add. (8428->149), mult. (13478->297), div. (726->12), fcn. (17045->13), ass. (0->128)
t270 = cos(pkin(6));
t272 = sin(qJ(2));
t344 = sin(qJ(1));
t305 = t344 * t272;
t295 = t270 * t305;
t299 = t344 * qJD(2);
t274 = cos(qJ(2));
t275 = cos(qJ(1));
t321 = t275 * t274;
t269 = sin(pkin(6));
t325 = t269 * t275;
t349 = -qJD(1) * t295 - t272 * t299 + (qJD(2) * t270 + qJD(1)) * t321 - qJD(4) * t325;
t268 = pkin(11) + qJ(4);
t266 = sin(t268);
t267 = cos(t268);
t304 = t344 * t274;
t322 = t275 * t272;
t287 = -t270 * t322 - t304;
t239 = -t266 * t287 + t267 * t325;
t327 = t269 * t272;
t249 = t266 * t327 - t270 * t267;
t228 = atan2(-t239, t249);
t215 = sin(t228);
t216 = cos(t228);
t206 = -t215 * t239 + t216 * t249;
t204 = 0.1e1 / t206 ^ 2;
t255 = -t295 + t321;
t306 = t269 * t344;
t286 = -t255 * t266 + t267 * t306;
t238 = t286 ^ 2;
t200 = t238 * t204 + 0.1e1;
t285 = -t270 * t304 - t322;
t232 = t287 * qJD(1) + t285 * qJD(2);
t245 = t255 * t267 + t266 * t306;
t303 = qJD(1) * t325;
t210 = t245 * qJD(4) + t232 * t266 - t267 * t303;
t337 = t210 * t204;
t237 = t239 ^ 2;
t247 = 0.1e1 / t249 ^ 2;
t227 = t237 * t247 + 0.1e1;
t221 = 0.1e1 / t227;
t300 = t344 * qJD(1);
t294 = t269 * t300;
t318 = qJD(4) * t267;
t212 = t349 * t266 - t267 * t294 - t287 * t318;
t250 = t270 * t266 + t267 * t327;
t319 = qJD(2) * t274;
t302 = t269 * t319;
t235 = t250 * qJD(4) + t266 * t302;
t246 = 0.1e1 / t249;
t331 = t239 * t247;
t291 = -t212 * t246 + t235 * t331;
t194 = t291 * t221;
t292 = -t215 * t249 - t216 * t239;
t189 = t292 * t194 - t215 * t212 + t216 * t235;
t203 = 0.1e1 / t206;
t205 = t203 * t204;
t342 = t189 * t205;
t316 = 0.2e1 * (-t238 * t342 - t286 * t337) / t200 ^ 2;
t348 = t235 * t247;
t307 = t270 * t321;
t252 = -t305 + t307;
t326 = t269 * t274;
t288 = -t246 * t252 + t326 * t331;
t347 = t266 * t288;
t213 = (qJD(4) * t287 + t294) * t266 + t349 * t267;
t273 = cos(qJ(5));
t271 = sin(qJ(5));
t329 = t285 * t271;
t226 = t245 * t273 - t329;
t218 = 0.1e1 / t226;
t219 = 0.1e1 / t226 ^ 2;
t346 = -0.2e1 * t239;
t345 = -0.2e1 * t286;
t211 = t286 * qJD(4) + t232 * t267 + t266 * t303;
t231 = -qJD(1) * t307 - t275 * t319 + (t270 * t299 + t300) * t272;
t201 = t226 * qJD(5) + t211 * t271 + t231 * t273;
t328 = t285 * t273;
t225 = t245 * t271 + t328;
t217 = t225 ^ 2;
t209 = t217 * t219 + 0.1e1;
t334 = t219 * t225;
t317 = qJD(5) * t225;
t202 = t211 * t273 - t231 * t271 - t317;
t339 = t202 * t218 * t219;
t341 = (t201 * t334 - t217 * t339) / t209 ^ 2;
t333 = t246 * t348;
t340 = (t212 * t331 - t237 * t333) / t227 ^ 2;
t338 = t204 * t286;
t336 = t215 * t286;
t335 = t216 * t286;
t332 = t239 * t246;
t330 = t285 * t266;
t324 = t271 * t218;
t323 = t273 * t225;
t320 = qJD(2) * t272;
t315 = -0.2e1 * t341;
t314 = 0.2e1 * t341;
t313 = -0.2e1 * t340;
t312 = t205 * t345;
t311 = t246 * t340;
t310 = t204 * t336;
t309 = t204 * t335;
t308 = t225 * t339;
t298 = 0.2e1 * t308;
t297 = t333 * t346;
t241 = -t266 * t325 - t267 * t287;
t293 = -qJD(5) * t267 * t285 + t232;
t224 = -t241 * t273 + t252 * t271;
t223 = -t241 * t271 - t252 * t273;
t290 = t219 * t323 - t324;
t289 = -t241 * t246 + t250 * t331;
t283 = -t215 + (t216 * t332 + t215) * t221;
t282 = -qJD(4) * t330 + qJD(5) * t255 + t231 * t267;
t236 = -t249 * qJD(4) + t267 * t302;
t233 = t285 * qJD(1) + t287 * qJD(2);
t230 = t255 * t271 + t267 * t328;
t229 = -t255 * t273 + t267 * t329;
t207 = 0.1e1 / t209;
t198 = 0.1e1 / t200;
t197 = t221 * t347;
t195 = t289 * t221;
t193 = t283 * t286;
t191 = (-t215 * t252 + t216 * t326) * t266 + t292 * t197;
t190 = t292 * t195 - t215 * t241 + t216 * t250;
t188 = t289 * t313 + (t250 * t297 - t213 * t246 + (t212 * t250 + t235 * t241 + t236 * t239) * t247) * t221;
t186 = t313 * t347 + (t288 * t318 + (t297 * t326 - t233 * t246 + (t235 * t252 + (t212 * t274 - t239 * t320) * t269) * t247) * t266) * t221;
t1 = [t311 * t345 + (-t210 * t246 - t286 * t348) * t221, t186, 0, t188, 0, 0; t239 * t203 * t316 + (-t212 * t203 + (t189 * t239 + t193 * t210) * t204) * t198 - (-t193 * t204 * t316 + (-0.2e1 * t193 * t342 + (-t194 * t221 * t332 + t313) * t310 + (t311 * t346 - t194 + (t194 - t291) * t221) * t309 - t283 * t337) * t198) * t286 (-t191 * t338 - t203 * t330) * t316 + (-t191 * t337 + (t231 * t266 + t285 * t318) * t203 + (t191 * t312 - t204 * t330) * t189 + (-t186 * t239 - t197 * t212 + (-t266 * t320 + t274 * t318) * t269 + (-t197 * t249 - t252 * t266) * t194) * t309 + (-t252 * t318 - t186 * t249 - t197 * t235 - t233 * t266 + (t197 * t239 - t266 * t326) * t194) * t310) * t198, 0 (-t190 * t338 - t203 * t245) * t316 + (t190 * t189 * t312 + t211 * t203 + (-t245 * t189 - t190 * t210 + (-t188 * t239 - t195 * t212 + t236 + (-t195 * t249 - t241) * t194) * t335 + (-t188 * t249 - t195 * t235 - t213 + (t195 * t239 - t250) * t194) * t336) * t204) * t198, 0, 0; (-t218 * t223 + t224 * t334) * t314 + ((t224 * qJD(5) - t213 * t271 - t233 * t273) * t218 + t224 * t298 + (-t223 * t202 - (-t223 * qJD(5) - t213 * t273 + t233 * t271) * t225 - t224 * t201) * t219) * t207 (-t218 * t229 + t230 * t334) * t314 + (t230 * t298 - t293 * t218 * t273 + t282 * t324 + (-t293 * t225 * t271 - t230 * t201 - t229 * t202 - t282 * t323) * t219) * t207, 0, -t290 * t286 * t315 + (t290 * t210 - ((-qJD(5) * t218 - 0.2e1 * t308) * t273 + (t201 * t273 + (t202 - t317) * t271) * t219) * t286) * t207, t315 + 0.2e1 * (t201 * t219 * t207 + (-t207 * t339 - t219 * t341) * t225) * t225, 0;];
JaD_rot  = t1;
