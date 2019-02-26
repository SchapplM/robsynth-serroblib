% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRR8_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:07
% EndTime: 2019-02-26 20:08:09
% DurationCPUTime: 1.16s
% Computational Cost: add. (5341->132), mult. (16984->251), div. (538->12), fcn. (21480->15), ass. (0->120)
t286 = cos(pkin(7));
t289 = sin(qJ(3));
t334 = t286 * t289;
t288 = sin(qJ(5));
t291 = cos(qJ(5));
t283 = sin(pkin(12));
t284 = sin(pkin(7));
t285 = sin(pkin(6));
t287 = cos(pkin(6));
t290 = sin(qJ(2));
t354 = cos(pkin(12));
t323 = t354 * t290;
t293 = cos(qJ(2));
t339 = t283 * t293;
t302 = t287 * t339 + t323;
t309 = t283 * t285 * t286 + t284 * t302;
t292 = cos(qJ(3));
t338 = t284 * t285;
t325 = t283 * t338;
t333 = t286 * t292;
t322 = t354 * t293;
t340 = t283 * t290;
t303 = t287 * t340 - t322;
t341 = t303 * t289;
t316 = -t292 * t325 + t302 * t333 - t341;
t237 = t316 * t288 + t309 * t291;
t301 = -t287 * t323 - t339;
t318 = t354 * t338;
t356 = -t289 * t301 + t292 * t318;
t274 = t287 * t322 - t340;
t270 = t274 * qJD(2);
t271 = t301 * qJD(2);
t344 = t274 * t292;
t227 = t271 * t334 + t270 * t292 + (t344 * t286 - t356) * qJD(3);
t342 = t301 * t292;
t253 = t274 * t334 - t289 * t318 - t342;
t250 = t253 ^ 2;
t330 = t290 * t292;
t331 = t289 * t293;
t305 = t286 * t331 + t330;
t337 = t284 * t287;
t265 = t305 * t285 + t289 * t337;
t262 = 0.1e1 / t265 ^ 2;
t242 = t250 * t262 + 0.1e1;
t345 = t253 * t262;
t329 = t292 * t293;
t332 = t289 * t290;
t304 = -t286 * t332 + t329;
t307 = t286 * t329 - t332;
t324 = qJD(3) * t337;
t248 = t292 * t324 + (t304 * qJD(2) + t307 * qJD(3)) * t285;
t261 = 0.1e1 / t265;
t346 = t248 * t261 * t262;
t355 = -0.2e1 * (t227 * t345 - t250 * t346) / t242 ^ 2;
t243 = atan2(-t253, t265);
t238 = sin(t243);
t239 = cos(t243);
t220 = -t238 * t253 + t239 * t265;
t217 = 0.1e1 / t220;
t233 = 0.1e1 / t237;
t218 = 0.1e1 / t220 ^ 2;
t234 = 0.1e1 / t237 ^ 2;
t240 = 0.1e1 / t242;
t210 = (-t227 * t261 + t248 * t345) * t240;
t315 = -t238 * t265 - t239 * t253;
t206 = t315 * t210 - t238 * t227 + t239 * t248;
t353 = t206 * t217 * t218;
t308 = -t286 * t302 + t325;
t255 = t308 * t289 - t292 * t303;
t352 = t218 * t255;
t272 = t302 * qJD(2);
t273 = t303 * qJD(2);
t228 = t255 * qJD(3) - t272 * t289 - t273 * t333;
t236 = t309 * t288 - t316 * t291;
t328 = qJD(5) * t236;
t335 = t284 * t291;
t222 = t228 * t288 - t273 * t335 - t328;
t351 = t222 * t233 * t234;
t336 = t284 * t288;
t221 = -t237 * qJD(5) + t228 * t291 + t273 * t336;
t232 = t236 ^ 2;
t225 = t232 * t234 + 0.1e1;
t349 = t234 * t236;
t350 = 0.1e1 / t225 ^ 2 * (-t221 * t349 - t232 * t351);
t348 = t238 * t255;
t347 = t239 * t255;
t251 = t255 ^ 2;
t216 = t218 * t251 + 0.1e1;
t229 = t273 * t334 - t272 * t292 + (t308 * t292 + t341) * qJD(3);
t327 = 0.2e1 * (t229 * t352 - t251 * t353) / t216 ^ 2;
t326 = 0.2e1 * t350;
t321 = 0.2e1 * t255 * t353;
t320 = 0.2e1 * t236 * t351;
t319 = -0.2e1 * t253 * t346;
t312 = t233 * t291 + t288 * t349;
t252 = -t274 * t333 + t356;
t264 = t307 * t285 + t292 * t337;
t311 = t252 * t261 + t264 * t345;
t257 = t301 * t334 + t344;
t268 = t304 * t285;
t310 = -t257 * t261 + t268 * t345;
t258 = -t289 * t302 - t303 * t333;
t245 = t258 * t288 - t303 * t335;
t244 = -t258 * t291 - t303 * t336;
t259 = -t292 * t302 + t303 * t334;
t306 = -t286 * t330 - t331;
t256 = (-t305 * qJD(2) + t306 * qJD(3)) * t285;
t247 = -t289 * t324 + (t306 * qJD(2) - t305 * qJD(3)) * t285;
t231 = t259 * qJD(3) - t272 * t333 + t273 * t289;
t230 = -t270 * t334 + t271 * t292 + (-t274 * t289 + t301 * t333) * qJD(3);
t226 = -t271 * t333 + t270 * t289 + (-t342 + (t274 * t286 - t318) * t289) * qJD(3);
t223 = 0.1e1 / t225;
t214 = 0.1e1 / t216;
t212 = t310 * t240;
t211 = t311 * t240;
t208 = t315 * t212 - t238 * t257 + t239 * t268;
t207 = t315 * t211 + t238 * t252 + t239 * t264;
t205 = t310 * t355 + (t268 * t319 - t230 * t261 + (t227 * t268 + t248 * t257 + t253 * t256) * t262) * t240;
t204 = t311 * t355 + (t264 * t319 + t226 * t261 + (t227 * t264 + t247 * t253 - t248 * t252) * t262) * t240;
t1 = [0, t205, t204, 0, 0, 0; 0 (t208 * t352 - t217 * t259) * t327 + ((-t258 * qJD(3) + t272 * t334 + t273 * t292) * t217 + t208 * t321 + (-t259 * t206 - t208 * t229 - (-t205 * t253 - t212 * t227 + t256 + (-t212 * t265 - t257) * t210) * t347 - (-t205 * t265 - t212 * t248 - t230 + (t212 * t253 - t268) * t210) * t348) * t218) * t214 (t207 * t352 + t316 * t217) * t327 + (t207 * t321 - t228 * t217 + (t316 * t206 - t207 * t229 - (-t204 * t253 - t211 * t227 + t247 + (-t211 * t265 + t252) * t210) * t347 - (-t204 * t265 - t211 * t248 + t226 + (t211 * t253 - t264) * t210) * t348) * t218) * t214, 0, 0, 0; 0 (-t233 * t244 + t245 * t349) * t326 + ((t245 * qJD(5) - t231 * t291 - t272 * t336) * t233 + t245 * t320 + (-t244 * t222 - (-t244 * qJD(5) + t231 * t288 - t272 * t335) * t236 + t245 * t221) * t234) * t223, t312 * t255 * t326 + (-t312 * t229 + ((qJD(5) * t233 + t320) * t288 + (t221 * t288 + (t222 - t328) * t291) * t234) * t255) * t223, 0, -0.2e1 * t350 + 0.2e1 * (-t221 * t234 * t223 + (-t223 * t351 - t234 * t350) * t236) * t236, 0;];
JaD_rot  = t1;
