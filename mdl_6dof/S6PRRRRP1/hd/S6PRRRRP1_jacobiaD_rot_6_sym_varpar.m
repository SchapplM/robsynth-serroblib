% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRRP1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:15:13
% EndTime: 2019-02-26 20:15:14
% DurationCPUTime: 1.14s
% Computational Cost: add. (9025->114), mult. (13312->231), div. (822->12), fcn. (17117->13), ass. (0->113)
t268 = sin(pkin(11));
t270 = cos(pkin(11));
t273 = sin(qJ(2));
t271 = cos(pkin(6));
t275 = cos(qJ(2));
t302 = t271 * t275;
t256 = -t268 * t273 + t270 * t302;
t252 = t256 * qJD(2);
t303 = t271 * t273;
t257 = t268 * t275 + t270 * t303;
t267 = qJ(3) + qJ(4);
t264 = sin(t267);
t266 = qJD(3) + qJD(4);
t269 = sin(pkin(6));
t306 = t269 * t270;
t291 = t264 * t306;
t265 = cos(t267);
t308 = t265 * t266;
t219 = t252 * t264 + t257 * t308 - t266 * t291;
t241 = t257 * t264 + t265 * t306;
t239 = t241 ^ 2;
t305 = t269 * t273;
t293 = t264 * t305;
t249 = -t271 * t265 + t293;
t247 = 0.1e1 / t249 ^ 2;
t233 = t239 * t247 + 0.1e1;
t231 = 0.1e1 / t233;
t300 = qJD(2) * t275;
t284 = t266 * t271 + t269 * t300;
t292 = t265 * t305;
t237 = t284 * t264 + t266 * t292;
t246 = 0.1e1 / t249;
t313 = t241 * t247;
t203 = (-t219 * t246 + t237 * t313) * t231;
t234 = atan2(-t241, t249);
t229 = sin(t234);
t230 = cos(t234);
t287 = -t229 * t249 - t230 * t241;
t199 = t287 * t203 - t229 * t219 + t230 * t237;
t213 = -t229 * t241 + t230 * t249;
t210 = 0.1e1 / t213;
t211 = 0.1e1 / t213 ^ 2;
t327 = t199 * t210 * t211;
t294 = t268 * t303;
t259 = t270 * t275 - t294;
t307 = t268 * t269;
t244 = t259 * t264 - t265 * t307;
t326 = 0.2e1 * t244 * t327;
t304 = t269 * t275;
t283 = -t246 * t256 + t304 * t313;
t325 = t264 * t283;
t314 = t237 * t246 * t247;
t324 = -0.2e1 * (t219 * t313 - t239 * t314) / t233 ^ 2;
t245 = t259 * t265 + t264 * t307;
t274 = cos(qJ(5));
t258 = t268 * t302 + t270 * t273;
t272 = sin(qJ(5));
t311 = t258 * t272;
t228 = t245 * t274 + t311;
t224 = 0.1e1 / t228;
t225 = 0.1e1 / t228 ^ 2;
t254 = t258 * qJD(2);
t289 = t266 * t307 - t254;
t309 = t264 * t266;
t222 = -t259 * t309 + t289 * t265;
t255 = -qJD(2) * t294 + t270 * t300;
t214 = t228 * qJD(5) + t222 * t272 - t255 * t274;
t310 = t258 * t274;
t227 = t245 * t272 - t310;
t223 = t227 ^ 2;
t218 = t223 * t225 + 0.1e1;
t318 = t225 * t227;
t299 = qJD(5) * t227;
t215 = t222 * t274 + t255 * t272 - t299;
t321 = t215 * t224 * t225;
t323 = (t214 * t318 - t223 * t321) / t218 ^ 2;
t322 = t211 * t244;
t221 = t259 * t308 + t289 * t264;
t320 = t221 * t211;
t319 = t224 * t272;
t317 = t227 * t274;
t316 = t229 * t244;
t315 = t230 * t244;
t312 = t258 * t264;
t301 = qJD(2) * t273;
t240 = t244 ^ 2;
t209 = t240 * t211 + 0.1e1;
t298 = 0.2e1 * (-t240 * t327 + t244 * t320) / t209 ^ 2;
t297 = -0.2e1 * t323;
t295 = t227 * t321;
t290 = -0.2e1 * t241 * t314;
t288 = qJD(5) * t258 * t265 - t254;
t286 = t225 * t317 - t319;
t243 = t257 * t265 - t291;
t250 = t271 * t264 + t292;
t285 = -t243 * t246 + t250 * t313;
t282 = qJD(5) * t259 - t255 * t265 + t258 * t309;
t253 = t257 * qJD(2);
t238 = t284 * t265 - t266 * t293;
t236 = t259 * t272 - t265 * t310;
t235 = -t259 * t274 - t265 * t311;
t220 = -t257 * t309 + (-t266 * t306 + t252) * t265;
t216 = 0.1e1 / t218;
t207 = 0.1e1 / t209;
t205 = t231 * t325;
t204 = t285 * t231;
t201 = (-t229 * t256 + t230 * t304) * t264 + t287 * t205;
t200 = t287 * t204 - t229 * t243 + t230 * t250;
t198 = t285 * t324 + (t250 * t290 - t220 * t246 + (t219 * t250 + t237 * t243 + t238 * t241) * t247) * t231;
t196 = t324 * t325 + (t283 * t308 + (t290 * t304 + t246 * t253 + (t237 * t256 + (t219 * t275 - t241 * t301) * t269) * t247) * t264) * t231;
t195 = t286 * t244 * t297 + (t286 * t221 + ((-qJD(5) * t224 - 0.2e1 * t295) * t274 + (t214 * t274 + (t215 - t299) * t272) * t225) * t244) * t216;
t194 = (t200 * t322 - t210 * t245) * t298 + (t200 * t326 + t222 * t210 + (-t245 * t199 - t200 * t221 - (-t198 * t241 - t204 * t219 + t238 + (-t204 * t249 - t243) * t203) * t315 - (-t198 * t249 - t204 * t237 - t220 + (t204 * t241 - t250) * t203) * t316) * t211) * t207;
t1 = [0, t196, t198, t198, 0, 0; 0 (t201 * t322 + t210 * t312) * t298 + ((-t255 * t264 - t258 * t308) * t210 + (-t320 + t326) * t201 + (t312 * t199 - (-t196 * t241 - t205 * t219 + (-t264 * t301 + t275 * t308) * t269 + (-t205 * t249 - t256 * t264) * t203) * t315 - (-t256 * t308 - t196 * t249 - t205 * t237 + t253 * t264 + (t205 * t241 - t264 * t304) * t203) * t316) * t211) * t207, t194, t194, 0, 0; 0, 0.2e1 * (-t224 * t235 + t236 * t318) * t323 + (0.2e1 * t236 * t295 - t288 * t224 * t274 + t282 * t319 + (-t288 * t227 * t272 - t236 * t214 - t235 * t215 - t282 * t317) * t225) * t216, t195, t195, t297 + 0.2e1 * (t214 * t225 * t216 + (-t216 * t321 - t225 * t323) * t227) * t227, 0;];
JaD_rot  = t1;
