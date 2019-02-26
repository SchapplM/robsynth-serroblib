% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPRRPR1_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_jacobiaD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:29
% EndTime: 2019-02-26 19:40:30
% DurationCPUTime: 1.52s
% Computational Cost: add. (8432->112), mult. (25028->236), div. (535->12), fcn. (32914->17), ass. (0->113)
t317 = sin(pkin(12));
t318 = sin(pkin(11));
t290 = t318 * t317;
t320 = cos(pkin(12));
t321 = cos(pkin(11));
t293 = t321 * t320;
t323 = cos(pkin(6));
t261 = -t290 * t323 + t293;
t267 = sin(qJ(3));
t269 = cos(qJ(3));
t291 = t318 * t320;
t292 = t321 * t317;
t282 = t291 * t323 + t292;
t264 = sin(pkin(6));
t297 = t264 * t318;
t319 = sin(pkin(7));
t322 = cos(pkin(7));
t325 = t282 * t322 - t319 * t297;
t245 = t261 * t267 + t325 * t269;
t260 = t292 * t323 + t291;
t281 = -t293 * t323 + t290;
t298 = t264 * t319;
t277 = -t281 * t322 - t298 * t321;
t244 = t260 * t269 + t267 * t277;
t266 = sin(qJ(4));
t268 = cos(qJ(4));
t276 = -t264 * t321 * t322 + t281 * t319;
t236 = t244 * t268 + t266 * t276;
t243 = -t260 * t267 + t269 * t277;
t239 = t243 * qJD(3);
t212 = qJD(4) * t236 + t239 * t266;
t234 = t244 * t266 - t268 * t276;
t232 = t234 ^ 2;
t294 = t322 * t320;
t295 = t323 * t319;
t257 = t264 * (t267 * t294 + t269 * t317) + t267 * t295;
t259 = -t298 * t320 + t322 * t323;
t250 = t257 * t266 - t259 * t268;
t248 = 0.1e1 / t250 ^ 2;
t226 = t232 * t248 + 0.1e1;
t224 = 0.1e1 / t226;
t251 = t257 * t268 + t259 * t266;
t256 = t269 * t295 + (-t267 * t317 + t269 * t294) * t264;
t252 = t256 * qJD(3);
t230 = qJD(4) * t251 + t252 * t266;
t247 = 0.1e1 / t250;
t308 = t234 * t248;
t196 = (-t212 * t247 + t230 * t308) * t224;
t227 = atan2(-t234, t250);
t222 = sin(t227);
t223 = cos(t227);
t289 = -t222 * t250 - t223 * t234;
t192 = t196 * t289 - t212 * t222 + t223 * t230;
t206 = -t222 * t234 + t223 * t250;
t203 = 0.1e1 / t206;
t204 = 0.1e1 / t206 ^ 2;
t328 = t192 * t203 * t204;
t285 = -t243 * t247 + t256 * t308;
t327 = t266 * t285;
t246 = t261 * t269 - t267 * t325;
t278 = t282 * t319 + t297 * t322;
t237 = t246 * t266 - t268 * t278;
t326 = 0.2e1 * t237 * t328;
t309 = t230 * t247 * t248;
t324 = -0.2e1 * (t212 * t308 - t232 * t309) / t226 ^ 2;
t238 = t246 * t268 + t266 * t278;
t263 = sin(pkin(13));
t265 = cos(pkin(13));
t221 = t238 * t265 + t245 * t263;
t217 = 0.1e1 / t221;
t218 = 0.1e1 / t221 ^ 2;
t316 = t204 * t237;
t241 = t245 * qJD(3);
t215 = -qJD(4) * t237 - t241 * t268;
t242 = t246 * qJD(3);
t211 = t215 * t265 + t242 * t263;
t315 = t211 * t217 * t218;
t314 = t217 * t263;
t220 = t238 * t263 - t245 * t265;
t313 = t218 * t220;
t312 = t220 * t265;
t311 = t222 * t237;
t310 = t223 * t237;
t307 = t245 * t266;
t306 = t245 * t268;
t303 = qJD(4) * t268;
t233 = t237 ^ 2;
t202 = t204 * t233 + 0.1e1;
t214 = qJD(4) * t238 - t241 * t266;
t302 = 0.2e1 * (t214 * t316 - t233 * t328) / t202 ^ 2;
t216 = t220 ^ 2;
t209 = t216 * t218 + 0.1e1;
t210 = t215 * t263 - t242 * t265;
t301 = 0.2e1 * (t210 * t313 - t216 * t315) / t209 ^ 2;
t299 = t220 * t315;
t296 = -0.2e1 * t234 * t309;
t286 = -t236 * t247 + t251 * t308;
t284 = qJD(4) * t307 - t242 * t268;
t253 = t257 * qJD(3);
t240 = t244 * qJD(3);
t231 = -qJD(4) * t250 + t252 * t268;
t229 = t246 * t263 - t265 * t306;
t228 = -t246 * t265 - t263 * t306;
t213 = -qJD(4) * t234 + t239 * t268;
t207 = 0.1e1 / t209;
t200 = 0.1e1 / t202;
t198 = t224 * t327;
t197 = t286 * t224;
t194 = (-t222 * t243 + t223 * t256) * t266 + t289 * t198;
t193 = t197 * t289 - t222 * t236 + t223 * t251;
t190 = t286 * t324 + (t251 * t296 - t213 * t247 + (t212 * t251 + t230 * t236 + t231 * t234) * t248) * t224;
t189 = t324 * t327 + (t285 * t303 + (t256 * t296 + t240 * t247 + (t212 * t256 + t230 * t243 - t234 * t253) * t248) * t266) * t224;
t1 = [0, 0, t189, t190, 0, 0; 0, 0 (t194 * t316 + t203 * t307) * t302 + ((-t242 * t266 - t245 * t303) * t203 + t194 * t326 + (-t194 * t214 + t307 * t192 - (t256 * t303 - t189 * t234 - t198 * t212 - t253 * t266 + (-t198 * t250 - t243 * t266) * t196) * t310 - (-t243 * t303 - t189 * t250 - t198 * t230 + t240 * t266 + (t198 * t234 - t256 * t266) * t196) * t311) * t204) * t200 (t193 * t316 - t203 * t238) * t302 + (t193 * t326 + t215 * t203 + (-t238 * t192 - t193 * t214 - (-t190 * t234 - t197 * t212 + t231 + (-t197 * t250 - t236) * t196) * t310 - (-t190 * t250 - t197 * t230 - t213 + (t197 * t234 - t251) * t196) * t311) * t204) * t200, 0, 0; 0, 0 (-t217 * t228 + t229 * t313) * t301 + ((t241 * t265 + t263 * t284) * t217 + 0.2e1 * t229 * t299 + (-t228 * t211 - (-t241 * t263 + t265 * t284) * t220 - t229 * t210) * t218) * t207 (-t218 * t312 + t314) * t237 * t301 + (-0.2e1 * t237 * t265 * t299 - t214 * t314 + (t214 * t312 + (t210 * t265 + t211 * t263) * t237) * t218) * t207, 0, 0;];
JaD_rot  = t1;
