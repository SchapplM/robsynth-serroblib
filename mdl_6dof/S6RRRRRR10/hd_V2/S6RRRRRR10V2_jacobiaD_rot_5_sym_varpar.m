% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR10V2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRR10V2_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10V2_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_jacobiaD_rot_5_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:56:32
% EndTime: 2019-04-11 14:56:34
% DurationCPUTime: 1.63s
% Computational Cost: add. (6566->153), mult. (9976->322), div. (1522->14), fcn. (12462->11), ass. (0->136)
t240 = qJ(2) + qJ(3);
t235 = cos(t240);
t242 = sin(qJ(4));
t330 = sin(qJ(1));
t276 = t330 * t242;
t244 = cos(qJ(4));
t245 = cos(qJ(1));
t302 = t245 * t244;
t219 = t235 * t302 + t276;
t241 = sin(qJ(5));
t243 = cos(qJ(5));
t234 = sin(t240);
t310 = t234 * t245;
t203 = t219 * t241 - t243 * t310;
t334 = 0.2e1 * t203;
t236 = qJD(2) + qJD(3);
t333 = -qJD(5) * t244 + t236;
t204 = t219 * t243 + t241 * t310;
t198 = 0.1e1 / t204;
t199 = 0.1e1 / t204 ^ 2;
t322 = t199 * t203;
t261 = t241 * t198 - t243 * t322;
t215 = t235 * t276 + t302;
t270 = t330 * qJD(4);
t264 = t242 * t270;
t298 = qJD(4) * t245;
t273 = t244 * t298;
t306 = t236 * t245;
t281 = t234 * t306;
t193 = t215 * qJD(1) - t235 * t273 + t242 * t281 - t264;
t275 = t330 * t244;
t303 = t245 * t242;
t218 = t235 * t303 - t275;
t231 = 0.1e1 / t234;
t232 = 0.1e1 / t234 ^ 2;
t238 = 0.1e1 / t242 ^ 2;
t299 = qJD(4) * t244;
t274 = t238 * t299;
t237 = 0.1e1 / t242;
t308 = t236 * t237;
t280 = t235 * t308;
t314 = t231 * t237;
t332 = (t231 * t274 + t232 * t280) * t218 + t193 * t314;
t312 = t234 * t242;
t209 = atan2(-t215, t312);
t206 = cos(t209);
t205 = sin(t209);
t320 = t205 * t215;
t189 = t206 * t312 - t320;
t186 = 0.1e1 / t189;
t187 = 0.1e1 / t189 ^ 2;
t331 = 0.2e1 * t218;
t213 = t215 ^ 2;
t313 = t232 * t238;
t210 = t213 * t313 + 0.1e1;
t207 = 0.1e1 / t210;
t309 = t235 * t236;
t258 = t234 * t299 + t242 * t309;
t283 = t215 * t313;
t271 = qJD(1) * t330;
t301 = qJD(1) * t245;
t195 = -t236 * t234 * t276 - t242 * t298 - t244 * t271 + (t242 * t301 + t244 * t270) * t235;
t286 = t195 * t314;
t177 = (t258 * t283 - t286) * t207;
t256 = -t177 * t215 + t258;
t172 = (-t177 * t312 - t195) * t205 + t256 * t206;
t188 = t186 * t187;
t329 = t172 * t188;
t255 = qJD(5) * t219 + t234 * t271 - t235 * t306;
t194 = (-qJD(4) * t235 + qJD(1)) * t303 + (-t235 * t271 + t270 - t281) * t244;
t263 = qJD(5) * t310 + t194;
t179 = t263 * t241 + t255 * t243;
t197 = t203 ^ 2;
t192 = t197 * t199 + 0.1e1;
t180 = -t255 * t241 + t263 * t243;
t200 = t198 * t199;
t326 = t180 * t200;
t328 = (t179 * t322 - t197 * t326) / t192 ^ 2;
t233 = t231 * t232;
t239 = t237 * t238;
t327 = (t195 * t283 + (-t232 * t239 * t299 - t233 * t238 * t309) * t213) / t210 ^ 2;
t325 = t187 * t218;
t190 = 0.1e1 / t192;
t324 = t190 * t199;
t323 = t193 * t187;
t311 = t234 * t244;
t212 = (t235 * t241 - t243 * t311) * t245;
t321 = t203 * t212;
t319 = t205 * t218;
t318 = t205 * t234;
t317 = t206 * t215;
t316 = t206 * t218;
t315 = t206 * t235;
t307 = t236 * t244;
t305 = t238 * t244;
t300 = qJD(4) * t242;
t214 = t218 ^ 2;
t184 = t214 * t187 + 0.1e1;
t296 = 0.2e1 * (-t214 * t329 - t218 * t323) / t184 ^ 2;
t295 = -0.2e1 * t328;
t294 = 0.2e1 * t328;
t293 = -0.2e1 * t327;
t292 = t188 * t331;
t291 = t199 * t328;
t290 = t231 * t327;
t289 = t179 * t324;
t288 = t187 * t319;
t285 = t203 * t326;
t284 = t215 * t314;
t282 = t232 * t235 * t237;
t278 = t330 * t234;
t277 = t330 * t235;
t259 = t215 * t282 + t330;
t185 = t259 * t207;
t272 = t330 - t185;
t269 = t186 * t296;
t268 = t187 * t296;
t266 = t237 * t290;
t265 = t234 * t275;
t196 = t219 * qJD(1) - t235 * t264 - t236 * t265 - t273;
t262 = -qJD(5) * t278 - t196;
t217 = t235 * t275 - t303;
t260 = t215 * t305 - t217 * t237;
t254 = -qJD(5) * t217 + t234 * t301 + t236 * t277;
t211 = (-t235 * t243 - t241 * t311) * t245;
t202 = -t217 * t243 - t241 * t278;
t182 = 0.1e1 / t184;
t181 = t260 * t231 * t207;
t176 = (-t205 + (t206 * t284 + t205) * t207) * t218;
t175 = -t185 * t317 + (t272 * t318 + t315) * t242;
t173 = t206 * t311 - t205 * t217 + (-t205 * t312 - t317) * t181;
t171 = t259 * t293 + (t195 * t282 + t301 + (-t231 * t308 + (-t232 * t274 - 0.2e1 * t233 * t280) * t235) * t215) * t207;
t169 = -0.2e1 * t260 * t290 + (-t260 * t232 * t309 + (t195 * t305 - t196 * t237 + (t217 * t305 + (-0.2e1 * t239 * t244 ^ 2 - t237) * t215) * qJD(4)) * t231) * t207;
t168 = (-t198 * t211 + t199 * t321) * t294 + (-t212 * t179 * t199 + (-t211 * t199 + 0.2e1 * t200 * t321) * t180 + ((t241 * t265 + t243 * t277) * t198 - (-t241 * t277 + t243 * t265) * t322) * qJD(1) + (((t241 * t300 + t243 * t333) * t198 - (-t241 * t333 + t243 * t300) * t322) * t234 + t261 * t235 * (qJD(5) - t307)) * t245) * t190;
t167 = t175 * t218 * t268 + (-(-t171 * t317 + (t177 * t320 - t195 * t206) * t185) * t325 + (t172 * t292 + t323) * t175 + (-t186 * t310 - (-t185 * t318 + t205 * t278 + t315) * t325) * t299) * t182 + (t269 * t310 + ((-t186 * t306 - (t272 * t236 - t177) * t288) * t235 + (t186 * t271 + (t245 * t172 - (-t171 + t301) * t319 - (t272 * t177 - t236) * t316) * t187) * t234) * t182) * t242;
t1 = [t207 * t332 + t266 * t331, t171, t171, t169, 0, 0; t215 * t269 + (-t195 * t186 + (t172 * t215 + t176 * t193) * t187) * t182 + (t176 * t268 + (0.2e1 * t176 * t329 + (t193 * t207 - t193 - (-t177 * t207 * t284 + t293) * t218) * t187 * t205 + (-(-0.2e1 * t215 * t266 - t177) * t325 + (-(t177 + t286) * t218 + t332 * t215) * t187 * t207) * t206) * t182) * t218, t167, t167 (t173 * t325 - t186 * t219) * t296 + (t173 * t323 + t194 * t186 + (t173 * t292 - t219 * t187) * t172 - (-t234 * t300 + t235 * t307 - t169 * t215 - t181 * t195 + (-t181 * t312 - t217) * t177) * t187 * t316 - (-t196 + (-t169 * t242 - t177 * t244) * t234 - t256 * t181) * t288) * t182, 0, 0; (t291 * t334 - t289) * t202 + (-t180 * t324 + t198 * t295) * (-t217 * t241 + t243 * t278) + ((t262 * t241 + t254 * t243) * t198 - (-t254 * t241 + t262 * t243) * t322 + 0.2e1 * t202 * t285) * t190, t168, t168, t261 * t218 * t294 + (t261 * t193 + ((-qJD(5) * t198 - 0.2e1 * t285) * t243 + (t179 * t243 + (-qJD(5) * t203 + t180) * t241) * t199) * t218) * t190, t295 + (t289 + (-t190 * t326 - t291) * t203) * t334, 0;];
JaD_rot  = t1;
