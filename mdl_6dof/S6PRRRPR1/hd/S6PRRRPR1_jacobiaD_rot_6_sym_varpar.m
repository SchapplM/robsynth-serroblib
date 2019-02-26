% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRPR1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:10:39
% EndTime: 2019-02-26 20:10:40
% DurationCPUTime: 1.19s
% Computational Cost: add. (13183->114), mult. (13312->231), div. (822->12), fcn. (17117->13), ass. (0->113)
t263 = sin(pkin(11));
t265 = cos(pkin(11));
t268 = sin(qJ(2));
t266 = cos(pkin(6));
t270 = cos(qJ(2));
t297 = t266 * t270;
t251 = -t263 * t268 + t265 * t297;
t247 = t251 * qJD(2);
t298 = t266 * t268;
t252 = t263 * t270 + t265 * t298;
t261 = qJ(3) + qJ(4) + pkin(12);
t259 = sin(t261);
t262 = qJD(3) + qJD(4);
t264 = sin(pkin(6));
t301 = t264 * t265;
t286 = t259 * t301;
t260 = cos(t261);
t303 = t260 * t262;
t214 = t247 * t259 + t252 * t303 - t262 * t286;
t236 = t252 * t259 + t260 * t301;
t234 = t236 ^ 2;
t300 = t264 * t268;
t288 = t259 * t300;
t244 = -t266 * t260 + t288;
t242 = 0.1e1 / t244 ^ 2;
t222 = t234 * t242 + 0.1e1;
t220 = 0.1e1 / t222;
t295 = qJD(2) * t270;
t279 = t262 * t266 + t264 * t295;
t287 = t260 * t300;
t232 = t279 * t259 + t262 * t287;
t241 = 0.1e1 / t244;
t308 = t236 * t242;
t198 = (-t214 * t241 + t232 * t308) * t220;
t223 = atan2(-t236, t244);
t218 = sin(t223);
t219 = cos(t223);
t282 = -t218 * t244 - t219 * t236;
t194 = t282 * t198 - t218 * t214 + t219 * t232;
t208 = -t218 * t236 + t219 * t244;
t205 = 0.1e1 / t208;
t206 = 0.1e1 / t208 ^ 2;
t322 = t194 * t205 * t206;
t289 = t263 * t298;
t254 = t265 * t270 - t289;
t302 = t263 * t264;
t239 = t254 * t259 - t260 * t302;
t321 = 0.2e1 * t239 * t322;
t299 = t264 * t270;
t278 = -t241 * t251 + t299 * t308;
t320 = t259 * t278;
t309 = t232 * t241 * t242;
t319 = -0.2e1 * (t214 * t308 - t234 * t309) / t222 ^ 2;
t240 = t254 * t260 + t259 * t302;
t269 = cos(qJ(6));
t253 = t263 * t297 + t265 * t268;
t267 = sin(qJ(6));
t306 = t253 * t267;
t229 = t240 * t269 + t306;
t225 = 0.1e1 / t229;
t226 = 0.1e1 / t229 ^ 2;
t249 = t253 * qJD(2);
t284 = t262 * t302 - t249;
t304 = t259 * t262;
t217 = -t254 * t304 + t284 * t260;
t250 = -qJD(2) * t289 + t265 * t295;
t209 = t229 * qJD(6) + t217 * t267 - t250 * t269;
t305 = t253 * t269;
t228 = t240 * t267 - t305;
t224 = t228 ^ 2;
t213 = t224 * t226 + 0.1e1;
t311 = t226 * t228;
t294 = qJD(6) * t228;
t210 = t217 * t269 + t250 * t267 - t294;
t316 = t210 * t225 * t226;
t318 = (t209 * t311 - t224 * t316) / t213 ^ 2;
t317 = t206 * t239;
t216 = t254 * t303 + t284 * t259;
t315 = t216 * t206;
t314 = t218 * t239;
t313 = t219 * t239;
t312 = t225 * t267;
t310 = t228 * t269;
t307 = t253 * t259;
t296 = qJD(2) * t268;
t235 = t239 ^ 2;
t204 = t235 * t206 + 0.1e1;
t293 = 0.2e1 * (-t235 * t322 + t239 * t315) / t204 ^ 2;
t292 = -0.2e1 * t318;
t290 = t228 * t316;
t285 = -0.2e1 * t236 * t309;
t283 = qJD(6) * t253 * t260 - t249;
t281 = t226 * t310 - t312;
t238 = t252 * t260 - t286;
t245 = t266 * t259 + t287;
t280 = -t238 * t241 + t245 * t308;
t277 = qJD(6) * t254 - t250 * t260 + t253 * t304;
t248 = t252 * qJD(2);
t233 = t279 * t260 - t262 * t288;
t231 = t254 * t267 - t260 * t305;
t230 = -t254 * t269 - t260 * t306;
t215 = -t252 * t304 + (-t262 * t301 + t247) * t260;
t211 = 0.1e1 / t213;
t202 = 0.1e1 / t204;
t200 = t220 * t320;
t199 = t280 * t220;
t196 = (-t218 * t251 + t219 * t299) * t259 + t282 * t200;
t195 = t282 * t199 - t218 * t238 + t219 * t245;
t192 = t280 * t319 + (t245 * t285 - t215 * t241 + (t214 * t245 + t232 * t238 + t233 * t236) * t242) * t220;
t191 = t319 * t320 + (t278 * t303 + (t285 * t299 + t241 * t248 + (t232 * t251 + (t214 * t270 - t236 * t296) * t264) * t242) * t259) * t220;
t190 = t281 * t239 * t292 + (t281 * t216 + ((-qJD(6) * t225 - 0.2e1 * t290) * t269 + (t209 * t269 + (t210 - t294) * t267) * t226) * t239) * t211;
t189 = (t195 * t317 - t205 * t240) * t293 + (t195 * t321 + t217 * t205 + (-t240 * t194 - t195 * t216 - (-t192 * t236 - t199 * t214 + t233 + (-t199 * t244 - t238) * t198) * t313 - (-t192 * t244 - t199 * t232 - t215 + (t199 * t236 - t245) * t198) * t314) * t206) * t202;
t1 = [0, t191, t192, t192, 0, 0; 0 (t196 * t317 + t205 * t307) * t293 + ((-t250 * t259 - t253 * t303) * t205 + (-t315 + t321) * t196 + (t307 * t194 - (-t191 * t236 - t200 * t214 + (-t259 * t296 + t270 * t303) * t264 + (-t200 * t244 - t251 * t259) * t198) * t313 - (-t251 * t303 - t191 * t244 - t200 * t232 + t248 * t259 + (t200 * t236 - t259 * t299) * t198) * t314) * t206) * t202, t189, t189, 0, 0; 0, 0.2e1 * (-t225 * t230 + t231 * t311) * t318 + (0.2e1 * t231 * t290 - t283 * t225 * t269 + t277 * t312 + (-t283 * t228 * t267 - t231 * t209 - t230 * t210 - t277 * t310) * t226) * t211, t190, t190, 0, t292 + 0.2e1 * (t209 * t226 * t211 + (-t211 * t316 - t226 * t318) * t228) * t228;];
JaD_rot  = t1;
