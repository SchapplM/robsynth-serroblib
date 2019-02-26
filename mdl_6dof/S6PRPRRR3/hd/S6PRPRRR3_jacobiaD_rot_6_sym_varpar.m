% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRR3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:54:49
% EndTime: 2019-02-26 19:54:50
% DurationCPUTime: 1.14s
% Computational Cost: add. (13183->114), mult. (13312->231), div. (822->12), fcn. (17117->13), ass. (0->113)
t260 = sin(pkin(11));
t262 = cos(pkin(11));
t265 = sin(qJ(2));
t263 = cos(pkin(6));
t267 = cos(qJ(2));
t294 = t263 * t267;
t248 = -t260 * t265 + t262 * t294;
t244 = t248 * qJD(2);
t295 = t263 * t265;
t249 = t260 * t267 + t262 * t295;
t258 = pkin(12) + qJ(4) + qJ(5);
t256 = sin(t258);
t259 = qJD(4) + qJD(5);
t261 = sin(pkin(6));
t298 = t261 * t262;
t283 = t256 * t298;
t257 = cos(t258);
t300 = t257 * t259;
t211 = t244 * t256 + t249 * t300 - t259 * t283;
t233 = t249 * t256 + t257 * t298;
t231 = t233 ^ 2;
t297 = t261 * t265;
t285 = t256 * t297;
t241 = -t263 * t257 + t285;
t239 = 0.1e1 / t241 ^ 2;
t219 = t231 * t239 + 0.1e1;
t217 = 0.1e1 / t219;
t292 = qJD(2) * t267;
t276 = t259 * t263 + t261 * t292;
t284 = t257 * t297;
t229 = t256 * t276 + t259 * t284;
t238 = 0.1e1 / t241;
t305 = t233 * t239;
t195 = (-t211 * t238 + t229 * t305) * t217;
t220 = atan2(-t233, t241);
t215 = sin(t220);
t216 = cos(t220);
t279 = -t215 * t241 - t216 * t233;
t191 = t195 * t279 - t215 * t211 + t216 * t229;
t205 = -t215 * t233 + t216 * t241;
t202 = 0.1e1 / t205;
t203 = 0.1e1 / t205 ^ 2;
t319 = t191 * t202 * t203;
t286 = t260 * t295;
t251 = t262 * t267 - t286;
t299 = t260 * t261;
t236 = t251 * t256 - t257 * t299;
t318 = 0.2e1 * t236 * t319;
t296 = t261 * t267;
t275 = -t238 * t248 + t296 * t305;
t317 = t256 * t275;
t306 = t229 * t238 * t239;
t316 = -0.2e1 * (t211 * t305 - t231 * t306) / t219 ^ 2;
t237 = t251 * t257 + t256 * t299;
t266 = cos(qJ(6));
t250 = t260 * t294 + t262 * t265;
t264 = sin(qJ(6));
t303 = t250 * t264;
t226 = t237 * t266 + t303;
t222 = 0.1e1 / t226;
t223 = 0.1e1 / t226 ^ 2;
t246 = t250 * qJD(2);
t281 = t259 * t299 - t246;
t301 = t256 * t259;
t214 = -t251 * t301 + t257 * t281;
t247 = -qJD(2) * t286 + t262 * t292;
t206 = qJD(6) * t226 + t214 * t264 - t247 * t266;
t302 = t250 * t266;
t225 = t237 * t264 - t302;
t221 = t225 ^ 2;
t210 = t221 * t223 + 0.1e1;
t308 = t223 * t225;
t291 = qJD(6) * t225;
t207 = t214 * t266 + t247 * t264 - t291;
t313 = t207 * t222 * t223;
t315 = (t206 * t308 - t221 * t313) / t210 ^ 2;
t314 = t203 * t236;
t213 = t251 * t300 + t256 * t281;
t312 = t213 * t203;
t311 = t215 * t236;
t310 = t216 * t236;
t309 = t222 * t264;
t307 = t225 * t266;
t304 = t250 * t256;
t293 = qJD(2) * t265;
t232 = t236 ^ 2;
t201 = t232 * t203 + 0.1e1;
t290 = 0.2e1 * (-t232 * t319 + t236 * t312) / t201 ^ 2;
t289 = -0.2e1 * t315;
t287 = t225 * t313;
t282 = -0.2e1 * t233 * t306;
t280 = qJD(6) * t250 * t257 - t246;
t278 = t223 * t307 - t309;
t235 = t249 * t257 - t283;
t242 = t263 * t256 + t284;
t277 = -t235 * t238 + t242 * t305;
t274 = qJD(6) * t251 - t247 * t257 + t250 * t301;
t245 = t249 * qJD(2);
t230 = t257 * t276 - t259 * t285;
t228 = t251 * t264 - t257 * t302;
t227 = -t251 * t266 - t257 * t303;
t212 = -t249 * t301 + (-t259 * t298 + t244) * t257;
t208 = 0.1e1 / t210;
t199 = 0.1e1 / t201;
t197 = t217 * t317;
t196 = t277 * t217;
t193 = (-t215 * t248 + t216 * t296) * t256 + t279 * t197;
t192 = t196 * t279 - t215 * t235 + t216 * t242;
t189 = t277 * t316 + (t242 * t282 - t212 * t238 + (t211 * t242 + t229 * t235 + t230 * t233) * t239) * t217;
t188 = t316 * t317 + (t275 * t300 + (t282 * t296 + t238 * t245 + (t229 * t248 + (t211 * t267 - t233 * t293) * t261) * t239) * t256) * t217;
t187 = t278 * t236 * t289 + (t278 * t213 + ((-qJD(6) * t222 - 0.2e1 * t287) * t266 + (t206 * t266 + (t207 - t291) * t264) * t223) * t236) * t208;
t186 = (t192 * t314 - t202 * t237) * t290 + (t192 * t318 + t214 * t202 + (-t237 * t191 - t192 * t213 - (-t189 * t233 - t196 * t211 + t230 + (-t196 * t241 - t235) * t195) * t310 - (-t189 * t241 - t196 * t229 - t212 + (t196 * t233 - t242) * t195) * t311) * t203) * t199;
t1 = [0, t188, 0, t189, t189, 0; 0 (t193 * t314 + t202 * t304) * t290 + ((-t247 * t256 - t250 * t300) * t202 + (-t312 + t318) * t193 + (t304 * t191 - (-t188 * t233 - t197 * t211 + (-t256 * t293 + t267 * t300) * t261 + (-t197 * t241 - t248 * t256) * t195) * t310 - (-t248 * t300 - t188 * t241 - t197 * t229 + t245 * t256 + (t197 * t233 - t256 * t296) * t195) * t311) * t203) * t199, 0, t186, t186, 0; 0, 0.2e1 * (-t222 * t227 + t228 * t308) * t315 + (0.2e1 * t228 * t287 - t280 * t222 * t266 + t274 * t309 + (-t225 * t264 * t280 - t228 * t206 - t227 * t207 - t274 * t307) * t223) * t208, 0, t187, t187, t289 + 0.2e1 * (t206 * t223 * t208 + (-t208 * t313 - t223 * t315) * t225) * t225;];
JaD_rot  = t1;
