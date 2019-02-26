% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR4
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
% Datum: 2019-02-26 19:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRR4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:55:23
% EndTime: 2019-02-26 19:55:24
% DurationCPUTime: 0.98s
% Computational Cost: add. (6552->112), mult. (9671->227), div. (577->12), fcn. (12382->13), ass. (0->111)
t254 = sin(pkin(11));
t256 = cos(pkin(11));
t258 = sin(qJ(2));
t257 = cos(pkin(6));
t259 = cos(qJ(2));
t286 = t257 * t259;
t239 = -t254 * t258 + t256 * t286;
t235 = t239 * qJD(2);
t287 = t257 * t258;
t240 = t254 * t259 + t256 * t287;
t251 = pkin(12) + qJ(4);
t247 = sin(t251);
t255 = sin(pkin(6));
t290 = t255 * t256;
t277 = t247 * t290;
t248 = cos(t251);
t283 = qJD(4) * t248;
t206 = -qJD(4) * t277 + t235 * t247 + t240 * t283;
t224 = t240 * t247 + t248 * t290;
t222 = t224 ^ 2;
t289 = t255 * t258;
t233 = t247 * t289 - t257 * t248;
t231 = 0.1e1 / t233 ^ 2;
t216 = t222 * t231 + 0.1e1;
t214 = 0.1e1 / t216;
t234 = t257 * t247 + t248 * t289;
t284 = qJD(2) * t259;
t276 = t255 * t284;
t220 = t234 * qJD(4) + t247 * t276;
t230 = 0.1e1 / t233;
t296 = t224 * t231;
t186 = (-t206 * t230 + t220 * t296) * t214;
t217 = atan2(-t224, t233);
t212 = sin(t217);
t213 = cos(t217);
t271 = -t212 * t233 - t213 * t224;
t182 = t271 * t186 - t212 * t206 + t213 * t220;
t196 = -t212 * t224 + t213 * t233;
t193 = 0.1e1 / t196;
t194 = 0.1e1 / t196 ^ 2;
t310 = t182 * t193 * t194;
t278 = t254 * t287;
t242 = t256 * t259 - t278;
t291 = t254 * t255;
t268 = -t242 * t247 + t248 * t291;
t309 = -0.2e1 * t268 * t310;
t288 = t255 * t259;
t267 = -t230 * t239 + t288 * t296;
t308 = t247 * t267;
t297 = t220 * t230 * t231;
t307 = -0.2e1 * (t206 * t296 - t222 * t297) / t216 ^ 2;
t228 = t242 * t248 + t247 * t291;
t253 = qJ(5) + qJ(6);
t250 = cos(t253);
t241 = t254 * t286 + t256 * t258;
t249 = sin(t253);
t294 = t241 * t249;
t211 = t228 * t250 + t294;
t203 = 0.1e1 / t211;
t204 = 0.1e1 / t211 ^ 2;
t238 = -qJD(2) * t278 + t256 * t284;
t252 = qJD(5) + qJD(6);
t273 = t228 * t252 - t238;
t237 = t241 * qJD(2);
t209 = t268 * qJD(4) - t237 * t248;
t292 = t241 * t252;
t274 = t209 + t292;
t197 = t274 * t249 + t273 * t250;
t293 = t241 * t250;
t210 = t228 * t249 - t293;
t202 = t210 ^ 2;
t201 = t202 * t204 + 0.1e1;
t302 = t204 * t210;
t198 = -t273 * t249 + t274 * t250;
t304 = t198 * t203 * t204;
t306 = (t197 * t302 - t202 * t304) / t201 ^ 2;
t305 = t194 * t268;
t303 = t203 * t249;
t208 = t228 * qJD(4) - t237 * t247;
t301 = t208 * t194;
t300 = t210 * t250;
t299 = t212 * t268;
t298 = t213 * t268;
t295 = t241 * t247;
t285 = qJD(2) * t258;
t223 = t268 ^ 2;
t192 = t223 * t194 + 0.1e1;
t282 = 0.2e1 * (-t223 * t310 - t268 * t301) / t192 ^ 2;
t281 = -0.2e1 * t306;
t279 = t210 * t304;
t275 = -0.2e1 * t224 * t297;
t272 = t248 * t292 - t237;
t270 = t204 * t300 - t303;
t226 = t240 * t248 - t277;
t269 = -t226 * t230 + t234 * t296;
t266 = qJD(4) * t295 - t238 * t248 + t242 * t252;
t236 = t240 * qJD(2);
t221 = -t233 * qJD(4) + t248 * t276;
t219 = t242 * t249 - t248 * t293;
t218 = -t242 * t250 - t248 * t294;
t207 = -t224 * qJD(4) + t235 * t248;
t199 = 0.1e1 / t201;
t189 = 0.1e1 / t192;
t188 = t214 * t308;
t187 = t269 * t214;
t184 = (-t212 * t239 + t213 * t288) * t247 + t271 * t188;
t183 = t271 * t187 - t212 * t226 + t213 * t234;
t181 = t269 * t307 + (t234 * t275 - t207 * t230 + (t206 * t234 + t220 * t226 + t221 * t224) * t231) * t214;
t179 = t307 * t308 + (t267 * t283 + (t275 * t288 + t230 * t236 + (t220 * t239 + (t206 * t259 - t224 * t285) * t255) * t231) * t247) * t214;
t178 = t281 + 0.2e1 * (t197 * t204 * t199 + (-t199 * t304 - t204 * t306) * t210) * t210;
t1 = [0, t179, 0, t181, 0, 0; 0 (-t184 * t305 + t193 * t295) * t282 + ((-t238 * t247 - t241 * t283) * t193 + (-t301 + t309) * t184 + (t295 * t182 + (-t179 * t224 - t188 * t206 + (-t247 * t285 + t259 * t283) * t255 + (-t188 * t233 - t239 * t247) * t186) * t298 + (-t239 * t283 - t179 * t233 - t188 * t220 + t236 * t247 + (t188 * t224 - t247 * t288) * t186) * t299) * t194) * t189, 0 (-t183 * t305 - t193 * t228) * t282 + (t183 * t309 + t209 * t193 + (-t228 * t182 - t183 * t208 + (-t181 * t224 - t187 * t206 + t221 + (-t187 * t233 - t226) * t186) * t298 + (-t181 * t233 - t187 * t220 - t207 + (t187 * t224 - t234) * t186) * t299) * t194) * t189, 0, 0; 0, 0.2e1 * (-t203 * t218 + t219 * t302) * t306 + (0.2e1 * t219 * t279 - t272 * t203 * t250 + t266 * t303 + (-t272 * t210 * t249 - t219 * t197 - t218 * t198 - t266 * t300) * t204) * t199, 0, -t270 * t268 * t281 + (t270 * t208 - ((-t203 * t252 - 0.2e1 * t279) * t250 + (t197 * t250 + (-t210 * t252 + t198) * t249) * t204) * t268) * t199, t178, t178;];
JaD_rot  = t1;
