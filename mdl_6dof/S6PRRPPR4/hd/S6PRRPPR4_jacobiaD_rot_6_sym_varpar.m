% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPPR4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:00:00
% EndTime: 2019-02-26 20:00:01
% DurationCPUTime: 1.13s
% Computational Cost: add. (3755->120), mult. (11172->251), div. (553->12), fcn. (14363->15), ass. (0->113)
t258 = sin(pkin(10));
t261 = cos(pkin(10));
t265 = sin(qJ(2));
t262 = cos(pkin(6));
t268 = cos(qJ(2));
t295 = t262 * t268;
t248 = -t258 * t265 + t261 * t295;
t241 = t248 * qJD(2);
t296 = t262 * t265;
t249 = t258 * t268 + t261 * t296;
t264 = sin(qJ(3));
t259 = sin(pkin(6));
t299 = t259 * t264;
t286 = t261 * t299;
t267 = cos(qJ(3));
t292 = qJD(3) * t267;
t219 = -qJD(3) * t286 + t241 * t264 + t249 * t292;
t298 = t259 * t267;
t235 = t249 * t264 + t261 * t298;
t233 = t235 ^ 2;
t277 = -t262 * t267 + t265 * t299;
t246 = 0.1e1 / t277 ^ 2;
t229 = t233 * t246 + 0.1e1;
t227 = 0.1e1 / t229;
t253 = -t262 * t264 - t265 * t298;
t293 = qJD(2) * t268;
t285 = t259 * t293;
t239 = t253 * qJD(3) - t264 * t285;
t245 = 0.1e1 / t277;
t304 = t235 * t246;
t190 = (-t219 * t245 - t239 * t304) * t227;
t230 = atan2(t235, -t277);
t225 = sin(t230);
t226 = cos(t230);
t282 = t225 * t277 + t226 * t235;
t185 = t282 * t190 + t225 * t219 + t226 * t239;
t203 = t225 * t235 - t226 * t277;
t200 = 0.1e1 / t203;
t201 = 0.1e1 / t203 ^ 2;
t315 = t185 * t200 * t201;
t287 = t258 * t296;
t251 = t261 * t268 - t287;
t237 = -t251 * t264 + t258 * t298;
t314 = 0.2e1 * t237 * t315;
t302 = t239 * t245 * t246;
t313 = (t219 * t304 + t233 * t302) / t229 ^ 2;
t297 = t259 * t268;
t288 = t235 * t297;
t275 = -t245 * t248 + t246 * t288;
t312 = t264 * t275;
t257 = sin(pkin(11));
t260 = cos(pkin(11));
t263 = sin(qJ(6));
t266 = cos(qJ(6));
t280 = t257 * t266 - t260 * t263;
t216 = t280 * t237;
t238 = t251 * t267 + t258 * t299;
t250 = t258 * t295 + t261 * t265;
t223 = t238 * t257 - t250 * t260;
t224 = t238 * t260 + t250 * t257;
t209 = t223 * t263 + t224 * t266;
t205 = 0.1e1 / t209;
t206 = 0.1e1 / t209 ^ 2;
t243 = t250 * qJD(2);
t222 = t237 * qJD(3) - t243 * t267;
t244 = -qJD(2) * t287 + t261 * t293;
t212 = t222 * t257 - t244 * t260;
t213 = t222 * t260 + t244 * t257;
t188 = t209 * qJD(6) - t212 * t266 + t213 * t263;
t283 = t223 * t266 - t224 * t263;
t204 = t283 ^ 2;
t193 = t204 * t206 + 0.1e1;
t308 = t206 * t283;
t189 = t283 * qJD(6) + t212 * t263 + t213 * t266;
t310 = t189 * t205 * t206;
t311 = (-t188 * t308 - t204 * t310) / t193 ^ 2;
t309 = t201 * t237;
t221 = -t238 * qJD(3) + t243 * t264;
t307 = t221 * t201;
t306 = t225 * t237;
t305 = t226 * t237;
t303 = t235 * t253;
t301 = t250 * t264;
t300 = t250 * t267;
t294 = qJD(2) * t265;
t234 = t237 ^ 2;
t199 = t234 * t201 + 0.1e1;
t291 = 0.2e1 * (-t234 * t315 + t237 * t307) / t199 ^ 2;
t290 = 0.2e1 * t311;
t284 = -0.2e1 * t283 * t310;
t231 = -t251 * t260 - t257 * t300;
t232 = t251 * t257 - t260 * t300;
t281 = t231 * t266 - t232 * t263;
t211 = t231 * t263 + t232 * t266;
t279 = t257 * t263 + t260 * t266;
t236 = t249 * t267 - t286;
t278 = t236 * t245 + t246 * t303;
t276 = qJD(3) * t301 - t244 * t267;
t217 = t279 * t237;
t242 = t249 * qJD(2);
t240 = t277 * qJD(3) - t267 * t285;
t220 = -t235 * qJD(3) + t241 * t267;
t215 = -t243 * t257 + t276 * t260;
t214 = t243 * t260 + t276 * t257;
t196 = 0.1e1 / t199;
t195 = t227 * t312;
t194 = t278 * t227;
t191 = 0.1e1 / t193;
t187 = (t225 * t248 - t226 * t297) * t264 + t282 * t195;
t186 = -t282 * t194 + t225 * t236 + t226 * t253;
t183 = 0.2e1 * t278 * t313 + (-0.2e1 * t302 * t303 - t220 * t245 + (-t219 * t253 - t235 * t240 - t236 * t239) * t246) * t227;
t181 = -0.2e1 * t312 * t313 + (t275 * t292 + (0.2e1 * t288 * t302 + t242 * t245 + (-t239 * t248 + (t219 * t268 - t235 * t294) * t259) * t246) * t264) * t227;
t1 = [0, t181, t183, 0, 0, 0; 0 (t187 * t309 - t200 * t301) * t291 + ((t244 * t264 + t250 * t292) * t200 + (-t307 + t314) * t187 + (-t301 * t185 - (t181 * t235 + t195 * t219 + (t264 * t294 - t268 * t292) * t259 + (t195 * t277 + t248 * t264) * t190) * t305 - (t248 * t292 + t181 * t277 - t195 * t239 - t242 * t264 + (-t195 * t235 + t264 * t297) * t190) * t306) * t201) * t196 (t186 * t309 + t200 * t238) * t291 + (t186 * t314 - t222 * t200 + (t238 * t185 - t186 * t221 - (t183 * t235 - t194 * t219 + t240 + (-t194 * t277 + t236) * t190) * t305 - (t183 * t277 + t194 * t239 + t220 + (t194 * t235 - t253) * t190) * t306) * t201) * t196, 0, 0, 0; 0 (t205 * t281 - t211 * t308) * t290 + ((t211 * qJD(6) - t214 * t266 + t215 * t263) * t205 + t211 * t284 + (t281 * t189 + (t281 * qJD(6) + t214 * t263 + t215 * t266) * t283 - t211 * t188) * t206) * t191 (t205 * t216 - t217 * t308) * t290 + ((qJD(6) * t217 - t280 * t221) * t205 + t217 * t284 + (t216 * t189 + (qJD(6) * t216 + t279 * t221) * t283 - t217 * t188) * t206) * t191, 0, 0, -0.2e1 * t311 - 0.2e1 * (t188 * t206 * t191 - (-t191 * t310 - t206 * t311) * t283) * t283;];
JaD_rot  = t1;
