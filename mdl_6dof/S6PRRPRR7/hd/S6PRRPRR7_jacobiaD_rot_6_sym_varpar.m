% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRR7_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:37
% EndTime: 2019-02-26 20:07:38
% DurationCPUTime: 1.00s
% Computational Cost: add. (3659->109), mult. (9671->221), div. (577->12), fcn. (12382->13), ass. (0->107)
t245 = cos(pkin(6));
t249 = cos(qJ(2));
t299 = cos(pkin(11));
t271 = t299 * t249;
t243 = sin(pkin(11));
t247 = sin(qJ(2));
t285 = t243 * t247;
t230 = t245 * t271 - t285;
t223 = t230 * qJD(2);
t248 = cos(qJ(3));
t246 = sin(qJ(3));
t272 = t299 * t247;
t284 = t243 * t249;
t258 = -t245 * t272 - t284;
t244 = sin(pkin(6));
t273 = t244 * t299;
t301 = t258 * t246 - t248 * t273;
t201 = t301 * qJD(3) + t223 * t248;
t216 = -t246 * t273 - t258 * t248;
t213 = t216 ^ 2;
t282 = t244 * t248;
t234 = t245 * t246 + t247 * t282;
t228 = 0.1e1 / t234 ^ 2;
t208 = t213 * t228 + 0.1e1;
t206 = 0.1e1 / t208;
t283 = t244 * t246;
t233 = t245 * t248 - t247 * t283;
t281 = t244 * t249;
t274 = qJD(2) * t281;
t221 = t233 * qJD(3) + t248 * t274;
t227 = 0.1e1 / t234;
t289 = t216 * t228;
t178 = (-t201 * t227 + t221 * t289) * t206;
t209 = atan2(-t216, t234);
t204 = sin(t209);
t205 = cos(t209);
t264 = -t204 * t234 - t205 * t216;
t174 = t264 * t178 - t201 * t204 + t205 * t221;
t190 = -t204 * t216 + t205 * t234;
t187 = 0.1e1 / t190;
t188 = 0.1e1 / t190 ^ 2;
t304 = t174 * t187 * t188;
t232 = -t245 * t285 + t271;
t219 = t232 * t248 + t243 * t283;
t303 = 0.2e1 * t219 * t304;
t260 = -t227 * t230 + t281 * t289;
t302 = t248 * t260;
t288 = t221 * t227 * t228;
t300 = -0.2e1 * (t201 * t289 - t213 * t288) / t208 ^ 2;
t242 = qJ(5) + qJ(6);
t239 = sin(t242);
t240 = cos(t242);
t259 = -t245 * t284 - t272;
t261 = -t232 * t246 + t243 * t282;
t199 = -t239 * t261 - t240 * t259;
t195 = 0.1e1 / t199;
t196 = 0.1e1 / t199 ^ 2;
t226 = t232 * qJD(2);
t241 = qJD(5) + qJD(6);
t267 = -t241 * t261 + t226;
t225 = t259 * qJD(2);
t202 = t219 * qJD(3) + t225 * t246;
t268 = -t241 * t259 - t202;
t185 = t267 * t239 + t268 * t240;
t198 = -t239 * t259 + t240 * t261;
t194 = t198 ^ 2;
t193 = t194 * t196 + 0.1e1;
t294 = t196 * t198;
t186 = -t268 * t239 + t267 * t240;
t297 = t186 * t195 * t196;
t298 = (t185 * t294 - t194 * t297) / t193 ^ 2;
t296 = t188 * t219;
t295 = t195 * t240;
t293 = t198 * t239;
t203 = t261 * qJD(3) + t225 * t248;
t292 = t203 * t188;
t291 = t204 * t219;
t290 = t205 * t219;
t287 = t259 * t246;
t286 = t259 * t248;
t280 = qJD(2) * t247;
t279 = qJD(3) * t246;
t214 = t219 ^ 2;
t184 = t188 * t214 + 0.1e1;
t278 = 0.2e1 * (-t214 * t304 + t219 * t292) / t184 ^ 2;
t277 = 0.2e1 * t298;
t270 = 0.2e1 * t198 * t297;
t269 = -0.2e1 * t216 * t288;
t265 = t241 * t287 + t225;
t263 = t196 * t293 + t295;
t262 = -t227 * t301 + t233 * t289;
t257 = -qJD(3) * t286 + t226 * t246 + t232 * t241;
t224 = t258 * qJD(2);
t220 = -t234 * qJD(3) - t246 * t274;
t211 = t232 * t240 + t239 * t287;
t210 = t232 * t239 - t240 * t287;
t200 = t216 * qJD(3) + t223 * t246;
t191 = 0.1e1 / t193;
t181 = 0.1e1 / t184;
t180 = t206 * t302;
t179 = t262 * t206;
t176 = (-t204 * t230 + t205 * t281) * t248 + t264 * t180;
t175 = t264 * t179 - t204 * t301 + t205 * t233;
t173 = t262 * t300 + (t233 * t269 + t200 * t227 + (t201 * t233 + t216 * t220 + t221 * t301) * t228) * t206;
t171 = t300 * t302 + (-t260 * t279 + (t269 * t281 - t224 * t227 + (t221 * t230 + (t201 * t249 - t216 * t280) * t244) * t228) * t248) * t206;
t170 = -0.2e1 * t298 + 0.2e1 * (t185 * t196 * t191 + (-t191 * t297 - t196 * t298) * t198) * t198;
t1 = [0, t171, t173, 0, 0, 0; 0 (t176 * t296 - t187 * t286) * t278 + ((-t226 * t248 - t259 * t279) * t187 + (-t292 + t303) * t176 + (-t286 * t174 - (-t171 * t216 - t180 * t201 + (-t248 * t280 - t249 * t279) * t244 + (-t180 * t234 - t230 * t248) * t178) * t290 - (t230 * t279 - t171 * t234 - t180 * t221 - t224 * t248 + (t180 * t216 - t248 * t281) * t178) * t291) * t188) * t181 (t175 * t296 - t187 * t261) * t278 + (t175 * t303 - t202 * t187 + (-t261 * t174 - t175 * t203 - (-t173 * t216 - t179 * t201 + t220 + (-t179 * t234 - t301) * t178) * t290 - (-t173 * t234 - t179 * t221 + t200 + (t179 * t216 - t233) * t178) * t291) * t188) * t181, 0, 0, 0; 0 (-t195 * t210 + t211 * t294) * t277 + (t211 * t270 + t265 * t195 * t239 + t257 * t295 + (-t265 * t198 * t240 - t211 * t185 - t210 * t186 + t257 * t293) * t196) * t191, t263 * t219 * t277 + (-t263 * t203 + ((t195 * t241 + t270) * t239 + (-t185 * t239 + (-t198 * t241 + t186) * t240) * t196) * t219) * t191, 0, t170, t170;];
JaD_rot  = t1;
