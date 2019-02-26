% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR13_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_jacobiaD_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:55:38
% EndTime: 2019-02-26 20:55:40
% DurationCPUTime: 1.08s
% Computational Cost: add. (3797->92), mult. (12220->188), div. (410->12), fcn. (15578->13), ass. (0->99)
t222 = sin(pkin(7));
t225 = cos(pkin(7));
t226 = cos(pkin(6));
t228 = cos(qJ(1));
t289 = sin(pkin(12));
t259 = t228 * t289;
t224 = cos(pkin(12));
t290 = sin(qJ(1));
t260 = t290 * t224;
t241 = t226 * t260 + t259;
t223 = sin(pkin(6));
t262 = t223 * t290;
t203 = t241 * t222 + t225 * t262;
t200 = 0.1e1 / t203 ^ 2;
t253 = t290 * t289;
t270 = t228 * t224;
t244 = t226 * t270 - t253;
t272 = t223 * t228;
t202 = t244 * t222 + t225 * t272;
t299 = t200 * t202;
t227 = sin(qJ(3));
t273 = t222 * t227;
t219 = t272 * t273;
t243 = -t226 * t259 - t260;
t271 = t225 * t227;
t291 = cos(qJ(3));
t298 = t243 * t291 - t244 * t271 + t219;
t263 = t222 * t291;
t252 = t263 * t272;
t261 = t225 * t291;
t297 = t244 * t261 - t252;
t184 = -t227 * t243 - t297;
t294 = (t224 * t261 - t289 * t227) * t223 + t226 * t263;
t177 = atan2(-t184, -t294);
t172 = sin(t177);
t173 = cos(t177);
t167 = -t172 * t184 - t173 * t294;
t165 = 0.1e1 / t167 ^ 2;
t213 = -t226 * t253 + t270;
t239 = t241 * t225;
t255 = t222 * t262;
t247 = t291 * t255;
t295 = -t213 * t227 - t291 * t239 + t247;
t182 = t295 ^ 2;
t163 = t182 * t165 + 0.1e1;
t190 = t213 * t291 + (-t239 + t255) * t227;
t208 = t243 * qJD(1);
t236 = qJD(1) * t225 * t244;
t168 = -qJD(1) * t252 + t190 * qJD(3) + t208 * t227 + t291 * t236;
t282 = t168 * t165;
t164 = 0.1e1 / t167;
t181 = t184 ^ 2;
t195 = 0.1e1 / t294 ^ 2;
t176 = t181 * t195 + 0.1e1;
t174 = 0.1e1 / t176;
t194 = 0.1e1 / t294;
t198 = t226 * t273 + (t224 * t271 + t289 * t291) * t223;
t192 = t198 * qJD(3);
t278 = t192 * t195;
t209 = t241 * qJD(1);
t210 = t213 * qJD(1);
t293 = qJD(1) * t247 + qJD(3) * t298 - t209 * t261 - t210 * t227;
t249 = t184 * t278 - t194 * t293;
t157 = t249 * t174;
t251 = t172 * t294 - t173 * t184;
t154 = t251 * t157 + t172 * t293 + t173 * t192;
t296 = t154 * t165;
t287 = t164 * t296;
t269 = 0.2e1 * (-t182 * t287 - t282 * t295) / t163 ^ 2;
t276 = qJD(1) * t299;
t254 = qJD(1) * t262;
t171 = (qJD(3) * t243 - t209 * t225 + t222 * t254) * t227 + t210 * t291 + t297 * qJD(3);
t199 = 0.1e1 / t203;
t292 = -0.2e1 * t295;
t277 = t194 * t278;
t285 = (-t184 * t195 * t293 + t181 * t277) / t176 ^ 2;
t169 = qJD(1) * t219 + t295 * qJD(3) + t208 * t291 - t227 * t236;
t183 = t190 ^ 2;
t180 = t183 * t200 + 0.1e1;
t275 = t199 * t276;
t279 = t190 * t200;
t284 = (t169 * t279 - t183 * t275) / t180 ^ 2;
t283 = t165 * t295;
t281 = t184 * t194;
t280 = t184 * t198;
t268 = -0.2e1 * t285;
t267 = 0.2e1 * t190 * t202;
t266 = t194 * t285;
t265 = t199 * t284;
t258 = t287 * t292;
t248 = -t194 * t298 + t195 * t280;
t245 = -t172 + (-t173 * t281 + t172) * t174;
t191 = t294 * qJD(3);
t178 = 0.1e1 / t180;
t161 = 0.1e1 / t163;
t158 = t248 * t174;
t155 = t251 * t158 + t172 * t298 + t173 * t198;
t153 = t248 * t268 + (0.2e1 * t277 * t280 + t171 * t194 + (t184 * t191 - t192 * t298 - t198 * t293) * t195) * t174;
t1 = [-t266 * t292 + (t168 * t194 - t278 * t295) * t174, 0, t153, 0, 0, 0; t184 * t164 * t269 + (t293 * t164 + t184 * t296 + (t245 * t168 - ((t157 * t174 * t281 + t268) * t172 + (0.2e1 * t184 * t266 - t157 + (t157 - t249) * t174) * t173) * t295) * t283) * t161 - (-t283 * t269 + (-t282 + t258) * t161) * t245 * t295, 0 (-t155 * t283 - t164 * t190) * t269 + (t155 * t258 + t169 * t164 + (-t190 * t154 - t155 * t168 - (-(-t153 * t184 + t158 * t293 + t191 + (t158 * t294 + t298) * t157) * t173 - (t153 * t294 - t158 * t192 - t171 + (t158 * t184 - t198) * t157) * t172) * t295) * t165) * t161, 0, 0, 0; t200 * t267 * t284 - 0.2e1 * t298 * t265 + (-t171 * t199 - t298 * t276 - (-t209 * t222 - t225 * t254) * t279 - t169 * t299 + t267 * t275) * t178, 0, t265 * t292 + (-t168 * t199 - t276 * t295) * t178, 0, 0, 0;];
JaD_rot  = t1;
