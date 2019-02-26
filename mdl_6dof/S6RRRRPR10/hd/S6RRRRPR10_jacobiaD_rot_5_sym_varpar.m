% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR10_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR10_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:35:48
% EndTime: 2019-02-26 22:35:50
% DurationCPUTime: 1.30s
% Computational Cost: add. (10248->127), mult. (14849->255), div. (934->12), fcn. (18926->11), ass. (0->114)
t233 = qJ(3) + qJ(4);
t230 = sin(t233);
t235 = cos(pkin(6));
t236 = sin(qJ(2));
t296 = cos(qJ(1));
t264 = t296 * t236;
t237 = sin(qJ(1));
t238 = cos(qJ(2));
t278 = t237 * t238;
t246 = -t235 * t264 - t278;
t231 = cos(t233);
t234 = sin(pkin(6));
t265 = t234 * t296;
t257 = t231 * t265;
t201 = -t230 * t246 + t257;
t282 = t234 * t236;
t268 = t230 * t282;
t211 = -t231 * t235 + t268;
t188 = atan2(-t201, t211);
t183 = sin(t188);
t184 = cos(t188);
t178 = -t183 * t201 + t184 * t211;
t176 = 0.1e1 / t178 ^ 2;
t263 = t296 * t238;
t279 = t237 * t236;
t266 = t235 * t279;
t219 = t263 - t266;
t281 = t234 * t237;
t206 = t219 * t230 - t231 * t281;
t199 = t206 ^ 2;
t174 = t176 * t199 + 0.1e1;
t247 = -t235 * t278 - t264;
t195 = qJD(1) * t246 + qJD(2) * t247;
t232 = qJD(3) + qJD(4);
t254 = t232 * t281 + t195;
t261 = t296 * qJD(1);
t255 = t234 * t261;
t283 = t231 * t232;
t179 = t219 * t283 + t230 * t254 - t231 * t255;
t291 = t179 * t176;
t198 = t201 ^ 2;
t209 = 0.1e1 / t211 ^ 2;
t187 = t198 * t209 + 0.1e1;
t185 = 0.1e1 / t187;
t260 = t296 * qJD(2);
t277 = qJD(2) * t236;
t197 = -qJD(1) * t266 - t237 * t277 + (t235 * t260 + t261) * t238;
t225 = t230 * t265;
t262 = qJD(1) * t281;
t181 = t197 * t230 - t225 * t232 - t231 * t262 - t246 * t283;
t280 = t234 * t238;
t249 = qJD(2) * t280 + t232 * t235;
t267 = t231 * t282;
t192 = t230 * t249 + t232 * t267;
t208 = 0.1e1 / t211;
t286 = t201 * t209;
t252 = -t181 * t208 + t192 * t286;
t167 = t252 * t185;
t253 = -t183 * t211 - t184 * t201;
t162 = t167 * t253 - t183 * t181 + t184 * t192;
t175 = 0.1e1 / t178;
t177 = t175 * t176;
t294 = t162 * t177;
t276 = 0.2e1 * (-t199 * t294 + t206 * t291) / t174 ^ 2;
t300 = t192 * t209;
t256 = t235 * t263;
t216 = t256 - t279;
t248 = -t208 * t216 + t280 * t286;
t299 = t230 * t248;
t182 = (t232 * t246 + t262) * t230 + t197 * t231 - t232 * t257;
t213 = 0.1e1 / t247;
t214 = 0.1e1 / t247 ^ 2;
t298 = -0.2e1 * t201;
t297 = 0.2e1 * t206;
t288 = t208 * t300;
t293 = (t181 * t286 - t198 * t288) / t187 ^ 2;
t292 = t176 * t206;
t290 = t183 * t206;
t289 = t184 * t206;
t287 = t201 * t208;
t207 = t219 * t231 + t230 * t281;
t285 = t207 * t214;
t284 = t247 * t230;
t275 = -0.2e1 * t293;
t180 = t254 * t231 + (-t219 * t232 + t255) * t230;
t200 = t207 ^ 2;
t191 = t200 * t214 + 0.1e1;
t194 = -qJD(1) * t256 - t238 * t260 + (qJD(2) * t235 + qJD(1)) * t279;
t215 = t213 * t214;
t274 = 0.2e1 * (-t194 * t200 * t215 + t180 * t285) / t191 ^ 2;
t273 = t177 * t297;
t272 = 0.2e1 * t207 * t215;
t271 = t208 * t293;
t270 = t176 * t290;
t269 = t176 * t289;
t259 = t288 * t298;
t203 = -t231 * t246 - t225;
t212 = t230 * t235 + t267;
t251 = -t203 * t208 + t212 * t286;
t245 = -t183 + (t184 * t287 + t183) * t185;
t196 = qJD(1) * t247 + qJD(2) * t246;
t193 = t231 * t249 - t232 * t268;
t189 = 0.1e1 / t191;
t172 = 0.1e1 / t174;
t171 = t185 * t299;
t169 = t251 * t185;
t166 = t245 * t206;
t165 = -t206 * t213 * t274 + (-t194 * t206 * t214 + t179 * t213) * t189;
t164 = (-t183 * t216 + t184 * t280) * t230 + t253 * t171;
t163 = t169 * t253 - t183 * t203 + t184 * t212;
t160 = t251 * t275 + (t212 * t259 - t182 * t208 + (t181 * t212 + t192 * t203 + t193 * t201) * t209) * t185;
t159 = t275 * t299 + (t248 * t283 + (t259 * t280 - t196 * t208 + (t192 * t216 + (t181 * t238 - t201 * t277) * t234) * t209) * t230) * t185;
t158 = (t163 * t292 - t175 * t207) * t276 + (t163 * t162 * t273 + t180 * t175 + (-t207 * t162 - t163 * t179 - (-t160 * t201 - t169 * t181 + t193 + (-t169 * t211 - t203) * t167) * t289 - (-t160 * t211 - t169 * t192 - t182 + (t169 * t201 - t212) * t167) * t290) * t176) * t172;
t1 = [t271 * t297 + (-t179 * t208 + t206 * t300) * t185, t159, t160, t160, 0, 0; t201 * t175 * t276 + (-t181 * t175 + (t162 * t201 - t166 * t179) * t176) * t172 + (t166 * t176 * t276 + (0.2e1 * t166 * t294 - (-t167 * t185 * t287 + t275) * t270 - (t271 * t298 - t167 + (t167 - t252) * t185) * t269 - t245 * t291) * t172) * t206 (t164 * t292 - t175 * t284) * t276 + (-t164 * t291 + (t194 * t230 + t247 * t283) * t175 + (t164 * t273 - t176 * t284) * t162 - (-t159 * t201 - t171 * t181 + (-t230 * t277 + t238 * t283) * t234 + (-t171 * t211 - t216 * t230) * t167) * t269 - (-t216 * t283 - t159 * t211 - t171 * t192 - t196 * t230 + (t171 * t201 - t230 * t280) * t167) * t270) * t172, t158, t158, 0, 0; (-t203 * t213 + t216 * t285) * t274 + (t182 * t213 + t216 * t194 * t272 + (-t180 * t216 - t194 * t203 - t196 * t207) * t214) * t189 (t213 * t231 * t247 + t219 * t285) * t274 + (t213 * t232 * t284 + (-t180 * t219 - t195 * t207) * t214 + (t219 * t272 + (t214 * t247 - t213) * t231) * t194) * t189, t165, t165, 0, 0;];
JaD_rot  = t1;
