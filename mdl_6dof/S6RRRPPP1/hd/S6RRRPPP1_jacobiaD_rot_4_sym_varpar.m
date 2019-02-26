% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPPP1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPPP1_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPP1_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_jacobiaD_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:02:48
% EndTime: 2019-02-26 22:02:50
% DurationCPUTime: 1.53s
% Computational Cost: add. (3960->153), mult. (12250->321), div. (690->12), fcn. (15383->13), ass. (0->139)
t218 = sin(qJ(2));
t219 = sin(qJ(1));
t221 = cos(qJ(2));
t272 = qJD(2) * t221;
t222 = cos(qJ(1));
t275 = qJD(1) * t222;
t232 = t218 * t275 + t219 * t272;
t216 = cos(pkin(6));
t217 = sin(qJ(3));
t214 = sin(pkin(6));
t288 = t214 * t221;
t200 = t216 * t218 + t217 * t288;
t220 = cos(qJ(3));
t269 = qJD(3) * t220;
t250 = t218 * t269;
t189 = t200 * qJD(2) + t214 * t250;
t285 = t217 * t218;
t199 = t214 * t285 - t216 * t221;
t197 = 0.1e1 / t199 ^ 2;
t309 = t189 * t197;
t271 = qJD(2) * t222;
t277 = qJD(1) * t219;
t286 = t216 * t217;
t308 = -t214 * t277 + t271 * t286;
t307 = -qJD(2) * t214 + t216 * t269;
t231 = t218 * t277 - t221 * t271;
t276 = qJD(1) * t221;
t242 = -qJD(3) + t276;
t268 = qJD(3) * t221;
t243 = -qJD(1) + t268;
t255 = t218 * t271;
t279 = t220 * t222;
t178 = -t243 * t279 + (t242 * t219 + t255) * t217;
t306 = -t178 * t216 + t231 * t214;
t273 = qJD(2) * t218;
t256 = t219 * t273;
t281 = t219 * t220;
t180 = t243 * t281 + (t242 * t222 - t256) * t217;
t305 = -t180 * t216 + t232 * t214;
t280 = t219 * t221;
t201 = t217 * t280 + t279;
t284 = t218 * t219;
t244 = -t201 * t214 - t216 * t284;
t177 = atan2(t244, t199);
t172 = sin(t177);
t173 = cos(t177);
t155 = t172 * t244 + t173 * t199;
t152 = 0.1e1 / t155;
t278 = t221 * t222;
t204 = t219 * t217 + t220 * t278;
t213 = sin(pkin(10));
t215 = cos(pkin(10));
t260 = t217 * t278;
t203 = -t260 + t281;
t282 = t218 * t222;
t236 = t203 * t216 + t214 * t282;
t171 = t204 * t215 + t236 * t213;
t165 = 0.1e1 / t171;
t196 = 0.1e1 / t199;
t153 = 0.1e1 / t155 ^ 2;
t166 = 0.1e1 / t171 ^ 2;
t192 = t203 * t214 - t216 * t282;
t187 = t192 ^ 2;
t151 = t153 * t187 + 0.1e1;
t161 = t178 * t214 + t216 * t231;
t299 = t161 * t153;
t186 = t244 ^ 2;
t176 = t186 * t197 + 0.1e1;
t174 = 0.1e1 / t176;
t162 = -t180 * t214 - t232 * t216;
t290 = t244 * t197;
t239 = t162 * t196 - t189 * t290;
t145 = t239 * t174;
t240 = -t172 * t199 + t173 * t244;
t140 = t240 * t145 + t162 * t172 + t173 * t189;
t303 = t140 * t152 * t153;
t304 = (-t187 * t303 + t192 * t299) / t151 ^ 2;
t292 = t196 * t309;
t302 = (t162 * t290 - t186 * t292) / t176 ^ 2;
t301 = t152 * t214;
t300 = t153 * t192;
t170 = t204 * t213 - t236 * t215;
t298 = t166 * t170;
t235 = t216 * t285 + t288;
t283 = t218 * t220;
t185 = (t235 * t213 - t215 * t283) * t222;
t297 = t166 * t185;
t296 = t172 * t192;
t295 = t173 * t192;
t291 = t244 * t196;
t289 = t213 * t216;
t287 = t215 * t216;
t270 = qJD(3) * t217;
t267 = -0.2e1 * t304;
t266 = 0.2e1 * t304;
t265 = -0.2e1 * t303;
t233 = t219 * t276 + t255;
t249 = t219 * t269;
t179 = qJD(3) * t260 - t217 * t275 + t233 * t220 - t249;
t156 = -t179 * t213 + t215 * t306;
t157 = -t179 * t215 - t213 * t306;
t164 = t170 ^ 2;
t160 = t164 * t166 + 0.1e1;
t167 = t165 * t166;
t264 = 0.2e1 * (-t157 * t164 * t167 + t156 * t298) / t160 ^ 2;
t263 = 0.2e1 * t302;
t262 = 0.2e1 * t167 * t170;
t261 = t214 * t284;
t253 = t220 * t271;
t248 = -0.2e1 * t196 * t302;
t247 = t157 * t262;
t246 = 0.2e1 * t244 * t292;
t245 = t192 * t265;
t194 = -t216 * t280 + t217 * t261;
t238 = -t194 * t196 + t200 * t290;
t237 = -t201 * t216 + t261;
t202 = t217 * t222 - t220 * t280;
t234 = -t196 * t202 + t283 * t290;
t230 = t218 * t270 - t220 * t272;
t229 = t172 + (t173 * t291 - t172) * t174;
t195 = t199 * t222;
t188 = t214 * t220 * t268 - t199 * qJD(2);
t184 = (-t213 * t283 - t235 * t215) * t222;
t183 = t203 * t215 - t204 * t289;
t182 = t203 * t213 + t204 * t287;
t181 = -t242 * t279 + (t243 * t217 + t220 * t273) * t219;
t169 = t202 * t215 - t237 * t213;
t168 = t202 * t213 + t237 * t215;
t163 = (-t221 * t275 + t256) * t216 + (t232 * t217 + t218 * t249) * t214;
t158 = 0.1e1 / t160;
t149 = 0.1e1 / t151;
t148 = t234 * t214 * t174;
t146 = t238 * t174;
t144 = t229 * t192;
t142 = (t172 * t202 + t173 * t283) * t214 - t240 * t148;
t141 = -t240 * t146 + t172 * t194 + t173 * t200;
t139 = t238 * t263 + (t200 * t246 + t163 * t196 + (-t162 * t200 - t188 * t244 - t189 * t194) * t197) * t174;
t137 = (t234 * t263 + (t246 * t283 + t181 * t196 + (-t162 * t283 - t189 * t202 + t230 * t244) * t197) * t174) * t214;
t1 = [t192 * t248 + (t161 * t196 - t192 * t309) * t174, t139, t137, 0, 0, 0; t244 * t152 * t267 + (t162 * t152 + (-t140 * t244 + t144 * t161) * t153) * t149 + ((t144 * t265 + t229 * t299) * t149 + (t144 * t267 + ((-t145 * t174 * t291 + t263) * t296 + (t244 * t248 + t145 + (-t145 + t239) * t174) * t295) * t149) * t153) * t192 (-t141 * t300 + t152 * t195) * t266 + (t141 * t245 - t233 * t152 * t216 + (t231 * t217 - t222 * t250) * t301 + (t195 * t140 + t141 * t161 + (t139 * t244 - t146 * t162 + t188 + (t146 * t199 + t194) * t145) * t295 + (-t139 * t199 + t146 * t189 + t163 + (t146 * t244 - t200) * t145) * t296) * t153) * t149 (-t142 * t300 - t204 * t301) * t266 + ((t240 * t137 - (-t155 * t145 + t162 * t173 - t172 * t189) * t148) * t300 + (t245 + t299) * t142 + (-t179 * t152 + (-t204 * t140 + (-t145 * t283 + t181) * t296 + (t145 * t202 - t230) * t295) * t153) * t214) * t149, 0, 0, 0; (-t165 * t168 + t169 * t298) * t264 + (t169 * t247 + (-t169 * t156 - t168 * t157 + (-t181 * t215 + t213 * t305) * t170) * t166 + (t181 * t213 + t215 * t305) * t165) * t158 (-t165 * t184 + t170 * t297) * t264 + (-t156 * t297 + (-t166 * t184 + t185 * t262) * t157 + ((-t213 * t253 - t215 * t308) * t165 - (t213 * t308 - t215 * t253) * t298) * t221 + (((t213 * t220 + t215 * t286) * t165 - (-t213 * t286 + t215 * t220) * t298) * t277 + ((t213 * t270 - t215 * t307) * t165 - (t213 * t307 + t215 * t270) * t298) * t222) * t218) * t158 (-t165 * t182 + t183 * t298) * t264 + ((t178 * t213 - t179 * t287) * t165 + t183 * t247 + (-t182 * t157 - (t178 * t215 + t179 * t289) * t170 - t183 * t156) * t166) * t158, 0, 0, 0;];
JaD_rot  = t1;
