% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRP4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRP4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:26:47
% EndTime: 2019-02-26 21:26:49
% DurationCPUTime: 1.53s
% Computational Cost: add. (4632->128), mult. (14196->269), div. (690->12), fcn. (17956->11), ass. (0->121)
t217 = cos(pkin(9));
t216 = sin(pkin(9));
t223 = cos(qJ(1));
t267 = t223 * t216;
t220 = sin(qJ(1));
t222 = cos(qJ(2));
t269 = t220 * t222;
t207 = t217 * t269 - t267;
t218 = sin(qJ(5));
t221 = cos(qJ(5));
t240 = t216 * t269 + t217 * t223;
t187 = t207 * t221 + t240 * t218;
t268 = t222 * t223;
t271 = t220 * t216;
t208 = t217 * t268 + t271;
t219 = sin(qJ(2));
t273 = t219 * t220;
t256 = qJD(2) * t273;
t195 = t208 * qJD(1) - t217 * t256;
t257 = t222 * t267;
t266 = qJD(1) * t220;
t237 = qJD(1) * t257 - t216 * t256 - t217 * t266;
t165 = t187 * qJD(5) + t195 * t218 - t237 * t221;
t296 = t207 * t218 - t240 * t221;
t303 = qJD(5) * t296 - t195 * t221 - t237 * t218;
t270 = t220 * t217;
t241 = t257 - t270;
t190 = t208 * t218 - t241 * t221;
t302 = -0.2e1 * t190;
t191 = t208 * t221 + t241 * t218;
t183 = 0.1e1 / t191 ^ 2;
t214 = t219 ^ 2;
t215 = t223 ^ 2;
t274 = t214 * t215;
t176 = t183 * t274 + 0.1e1;
t264 = qJD(2) * t222;
t265 = qJD(1) * t223;
t263 = qJD(2) * t223;
t255 = t219 * t263;
t194 = -t207 * qJD(1) - t217 * t255;
t230 = -t240 * qJD(1) - t216 * t255;
t164 = -t190 * qJD(5) + t194 * t221 + t230 * t218;
t182 = 0.1e1 / t191;
t288 = t164 * t182 * t183;
t301 = (-t274 * t288 + (-t214 * t220 * t265 + t215 * t219 * t264) * t183) / t176 ^ 2;
t247 = t216 * t221 - t217 * t218;
t206 = t247 * t222;
t232 = qJD(2) * t206;
t246 = t216 * t218 + t217 * t221;
t234 = qJD(5) * t246;
t178 = t219 * t234 - t232;
t299 = t247 * t219;
t202 = 0.1e1 / t299 ^ 2;
t300 = t178 * t202;
t173 = atan2(-t296, -t299);
t168 = sin(t173);
t169 = cos(t173);
t162 = -t168 * t296 - t169 * t299;
t159 = 0.1e1 / t162;
t201 = 0.1e1 / t299;
t160 = 0.1e1 / t162 ^ 2;
t295 = -0.2e1 * t296;
t294 = 0.2e1 * t190;
t181 = t190 ^ 2;
t158 = t160 * t181 + 0.1e1;
t163 = t191 * qJD(5) + t194 * t218 - t230 * t221;
t289 = t160 * t190;
t180 = t296 ^ 2;
t172 = t180 * t202 + 0.1e1;
t170 = 0.1e1 / t172;
t279 = t296 * t202;
t245 = t165 * t201 + t178 * t279;
t151 = t245 * t170;
t248 = t168 * t299 - t169 * t296;
t147 = t248 * t151 - t165 * t168 + t169 * t178;
t292 = t147 * t159 * t160;
t293 = (t163 * t289 - t181 * t292) / t158 ^ 2;
t283 = t201 * t300;
t291 = (t165 * t279 + t180 * t283) / t172 ^ 2;
t156 = 0.1e1 / t158;
t290 = t156 * t160;
t287 = t168 * t190;
t286 = t169 * t190;
t174 = 0.1e1 / t176;
t284 = t174 * t183;
t282 = t183 * t219;
t281 = t183 * t222;
t280 = t296 * t201;
t272 = t219 * t223;
t262 = 0.2e1 * t293;
t261 = -0.2e1 * t291;
t260 = 0.2e1 * t301;
t259 = t201 * t291;
t258 = t183 * t272;
t254 = t222 * t263;
t253 = t292 * t294;
t252 = t283 * t295;
t251 = t272 * t288;
t250 = t273 * t284;
t249 = t258 * t301;
t238 = t246 * t223;
t198 = t219 * t238;
t244 = -t182 * t222 - t198 * t282;
t205 = t246 * t219;
t243 = t187 * t201 + t205 * t279;
t196 = (t218 * t270 - t221 * t271) * t219;
t242 = -t196 * t201 - t206 * t279;
t236 = qJD(1) * t247;
t235 = qJD(2) * t247;
t233 = -t168 + (-t169 * t280 + t168) * t170;
t197 = t223 * t299;
t179 = qJD(5) * t299 + t246 * t264;
t177 = t219 * t235 + t222 * t234;
t167 = -t220 * t232 + (t220 * t234 - t223 * t236) * t219;
t155 = t242 * t170;
t153 = t243 * t170;
t149 = t248 * t155 + t168 * t196 - t169 * t206;
t148 = t248 * t153 - t168 * t187 + t169 * t205;
t146 = t242 * t261 + (t206 * t252 - t167 * t201 + (-t165 * t206 + t177 * t296 - t178 * t196) * t202) * t170;
t145 = t243 * t261 + (-t205 * t252 - t303 * t201 + (t165 * t205 + t178 * t187 + t179 * t296) * t202) * t170;
t1 = [-t259 * t294 + (t163 * t201 + t190 * t300) * t170, t146, 0, 0, t145, 0; -(-t147 * t290 - 0.2e1 * t159 * t293) * t296 + (-t165 * t159 - (t233 * t163 + ((t151 * t170 * t280 + t261) * t168 + (-t259 * t295 - t151 + (t151 - t245) * t170) * t169) * t190) * t289) * t156 + (t156 * t253 - t163 * t290 + t289 * t262) * t233 * t190 (t149 * t289 - t159 * t197) * t262 + (t149 * t253 + (-t197 * t147 - t149 * t163 - (-t146 * t296 - t155 * t165 + t177 + (t155 * t299 + t196) * t151) * t286 - (t146 * t299 - t155 * t178 + t167 + (t155 * t296 + t206) * t151) * t287) * t160 + (t235 * t268 + (-qJD(5) * t238 - t220 * t236) * t219) * t159) * t156, 0, 0 (t148 * t289 - t159 * t191) * t262 + (t148 * t253 + t164 * t159 + (-t191 * t147 - t148 * t163 - (-t145 * t296 - t153 * t165 + t179 + (t153 * t299 - t187) * t151) * t286 - (t145 * t299 - t153 * t178 + t303 + (t153 * t296 - t205) * t151) * t287) * t160) * t156, 0; t164 * t250 + t182 * t260 * t273 - (qJD(1) * t250 - t254 * t284 + 0.2e1 * t249) * t187 + (-t303 * t258 + (-t219 * t265 - t220 * t264) * t182 - 0.2e1 * t187 * t251) * t174, t244 * t223 * t260 + (t244 * t266 + ((qJD(2) * t198 - t164) * t281 + (-0.2e1 * t198 * t288 + (t247 * t223 * qJD(5) - t246 * t266) * t282 + (t238 * t281 - t182) * qJD(2)) * t219) * t223) * t174, 0, 0, t249 * t302 + (t251 * t302 + (t163 * t272 + (-t219 * t266 + t254) * t190) * t183) * t174, 0;];
JaD_rot  = t1;
