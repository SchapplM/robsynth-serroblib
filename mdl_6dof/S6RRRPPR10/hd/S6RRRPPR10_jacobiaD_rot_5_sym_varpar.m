% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPPR10_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR10_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:08:47
% EndTime: 2019-02-26 22:08:48
% DurationCPUTime: 1.25s
% Computational Cost: add. (4128->140), mult. (12381->286), div. (702->12), fcn. (15714->13), ass. (0->123)
t225 = sin(qJ(1));
t222 = cos(pkin(6));
t241 = qJD(2) * t222 + qJD(1);
t224 = sin(qJ(2));
t266 = t225 * t224;
t250 = t222 * t266;
t261 = qJD(2) * t224;
t227 = cos(qJ(2));
t228 = cos(qJ(1));
t263 = t228 * t227;
t220 = sin(pkin(6));
t269 = t220 * t228;
t292 = -qJD(1) * t250 - qJD(3) * t269 - t225 * t261 + t241 * t263;
t264 = t228 * t224;
t265 = t225 * t227;
t207 = t222 * t264 + t265;
t226 = cos(qJ(3));
t223 = sin(qJ(3));
t267 = t223 * t228;
t193 = t207 * t226 - t220 * t267;
t271 = t220 * t226;
t205 = t222 * t223 + t224 * t271;
t181 = atan2(-t193, t205);
t176 = sin(t181);
t177 = cos(t181);
t157 = -t176 * t193 + t177 * t205;
t155 = 0.1e1 / t157 ^ 2;
t209 = -t250 + t263;
t272 = t220 * t223;
t198 = t209 * t226 + t225 * t272;
t189 = t198 ^ 2;
t153 = t189 * t155 + 0.1e1;
t208 = t222 * t265 + t264;
t185 = -qJD(1) * t207 - qJD(2) * t208;
t258 = qJD(3) * t226;
t259 = qJD(3) * t223;
t164 = t185 * t226 - t209 * t259 + (qJD(1) * t267 + t225 * t258) * t220;
t283 = t155 * t198;
t188 = t193 ^ 2;
t202 = 0.1e1 / t205 ^ 2;
t180 = t188 * t202 + 0.1e1;
t178 = 0.1e1 / t180;
t242 = t292 * t226;
t262 = qJD(1) * t220;
t248 = t225 * t262;
t166 = -t207 * t259 + t223 * t248 + t242;
t204 = t222 * t226 - t224 * t272;
t260 = qJD(2) * t227;
t246 = t220 * t260;
t191 = qJD(3) * t204 + t226 * t246;
t201 = 0.1e1 / t205;
t275 = t193 * t202;
t239 = -t166 * t201 + t191 * t275;
t147 = t239 * t178;
t240 = -t176 * t205 - t177 * t193;
t142 = t147 * t240 - t176 * t166 + t177 * t191;
t154 = 0.1e1 / t157;
t156 = t154 * t155;
t286 = t142 * t156;
t257 = 0.2e1 * (t164 * t283 - t189 * t286) / t153 ^ 2;
t291 = t191 * t202;
t249 = t222 * t263;
t206 = t249 - t266;
t270 = t220 * t227;
t236 = -t201 * t206 + t270 * t275;
t290 = t226 * t236;
t197 = t209 * t223 - t225 * t271;
t219 = sin(pkin(11));
t221 = cos(pkin(11));
t175 = t197 * t219 + t208 * t221;
t169 = 0.1e1 / t175;
t170 = 0.1e1 / t175 ^ 2;
t289 = -0.2e1 * t193;
t288 = 0.2e1 * t198;
t277 = t201 * t291;
t285 = (t166 * t275 - t188 * t277) / t180 ^ 2;
t284 = t155 * t164;
t247 = t226 * t262;
t163 = qJD(3) * t198 + t185 * t223 - t228 * t247;
t184 = -qJD(1) * t249 - t228 * t260 + t241 * t266;
t159 = t163 * t219 - t184 * t221;
t282 = t159 * t169 * t170;
t174 = -t197 * t221 + t208 * t219;
t281 = t170 * t174;
t280 = t174 * t219;
t279 = t176 * t198;
t278 = t177 * t198;
t276 = t193 * t201;
t274 = t208 * t223;
t273 = t208 * t226;
t268 = t221 * t169;
t158 = -t163 * t221 - t184 * t219;
t168 = t174 ^ 2;
t162 = t168 * t170 + 0.1e1;
t256 = 0.2e1 * (t158 * t281 - t168 * t282) / t162 ^ 2;
t255 = -0.2e1 * t285;
t254 = t156 * t288;
t253 = t201 * t285;
t252 = t155 * t279;
t251 = t155 * t278;
t244 = 0.2e1 * t174 * t282;
t243 = t277 * t289;
t192 = t207 * t223 + t226 * t269;
t238 = t192 * t201 + t204 * t275;
t237 = -t184 * t223 + t208 * t258;
t235 = -t176 + (t177 * t276 + t176) * t178;
t165 = t207 * t258 + t292 * t223 - t225 * t247;
t190 = -qJD(3) * t205 - t223 * t246;
t186 = -qJD(1) * t208 - qJD(2) * t207;
t183 = t209 * t221 - t219 * t274;
t182 = t209 * t219 + t221 * t274;
t173 = -t192 * t219 + t206 * t221;
t172 = t192 * t221 + t206 * t219;
t160 = 0.1e1 / t162;
t151 = 0.1e1 / t153;
t150 = t178 * t290;
t149 = t238 * t178;
t146 = t235 * t198;
t144 = (-t176 * t206 + t177 * t270) * t226 + t240 * t150;
t143 = t149 * t240 + t176 * t192 + t177 * t204;
t141 = t238 * t255 + (t204 * t243 + t165 * t201 + (t166 * t204 + t190 * t193 - t191 * t192) * t202) * t178;
t139 = t255 * t290 + (-t236 * t259 + (t243 * t270 - t186 * t201 + (t191 * t206 + (t166 * t227 - t193 * t261) * t220) * t202) * t226) * t178;
t1 = [t253 * t288 + (-t164 * t201 + t198 * t291) * t178, t139, t141, 0, 0, 0; t193 * t154 * t257 + (((qJD(3) * t207 - t248) * t223 - t242) * t154 + (t142 * t193 - t146 * t164) * t155) * t151 + (t146 * t155 * t257 + (0.2e1 * t146 * t286 - (-t147 * t178 * t276 + t255) * t252 - (t253 * t289 - t147 + (t147 - t239) * t178) * t251 - t235 * t284) * t151) * t198 (t144 * t283 + t154 * t273) * t257 + (-t144 * t284 + (t184 * t226 + t208 * t259) * t154 + (t144 * t254 + t155 * t273) * t142 - (-t139 * t193 - t150 * t166 + (-t226 * t261 - t227 * t259) * t220 + (-t150 * t205 - t206 * t226) * t147) * t251 - (t206 * t259 - t139 * t205 - t150 * t191 - t186 * t226 + (t150 * t193 - t226 * t270) * t147) * t252) * t151 (t143 * t283 + t154 * t197) * t257 + (t143 * t142 * t254 - t163 * t154 + (t197 * t142 - t143 * t164 - (-t141 * t193 - t149 * t166 + t190 + (-t149 * t205 + t192) * t147) * t278 - (-t141 * t205 - t149 * t191 + t165 + (t149 * t193 - t204) * t147) * t279) * t155) * t151, 0, 0, 0; (-t169 * t172 + t173 * t281) * t256 + ((t165 * t221 + t186 * t219) * t169 + t173 * t244 + (-t172 * t159 - (-t165 * t219 + t186 * t221) * t174 - t173 * t158) * t170) * t160 (-t169 * t182 + t183 * t281) * t256 + ((t185 * t219 + t221 * t237) * t169 + t183 * t244 + (-t182 * t159 - (t185 * t221 - t219 * t237) * t174 - t183 * t158) * t170) * t160 (t170 * t280 + t268) * t198 * t256 + (t198 * t219 * t244 - t164 * t268 + (-t164 * t280 + (-t158 * t219 + t159 * t221) * t198) * t170) * t160, 0, 0, 0;];
JaD_rot  = t1;
