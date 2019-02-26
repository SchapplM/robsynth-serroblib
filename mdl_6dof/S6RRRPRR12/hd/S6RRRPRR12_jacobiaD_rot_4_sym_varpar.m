% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR12_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR12_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_jacobiaD_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:22:26
% EndTime: 2019-02-26 22:22:27
% DurationCPUTime: 1.30s
% Computational Cost: add. (4128->138), mult. (12381->282), div. (702->12), fcn. (15714->13), ass. (0->120)
t224 = cos(pkin(6));
t226 = sin(qJ(2));
t292 = sin(qJ(1));
t257 = t292 * t226;
t247 = t224 * t257;
t252 = qJD(2) * t292;
t228 = cos(qJ(2));
t229 = cos(qJ(1));
t271 = t229 * t228;
t222 = sin(pkin(6));
t273 = t222 * t229;
t297 = -qJD(1) * t247 - t226 * t252 + (qJD(2) * t224 + qJD(1)) * t271 - qJD(3) * t273;
t210 = -t247 + t271;
t225 = sin(qJ(3));
t227 = cos(qJ(3));
t258 = t222 * t292;
t201 = t210 * t227 + t225 * t258;
t221 = sin(pkin(12));
t223 = cos(pkin(12));
t256 = t292 * t228;
t272 = t229 * t226;
t238 = -t224 * t256 - t272;
t177 = t201 * t221 + t223 * t238;
t240 = -t224 * t272 - t256;
t188 = qJD(1) * t240 + qJD(2) * t238;
t239 = -t210 * t225 + t227 * t258;
t255 = qJD(1) * t273;
t167 = qJD(3) * t239 + t188 * t227 + t225 * t255;
t251 = t292 * qJD(1);
t259 = t224 * t271;
t269 = qJD(2) * t228;
t187 = -qJD(1) * t259 - t229 * t269 + (t224 * t252 + t251) * t226;
t162 = t167 * t223 - t187 * t221;
t178 = t201 * t223 - t221 * t238;
t172 = 0.1e1 / t178;
t173 = 0.1e1 / t178 ^ 2;
t286 = t162 * t172 * t173;
t250 = 0.2e1 * t177 * t286;
t195 = -t225 * t240 + t227 * t273;
t275 = t222 * t226;
t205 = -t224 * t227 + t225 * t275;
t184 = atan2(-t195, t205);
t179 = sin(t184);
t180 = cos(t184);
t160 = -t179 * t195 + t180 * t205;
t158 = 0.1e1 / t160 ^ 2;
t192 = t239 ^ 2;
t156 = t192 * t158 + 0.1e1;
t166 = qJD(3) * t201 + t188 * t225 - t227 * t255;
t287 = t158 * t239;
t191 = t195 ^ 2;
t203 = 0.1e1 / t205 ^ 2;
t183 = t191 * t203 + 0.1e1;
t181 = 0.1e1 / t183;
t246 = t222 * t251;
t268 = qJD(3) * t227;
t168 = t297 * t225 - t227 * t246 - t240 * t268;
t206 = t224 * t225 + t227 * t275;
t254 = t222 * t269;
t193 = qJD(3) * t206 + t225 * t254;
t202 = 0.1e1 / t205;
t279 = t195 * t203;
t244 = -t168 * t202 + t193 * t279;
t150 = t244 * t181;
t245 = -t179 * t205 - t180 * t195;
t145 = t150 * t245 - t179 * t168 + t180 * t193;
t157 = 0.1e1 / t160;
t159 = t157 * t158;
t290 = t145 * t159;
t267 = 0.2e1 * (-t166 * t287 - t192 * t290) / t156 ^ 2;
t296 = t193 * t203;
t207 = -t257 + t259;
t274 = t222 * t228;
t241 = -t202 * t207 + t274 * t279;
t295 = t225 * t241;
t169 = t225 * (qJD(3) * t240 + t246) + t297 * t227;
t294 = -0.2e1 * t195;
t293 = -0.2e1 * t239;
t281 = t202 * t296;
t289 = (t168 * t279 - t191 * t281) / t183 ^ 2;
t288 = t158 * t166;
t285 = t173 * t177;
t284 = t177 * t223;
t283 = t179 * t239;
t282 = t180 * t239;
t280 = t195 * t202;
t278 = t238 * t225;
t277 = t238 * t227;
t276 = t221 * t172;
t270 = qJD(2) * t226;
t161 = t167 * t221 + t187 * t223;
t171 = t177 ^ 2;
t165 = t171 * t173 + 0.1e1;
t266 = 0.2e1 * (t161 * t285 - t171 * t286) / t165 ^ 2;
t265 = -0.2e1 * t289;
t264 = t159 * t293;
t263 = t202 * t289;
t262 = t158 * t283;
t261 = t158 * t282;
t249 = t281 * t294;
t197 = -t225 * t273 - t227 * t240;
t243 = -t197 * t202 + t206 * t279;
t242 = -qJD(3) * t278 + t187 * t227;
t236 = -t179 + (t180 * t280 + t179) * t181;
t194 = -qJD(3) * t205 + t227 * t254;
t189 = qJD(1) * t238 + qJD(2) * t240;
t186 = t210 * t221 + t223 * t277;
t185 = -t210 * t223 + t221 * t277;
t176 = -t197 * t223 + t207 * t221;
t175 = -t197 * t221 - t207 * t223;
t163 = 0.1e1 / t165;
t154 = 0.1e1 / t156;
t153 = t181 * t295;
t152 = t243 * t181;
t149 = t236 * t239;
t147 = (-t179 * t207 + t180 * t274) * t225 + t245 * t153;
t146 = t152 * t245 - t179 * t197 + t180 * t206;
t144 = t243 * t265 + (t206 * t249 - t169 * t202 + (t168 * t206 + t193 * t197 + t194 * t195) * t203) * t181;
t142 = t265 * t295 + (t241 * t268 + (t249 * t274 - t189 * t202 + (t193 * t207 + (t168 * t228 - t195 * t270) * t222) * t203) * t225) * t181;
t1 = [t263 * t293 + (-t166 * t202 - t239 * t296) * t181, t142, t144, 0, 0, 0; t195 * t157 * t267 + (-t168 * t157 + (t145 * t195 + t149 * t166) * t158) * t154 - (-t149 * t158 * t267 + (-0.2e1 * t149 * t290 + (-t150 * t181 * t280 + t265) * t262 + (t263 * t294 - t150 + (t150 - t244) * t181) * t261 - t236 * t288) * t154) * t239 (-t147 * t287 - t157 * t278) * t267 + (-t147 * t288 + (t187 * t225 + t238 * t268) * t157 + (t147 * t264 - t158 * t278) * t145 + (-t142 * t195 - t153 * t168 + (-t225 * t270 + t228 * t268) * t222 + (-t153 * t205 - t207 * t225) * t150) * t261 + (-t207 * t268 - t142 * t205 - t153 * t193 - t189 * t225 + (t153 * t195 - t225 * t274) * t150) * t262) * t154 (-t146 * t287 - t157 * t201) * t267 + (t146 * t145 * t264 + t167 * t157 + (-t201 * t145 - t146 * t166 + (-t144 * t195 - t152 * t168 + t194 + (-t152 * t205 - t197) * t150) * t282 + (-t144 * t205 - t152 * t193 - t169 + (t152 * t195 - t206) * t150) * t283) * t158) * t154, 0, 0, 0; (-t172 * t175 + t176 * t285) * t266 + ((-t169 * t221 - t189 * t223) * t172 + t176 * t250 + (-t175 * t162 - (-t169 * t223 + t189 * t221) * t177 - t176 * t161) * t173) * t163 (-t172 * t185 + t186 * t285) * t266 + ((-t188 * t223 + t221 * t242) * t172 + t186 * t250 + (-t185 * t162 - (t188 * t221 + t223 * t242) * t177 - t186 * t161) * t173) * t163 -(-t173 * t284 + t276) * t239 * t266 + (t239 * t223 * t250 - t166 * t276 + (t166 * t284 - (t161 * t223 + t162 * t221) * t239) * t173) * t163, 0, 0, 0;];
JaD_rot  = t1;
