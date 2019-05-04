% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR14V3
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR14V3_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14V3_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobiaD_rot_5_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:12:04
% EndTime: 2019-04-12 15:12:06
% DurationCPUTime: 1.42s
% Computational Cost: add. (2288->151), mult. (7495->320), div. (1129->14), fcn. (9361->11), ass. (0->133)
t200 = cos(qJ(2));
t196 = sin(qJ(4));
t287 = sin(qJ(1));
t237 = t287 * t196;
t199 = cos(qJ(4));
t201 = cos(qJ(1));
t263 = t201 * t199;
t177 = t200 * t263 + t237;
t195 = sin(qJ(5));
t198 = cos(qJ(5));
t197 = sin(qJ(2));
t265 = t197 * t201;
t161 = t177 * t195 - t198 * t265;
t291 = 0.2e1 * t161;
t162 = t177 * t198 + t195 * t265;
t156 = 0.1e1 / t162;
t157 = 0.1e1 / t162 ^ 2;
t278 = t157 * t161;
t217 = t156 * t195 - t198 * t278;
t290 = -qJD(5) * t199 + qJD(2);
t173 = t200 * t237 + t263;
t228 = t287 * qJD(4);
t220 = t196 * t228;
t257 = qJD(4) * t201;
t232 = t199 * t257;
t260 = qJD(2) * t201;
t234 = t197 * t260;
t151 = qJD(1) * t173 + t196 * t234 - t200 * t232 - t220;
t236 = t287 * t199;
t264 = t201 * t196;
t176 = t200 * t264 - t236;
t189 = 0.1e1 / t196;
t190 = 0.1e1 / t196 ^ 2;
t192 = 0.1e1 / t197;
t193 = 0.1e1 / t197 ^ 2;
t261 = qJD(2) * t200;
t235 = t193 * t261;
t258 = qJD(4) * t199;
t270 = t189 * t192;
t289 = t176 * (t190 * t192 * t258 + t189 * t235) + t151 * t270;
t267 = t197 * t196;
t167 = atan2(-t173, t267);
t164 = cos(t167);
t163 = sin(t167);
t276 = t163 * t173;
t147 = t164 * t267 - t276;
t144 = 0.1e1 / t147;
t145 = 0.1e1 / t147 ^ 2;
t288 = 0.2e1 * t176;
t171 = t173 ^ 2;
t269 = t190 * t193;
t168 = t171 * t269 + 0.1e1;
t165 = 0.1e1 / t168;
t214 = t196 * t261 + t197 * t258;
t241 = t173 * t269;
t229 = qJD(2) * t287;
t221 = t197 * t229;
t230 = qJD(1) * t287;
t262 = qJD(1) * t201;
t153 = (t200 * t228 - t230) * t199 + (t200 * t262 - t221 - t257) * t196;
t244 = t153 * t270;
t135 = (t214 * t241 - t244) * t165;
t212 = -t135 * t173 + t214;
t130 = (-t135 * t267 - t153) * t163 + t212 * t164;
t146 = t144 * t145;
t286 = t130 * t146;
t211 = qJD(5) * t177 + t197 * t230 - t200 * t260;
t152 = (-qJD(4) * t200 + qJD(1)) * t264 + (-t200 * t230 + t228 - t234) * t199;
t219 = qJD(5) * t265 + t152;
t137 = t195 * t219 + t198 * t211;
t155 = t161 ^ 2;
t150 = t155 * t157 + 0.1e1;
t138 = -t195 * t211 + t198 * t219;
t158 = t156 * t157;
t283 = t138 * t158;
t285 = (t137 * t278 - t155 * t283) / t150 ^ 2;
t191 = t189 * t190;
t194 = t192 * t193;
t233 = t193 * t258;
t284 = (t153 * t241 + (-t190 * t194 * t261 - t191 * t233) * t171) / t168 ^ 2;
t282 = t145 * t176;
t148 = 0.1e1 / t150;
t281 = t148 * t157;
t280 = t151 * t145;
t266 = t197 * t199;
t170 = (t195 * t200 - t198 * t266) * t201;
t277 = t161 * t170;
t275 = t163 * t176;
t274 = t163 * t197;
t273 = t164 * t173;
t272 = t164 * t176;
t271 = t164 * t200;
t268 = t190 * t199;
t259 = qJD(4) * t196;
t172 = t176 ^ 2;
t142 = t145 * t172 + 0.1e1;
t255 = 0.2e1 * (-t172 * t286 - t176 * t280) / t142 ^ 2;
t254 = -0.2e1 * t285;
t253 = 0.2e1 * t285;
t252 = -0.2e1 * t284;
t251 = t146 * t288;
t250 = t157 * t285;
t249 = t192 * t284;
t248 = t137 * t281;
t247 = t161 * t283;
t246 = t145 * t275;
t242 = t173 * t270;
t240 = t189 * t193 * t200;
t239 = t197 * t287;
t238 = t200 * t287;
t215 = t173 * t240 + t287;
t143 = t215 * t165;
t231 = t287 - t143;
t227 = t144 * t255;
t226 = t145 * t255;
t224 = t189 * t249;
t223 = t195 * t239;
t222 = t198 * t239;
t154 = qJD(1) * t177 - t199 * t221 - t200 * t220 - t232;
t218 = -qJD(5) * t239 - t154;
t175 = t200 * t236 - t264;
t216 = t173 * t268 - t175 * t189;
t210 = -qJD(5) * t175 + t197 * t262 + t200 * t229;
t169 = (-t195 * t266 - t198 * t200) * t201;
t160 = -t175 * t198 - t223;
t140 = 0.1e1 / t142;
t139 = t216 * t192 * t165;
t134 = (-t163 + (t164 * t242 + t163) * t165) * t176;
t133 = -t143 * t273 + (t231 * t274 + t271) * t196;
t131 = t164 * t266 - t163 * t175 + (-t163 * t267 - t273) * t139;
t129 = t215 * t252 + (t153 * t240 + t262 + (-t190 * t200 * t233 + (-0.2e1 * t194 * t200 ^ 2 - t192) * t189 * qJD(2)) * t173) * t165;
t127 = -0.2e1 * t216 * t249 + (-t216 * t235 + (t153 * t268 - t154 * t189 + (t175 * t268 + (-0.2e1 * t191 * t199 ^ 2 - t189) * t173) * qJD(4)) * t192) * t165;
t1 = [t289 * t165 + t224 * t288, t129, 0, t127, 0, 0; t173 * t227 + (-t153 * t144 + (t130 * t173 + t134 * t151) * t145) * t140 + (t134 * t226 + (0.2e1 * t134 * t286 + (t151 * t165 - t151 - (-t135 * t165 * t242 + t252) * t176) * t145 * t163 + (-(-0.2e1 * t173 * t224 - t135) * t282 + (-(t135 + t244) * t176 + t289 * t173) * t145 * t165) * t164) * t140) * t176, t133 * t176 * t226 + (-(-t129 * t273 + (t135 * t276 - t153 * t164) * t143) * t282 + (t130 * t251 + t280) * t133 + (-t144 * t265 - (-t143 * t274 + t163 * t239 + t271) * t282) * t258) * t140 + (t227 * t265 + ((-t144 * t260 - (qJD(2) * t231 - t135) * t246) * t200 + (t144 * t230 + (t201 * t130 - (-t129 + t262) * t275 - (t135 * t231 - qJD(2)) * t272) * t145) * t197) * t140) * t196, 0 (t131 * t282 - t144 * t177) * t255 + (t131 * t280 + t152 * t144 + (t131 * t251 - t145 * t177) * t130 - (t199 * t261 - t197 * t259 - t127 * t173 - t139 * t153 + (-t139 * t267 - t175) * t135) * t145 * t272 - (-t154 + (-t127 * t196 - t135 * t199) * t197 - t212 * t139) * t246) * t140, 0, 0; (t250 * t291 - t248) * t160 + (-t138 * t281 + t156 * t254) * (-t175 * t195 + t222) + ((t195 * t218 + t198 * t210) * t156 - (-t195 * t210 + t198 * t218) * t278 + 0.2e1 * t160 * t247) * t148 (-t156 * t169 + t157 * t277) * t253 + (-t170 * t137 * t157 + (-t157 * t169 + 0.2e1 * t158 * t277) * t138 + ((t198 * t238 + t199 * t223) * t156 - (-t195 * t238 + t199 * t222) * t278) * qJD(1) + (((t195 * t259 + t290 * t198) * t156 - (-t290 * t195 + t198 * t259) * t278) * t197 + t217 * t200 * (-qJD(2) * t199 + qJD(5))) * t201) * t148, 0, t217 * t176 * t253 + (t217 * t151 + ((-qJD(5) * t156 - 0.2e1 * t247) * t198 + (t137 * t198 + (-qJD(5) * t161 + t138) * t195) * t157) * t176) * t148, t254 + (t248 + (-t148 * t283 - t250) * t161) * t291, 0;];
JaD_rot  = t1;
