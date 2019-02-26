% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRP7_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP7_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:11:19
% EndTime: 2019-02-26 21:11:21
% DurationCPUTime: 1.19s
% Computational Cost: add. (9776->126), mult. (8378->269), div. (1558->15), fcn. (10537->9), ass. (0->118)
t200 = pkin(10) + qJ(3);
t197 = cos(t200);
t203 = qJ(4) + qJ(5);
t198 = sin(t203);
t275 = sin(qJ(1));
t234 = t275 * t198;
t199 = cos(t203);
t204 = cos(qJ(1));
t254 = t204 * t199;
t178 = t197 * t254 + t234;
t172 = 0.1e1 / t178 ^ 2;
t196 = sin(t200);
t189 = t196 ^ 2;
t202 = t204 ^ 2;
t263 = t189 * t202;
t242 = t172 * t263;
t168 = 0.1e1 + t242;
t201 = qJD(4) + qJD(5);
t227 = qJD(1) * t275;
t251 = qJD(3) * t204;
t229 = t196 * t251;
t214 = t197 * t227 + t229;
t232 = t275 * t201;
t255 = t204 * t198;
t157 = (-t197 * t201 + qJD(1)) * t255 + (t232 - t214) * t199;
t171 = 0.1e1 / t178;
t270 = t157 * t171 * t172;
t222 = t263 * t270;
t230 = qJD(3) * t196 * t202;
t279 = (-t222 + (-t189 * t204 * t227 + t197 * t230) * t172) / t168 ^ 2;
t194 = 0.1e1 / t198 ^ 2;
t256 = t199 * t201;
t238 = t194 * t256;
t258 = t196 * t204;
t174 = t197 * t234 + t254;
t221 = t198 * t232;
t236 = t201 * t254;
t156 = t174 * qJD(1) - t197 * t236 + t198 * t229 - t221;
t233 = t275 * t199;
t177 = t197 * t255 - t233;
t190 = 0.1e1 / t196;
t193 = 0.1e1 / t198;
t191 = 0.1e1 / t196 ^ 2;
t252 = qJD(3) * t197;
t231 = t191 * t252;
t262 = t190 * t193;
t278 = (t190 * t238 + t193 * t231) * t177 + t156 * t262;
t259 = t196 * t198;
t164 = atan2(-t174, t259);
t161 = cos(t164);
t160 = sin(t164);
t269 = t160 * t174;
t155 = t161 * t259 - t269;
t152 = 0.1e1 / t155;
t153 = 0.1e1 / t155 ^ 2;
t277 = -0.2e1 * t174;
t276 = 0.2e1 * t177;
t169 = t174 ^ 2;
t261 = t191 * t194;
t165 = t169 * t261 + 0.1e1;
t162 = 0.1e1 / t165;
t215 = t196 * t256 + t198 * t252;
t240 = t174 * t261;
t235 = t196 * t275;
t219 = qJD(3) * t235;
t220 = t199 * t227;
t253 = qJD(1) * t204;
t158 = -t198 * t219 - t201 * t255 - t220 + (t198 * t253 + t199 * t232) * t197;
t243 = t158 * t262;
t144 = (t215 * t240 - t243) * t162;
t212 = -t144 * t174 + t215;
t139 = (-t144 * t259 - t158) * t160 + t212 * t161;
t154 = t152 * t153;
t274 = t139 * t154;
t192 = t190 / t189;
t237 = t193 * t238;
t273 = (t158 * t240 + (-t192 * t194 * t252 - t191 * t237) * t169) / t165 ^ 2;
t272 = t153 * t177;
t271 = t156 * t153;
t268 = t160 * t177;
t267 = t160 * t196;
t266 = t161 * t174;
t265 = t161 * t177;
t264 = t161 * t197;
t260 = t191 * t197;
t257 = t198 * t201;
t170 = t177 ^ 2;
t150 = t153 * t170 + 0.1e1;
t250 = 0.2e1 * (-t170 * t274 - t177 * t271) / t150 ^ 2;
t249 = -0.2e1 * t273;
t248 = 0.2e1 * t279;
t247 = t154 * t276;
t246 = t190 * t273;
t245 = t153 * t268;
t241 = t174 * t262;
t239 = t193 * t260;
t217 = t174 * t239 + t275;
t151 = t217 * t162;
t228 = t275 - t151;
t226 = t152 * t250;
t225 = t153 * t250;
t224 = t258 * t276;
t223 = t193 * t246;
t176 = t197 * t233 - t255;
t218 = t174 * t194 * t199 - t176 * t193;
t216 = t172 * t176 * t204 - t275 * t171;
t166 = 0.1e1 / t168;
t159 = t178 * qJD(1) - t197 * t221 - t199 * t219 - t236;
t148 = 0.1e1 / t150;
t147 = t218 * t190 * t162;
t143 = (-t160 + (t161 * t241 + t160) * t162) * t177;
t142 = -t151 * t266 + (t228 * t267 + t264) * t198;
t141 = t161 * t196 * t199 - t160 * t176 + (-t160 * t259 - t266) * t147;
t140 = t172 * t224 * t279 + (t224 * t270 + (t156 * t258 + (t196 * t227 - t197 * t251) * t177) * t172) * t166;
t138 = t217 * t249 + (t158 * t239 + t253 + (-t238 * t260 + (-0.2e1 * t192 * t197 ^ 2 - t190) * t193 * qJD(3)) * t174) * t162;
t136 = -0.2e1 * t218 * t246 + (-t218 * t231 + ((-t174 * t201 - t159) * t193 + (t237 * t277 + (t176 * t201 + t158) * t194) * t199) * t190) * t162;
t135 = (t141 * t272 - t152 * t178) * t250 + (t141 * t271 + t157 * t152 + (t141 * t247 - t153 * t178) * t139 - (t199 * t252 - t196 * t257 - t136 * t174 - t147 * t158 + (-t147 * t259 - t176) * t144) * t153 * t265 - (-t159 + (-t136 * t198 - t144 * t199) * t196 - t212 * t147) * t245) * t148;
t1 = [t162 * t278 + t223 * t276, 0, t138, t136, t136, 0; t174 * t226 + (-t158 * t152 + (t139 * t174 + t143 * t156) * t153) * t148 + (t143 * t225 + (0.2e1 * t143 * t274 + (t156 * t162 - t156 - (-t144 * t162 * t241 + t249) * t177) * t153 * t160 + (-(t223 * t277 - t144) * t272 + (-(t144 + t243) * t177 + t278 * t174) * t153 * t162) * t161) * t148) * t177, 0, t142 * t177 * t225 + (-(-t138 * t266 + (t144 * t269 - t158 * t161) * t151) * t272 + (-t152 * t258 - (-t151 * t267 + t160 * t235 + t264) * t272) * t256 + (t139 * t247 + t271) * t142) * t148 + (t226 * t258 + ((-t152 * t251 - (t228 * qJD(3) - t144) * t245) * t197 + (t152 * t227 + (t204 * t139 - (-t138 + t253) * t268 - (t228 * t144 - qJD(3)) * t265) * t153) * t196) * t148) * t198, t135, t135, 0; t216 * t196 * t248 + (-t216 * t252 + ((qJD(1) * t171 + 0.2e1 * t176 * t270) * t204 + (-t275 * t157 - t159 * t204 + t176 * t227) * t172) * t196) * t166, 0 (t171 * t197 * t204 + t199 * t242) * t248 + (0.2e1 * t199 * t222 + t214 * t171 + ((t157 * t204 - 0.2e1 * t199 * t230) * t197 + (t202 * t257 + 0.2e1 * t204 * t220) * t189) * t172) * t166, t140, t140, 0;];
JaD_rot  = t1;
