% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPP3_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP3_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_jacobiaD_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:26:20
% EndTime: 2019-02-26 22:26:23
% DurationCPUTime: 1.10s
% Computational Cost: add. (5428->122), mult. (8127->262), div. (1567->14), fcn. (10302->9), ass. (0->115)
t200 = qJ(2) + qJ(3);
t191 = sin(t200);
t204 = cos(qJ(1));
t277 = t191 * t204;
t187 = 0.1e1 / t191;
t188 = 0.1e1 / t191 ^ 2;
t189 = t187 * t188;
t192 = cos(t200);
t193 = qJD(2) + qJD(3);
t276 = t193 * (0.2e1 * t189 * t192 ^ 2 + t187);
t203 = cos(qJ(4));
t252 = t203 * t204;
t201 = sin(qJ(4));
t202 = sin(qJ(1));
t254 = t202 * t201;
t172 = t192 * t254 + t252;
t244 = qJD(4) * t204;
t224 = t203 * t244;
t246 = qJD(4) * t201;
t225 = t202 * t246;
t228 = t193 * t277;
t156 = qJD(1) * t172 - t192 * t224 + t201 * t228 - t225;
t250 = t204 * t201;
t253 = t202 * t203;
t175 = t192 * t250 - t253;
t195 = 0.1e1 / t201 ^ 2;
t245 = qJD(4) * t203;
t226 = t195 * t245;
t194 = 0.1e1 / t201;
t262 = t188 * t192;
t232 = t194 * t262;
t263 = t187 * t194;
t275 = t175 * (t187 * t226 + t193 * t232) + t156 * t263;
t259 = t191 * t201;
t166 = atan2(-t172, t259);
t161 = cos(t166);
t160 = sin(t166);
t269 = t160 * t172;
t155 = t161 * t259 - t269;
t152 = 0.1e1 / t155;
t197 = 0.1e1 / t204;
t153 = 0.1e1 / t155 ^ 2;
t198 = 0.1e1 / t204 ^ 2;
t274 = 0.2e1 * t175;
t169 = t172 ^ 2;
t261 = t188 * t195;
t167 = t169 * t261 + 0.1e1;
t162 = 0.1e1 / t167;
t257 = t192 * t193;
t215 = t191 * t245 + t201 * t257;
t234 = t172 * t261;
t258 = t191 * t202;
t229 = t193 * t258;
t247 = qJD(1) * t204;
t248 = qJD(1) * t202;
t158 = t202 * t245 * t192 - t203 * t248 + (t247 * t192 - t229 - t244) * t201;
t236 = t158 * t263;
t144 = (t215 * t234 - t236) * t162;
t213 = -t144 * t172 + t215;
t140 = (-t144 * t259 - t158) * t160 + t213 * t161;
t154 = t152 * t153;
t273 = t140 * t154;
t196 = t194 * t195;
t230 = t189 * t257;
t272 = (t158 * t234 + (-t188 * t196 * t245 - t195 * t230) * t169) / t167 ^ 2;
t271 = t153 * t175;
t270 = t156 * t153;
t268 = t160 * t175;
t267 = t160 * t191;
t266 = t161 * t172;
t265 = t161 * t175;
t264 = t161 * t192;
t260 = t188 * t198;
t256 = t195 * t203;
t255 = t198 * t202;
t251 = t204 * t152;
t218 = t172 * t232 + t202;
t151 = t218 * t162;
t249 = -t151 + t202;
t170 = t175 ^ 2;
t150 = t153 * t170 + 0.1e1;
t243 = 0.2e1 * (-t170 * t273 - t175 * t270) / t150 ^ 2;
t242 = -0.2e1 * t272;
t157 = (-qJD(4) * t192 + qJD(1)) * t250 + (-t228 + (-qJD(1) * t192 + qJD(4)) * t202) * t203;
t176 = t192 * t252 + t254;
t171 = t176 ^ 2;
t168 = t171 * t260 + 0.1e1;
t199 = t197 * t198;
t241 = 0.2e1 * (t176 * t157 * t260 + (t188 * t199 * t248 - t198 * t230) * t171) / t168 ^ 2;
t240 = t154 * t274;
t239 = t187 * t272;
t238 = t153 * t268;
t235 = t172 * t263;
t233 = t188 * t257;
t231 = t197 * t262;
t227 = t198 * t248;
t223 = t152 * t243;
t222 = t153 * t243;
t221 = t187 * t241;
t219 = t194 * t239;
t174 = t192 * t253 - t250;
t217 = t172 * t256 - t174 * t194;
t216 = t174 * t197 - t176 * t255;
t164 = 0.1e1 / t168;
t159 = qJD(1) * t176 - t192 * t225 - t203 * t229 - t224;
t148 = 0.1e1 / t150;
t147 = t217 * t187 * t162;
t143 = (-t160 + (t161 * t235 + t160) * t162) * t175;
t142 = -t151 * t266 + (t249 * t267 + t264) * t201;
t141 = t161 * t191 * t203 - t160 * t174 + (-t160 * t259 - t266) * t147;
t139 = t218 * t242 + (t158 * t232 + t247 + (-t194 * t276 - t226 * t262) * t172) * t162;
t138 = (t176 * t231 + t203) * t241 + (-t157 * t231 + t246 + (t197 * t276 - t227 * t262) * t176) * t164;
t136 = -0.2e1 * t217 * t239 + (-t217 * t233 + (t158 * t256 - t159 * t194 + (t174 * t256 + (-0.2e1 * t196 * t203 ^ 2 - t194) * t172) * qJD(4)) * t187) * t162;
t135 = t142 * t175 * t222 + (-(-t139 * t266 + (t144 * t269 - t158 * t161) * t151) * t271 + (t140 * t240 + t270) * t142 + (-t191 * t251 - (-t151 * t267 + t160 * t258 + t264) * t271) * t245) * t148 + (t223 * t277 + ((-t193 * t251 - (t249 * t193 - t144) * t238) * t192 + (t152 * t248 + (t204 * t140 - (-t139 + t247) * t268 - (t249 * t144 - t193) * t265) * t153) * t191) * t148) * t201;
t1 = [t162 * t275 + t219 * t274, t139, t139, t136, 0, 0; t172 * t223 + (-t158 * t152 + (t140 * t172 + t143 * t156) * t153) * t148 + (t143 * t222 + (0.2e1 * t143 * t273 + (t156 * t162 - t156 - (-t144 * t162 * t235 + t242) * t175) * t153 * t160 + (-(-0.2e1 * t172 * t219 - t144) * t271 + (-(t144 + t236) * t175 + t275 * t172) * t153 * t162) * t161) * t148) * t175, t135, t135 (t141 * t271 - t152 * t176) * t243 + (t141 * t270 + t157 * t152 + (t141 * t240 - t153 * t176) * t140 - (-t191 * t246 + t203 * t257 - t136 * t172 - t147 * t158 + (-t147 * t259 - t174) * t144) * t153 * t265 - (-t159 + (-t136 * t201 - t144 * t203) * t191 - t213 * t147) * t238) * t148, 0, 0; t216 * t221 + (t216 * t233 + (t157 * t255 - t159 * t197 + (-t174 * t255 + (0.2e1 * t199 * t202 ^ 2 + t197) * t176) * qJD(1)) * t187) * t164, t138, t138, t175 * t197 * t221 + (t156 * t187 * t197 + (-t187 * t227 + t193 * t231) * t175) * t164, 0, 0;];
JaD_rot  = t1;
