% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% JaD_rot [3x7]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S7RRRRRRR1_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobiaD_rot_4_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_jacobiaD_rot_4_sym_varpar: qJD has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobiaD_rot_4_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:54:20
% EndTime: 2019-02-26 22:54:22
% DurationCPUTime: 1.21s
% Computational Cost: add. (2534->146), mult. (7495->313), div. (1129->14), fcn. (9361->11), ass. (0->130)
t197 = cos(qJ(3));
t193 = sin(qJ(3));
t199 = cos(qJ(1));
t252 = t199 * t193;
t195 = sin(qJ(1));
t198 = cos(qJ(2));
t254 = t195 * t198;
t173 = t197 * t254 + t252;
t174 = t195 * t197 + t198 * t252;
t194 = sin(qJ(2));
t257 = t194 * t195;
t227 = qJD(2) * t257;
t152 = t174 * qJD(1) + t173 * qJD(3) - t193 * t227;
t230 = t193 * t254;
t251 = t199 * t197;
t172 = t230 - t251;
t170 = t172 ^ 2;
t187 = 0.1e1 / t193 ^ 2;
t190 = 0.1e1 / t194 ^ 2;
t261 = t187 * t190;
t167 = t170 * t261 + 0.1e1;
t186 = 0.1e1 / t193;
t188 = t186 * t187;
t189 = 0.1e1 / t194;
t191 = t189 * t190;
t232 = t172 * t261;
t244 = qJD(3) * t197;
t247 = qJD(2) * t198;
t281 = (t152 * t232 + (-t187 * t191 * t247 - t188 * t190 * t244) * t170) / t167 ^ 2;
t260 = t187 * t197;
t216 = t172 * t260 - t173 * t186;
t280 = t189 * t216;
t255 = t194 * t199;
t176 = -t195 * t193 + t198 * t251;
t192 = sin(qJ(4));
t196 = cos(qJ(4));
t161 = t176 * t196 + t192 * t255;
t156 = 0.1e1 / t161 ^ 2;
t160 = t176 * t192 - t196 * t255;
t253 = t196 * t160;
t155 = 0.1e1 / t161;
t259 = t192 * t155;
t217 = -t156 * t253 + t259;
t279 = -qJD(4) * t197 + qJD(2);
t278 = -0.2e1 * t281;
t258 = t194 * t193;
t164 = atan2(t172, -t258);
t163 = cos(t164);
t162 = sin(t164);
t268 = t162 * t172;
t146 = -t163 * t258 + t268;
t142 = 0.1e1 / t146;
t143 = 0.1e1 / t146 ^ 2;
t171 = t174 ^ 2;
t141 = t171 * t143 + 0.1e1;
t223 = qJD(1) * t198 + qJD(3);
t246 = qJD(2) * t199;
t212 = t194 * t246 + t223 * t195;
t243 = qJD(3) * t198;
t226 = t197 * t243;
t248 = qJD(1) * t199;
t150 = t212 * t193 - t197 * t248 - t199 * t226;
t271 = t150 * t143;
t165 = 0.1e1 / t167;
t213 = -t193 * t247 - t194 * t244;
t262 = t186 * t189;
t235 = t152 * t262;
t134 = (-t213 * t232 - t235) * t165;
t211 = -t134 * t172 - t213;
t129 = (t134 * t258 + t152) * t162 - t211 * t163;
t144 = t142 * t143;
t276 = t129 * t144;
t277 = (-t171 * t276 - t174 * t271) / t141 ^ 2;
t249 = qJD(1) * t195;
t209 = qJD(4) * t176 + t194 * t249 - t198 * t246;
t151 = (qJD(1) + t243) * t252 + t212 * t197;
t242 = qJD(4) * t194;
t220 = t199 * t242 - t151;
t136 = t220 * t192 + t209 * t196;
t154 = t160 ^ 2;
t149 = t154 * t156 + 0.1e1;
t270 = t156 * t160;
t137 = -t209 * t192 + t220 * t196;
t157 = t155 * t156;
t273 = t137 * t157;
t275 = (t136 * t270 - t154 * t273) / t149 ^ 2;
t274 = t136 * t156;
t272 = t143 * t174;
t256 = t194 * t197;
t215 = -t192 * t198 + t196 * t256;
t169 = t215 * t199;
t269 = t160 * t169;
t267 = t162 * t174;
t266 = t162 * t194;
t265 = t163 * t172;
t264 = t163 * t174;
t263 = t163 * t198;
t231 = t186 * t190 * t198;
t218 = t172 * t231 + t195;
t145 = t218 * t165;
t250 = t145 - t195;
t245 = qJD(3) * t193;
t240 = -0.2e1 * t277;
t239 = 0.2e1 * t275;
t238 = 0.2e1 * t281;
t237 = -0.2e1 * t144 * t174;
t236 = t143 * t267;
t234 = t160 * t273;
t233 = t172 * t262;
t228 = t190 * t247;
t225 = t142 * t240;
t224 = t143 * t240;
t221 = t238 * t262;
t153 = -qJD(3) * t230 - t193 * t249 - t197 * t227 + t223 * t251;
t219 = -t195 * t242 - t153;
t214 = t192 * t256 + t196 * t198;
t210 = -qJD(4) * t173 + t194 * t248 + t195 * t247;
t208 = t150 * t262 + (t187 * t189 * t244 + t186 * t228) * t174;
t168 = t214 * t199;
t159 = -t173 * t196 - t192 * t257;
t158 = -t173 * t192 + t196 * t257;
t147 = 0.1e1 / t149;
t139 = 0.1e1 / t141;
t138 = t165 * t280;
t133 = (t162 + (-t163 * t233 - t162) * t165) * t174;
t132 = t145 * t265 + (t250 * t266 - t263) * t193;
t130 = -t163 * t256 + t162 * t173 + (t162 * t258 + t265) * t138;
t128 = t218 * t278 + (t152 * t231 + t248 + (-t226 * t261 + (-0.2e1 * t191 * t198 ^ 2 - t189) * t186 * qJD(2)) * t172) * t165;
t126 = t278 * t280 + (-t216 * t228 + (t152 * t260 - t153 * t186 + (t173 * t260 + (-0.2e1 * t188 * t197 ^ 2 - t186) * t172) * qJD(3)) * t189) * t165;
t1 = [t208 * t165 + t174 * t221, t128, t126, 0, 0, 0, 0; t172 * t225 + (t152 * t142 + (-t129 * t172 - t133 * t150) * t143) * t139 + (t133 * t224 + (-0.2e1 * t133 * t276 + (t150 * t165 - t150 + (t134 * t165 * t233 + t238) * t174) * t143 * t162 + ((t172 * t221 + t134) * t272 + ((-t134 - t235) * t174 + t208 * t172) * t143 * t165) * t163) * t139) * t174, t132 * t174 * t224 + ((t128 * t265 + (-t134 * t268 + t152 * t163) * t145) * t272 + (t129 * t237 - t271) * t132 + (t142 * t255 + (t145 * t266 - t162 * t257 - t263) * t272) * t244) * t139 + (t225 * t255 + ((t142 * t246 + (t250 * qJD(2) + t134) * t236) * t198 + (-t142 * t249 + (-t199 * t129 + (t128 - t248) * t267 + (t250 * t134 + qJD(2)) * t264) * t143) * t194) * t139) * t193, 0.2e1 * (-t130 * t272 + t142 * t176) * t277 + (-t130 * t271 + t151 * t142 + (t130 * t237 + t176 * t143) * t129 + (-t197 * t247 + t194 * t245 + t126 * t172 + t138 * t152 + (t138 * t258 + t173) * t134) * t143 * t264 + (t153 + (t126 * t193 + t134 * t197) * t194 + t211 * t138) * t236) * t139, 0, 0, 0, 0; (-t155 * t158 + t159 * t270) * t239 + (0.2e1 * t159 * t234 + t219 * t259 + t210 * t155 * t196 + (t210 * t160 * t192 - t159 * t136 - t158 * t137 - t219 * t253) * t156) * t147 (t155 * t168 - t156 * t269) * t239 + (t169 * t274 + (t168 * t156 - 0.2e1 * t157 * t269) * t137 + (t214 * t155 - t215 * t270) * t249 + (((t192 * t245 + t279 * t196) * t155 - (-t279 * t192 + t196 * t245) * t270) * t194 + t217 * t198 * (-qJD(2) * t197 + qJD(4))) * t199) * t147, t217 * t174 * t239 + (t217 * t150 + ((-qJD(4) * t155 - 0.2e1 * t234) * t196 + (t136 * t196 + (-qJD(4) * t160 + t137) * t192) * t156) * t174) * t147, -0.2e1 * t275 + 0.2e1 * (t147 * t274 + (-t147 * t273 - t156 * t275) * t160) * t160, 0, 0, 0;];
JaD_rot  = t1;
