% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRP3_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:51:34
% EndTime: 2019-02-26 19:51:35
% DurationCPUTime: 0.98s
% Computational Cost: add. (5804->110), mult. (9085->227), div. (559->12), fcn. (11668->13), ass. (0->106)
t220 = sin(pkin(10));
t222 = cos(pkin(10));
t225 = sin(qJ(2));
t223 = cos(pkin(6));
t227 = cos(qJ(2));
t253 = t223 * t227;
t209 = -t220 * t225 + t222 * t253;
t205 = t209 * qJD(2);
t254 = t223 * t225;
t210 = t220 * t227 + t222 * t254;
t219 = pkin(11) + qJ(4);
t217 = sin(t219);
t221 = sin(pkin(6));
t257 = t221 * t222;
t243 = t217 * t257;
t218 = cos(t219);
t250 = qJD(4) * t218;
t172 = -qJD(4) * t243 + t205 * t217 + t210 * t250;
t194 = t210 * t217 + t218 * t257;
t192 = t194 ^ 2;
t256 = t221 * t225;
t202 = t217 * t256 - t218 * t223;
t200 = 0.1e1 / t202 ^ 2;
t186 = t192 * t200 + 0.1e1;
t184 = 0.1e1 / t186;
t203 = t217 * t223 + t218 * t256;
t251 = qJD(2) * t227;
t242 = t221 * t251;
t190 = t203 * qJD(4) + t217 * t242;
t199 = 0.1e1 / t202;
t262 = t194 * t200;
t156 = (-t172 * t199 + t190 * t262) * t184;
t187 = atan2(-t194, t202);
t180 = sin(t187);
t181 = cos(t187);
t239 = -t180 * t202 - t181 * t194;
t152 = t239 * t156 - t172 * t180 + t181 * t190;
t166 = -t180 * t194 + t181 * t202;
t163 = 0.1e1 / t166;
t164 = 0.1e1 / t166 ^ 2;
t276 = t152 * t163 * t164;
t244 = t220 * t254;
t212 = t222 * t227 - t244;
t258 = t220 * t221;
t236 = -t212 * t217 + t218 * t258;
t275 = -0.2e1 * t236 * t276;
t255 = t221 * t227;
t235 = -t199 * t209 + t255 * t262;
t274 = t217 * t235;
t263 = t190 * t199 * t200;
t273 = -0.2e1 * (t172 * t262 - t192 * t263) / t186 ^ 2;
t198 = t212 * t218 + t217 * t258;
t226 = cos(qJ(5));
t211 = t220 * t253 + t222 * t225;
t224 = sin(qJ(5));
t260 = t211 * t224;
t183 = t198 * t226 + t260;
t177 = 0.1e1 / t183;
t178 = 0.1e1 / t183 ^ 2;
t207 = t211 * qJD(2);
t175 = t236 * qJD(4) - t207 * t218;
t208 = -qJD(2) * t244 + t222 * t251;
t167 = t183 * qJD(5) + t175 * t224 - t208 * t226;
t259 = t211 * t226;
t182 = t198 * t224 - t259;
t176 = t182 ^ 2;
t171 = t176 * t178 + 0.1e1;
t267 = t178 * t182;
t249 = qJD(5) * t182;
t168 = t175 * t226 + t208 * t224 - t249;
t270 = t168 * t177 * t178;
t272 = (t167 * t267 - t176 * t270) / t171 ^ 2;
t271 = t164 * t236;
t174 = t198 * qJD(4) - t207 * t217;
t269 = t174 * t164;
t268 = t177 * t224;
t266 = t180 * t236;
t265 = t181 * t236;
t264 = t182 * t226;
t261 = t211 * t217;
t252 = qJD(2) * t225;
t193 = t236 ^ 2;
t162 = t164 * t193 + 0.1e1;
t248 = 0.2e1 * (-t193 * t276 - t236 * t269) / t162 ^ 2;
t247 = -0.2e1 * t272;
t245 = t182 * t270;
t241 = -0.2e1 * t194 * t263;
t240 = qJD(5) * t211 * t218 - t207;
t238 = t178 * t264 - t268;
t196 = t210 * t218 - t243;
t237 = -t196 * t199 + t203 * t262;
t234 = qJD(4) * t261 + qJD(5) * t212 - t208 * t218;
t206 = t210 * qJD(2);
t191 = -t202 * qJD(4) + t218 * t242;
t189 = t212 * t224 - t218 * t259;
t188 = -t212 * t226 - t218 * t260;
t173 = -t194 * qJD(4) + t205 * t218;
t169 = 0.1e1 / t171;
t159 = 0.1e1 / t162;
t158 = t184 * t274;
t157 = t237 * t184;
t154 = (-t180 * t209 + t181 * t255) * t217 + t239 * t158;
t153 = t239 * t157 - t180 * t196 + t181 * t203;
t151 = t237 * t273 + (t203 * t241 - t173 * t199 + (t172 * t203 + t190 * t196 + t191 * t194) * t200) * t184;
t149 = t273 * t274 + (t235 * t250 + (t241 * t255 + t199 * t206 + (t190 * t209 + (t172 * t227 - t194 * t252) * t221) * t200) * t217) * t184;
t1 = [0, t149, 0, t151, 0, 0; 0 (-t154 * t271 + t163 * t261) * t248 + ((-t208 * t217 - t211 * t250) * t163 + (-t269 + t275) * t154 + (t261 * t152 + (-t149 * t194 - t158 * t172 + (-t217 * t252 + t227 * t250) * t221 + (-t158 * t202 - t209 * t217) * t156) * t265 + (-t209 * t250 - t149 * t202 - t158 * t190 + t206 * t217 + (t158 * t194 - t217 * t255) * t156) * t266) * t164) * t159, 0 (-t153 * t271 - t163 * t198) * t248 + (t153 * t275 + t175 * t163 + (-t198 * t152 - t153 * t174 + (-t151 * t194 - t157 * t172 + t191 + (-t157 * t202 - t196) * t156) * t265 + (-t151 * t202 - t157 * t190 - t173 + (t157 * t194 - t203) * t156) * t266) * t164) * t159, 0, 0; 0, 0.2e1 * (-t177 * t188 + t189 * t267) * t272 + (0.2e1 * t189 * t245 - t240 * t177 * t226 + t234 * t268 + (-t240 * t182 * t224 - t189 * t167 - t188 * t168 - t234 * t264) * t178) * t169, 0, -t238 * t236 * t247 + (t238 * t174 - ((-qJD(5) * t177 - 0.2e1 * t245) * t226 + (t167 * t226 + (t168 - t249) * t224) * t178) * t236) * t169, t247 + 0.2e1 * (t167 * t169 * t178 + (-t169 * t270 - t178 * t272) * t182) * t182, 0;];
JaD_rot  = t1;
