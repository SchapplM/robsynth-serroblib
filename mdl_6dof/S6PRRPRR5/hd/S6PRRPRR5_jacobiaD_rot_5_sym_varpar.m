% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRR5_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR5_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:06:37
% EndTime: 2019-02-26 20:06:38
% DurationCPUTime: 0.87s
% Computational Cost: add. (3313->110), mult. (9085->226), div. (559->12), fcn. (11668->13), ass. (0->103)
t218 = sin(pkin(11));
t220 = cos(pkin(11));
t223 = sin(qJ(2));
t221 = cos(pkin(6));
t225 = cos(qJ(2));
t251 = t221 * t225;
t205 = -t218 * t223 + t220 * t251;
t198 = t205 * qJD(2);
t252 = t221 * t223;
t206 = t218 * t225 + t220 * t252;
t222 = sin(qJ(3));
t219 = sin(pkin(6));
t255 = t219 * t222;
t241 = t220 * t255;
t224 = cos(qJ(3));
t248 = qJD(3) * t224;
t176 = -qJD(3) * t241 + t198 * t222 + t206 * t248;
t254 = t219 * t224;
t190 = t206 * t222 + t220 * t254;
t188 = t190 ^ 2;
t209 = -t221 * t224 + t223 * t255;
t203 = 0.1e1 / t209 ^ 2;
t184 = t188 * t203 + 0.1e1;
t182 = 0.1e1 / t184;
t210 = t221 * t222 + t223 * t254;
t249 = qJD(2) * t225;
t240 = t219 * t249;
t195 = qJD(3) * t210 + t222 * t240;
t202 = 0.1e1 / t209;
t259 = t190 * t203;
t154 = (-t176 * t202 + t195 * t259) * t182;
t185 = atan2(-t190, t209);
t180 = sin(t185);
t181 = cos(t185);
t237 = -t180 * t209 - t181 * t190;
t150 = t154 * t237 - t176 * t180 + t181 * t195;
t166 = -t180 * t190 + t181 * t209;
t163 = 0.1e1 / t166;
t164 = 0.1e1 / t166 ^ 2;
t271 = t150 * t163 * t164;
t253 = t219 * t225;
t233 = -t202 * t205 + t253 * t259;
t270 = t222 * t233;
t242 = t218 * t252;
t208 = t220 * t225 - t242;
t234 = -t208 * t222 + t218 * t254;
t269 = -0.2e1 * t234 * t271;
t258 = t195 * t202 * t203;
t268 = -0.2e1 * (t176 * t259 - t188 * t258) / t184 ^ 2;
t194 = t208 * t224 + t218 * t255;
t207 = t218 * t251 + t220 * t223;
t217 = pkin(12) + qJ(5);
t215 = sin(t217);
t216 = cos(t217);
t175 = t194 * t216 + t207 * t215;
t171 = 0.1e1 / t175;
t172 = 0.1e1 / t175 ^ 2;
t200 = t207 * qJD(2);
t179 = qJD(3) * t234 - t200 * t224;
t201 = -qJD(2) * t242 + t220 * t249;
t161 = qJD(5) * t175 + t179 * t215 - t201 * t216;
t174 = t194 * t215 - t207 * t216;
t170 = t174 ^ 2;
t169 = t170 * t172 + 0.1e1;
t263 = t172 * t174;
t247 = qJD(5) * t174;
t162 = t179 * t216 + t201 * t215 - t247;
t266 = t162 * t171 * t172;
t267 = (t161 * t263 - t170 * t266) / t169 ^ 2;
t265 = t164 * t234;
t264 = t171 * t215;
t262 = t174 * t216;
t261 = t180 * t234;
t260 = t181 * t234;
t257 = t207 * t222;
t256 = t207 * t224;
t250 = qJD(2) * t223;
t189 = t234 ^ 2;
t160 = t164 * t189 + 0.1e1;
t178 = qJD(3) * t194 - t200 * t222;
t246 = 0.2e1 * (-t178 * t265 - t189 * t271) / t160 ^ 2;
t245 = -0.2e1 * t267;
t243 = t174 * t266;
t239 = -0.2e1 * t190 * t258;
t238 = qJD(5) * t256 - t200;
t236 = t172 * t262 - t264;
t192 = t206 * t224 - t241;
t235 = -t192 * t202 + t210 * t259;
t232 = qJD(3) * t257 + qJD(5) * t208 - t201 * t224;
t199 = t206 * qJD(2);
t196 = -qJD(3) * t209 + t224 * t240;
t187 = t208 * t215 - t216 * t256;
t186 = -t208 * t216 - t215 * t256;
t177 = -qJD(3) * t190 + t198 * t224;
t167 = 0.1e1 / t169;
t157 = 0.1e1 / t160;
t156 = t182 * t270;
t155 = t235 * t182;
t152 = (-t180 * t205 + t181 * t253) * t222 + t237 * t156;
t151 = t155 * t237 - t180 * t192 + t181 * t210;
t149 = t235 * t268 + (t210 * t239 - t177 * t202 + (t176 * t210 + t190 * t196 + t192 * t195) * t203) * t182;
t147 = t268 * t270 + (t233 * t248 + (t239 * t253 + t199 * t202 + (t195 * t205 + (t176 * t225 - t190 * t250) * t219) * t203) * t222) * t182;
t1 = [0, t147, t149, 0, 0, 0; 0 (-t152 * t265 + t163 * t257) * t246 + ((-t201 * t222 - t207 * t248) * t163 + t152 * t269 + (-t152 * t178 + t257 * t150 + (-t147 * t190 - t156 * t176 + (-t222 * t250 + t225 * t248) * t219 + (-t156 * t209 - t205 * t222) * t154) * t260 + (-t205 * t248 - t147 * t209 - t156 * t195 + t199 * t222 + (t156 * t190 - t222 * t253) * t154) * t261) * t164) * t157 (-t151 * t265 - t163 * t194) * t246 + (t151 * t269 + t179 * t163 + (-t194 * t150 - t151 * t178 + (-t149 * t190 - t155 * t176 + t196 + (-t155 * t209 - t192) * t154) * t260 + (-t149 * t209 - t155 * t195 - t177 + (t155 * t190 - t210) * t154) * t261) * t164) * t157, 0, 0, 0; 0, 0.2e1 * (-t171 * t186 + t187 * t263) * t267 + (0.2e1 * t187 * t243 - t238 * t171 * t216 + t232 * t264 + (-t174 * t215 * t238 - t187 * t161 - t186 * t162 - t232 * t262) * t172) * t167, -t236 * t234 * t245 + (t236 * t178 - ((-qJD(5) * t171 - 0.2e1 * t243) * t216 + (t161 * t216 + (t162 - t247) * t215) * t172) * t234) * t167, 0, t245 + 0.2e1 * (t161 * t167 * t172 + (-t167 * t266 - t172 * t267) * t174) * t174, 0;];
JaD_rot  = t1;
