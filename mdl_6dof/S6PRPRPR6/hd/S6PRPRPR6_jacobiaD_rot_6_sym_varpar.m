% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRPR6_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR6_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:49:13
% EndTime: 2019-02-26 19:49:13
% DurationCPUTime: 0.90s
% Computational Cost: add. (3313->109), mult. (9085->231), div. (559->12), fcn. (11668->13), ass. (0->105)
t217 = sin(pkin(10));
t219 = cos(pkin(10));
t224 = cos(qJ(2));
t220 = cos(pkin(6));
t222 = sin(qJ(2));
t250 = t220 * t222;
t208 = t217 * t224 + t219 * t250;
t202 = t208 * qJD(2);
t218 = sin(pkin(6));
t223 = cos(qJ(4));
t221 = sin(qJ(4));
t249 = t220 * t224;
t234 = -t217 * t222 + t219 * t249;
t232 = t234 * t221;
t179 = qJD(4) * t232 + (qJD(4) * t218 * t219 + t202) * t223;
t254 = t218 * t221;
t192 = t219 * t254 - t234 * t223;
t189 = t192 ^ 2;
t251 = t218 * t224;
t211 = t220 * t221 + t223 * t251;
t206 = 0.1e1 / t211 ^ 2;
t184 = t189 * t206 + 0.1e1;
t182 = 0.1e1 / t184;
t212 = t220 * t223 - t221 * t251;
t253 = t218 * t222;
t239 = qJD(2) * t253;
t195 = t212 * qJD(4) - t223 * t239;
t205 = 0.1e1 / t211;
t259 = t192 * t206;
t154 = (t179 * t205 - t195 * t259) * t182;
t185 = atan2(t192, t211);
t180 = sin(t185);
t181 = cos(t185);
t237 = -t180 * t211 + t181 * t192;
t150 = t237 * t154 + t180 * t179 + t181 * t195;
t166 = t180 * t192 + t181 * t211;
t163 = 0.1e1 / t166;
t164 = 0.1e1 / t166 ^ 2;
t272 = t150 * t163 * t164;
t209 = t217 * t249 + t219 * t222;
t190 = -t209 * t223 + t217 * t254;
t271 = 0.2e1 * t190 * t272;
t257 = t195 * t205 * t206;
t270 = (t179 * t259 - t189 * t257) / t184 ^ 2;
t241 = t192 * t253;
t233 = t205 * t208 + t206 * t241;
t269 = t223 * t233;
t252 = t218 * t223;
t191 = t209 * t221 + t217 * t252;
t240 = t217 * t250;
t210 = t219 * t224 - t240;
t216 = pkin(11) + qJ(6);
t214 = sin(t216);
t215 = cos(t216);
t175 = t191 * t215 + t210 * t214;
t171 = 0.1e1 / t175;
t172 = 0.1e1 / t175 ^ 2;
t248 = qJD(2) * t224;
t204 = -qJD(2) * t240 + t219 * t248;
t176 = -t190 * qJD(4) + t204 * t221;
t203 = t209 * qJD(2);
t161 = t175 * qJD(6) + t176 * t214 + t203 * t215;
t174 = t191 * t214 - t210 * t215;
t170 = t174 ^ 2;
t169 = t170 * t172 + 0.1e1;
t264 = t172 * t174;
t246 = qJD(6) * t174;
t162 = t176 * t215 - t203 * t214 - t246;
t267 = t162 * t171 * t172;
t268 = (t161 * t264 - t170 * t267) / t169 ^ 2;
t266 = t164 * t190;
t265 = t171 * t214;
t263 = t174 * t215;
t177 = t191 * qJD(4) - t204 * t223;
t262 = t177 * t164;
t261 = t180 * t190;
t260 = t181 * t190;
t258 = t192 * t212;
t256 = t210 * t221;
t255 = t210 * t223;
t247 = qJD(4) * t221;
t188 = t190 ^ 2;
t160 = t164 * t188 + 0.1e1;
t245 = 0.2e1 * (-t188 * t272 + t190 * t262) / t160 ^ 2;
t244 = -0.2e1 * t268;
t242 = t174 * t267;
t238 = qJD(6) * t256 + t204;
t236 = t172 * t263 - t265;
t193 = t219 * t252 + t232;
t235 = -t193 * t205 + t206 * t258;
t231 = qJD(4) * t255 - qJD(6) * t209 - t203 * t221;
t201 = t234 * qJD(2);
t194 = -t211 * qJD(4) + t221 * t239;
t187 = -t209 * t214 + t215 * t256;
t186 = t209 * t215 + t214 * t256;
t178 = t192 * qJD(4) + t202 * t221;
t167 = 0.1e1 / t169;
t157 = 0.1e1 / t160;
t156 = t182 * t269;
t155 = t235 * t182;
t152 = (t180 * t208 - t181 * t253) * t223 + t237 * t156;
t151 = -t237 * t155 + t180 * t193 + t181 * t212;
t149 = 0.2e1 * t235 * t270 + (0.2e1 * t257 * t258 - t178 * t205 + (-t179 * t212 - t192 * t194 - t193 * t195) * t206) * t182;
t147 = -0.2e1 * t269 * t270 + (-t233 * t247 + (-0.2e1 * t241 * t257 + t201 * t205 + (-t195 * t208 + (t179 * t222 + t192 * t248) * t218) * t206) * t223) * t182;
t1 = [0, t147, 0, t149, 0, 0; 0 (t152 * t266 + t163 * t255) * t245 + ((t203 * t223 + t210 * t247) * t163 + (-t262 + t271) * t152 + (t255 * t150 - (t147 * t192 + t156 * t179 + (t222 * t247 - t223 * t248) * t218 + (-t156 * t211 + t208 * t223) * t154) * t260 - (-t208 * t247 - t147 * t211 - t156 * t195 + t201 * t223 + (-t156 * t192 + t222 * t252) * t154) * t261) * t164) * t157, 0 (t151 * t266 - t163 * t191) * t245 + (t151 * t271 + t176 * t163 + (-t191 * t150 - t151 * t177 - (t149 * t192 - t155 * t179 + t194 + (t155 * t211 + t193) * t154) * t260 - (-t149 * t211 + t155 * t195 - t178 + (t155 * t192 - t212) * t154) * t261) * t164) * t157, 0, 0; 0, 0.2e1 * (-t171 * t186 + t187 * t264) * t268 + (0.2e1 * t187 * t242 + t238 * t171 * t215 + t231 * t265 + (t238 * t174 * t214 - t187 * t161 - t186 * t162 - t231 * t263) * t172) * t167, 0, t236 * t190 * t244 + (t236 * t177 + ((-qJD(6) * t171 - 0.2e1 * t242) * t215 + (t161 * t215 + (t162 - t246) * t214) * t172) * t190) * t167, 0, t244 + 0.2e1 * (t161 * t167 * t172 + (-t167 * t267 - t172 * t268) * t174) * t174;];
JaD_rot  = t1;
