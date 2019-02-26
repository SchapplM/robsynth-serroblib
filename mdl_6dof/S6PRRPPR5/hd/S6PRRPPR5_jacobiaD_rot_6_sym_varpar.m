% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPPR5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:00:37
% EndTime: 2019-02-26 20:00:38
% DurationCPUTime: 0.95s
% Computational Cost: add. (3313->108), mult. (9085->220), div. (559->12), fcn. (11668->13), ass. (0->103)
t220 = cos(pkin(6));
t224 = cos(qJ(2));
t272 = cos(pkin(10));
t245 = t272 * t224;
t218 = sin(pkin(10));
t222 = sin(qJ(2));
t259 = t218 * t222;
t206 = t220 * t245 - t259;
t199 = t206 * qJD(2);
t223 = cos(qJ(3));
t221 = sin(qJ(3));
t246 = t272 * t222;
t258 = t218 * t224;
t233 = -t220 * t246 - t258;
t219 = sin(pkin(6));
t247 = t219 * t272;
t274 = t233 * t221 - t223 * t247;
t177 = t274 * qJD(3) + t199 * t223;
t192 = -t221 * t247 - t223 * t233;
t189 = t192 ^ 2;
t256 = t219 * t223;
t210 = t220 * t221 + t222 * t256;
t204 = 0.1e1 / t210 ^ 2;
t184 = t189 * t204 + 0.1e1;
t182 = 0.1e1 / t184;
t257 = t219 * t221;
t209 = t220 * t223 - t222 * t257;
t255 = t219 * t224;
t248 = qJD(2) * t255;
t197 = qJD(3) * t209 + t223 * t248;
t203 = 0.1e1 / t210;
t263 = t192 * t204;
t154 = (-t177 * t203 + t197 * t263) * t182;
t185 = atan2(-t192, t210);
t180 = sin(t185);
t181 = cos(t185);
t240 = -t180 * t210 - t181 * t192;
t150 = t154 * t240 - t180 * t177 + t181 * t197;
t166 = -t180 * t192 + t181 * t210;
t163 = 0.1e1 / t166;
t164 = 0.1e1 / t166 ^ 2;
t278 = t150 * t163 * t164;
t217 = pkin(11) + qJ(6);
t215 = sin(t217);
t216 = cos(t217);
t234 = -t220 * t258 - t246;
t208 = -t220 * t259 + t245;
t236 = -t208 * t221 + t218 * t256;
t239 = t215 * t234 - t216 * t236;
t277 = qJD(6) * t239;
t195 = t208 * t223 + t218 * t257;
t276 = 0.2e1 * t195 * t278;
t235 = -t203 * t206 + t255 * t263;
t275 = t223 * t235;
t262 = t197 * t203 * t204;
t273 = -0.2e1 * (t177 * t263 - t189 * t262) / t184 ^ 2;
t175 = -t215 * t236 - t216 * t234;
t171 = 0.1e1 / t175;
t172 = 0.1e1 / t175 ^ 2;
t201 = t234 * qJD(2);
t178 = qJD(3) * t195 + t201 * t221;
t202 = t208 * qJD(2);
t162 = t178 * t215 + t202 * t216 + t277;
t271 = t162 * t171 * t172;
t270 = t164 * t195;
t161 = qJD(6) * t175 - t178 * t216 + t202 * t215;
t170 = t239 ^ 2;
t169 = t170 * t172 + 0.1e1;
t267 = t172 * t239;
t269 = 0.1e1 / t169 ^ 2 * (-t161 * t267 - t170 * t271);
t268 = t171 * t216;
t266 = t239 * t215;
t265 = t180 * t195;
t264 = t181 * t195;
t261 = t234 * t221;
t260 = t234 * t223;
t254 = qJD(2) * t222;
t253 = qJD(3) * t221;
t190 = t195 ^ 2;
t160 = t164 * t190 + 0.1e1;
t179 = qJD(3) * t236 + t201 * t223;
t252 = 0.2e1 * (t179 * t270 - t190 * t278) / t160 ^ 2;
t250 = 0.2e1 * t269;
t244 = -0.2e1 * t239 * t271;
t243 = -0.2e1 * t192 * t262;
t241 = qJD(6) * t261 + t201;
t238 = -t172 * t266 + t268;
t237 = -t203 * t274 + t209 * t263;
t232 = -qJD(3) * t260 + qJD(6) * t208 + t202 * t221;
t200 = t233 * qJD(2);
t196 = -qJD(3) * t210 - t221 * t248;
t187 = t208 * t216 + t215 * t261;
t186 = t208 * t215 - t216 * t261;
t176 = qJD(3) * t192 + t199 * t221;
t167 = 0.1e1 / t169;
t157 = 0.1e1 / t160;
t156 = t182 * t275;
t155 = t237 * t182;
t152 = (-t180 * t206 + t181 * t255) * t223 + t240 * t156;
t151 = t155 * t240 - t180 * t274 + t181 * t209;
t149 = t237 * t273 + (t209 * t243 + t176 * t203 + (t177 * t209 + t192 * t196 + t197 * t274) * t204) * t182;
t147 = t273 * t275 + (-t235 * t253 + (t243 * t255 - t200 * t203 + (t197 * t206 + (t177 * t224 - t192 * t254) * t219) * t204) * t223) * t182;
t1 = [0, t147, t149, 0, 0, 0; 0 (t152 * t270 - t163 * t260) * t252 + ((-t202 * t223 - t234 * t253) * t163 + t152 * t276 + (-t152 * t179 - t260 * t150 - (-t147 * t192 - t156 * t177 + (-t223 * t254 - t224 * t253) * t219 + (-t156 * t210 - t206 * t223) * t154) * t264 - (t206 * t253 - t147 * t210 - t156 * t197 - t200 * t223 + (t156 * t192 - t223 * t255) * t154) * t265) * t164) * t157 (t151 * t270 - t163 * t236) * t252 + (t151 * t276 - t178 * t163 + (-t236 * t150 - t151 * t179 - (-t149 * t192 - t155 * t177 + t196 + (-t155 * t210 - t274) * t154) * t264 - (-t149 * t210 - t155 * t197 + t176 + (t155 * t192 - t209) * t154) * t265) * t164) * t157, 0, 0, 0; 0 (-t171 * t186 - t187 * t267) * t250 + (t187 * t244 + t241 * t171 * t215 + t232 * t268 + (t216 * t239 * t241 - t187 * t161 - t186 * t162 - t232 * t266) * t172) * t167, t238 * t195 * t250 + (-t238 * t179 + ((qJD(6) * t171 + t244) * t215 + (-t161 * t215 + (t162 + t277) * t216) * t172) * t195) * t167, 0, 0, -0.2e1 * t269 - 0.2e1 * (t161 * t172 * t167 - (-t167 * t271 - t172 * t269) * t239) * t239;];
JaD_rot  = t1;
