% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR4_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:11
% EndTime: 2019-02-26 21:30:12
% DurationCPUTime: 0.91s
% Computational Cost: add. (3306->98), mult. (9869->206), div. (448->12), fcn. (12731->13), ass. (0->97)
t223 = cos(qJ(1));
t221 = sin(qJ(1));
t220 = sin(qJ(2));
t272 = cos(pkin(11));
t273 = cos(pkin(6));
t240 = t273 * t272;
t217 = sin(pkin(11));
t244 = t217 * t273;
t274 = cos(qJ(2));
t231 = -t220 * t240 - t274 * t244;
t230 = t221 * t231;
t233 = -t220 * t217 + t274 * t272;
t278 = t223 * t233 + t230;
t211 = -t220 * t244 + t274 * t240;
t203 = t211 * qJD(2);
t234 = t274 * t217 + t220 * t272;
t232 = t234 * qJD(2);
t277 = t221 * t203 + t223 * t232;
t218 = sin(pkin(6));
t209 = t233 * t218;
t202 = qJD(2) * t209;
t210 = t234 * t218;
t206 = 0.1e1 / t210 ^ 2;
t260 = t202 * t206;
t193 = t221 * t233 - t223 * t231;
t180 = atan2(-t193, t210);
t175 = sin(t180);
t176 = cos(t180);
t190 = t193 ^ 2;
t179 = t190 * t206 + 0.1e1;
t177 = 0.1e1 / t179;
t205 = 0.1e1 / t210;
t262 = t193 * t205;
t275 = (t176 * t262 + t175) * t177 - t175;
t164 = -t175 * t193 + t176 * t210;
t161 = 0.1e1 / t164;
t196 = -t221 * t211 - t223 * t234;
t219 = sin(qJ(5));
t222 = cos(qJ(5));
t258 = t218 * t221;
t186 = -t196 * t219 + t222 * t258;
t182 = 0.1e1 / t186;
t162 = 0.1e1 / t164 ^ 2;
t183 = 0.1e1 / t186 ^ 2;
t253 = qJD(1) * t223;
t174 = qJD(1) * t230 + t223 * t203 - t221 * t232 + t233 * t253;
t238 = -t174 * t205 + t193 * t260;
t155 = t238 * t177;
t239 = -t175 * t210 - t176 * t193;
t151 = t239 * t155 - t175 * t174 + t176 * t202;
t271 = t151 * t161 * t162;
t192 = t223 * t211 - t221 * t234;
t261 = t193 * t209;
t236 = -t192 * t205 + t206 * t261;
t156 = t236 * t177;
t152 = t239 * t156 - t175 * t192 + t176 * t209;
t270 = t152 * t278;
t204 = t231 * qJD(2);
t212 = t233 * qJD(2);
t171 = t192 * qJD(1) + t221 * t204 + t223 * t212;
t245 = t218 * t253;
t165 = t186 * qJD(5) - t171 * t222 + t219 * t245;
t185 = t196 * t222 + t219 * t258;
t181 = t185 ^ 2;
t169 = t181 * t183 + 0.1e1;
t263 = t183 * t185;
t252 = qJD(5) * t185;
t166 = t171 * t219 + t222 * t245 - t252;
t267 = t166 * t182 * t183;
t269 = (t165 * t263 - t181 * t267) / t169 ^ 2;
t259 = t205 * t260;
t268 = (t193 * t206 * t174 - t190 * t259) / t179 ^ 2;
t172 = t193 * qJD(1) + t277;
t266 = t172 * t162;
t265 = t175 * t278;
t264 = t176 * t278;
t257 = t218 * t223;
t254 = qJD(1) * t221;
t191 = t278 ^ 2;
t160 = t191 * t162 + 0.1e1;
t251 = 0.2e1 * (-t191 * t271 - t266 * t278) / t160 ^ 2;
t250 = 0.2e1 * t271;
t249 = 0.2e1 * t269;
t248 = -0.2e1 * t268;
t247 = t205 * t268;
t246 = t218 * t254;
t243 = 0.2e1 * t185 * t267;
t237 = t222 * t182 + t219 * t263;
t188 = t192 * t219 + t222 * t257;
t187 = -t192 * t222 + t219 * t257;
t201 = t218 * t232;
t173 = t196 * qJD(1) + t223 * t204 - t221 * t212;
t167 = 0.1e1 / t169;
t158 = 0.1e1 / t160;
t154 = t275 * t278;
t150 = t236 * t248 + (-0.2e1 * t259 * t261 - t173 * t205 + (t174 * t209 + t192 * t202 - t193 * t201) * t206) * t177;
t1 = [0.2e1 * t278 * t247 + (t172 * t205 + t260 * t278) * t177, t150, 0, 0, 0, 0; t193 * t161 * t251 + (-t174 * t161 + (t151 * t193 + t154 * t172) * t162) * t158 + ((t154 * t250 + t275 * t266) * t158 + (t154 * t251 + (-(-t155 * t177 * t262 + t248) * t265 - (-0.2e1 * t193 * t247 - t155 + (t155 - t238) * t177) * t264) * t158) * t162) * t278 (-t161 * t196 + t162 * t270) * t251 + (t250 * t270 - t171 * t161 + (-t196 * t151 + t152 * t172 - (-t150 * t193 - t156 * t174 - t201 + (-t156 * t210 - t192) * t155) * t264 - (-t150 * t210 - t156 * t202 - t173 + (t156 * t193 - t209) * t155) * t265) * t162) * t158, 0, 0, 0, 0; (-t182 * t187 + t188 * t263) * t249 + ((t188 * qJD(5) - t173 * t222 - t219 * t246) * t182 + t188 * t243 + (-t187 * t166 - (-t187 * qJD(5) + t173 * t219 - t222 * t246) * t185 - t188 * t165) * t183) * t167, t237 * t278 * t249 + (-t237 * (t231 * t253 - t233 * t254 - t277) + ((qJD(5) * t182 + t243) * t219 + (-t165 * t219 + (t166 - t252) * t222) * t183) * t278) * t167, 0, 0, -0.2e1 * t269 + 0.2e1 * (t165 * t183 * t167 + (-t167 * t267 - t183 * t269) * t185) * t185, 0;];
JaD_rot  = t1;
