% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRPR4_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR4_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:12:19
% EndTime: 2019-02-26 20:12:20
% DurationCPUTime: 0.97s
% Computational Cost: add. (3313->110), mult. (9085->226), div. (559->12), fcn. (11668->13), ass. (0->104)
t219 = sin(pkin(11));
t221 = cos(pkin(11));
t224 = sin(qJ(2));
t222 = cos(pkin(6));
t226 = cos(qJ(2));
t252 = t222 * t226;
t206 = -t219 * t224 + t221 * t252;
t199 = t206 * qJD(2);
t253 = t222 * t224;
t207 = t219 * t226 + t221 * t253;
t223 = sin(qJ(3));
t220 = sin(pkin(6));
t256 = t220 * t223;
t242 = t221 * t256;
t225 = cos(qJ(3));
t249 = qJD(3) * t225;
t177 = -qJD(3) * t242 + t199 * t223 + t207 * t249;
t255 = t220 * t225;
t191 = t207 * t223 + t221 * t255;
t189 = t191 ^ 2;
t210 = -t222 * t225 + t224 * t256;
t204 = 0.1e1 / t210 ^ 2;
t185 = t189 * t204 + 0.1e1;
t183 = 0.1e1 / t185;
t211 = t222 * t223 + t224 * t255;
t250 = qJD(2) * t226;
t241 = t220 * t250;
t196 = t211 * qJD(3) + t223 * t241;
t203 = 0.1e1 / t210;
t260 = t191 * t204;
t155 = (-t177 * t203 + t196 * t260) * t183;
t186 = atan2(-t191, t210);
t181 = sin(t186);
t182 = cos(t186);
t238 = -t181 * t210 - t182 * t191;
t151 = t238 * t155 - t181 * t177 + t182 * t196;
t167 = -t181 * t191 + t182 * t210;
t164 = 0.1e1 / t167;
t165 = 0.1e1 / t167 ^ 2;
t273 = t151 * t164 * t165;
t243 = t219 * t253;
t209 = t221 * t226 - t243;
t235 = -t209 * t223 + t219 * t255;
t272 = -0.2e1 * t235 * t273;
t254 = t220 * t226;
t234 = -t203 * t206 + t254 * t260;
t271 = t223 * t234;
t259 = t196 * t203 * t204;
t270 = -0.2e1 * (t177 * t260 - t189 * t259) / t185 ^ 2;
t195 = t209 * t225 + t219 * t256;
t208 = t219 * t252 + t221 * t224;
t218 = qJ(4) + pkin(12);
t216 = sin(t218);
t217 = cos(t218);
t176 = t195 * t217 + t208 * t216;
t172 = 0.1e1 / t176;
t173 = 0.1e1 / t176 ^ 2;
t201 = t208 * qJD(2);
t180 = t235 * qJD(3) - t201 * t225;
t202 = -qJD(2) * t243 + t221 * t250;
t162 = t176 * qJD(4) + t180 * t216 - t202 * t217;
t175 = t195 * t216 - t208 * t217;
t171 = t175 ^ 2;
t170 = t171 * t173 + 0.1e1;
t265 = t173 * t175;
t248 = qJD(4) * t175;
t163 = t180 * t217 + t202 * t216 - t248;
t268 = t163 * t172 * t173;
t269 = (t162 * t265 - t171 * t268) / t170 ^ 2;
t267 = t165 * t235;
t266 = t172 * t216;
t264 = t175 * t217;
t179 = t195 * qJD(3) - t201 * t223;
t263 = t179 * t165;
t262 = t181 * t235;
t261 = t182 * t235;
t258 = t208 * t223;
t257 = t208 * t225;
t251 = qJD(2) * t224;
t190 = t235 ^ 2;
t161 = t165 * t190 + 0.1e1;
t247 = 0.2e1 * (-t190 * t273 - t235 * t263) / t161 ^ 2;
t246 = -0.2e1 * t269;
t244 = t175 * t268;
t240 = -0.2e1 * t191 * t259;
t239 = qJD(4) * t257 - t201;
t237 = t173 * t264 - t266;
t193 = t207 * t225 - t242;
t236 = -t193 * t203 + t211 * t260;
t233 = qJD(3) * t258 + qJD(4) * t209 - t202 * t225;
t200 = t207 * qJD(2);
t197 = -t210 * qJD(3) + t225 * t241;
t188 = t209 * t216 - t217 * t257;
t187 = -t209 * t217 - t216 * t257;
t178 = -t191 * qJD(3) + t199 * t225;
t168 = 0.1e1 / t170;
t158 = 0.1e1 / t161;
t157 = t183 * t271;
t156 = t236 * t183;
t153 = (-t181 * t206 + t182 * t254) * t223 + t238 * t157;
t152 = t238 * t156 - t181 * t193 + t182 * t211;
t150 = t236 * t270 + (t211 * t240 - t178 * t203 + (t177 * t211 + t191 * t197 + t193 * t196) * t204) * t183;
t148 = t270 * t271 + (t234 * t249 + (t240 * t254 + t200 * t203 + (t196 * t206 + (t177 * t226 - t191 * t251) * t220) * t204) * t223) * t183;
t1 = [0, t148, t150, 0, 0, 0; 0 (-t153 * t267 + t164 * t258) * t247 + ((-t202 * t223 - t208 * t249) * t164 + (-t263 + t272) * t153 + (t258 * t151 + (-t148 * t191 - t157 * t177 + (-t223 * t251 + t226 * t249) * t220 + (-t157 * t210 - t206 * t223) * t155) * t261 + (-t206 * t249 - t148 * t210 - t157 * t196 + t200 * t223 + (t157 * t191 - t223 * t254) * t155) * t262) * t165) * t158 (-t152 * t267 - t164 * t195) * t247 + (t152 * t272 + t180 * t164 + (-t195 * t151 - t152 * t179 + (-t150 * t191 - t156 * t177 + t197 + (-t156 * t210 - t193) * t155) * t261 + (-t150 * t210 - t156 * t196 - t178 + (t156 * t191 - t211) * t155) * t262) * t165) * t158, 0, 0, 0; 0, 0.2e1 * (-t172 * t187 + t188 * t265) * t269 + (0.2e1 * t188 * t244 - t239 * t172 * t217 + t233 * t266 + (-t239 * t175 * t216 - t188 * t162 - t187 * t163 - t233 * t264) * t173) * t168, -t237 * t235 * t246 + (t237 * t179 - ((-qJD(4) * t172 - 0.2e1 * t244) * t217 + (t162 * t217 + (t163 - t248) * t216) * t173) * t235) * t168, t246 + 0.2e1 * (t162 * t168 * t173 + (-t168 * t268 - t173 * t269) * t175) * t175, 0, 0;];
JaD_rot  = t1;
