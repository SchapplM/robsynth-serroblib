% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRR8_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_jacobiaD_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:17
% EndTime: 2019-02-26 20:08:18
% DurationCPUTime: 0.95s
% Computational Cost: add. (4346->110), mult. (14031->217), div. (500->12), fcn. (17771->13), ass. (0->108)
t216 = sin(pkin(7));
t218 = cos(pkin(7));
t220 = sin(qJ(2));
t222 = cos(qJ(2));
t279 = sin(pkin(12));
t281 = cos(pkin(6));
t249 = t281 * t279;
t280 = cos(pkin(12));
t235 = t280 * t220 + t222 * t249;
t217 = sin(pkin(6));
t258 = t217 * t279;
t199 = t216 * t235 + t218 * t258;
t195 = 0.1e1 / t199 ^ 2;
t236 = t220 * t249 - t280 * t222;
t207 = t236 * qJD(2);
t285 = t207 * t195 * t216;
t250 = t281 * t280;
t238 = -t279 * t220 + t222 * t250;
t205 = t238 * qJD(2);
t253 = t280 * t217 * t216;
t284 = -qJD(3) * t253 + t205;
t219 = sin(qJ(3));
t237 = -t220 * t250 - t279 * t222;
t232 = t237 * qJD(2);
t221 = cos(qJ(3));
t233 = t238 * t221;
t283 = qJD(3) * t233 + t219 * t232;
t234 = t238 * t219;
t229 = t218 * t234 - t221 * t237;
t267 = t218 * t221;
t163 = t229 * qJD(3) + t219 * t284 - t232 * t267;
t270 = t237 * t219;
t182 = -t218 * t233 + t221 * t253 - t270;
t179 = t182 ^ 2;
t263 = t221 * t222;
t266 = t219 * t220;
t243 = t218 * t263 - t266;
t257 = t281 * t216;
t197 = -t243 * t217 - t221 * t257;
t192 = 0.1e1 / t197 ^ 2;
t174 = t179 * t192 + 0.1e1;
t273 = t182 * t192;
t264 = t220 * t221;
t265 = t219 * t222;
t241 = t218 * t265 + t264;
t242 = t218 * t264 + t265;
t251 = qJD(3) * t257;
t177 = t219 * t251 + (t242 * qJD(2) + t241 * qJD(3)) * t217;
t191 = 0.1e1 / t197;
t274 = t177 * t191 * t192;
t282 = -0.2e1 * (t163 * t273 - t179 * t274) / t174 ^ 2;
t175 = atan2(-t182, t197);
t168 = sin(t175);
t169 = cos(t175);
t162 = -t168 * t182 + t169 * t197;
t159 = 0.1e1 / t162;
t194 = 0.1e1 / t199;
t160 = 0.1e1 / t162 ^ 2;
t170 = 0.1e1 / t174;
t151 = (-t163 * t191 + t177 * t273) * t170;
t248 = -t168 * t197 - t169 * t182;
t148 = t248 * t151 - t168 * t163 + t169 * t177;
t278 = t148 * t159 * t160;
t252 = t216 * t258;
t269 = t236 * t219;
t185 = -t221 * t252 + t235 * t267 - t269;
t277 = t160 * t185;
t276 = t168 * t185;
t275 = t169 * t185;
t239 = -t218 * t235 + t252;
t186 = t239 * t219 - t221 * t236;
t272 = t186 * t195;
t268 = t218 * t219;
t180 = t185 ^ 2;
t157 = t180 * t160 + 0.1e1;
t206 = t235 * qJD(2);
t165 = qJD(3) * t186 - t206 * t219 - t207 * t267;
t262 = 0.2e1 * (t165 * t277 - t180 * t278) / t157 ^ 2;
t166 = t207 * t268 - t206 * t221 + (t239 * t221 + t269) * qJD(3);
t181 = t186 ^ 2;
t176 = t181 * t195 + 0.1e1;
t260 = t194 * t285;
t261 = 0.2e1 * (t166 * t272 + t181 * t260) / t176 ^ 2;
t259 = qJD(3) * t270;
t255 = -0.2e1 * t182 * t274;
t254 = 0.2e1 * t185 * t278;
t184 = -t219 * t253 + t229;
t198 = t241 * t217 + t219 * t257;
t245 = -t184 * t191 + t198 * t273;
t188 = -t237 * t267 + t234;
t204 = t242 * t217;
t244 = -t188 * t191 + t204 * t273;
t189 = -t219 * t235 - t236 * t267;
t190 = -t221 * t235 + t236 * t268;
t240 = -t218 * t266 + t263;
t187 = (t243 * qJD(2) + t240 * qJD(3)) * t217;
t178 = t221 * t251 + (t240 * qJD(2) + t243 * qJD(3)) * t217;
t172 = 0.1e1 / t176;
t167 = t205 * t267 + t218 * t259 + t283;
t164 = t218 * t283 + t221 * t284 + t259;
t155 = 0.1e1 / t157;
t153 = t244 * t170;
t152 = t245 * t170;
t150 = t248 * t153 - t168 * t188 + t169 * t204;
t149 = t248 * t152 - t168 * t184 + t169 * t198;
t147 = t244 * t282 + (t204 * t255 - t167 * t191 + (t163 * t204 + t177 * t188 + t182 * t187) * t192) * t170;
t146 = t245 * t282 + (t198 * t255 - t164 * t191 + (t163 * t198 + t177 * t184 + t178 * t182) * t192) * t170;
t1 = [0, t147, t146, 0, 0, 0; 0 (t150 * t277 - t159 * t189) * t262 + ((t190 * qJD(3) - t206 * t267 + t207 * t219) * t159 + t150 * t254 + (-t189 * t148 - t150 * t165 - (-t147 * t182 - t153 * t163 + t187 + (-t153 * t197 - t188) * t151) * t275 - (-t147 * t197 - t153 * t177 - t167 + (t153 * t182 - t204) * t151) * t276) * t160) * t155 (t149 * t277 - t159 * t186) * t262 + (t149 * t254 + t166 * t159 + (-t186 * t148 - t149 * t165 - (-t146 * t182 - t152 * t163 + t178 + (-t152 * t197 - t184) * t151) * t275 - (-t146 * t197 - t152 * t177 - t164 + (t152 * t182 - t198) * t151) * t276) * t160) * t155, 0, 0, 0; 0 (-t216 * t236 * t272 - t190 * t194) * t261 + ((-t189 * qJD(3) + t206 * t268 + t207 * t221) * t194 + (0.2e1 * t236 * t186 * t260 + (t166 * t236 + t186 * t206 + t190 * t207) * t195) * t216) * t172, t185 * t194 * t261 + (-t165 * t194 - t185 * t285) * t172, 0, 0, 0;];
JaD_rot  = t1;
