% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPPR1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:58:14
% EndTime: 2019-02-26 19:58:14
% DurationCPUTime: 0.94s
% Computational Cost: add. (6115->111), mult. (9085->227), div. (559->12), fcn. (11668->13), ass. (0->106)
t233 = sin(pkin(10));
t235 = cos(pkin(10));
t237 = sin(qJ(2));
t236 = cos(pkin(6));
t238 = cos(qJ(2));
t264 = t236 * t238;
t219 = -t233 * t237 + t235 * t264;
t215 = t219 * qJD(2);
t265 = t236 * t237;
t220 = t233 * t238 + t235 * t265;
t232 = qJ(3) + pkin(11);
t228 = sin(t232);
t234 = sin(pkin(6));
t268 = t234 * t235;
t254 = t228 * t268;
t230 = cos(t232);
t261 = qJD(3) * t230;
t188 = -qJD(3) * t254 + t215 * t228 + t220 * t261;
t204 = t220 * t228 + t230 * t268;
t202 = t204 ^ 2;
t267 = t234 * t237;
t213 = t228 * t267 - t236 * t230;
t211 = 0.1e1 / t213 ^ 2;
t196 = t202 * t211 + 0.1e1;
t194 = 0.1e1 / t196;
t214 = t236 * t228 + t230 * t267;
t262 = qJD(2) * t238;
t253 = t234 * t262;
t200 = t214 * qJD(3) + t228 * t253;
t210 = 0.1e1 / t213;
t272 = t204 * t211;
t166 = (-t188 * t210 + t200 * t272) * t194;
t197 = atan2(-t204, t213);
t192 = sin(t197);
t193 = cos(t197);
t250 = -t192 * t213 - t193 * t204;
t162 = t250 * t166 - t188 * t192 + t193 * t200;
t176 = -t192 * t204 + t193 * t213;
t173 = 0.1e1 / t176;
t174 = 0.1e1 / t176 ^ 2;
t286 = t162 * t173 * t174;
t255 = t233 * t265;
t222 = t235 * t238 - t255;
t269 = t233 * t234;
t247 = -t222 * t228 + t230 * t269;
t285 = -0.2e1 * t247 * t286;
t266 = t234 * t238;
t246 = -t210 * t219 + t266 * t272;
t284 = t228 * t246;
t273 = t200 * t210 * t211;
t283 = -0.2e1 * (t188 * t272 - t202 * t273) / t196 ^ 2;
t208 = t222 * t230 + t228 * t269;
t221 = t233 * t264 + t235 * t237;
t231 = pkin(12) + qJ(6);
t227 = sin(t231);
t229 = cos(t231);
t187 = t208 * t229 + t221 * t227;
t183 = 0.1e1 / t187;
t184 = 0.1e1 / t187 ^ 2;
t217 = t221 * qJD(2);
t191 = t247 * qJD(3) - t217 * t230;
t218 = -qJD(2) * t255 + t235 * t262;
t177 = t187 * qJD(6) + t191 * t227 - t218 * t229;
t186 = t208 * t227 - t221 * t229;
t182 = t186 ^ 2;
t181 = t182 * t184 + 0.1e1;
t278 = t184 * t186;
t260 = qJD(6) * t186;
t178 = t191 * t229 + t218 * t227 - t260;
t280 = t178 * t183 * t184;
t282 = (t177 * t278 - t182 * t280) / t181 ^ 2;
t281 = t174 * t247;
t279 = t183 * t227;
t277 = t186 * t229;
t190 = t208 * qJD(3) - t217 * t228;
t276 = t190 * t174;
t275 = t192 * t247;
t274 = t193 * t247;
t271 = t221 * t228;
t270 = t221 * t230;
t263 = qJD(2) * t237;
t203 = t247 ^ 2;
t172 = t174 * t203 + 0.1e1;
t259 = 0.2e1 * (-t203 * t286 - t247 * t276) / t172 ^ 2;
t258 = -0.2e1 * t282;
t256 = t186 * t280;
t252 = -0.2e1 * t204 * t273;
t251 = qJD(6) * t270 - t217;
t249 = t184 * t277 - t279;
t206 = t220 * t230 - t254;
t248 = -t206 * t210 + t214 * t272;
t245 = qJD(3) * t271 + qJD(6) * t222 - t218 * t230;
t216 = t220 * qJD(2);
t201 = -t213 * qJD(3) + t230 * t253;
t199 = t222 * t227 - t229 * t270;
t198 = -t222 * t229 - t227 * t270;
t189 = -t204 * qJD(3) + t215 * t230;
t179 = 0.1e1 / t181;
t169 = 0.1e1 / t172;
t168 = t194 * t284;
t167 = t248 * t194;
t164 = (-t192 * t219 + t193 * t266) * t228 + t250 * t168;
t163 = t250 * t167 - t192 * t206 + t193 * t214;
t161 = t248 * t283 + (t214 * t252 - t189 * t210 + (t188 * t214 + t200 * t206 + t201 * t204) * t211) * t194;
t159 = t283 * t284 + (t246 * t261 + (t252 * t266 + t210 * t216 + (t200 * t219 + (t188 * t238 - t204 * t263) * t234) * t211) * t228) * t194;
t1 = [0, t159, t161, 0, 0, 0; 0 (-t164 * t281 + t173 * t271) * t259 + ((-t218 * t228 - t221 * t261) * t173 + (-t276 + t285) * t164 + (t271 * t162 + (-t159 * t204 - t168 * t188 + (-t228 * t263 + t238 * t261) * t234 + (-t168 * t213 - t219 * t228) * t166) * t274 + (-t219 * t261 - t159 * t213 - t168 * t200 + t216 * t228 + (t168 * t204 - t228 * t266) * t166) * t275) * t174) * t169 (-t163 * t281 - t173 * t208) * t259 + (t163 * t285 + t191 * t173 + (-t208 * t162 - t163 * t190 + (-t161 * t204 - t167 * t188 + t201 + (-t167 * t213 - t206) * t166) * t274 + (-t161 * t213 - t167 * t200 - t189 + (t167 * t204 - t214) * t166) * t275) * t174) * t169, 0, 0, 0; 0, 0.2e1 * (-t183 * t198 + t199 * t278) * t282 + (0.2e1 * t199 * t256 - t251 * t183 * t229 + t245 * t279 + (-t251 * t186 * t227 - t199 * t177 - t198 * t178 - t245 * t277) * t184) * t179, -t249 * t247 * t258 + (t249 * t190 - ((-qJD(6) * t183 - 0.2e1 * t256) * t229 + (t177 * t229 + (t178 - t260) * t227) * t184) * t247) * t179, 0, 0, t258 + 0.2e1 * (t177 * t184 * t179 + (-t179 * t280 - t184 * t282) * t186) * t186;];
JaD_rot  = t1;
