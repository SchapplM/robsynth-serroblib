% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRP1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:01:11
% EndTime: 2019-02-26 20:01:12
% DurationCPUTime: 0.92s
% Computational Cost: add. (5804->110), mult. (9085->227), div. (559->12), fcn. (11668->13), ass. (0->106)
t231 = sin(pkin(10));
t233 = cos(pkin(10));
t236 = sin(qJ(2));
t234 = cos(pkin(6));
t238 = cos(qJ(2));
t264 = t234 * t238;
t220 = -t231 * t236 + t233 * t264;
t216 = t220 * qJD(2);
t265 = t234 * t236;
t221 = t231 * t238 + t233 * t265;
t230 = qJ(3) + pkin(11);
t228 = sin(t230);
t232 = sin(pkin(6));
t268 = t232 * t233;
t254 = t228 * t268;
t229 = cos(t230);
t261 = qJD(3) * t229;
t183 = -qJD(3) * t254 + t216 * t228 + t221 * t261;
t205 = t221 * t228 + t229 * t268;
t203 = t205 ^ 2;
t267 = t232 * t236;
t213 = t228 * t267 - t229 * t234;
t211 = 0.1e1 / t213 ^ 2;
t197 = t203 * t211 + 0.1e1;
t195 = 0.1e1 / t197;
t214 = t228 * t234 + t229 * t267;
t262 = qJD(2) * t238;
t253 = t232 * t262;
t201 = t214 * qJD(3) + t228 * t253;
t210 = 0.1e1 / t213;
t273 = t205 * t211;
t167 = (-t183 * t210 + t201 * t273) * t195;
t198 = atan2(-t205, t213);
t191 = sin(t198);
t192 = cos(t198);
t250 = -t191 * t213 - t192 * t205;
t163 = t250 * t167 - t183 * t191 + t192 * t201;
t177 = -t191 * t205 + t192 * t213;
t174 = 0.1e1 / t177;
t175 = 0.1e1 / t177 ^ 2;
t287 = t163 * t174 * t175;
t255 = t231 * t265;
t223 = t233 * t238 - t255;
t269 = t231 * t232;
t247 = -t223 * t228 + t229 * t269;
t286 = -0.2e1 * t247 * t287;
t266 = t232 * t238;
t246 = -t210 * t220 + t266 * t273;
t285 = t228 * t246;
t274 = t201 * t210 * t211;
t284 = -0.2e1 * (t183 * t273 - t203 * t274) / t197 ^ 2;
t209 = t223 * t229 + t228 * t269;
t237 = cos(qJ(5));
t222 = t231 * t264 + t233 * t236;
t235 = sin(qJ(5));
t271 = t222 * t235;
t194 = t209 * t237 + t271;
t188 = 0.1e1 / t194;
t189 = 0.1e1 / t194 ^ 2;
t218 = t222 * qJD(2);
t186 = t247 * qJD(3) - t218 * t229;
t219 = -qJD(2) * t255 + t233 * t262;
t178 = t194 * qJD(5) + t186 * t235 - t219 * t237;
t270 = t222 * t237;
t193 = t209 * t235 - t270;
t187 = t193 ^ 2;
t182 = t187 * t189 + 0.1e1;
t278 = t189 * t193;
t260 = qJD(5) * t193;
t179 = t186 * t237 + t219 * t235 - t260;
t281 = t179 * t188 * t189;
t283 = (t178 * t278 - t187 * t281) / t182 ^ 2;
t282 = t175 * t247;
t185 = t209 * qJD(3) - t218 * t228;
t280 = t185 * t175;
t279 = t188 * t235;
t277 = t191 * t247;
t276 = t192 * t247;
t275 = t193 * t237;
t272 = t222 * t228;
t263 = qJD(2) * t236;
t204 = t247 ^ 2;
t173 = t175 * t204 + 0.1e1;
t259 = 0.2e1 * (-t204 * t287 - t247 * t280) / t173 ^ 2;
t258 = -0.2e1 * t283;
t256 = t193 * t281;
t252 = -0.2e1 * t205 * t274;
t251 = qJD(5) * t222 * t229 - t218;
t249 = t189 * t275 - t279;
t207 = t221 * t229 - t254;
t248 = -t207 * t210 + t214 * t273;
t245 = qJD(3) * t272 + qJD(5) * t223 - t219 * t229;
t217 = t221 * qJD(2);
t202 = -t213 * qJD(3) + t229 * t253;
t200 = t223 * t235 - t229 * t270;
t199 = -t223 * t237 - t229 * t271;
t184 = -t205 * qJD(3) + t216 * t229;
t180 = 0.1e1 / t182;
t170 = 0.1e1 / t173;
t169 = t195 * t285;
t168 = t248 * t195;
t165 = (-t191 * t220 + t192 * t266) * t228 + t250 * t169;
t164 = t250 * t168 - t191 * t207 + t192 * t214;
t162 = t248 * t284 + (t214 * t252 - t184 * t210 + (t183 * t214 + t201 * t207 + t202 * t205) * t211) * t195;
t160 = t284 * t285 + (t246 * t261 + (t252 * t266 + t210 * t217 + (t201 * t220 + (t183 * t238 - t205 * t263) * t232) * t211) * t228) * t195;
t1 = [0, t160, t162, 0, 0, 0; 0 (-t165 * t282 + t174 * t272) * t259 + ((-t219 * t228 - t222 * t261) * t174 + (-t280 + t286) * t165 + (t272 * t163 + (-t160 * t205 - t169 * t183 + (-t228 * t263 + t238 * t261) * t232 + (-t169 * t213 - t220 * t228) * t167) * t276 + (-t220 * t261 - t160 * t213 - t169 * t201 + t217 * t228 + (t169 * t205 - t228 * t266) * t167) * t277) * t175) * t170 (-t164 * t282 - t174 * t209) * t259 + (t164 * t286 + t186 * t174 + (-t209 * t163 - t164 * t185 + (-t162 * t205 - t168 * t183 + t202 + (-t168 * t213 - t207) * t167) * t276 + (-t162 * t213 - t168 * t201 - t184 + (t168 * t205 - t214) * t167) * t277) * t175) * t170, 0, 0, 0; 0, 0.2e1 * (-t188 * t199 + t200 * t278) * t283 + (0.2e1 * t200 * t256 - t251 * t188 * t237 + t245 * t279 + (-t251 * t193 * t235 - t200 * t178 - t199 * t179 - t245 * t275) * t189) * t180, -t249 * t247 * t258 + (t249 * t185 - ((-qJD(5) * t188 - 0.2e1 * t256) * t237 + (t178 * t237 + (t179 - t260) * t235) * t189) * t247) * t180, 0, t258 + 0.2e1 * (t178 * t189 * t180 + (-t180 * t281 - t189 * t283) * t193) * t193, 0;];
JaD_rot  = t1;
