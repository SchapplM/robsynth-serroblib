% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR3_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:29:34
% EndTime: 2019-02-26 21:29:35
% DurationCPUTime: 0.90s
% Computational Cost: add. (3615->96), mult. (9869->202), div. (448->12), fcn. (12731->13), ass. (0->97)
t226 = pkin(12) + qJ(5);
t224 = sin(t226);
t225 = cos(t226);
t281 = sin(pkin(11));
t283 = cos(pkin(6));
t251 = t283 * t281;
t282 = cos(pkin(11));
t252 = t283 * t282;
t284 = sin(qJ(2));
t285 = cos(qJ(2));
t215 = t251 * t285 + t252 * t284;
t217 = t281 * t284 - t282 * t285;
t228 = sin(qJ(1));
t229 = cos(qJ(1));
t248 = t228 * t215 + t229 * t217;
t227 = sin(pkin(6));
t266 = t227 * t228;
t243 = t224 * t248 + t225 * t266;
t288 = qJD(5) * t243;
t240 = t285 * t281 + t284 * t282;
t214 = t240 * t227;
t205 = qJD(2) * t214;
t213 = t217 * t227;
t210 = 0.1e1 / t213 ^ 2;
t268 = t205 * t210;
t287 = t284 * t251 - t285 * t252;
t239 = t217 * qJD(2);
t197 = -t228 * t240 - t229 * t287;
t185 = atan2(t197, t213);
t180 = sin(t185);
t181 = cos(t185);
t194 = t197 ^ 2;
t184 = t194 * t210 + 0.1e1;
t182 = 0.1e1 / t184;
t209 = 0.1e1 / t213;
t270 = t197 * t209;
t286 = t182 * (t181 * t270 - t180) + t180;
t169 = t180 * t197 + t181 * t213;
t166 = 0.1e1 / t169;
t193 = t224 * t266 - t225 * t248;
t187 = 0.1e1 / t193;
t167 = 0.1e1 / t169 ^ 2;
t188 = 0.1e1 / t193 ^ 2;
t236 = t228 * t287;
t200 = -t229 * t240 + t236;
t195 = t200 ^ 2;
t165 = t195 * t167 + 0.1e1;
t208 = t215 * qJD(2);
t175 = qJD(1) * t197 - t228 * t208 - t229 * t239;
t274 = t175 * t167;
t264 = qJD(1) * t229;
t178 = qJD(1) * t236 - t229 * t208 + t228 * t239 - t240 * t264;
t247 = t178 * t209 - t197 * t268;
t160 = t247 * t182;
t250 = -t180 * t213 + t181 * t197;
t156 = t160 * t250 + t180 * t178 + t181 * t205;
t279 = t156 * t166 * t167;
t280 = (-t195 * t279 - t200 * t274) / t165 ^ 2;
t249 = -t229 * t215 + t228 * t217;
t269 = t197 * t214;
t245 = -t209 * t249 + t210 * t269;
t161 = t245 * t182;
t157 = -t161 * t250 + t180 * t249 + t181 * t214;
t278 = t157 * t200;
t207 = t287 * qJD(2);
t216 = t240 * qJD(2);
t176 = qJD(1) * t249 + t228 * t207 - t229 * t216;
t258 = t227 * t264;
t170 = qJD(5) * t193 + t176 * t224 - t225 * t258;
t186 = t243 ^ 2;
t174 = t186 * t188 + 0.1e1;
t271 = t188 * t243;
t171 = t176 * t225 + t224 * t258 + t288;
t275 = t171 * t187 * t188;
t277 = (-t170 * t271 - t186 * t275) / t174 ^ 2;
t267 = t209 * t268;
t276 = (t197 * t210 * t178 - t194 * t267) / t184 ^ 2;
t273 = t180 * t200;
t272 = t181 * t200;
t265 = t227 * t229;
t263 = -0.2e1 * t280;
t262 = -0.2e1 * t279;
t261 = 0.2e1 * t277;
t260 = 0.2e1 * t276;
t259 = qJD(1) * t266;
t257 = -0.2e1 * t209 * t276;
t256 = -0.2e1 * t243 * t275;
t246 = -t224 * t187 - t225 * t271;
t244 = -t224 * t249 + t225 * t265;
t191 = t224 * t265 + t225 * t249;
t238 = qJD(1) * t248 + t229 * t207 + t228 * t216;
t206 = t227 * t239;
t172 = 0.1e1 / t174;
t163 = 0.1e1 / t165;
t159 = t286 * t200;
t155 = t245 * t260 + (0.2e1 * t267 * t269 + t238 * t209 + (-t178 * t214 + t197 * t206 - t205 * t249) * t210) * t182;
t1 = [t200 * t257 + (-t175 * t209 - t200 * t268) * t182, t155, 0, 0, 0, 0; t197 * t166 * t263 + (t178 * t166 + (-t156 * t197 - t159 * t175) * t167) * t163 + ((t159 * t262 - t286 * t274) * t163 + (t159 * t263 + ((-t160 * t182 * t270 + t260) * t273 + (t197 * t257 + t160 + (-t160 + t247) * t182) * t272) * t163) * t167) * t200, 0.2e1 * (t166 * t248 - t167 * t278) * t280 + (t176 * t166 + t262 * t278 + (t248 * t156 - t157 * t175 + (t155 * t197 - t161 * t178 - t206 + (t161 * t213 + t249) * t160) * t272 + (-t155 * t213 + t161 * t205 + t238 + (t161 * t197 - t214) * t160) * t273) * t167) * t163, 0, 0, 0, 0; (t187 * t244 - t191 * t271) * t261 + ((qJD(5) * t191 + t224 * t238 + t225 * t259) * t187 + t191 * t256 + (t244 * t171 + (qJD(5) * t244 - t224 * t259 + t225 * t238) * t243 - t191 * t170) * t188) * t172, t246 * t200 * t261 + (t246 * t175 + ((qJD(5) * t187 + t256) * t225 + (-t170 * t225 + (-t171 - t288) * t224) * t188) * t200) * t172, 0, 0, -0.2e1 * t277 - 0.2e1 * (t170 * t188 * t172 - (-t172 * t275 - t188 * t277) * t243) * t243, 0;];
JaD_rot  = t1;
