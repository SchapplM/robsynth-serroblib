% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPPR2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:58:46
% EndTime: 2019-02-26 19:58:47
% DurationCPUTime: 1.06s
% Computational Cost: add. (5804->108), mult. (9085->220), div. (559->12), fcn. (11668->13), ass. (0->105)
t228 = cos(pkin(6));
t232 = cos(qJ(2));
t282 = cos(pkin(10));
t253 = t282 * t232;
t226 = sin(pkin(10));
t230 = sin(qJ(2));
t266 = t226 * t230;
t216 = t228 * t253 - t266;
t212 = t216 * qJD(2);
t225 = qJ(3) + pkin(11);
t224 = cos(t225);
t223 = sin(t225);
t254 = t282 * t230;
t265 = t226 * t232;
t241 = -t228 * t254 - t265;
t227 = sin(pkin(6));
t255 = t227 * t282;
t284 = t241 * t223 - t224 * t255;
t179 = qJD(3) * t284 + t212 * t224;
t202 = -t223 * t255 - t241 * t224;
t199 = t202 ^ 2;
t264 = t227 * t230;
t211 = t223 * t228 + t224 * t264;
t208 = 0.1e1 / t211 ^ 2;
t192 = t199 * t208 + 0.1e1;
t190 = 0.1e1 / t192;
t210 = -t223 * t264 + t224 * t228;
t263 = t227 * t232;
t256 = qJD(2) * t263;
t198 = t210 * qJD(3) + t224 * t256;
t207 = 0.1e1 / t211;
t271 = t202 * t208;
t162 = (-t179 * t207 + t198 * t271) * t190;
t193 = atan2(-t202, t211);
t186 = sin(t193);
t187 = cos(t193);
t248 = -t186 * t211 - t187 * t202;
t158 = t248 * t162 - t186 * t179 + t187 * t198;
t172 = -t186 * t202 + t187 * t211;
t169 = 0.1e1 / t172;
t170 = 0.1e1 / t172 ^ 2;
t288 = t158 * t169 * t170;
t231 = cos(qJ(6));
t218 = -t228 * t266 + t253;
t267 = t226 * t227;
t244 = -t218 * t223 + t224 * t267;
t229 = sin(qJ(6));
t242 = -t228 * t265 - t254;
t269 = t242 * t229;
t247 = -t231 * t244 + t269;
t287 = t247 * qJD(6);
t205 = t218 * t224 + t223 * t267;
t286 = 0.2e1 * t205 * t288;
t243 = -t207 * t216 + t263 * t271;
t285 = t224 * t243;
t272 = t198 * t207 * t208;
t283 = -0.2e1 * (t179 * t271 - t199 * t272) / t192 ^ 2;
t268 = t242 * t231;
t189 = -t229 * t244 - t268;
t183 = 0.1e1 / t189;
t184 = 0.1e1 / t189 ^ 2;
t214 = t242 * qJD(2);
t180 = t205 * qJD(3) + t214 * t223;
t215 = t218 * qJD(2);
t173 = t189 * qJD(6) - t180 * t231 + t215 * t229;
t182 = t247 ^ 2;
t177 = t182 * t184 + 0.1e1;
t276 = t184 * t247;
t174 = t180 * t229 + t215 * t231 + t287;
t279 = t174 * t183 * t184;
t281 = (-t173 * t276 - t182 * t279) / t177 ^ 2;
t280 = t170 * t205;
t181 = t244 * qJD(3) + t214 * t224;
t278 = t181 * t170;
t277 = t183 * t231;
t275 = t186 * t205;
t274 = t187 * t205;
t273 = t247 * t229;
t270 = t242 * t224;
t262 = qJD(2) * t230;
t261 = qJD(3) * t223;
t200 = t205 ^ 2;
t168 = t170 * t200 + 0.1e1;
t260 = 0.2e1 * (-t200 * t288 + t205 * t278) / t168 ^ 2;
t259 = 0.2e1 * t281;
t252 = -0.2e1 * t247 * t279;
t251 = -0.2e1 * t202 * t272;
t249 = qJD(6) * t223 * t242 + t214;
t246 = -t184 * t273 + t277;
t245 = -t207 * t284 + t210 * t271;
t240 = -qJD(3) * t270 + qJD(6) * t218 + t215 * t223;
t213 = t241 * qJD(2);
t197 = -t211 * qJD(3) - t223 * t256;
t195 = t218 * t231 + t223 * t269;
t194 = t218 * t229 - t223 * t268;
t178 = t202 * qJD(3) + t212 * t223;
t175 = 0.1e1 / t177;
t165 = 0.1e1 / t168;
t164 = t190 * t285;
t163 = t245 * t190;
t160 = (-t186 * t216 + t187 * t263) * t224 + t248 * t164;
t159 = t248 * t163 - t186 * t284 + t187 * t210;
t157 = t245 * t283 + (t210 * t251 + t178 * t207 + (t179 * t210 + t197 * t202 + t198 * t284) * t208) * t190;
t155 = t283 * t285 + (-t243 * t261 + (t251 * t263 - t207 * t213 + (t198 * t216 + (t179 * t232 - t202 * t262) * t227) * t208) * t224) * t190;
t1 = [0, t155, t157, 0, 0, 0; 0 (t160 * t280 - t169 * t270) * t260 + ((-t215 * t224 - t242 * t261) * t169 + (-t278 + t286) * t160 + (-t270 * t158 - (-t155 * t202 - t164 * t179 + (-t224 * t262 - t232 * t261) * t227 + (-t164 * t211 - t216 * t224) * t162) * t274 - (t216 * t261 - t155 * t211 - t164 * t198 - t213 * t224 + (t164 * t202 - t224 * t263) * t162) * t275) * t170) * t165 (t159 * t280 - t169 * t244) * t260 + (t159 * t286 - t180 * t169 + (-t244 * t158 - t159 * t181 - (-t157 * t202 - t163 * t179 + t197 + (-t163 * t211 - t284) * t162) * t274 - (-t157 * t211 - t163 * t198 + t178 + (t163 * t202 - t210) * t162) * t275) * t170) * t165, 0, 0, 0; 0 (-t183 * t194 - t195 * t276) * t259 + (t195 * t252 + t249 * t183 * t229 + t240 * t277 + (t231 * t247 * t249 - t195 * t173 - t194 * t174 - t240 * t273) * t184) * t175, t246 * t205 * t259 + (-t246 * t181 + ((qJD(6) * t183 + t252) * t229 + (-t173 * t229 + (t174 + t287) * t231) * t184) * t205) * t175, 0, 0, -0.2e1 * t281 - 0.2e1 * (t173 * t175 * t184 - (-t175 * t279 - t184 * t281) * t247) * t247;];
JaD_rot  = t1;
