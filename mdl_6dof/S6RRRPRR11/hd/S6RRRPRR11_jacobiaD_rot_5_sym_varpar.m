% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR11_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR11_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:21:50
% EndTime: 2019-02-26 22:21:51
% DurationCPUTime: 1.03s
% Computational Cost: add. (2350->112), mult. (7168->247), div. (687->14), fcn. (9017->13), ass. (0->111)
t301 = qJD(3) - qJD(5);
t232 = sin(qJ(2));
t233 = sin(qJ(1));
t236 = cos(qJ(2));
t237 = cos(qJ(1));
t295 = cos(pkin(6));
t260 = t237 * t295;
t212 = t232 * t260 + t233 * t236;
t261 = t233 * t295;
t213 = t237 * t232 + t236 * t261;
t196 = t213 * qJD(1) + t212 * qJD(2);
t256 = t236 * t260;
t279 = t233 * t232;
t211 = -t256 + t279;
t209 = t211 ^ 2;
t229 = sin(pkin(6));
t225 = 0.1e1 / t229 ^ 2;
t227 = 0.1e1 / t236 ^ 2;
t208 = t209 * t225 * t227 + 0.1e1;
t226 = 0.1e1 / t236;
t228 = t226 * t227;
t276 = qJD(2) * t232;
t289 = (t196 * t211 * t227 + t209 * t228 * t276) * t225 / t208 ^ 2;
t300 = -0.2e1 * t289;
t224 = 0.1e1 / t229;
t299 = t211 * t224;
t235 = cos(qJ(3));
t298 = t301 * t235;
t231 = sin(qJ(3));
t297 = t301 * t231;
t281 = t229 * t236;
t207 = atan2(t211, t281);
t203 = sin(t207);
t204 = cos(t207);
t205 = 0.1e1 / t208;
t265 = t226 * t299;
t296 = (t204 * t265 - t203) * t205 + t203;
t195 = t212 * qJD(1) + t213 * qJD(2);
t257 = t232 * t261;
t278 = t237 * t236;
t215 = -t257 + t278;
t282 = t229 * t233;
t202 = t215 * t235 + t231 * t282;
t277 = qJD(1) * t229;
t263 = t237 * t277;
t174 = t202 * qJD(3) - t195 * t231 - t235 * t263;
t245 = -t215 * t231 + t235 * t282;
t175 = t245 * qJD(3) - t195 * t235 + t231 * t263;
t230 = sin(qJ(5));
t234 = cos(qJ(5));
t253 = -t202 * t230 - t234 * t245;
t163 = t253 * qJD(5) + t174 * t230 + t175 * t234;
t181 = t203 * t211 + t204 * t281;
t178 = 0.1e1 / t181;
t191 = t202 * t234 - t230 * t245;
t183 = 0.1e1 / t191;
t179 = 0.1e1 / t181 ^ 2;
t184 = 0.1e1 / t191 ^ 2;
t162 = t191 * qJD(5) - t174 * t234 + t175 * t230;
t182 = t253 ^ 2;
t167 = t182 * t184 + 0.1e1;
t288 = t184 * t253;
t185 = t183 * t184;
t290 = t163 * t185;
t294 = (-t162 * t288 - t182 * t290) / t167 ^ 2;
t210 = t213 ^ 2;
t172 = t210 * t179 + 0.1e1;
t250 = qJD(2) * t295 + qJD(1);
t275 = qJD(2) * t236;
t194 = -qJD(1) * t256 - t237 * t275 + t250 * t279;
t286 = t194 * t179;
t262 = t227 * t276;
t244 = (t196 * t226 + t211 * t262) * t224;
t168 = t205 * t244;
t248 = -t203 * t281 + t204 * t211;
t266 = t204 * t229 * t232;
t160 = -qJD(2) * t266 + t248 * t168 + t203 * t196;
t292 = t160 * t178 * t179;
t293 = (-t210 * t292 - t213 * t286) / t172 ^ 2;
t283 = t227 * t232;
t247 = t211 * t283 + t212 * t226;
t169 = t247 * t224 * t205;
t161 = t248 * t169 + t203 * t212 - t266;
t291 = t161 * t213;
t251 = t230 * t231 + t234 * t235;
t193 = t251 * t213;
t287 = t184 * t193;
t285 = t203 * t213;
t284 = t204 * t213;
t280 = t229 * t237;
t270 = 0.2e1 * t294;
t269 = -0.2e1 * t293;
t268 = -0.2e1 * t292;
t267 = -0.2e1 * t185 * t253;
t264 = t233 * t277;
t259 = t226 * t300;
t258 = t163 * t267;
t200 = -t212 * t235 + t231 * t280;
t246 = t212 * t231 + t235 * t280;
t254 = -t200 * t230 - t234 * t246;
t187 = t200 * t234 - t230 * t246;
t252 = t230 * t235 - t231 * t234;
t197 = -qJD(1) * t257 - t233 * t276 + t250 * t278;
t192 = t252 * t213;
t177 = t246 * qJD(3) - t197 * t235 - t231 * t264;
t176 = t200 * qJD(3) - t197 * t231 + t235 * t264;
t170 = 0.1e1 / t172;
t165 = 0.1e1 / t167;
t164 = t296 * t213;
t159 = (t247 * t300 + (t196 * t283 + t197 * t226 + (t212 * t283 + (0.2e1 * t228 * t232 ^ 2 + t226) * t211) * qJD(2)) * t205) * t224;
t1 = [(t213 * t259 + (-t194 * t226 + t213 * t262) * t205) * t224, t159, 0, 0, 0, 0; t211 * t178 * t269 + (t196 * t178 + (-t160 * t211 - t164 * t194) * t179) * t170 + ((t164 * t268 - t296 * t286) * t170 + (t164 * t269 + ((-t168 * t205 * t265 + 0.2e1 * t289) * t285 + (t259 * t299 + t168 + (-t168 + t244) * t205) * t284) * t170) * t179) * t213, 0.2e1 * (t178 * t215 - t179 * t291) * t293 + (t268 * t291 + t195 * t178 + (t215 * t160 - t161 * t194 + (-t229 * t275 + t159 * t211 + t169 * t196 + (-t169 * t281 + t212) * t168) * t284 + (-t168 * t169 * t211 + t197 + (-t159 * t236 + (qJD(2) * t169 + t168) * t232) * t229) * t285) * t179) * t170, 0, 0, 0, 0; (t183 * t254 - t187 * t288) * t270 + ((t187 * qJD(5) - t176 * t234 + t177 * t230) * t183 + t187 * t258 + (t254 * t163 + (t254 * qJD(5) + t176 * t230 + t177 * t234) * t253 - t187 * t162) * t184) * t165 (t183 * t192 + t253 * t287) * t270 + (t162 * t287 + (t192 * t184 - t193 * t267) * t163 + (t252 * t183 + t251 * t288) * t194 + ((t298 * t183 + t297 * t288) * t234 + (t297 * t183 - t298 * t288) * t230) * t213) * t165 (t183 * t191 + t253 * t288) * t270 + (-t163 * t183 - t253 * t258 + (0.2e1 * t253 * t162 + t191 * t163) * t184) * t165, 0, -0.2e1 * t294 - 0.2e1 * (t162 * t184 * t165 - (-t165 * t290 - t184 * t294) * t253) * t253, 0;];
JaD_rot  = t1;
