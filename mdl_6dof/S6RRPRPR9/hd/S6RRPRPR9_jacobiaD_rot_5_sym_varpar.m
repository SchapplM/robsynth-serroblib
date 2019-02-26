% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR9_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR9_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:42:33
% EndTime: 2019-02-26 21:42:34
% DurationCPUTime: 1.32s
% Computational Cost: add. (7870->139), mult. (12381->281), div. (702->12), fcn. (15714->13), ass. (0->122)
t235 = cos(pkin(6));
t236 = sin(qJ(2));
t302 = sin(qJ(1));
t266 = t302 * t236;
t256 = t235 * t266;
t261 = qJD(2) * t302;
t237 = cos(qJ(2));
t238 = cos(qJ(1));
t280 = t238 * t237;
t233 = sin(pkin(6));
t282 = t233 * t238;
t307 = -qJD(1) * t256 - t236 * t261 + (qJD(2) * t235 + qJD(1)) * t280 - qJD(4) * t282;
t218 = -t256 + t280;
t231 = pkin(11) + qJ(4);
t229 = sin(t231);
t230 = cos(t231);
t267 = t233 * t302;
t209 = t218 * t230 + t229 * t267;
t232 = sin(pkin(12));
t234 = cos(pkin(12));
t265 = t302 * t237;
t281 = t238 * t236;
t247 = -t235 * t265 - t281;
t285 = t247 * t234;
t187 = t209 * t232 + t285;
t249 = -t235 * t281 - t265;
t196 = t249 * qJD(1) + t247 * qJD(2);
t248 = -t218 * t229 + t230 * t267;
t264 = qJD(1) * t282;
t175 = t248 * qJD(4) + t196 * t230 + t229 * t264;
t260 = t302 * qJD(1);
t268 = t235 * t280;
t278 = qJD(2) * t237;
t195 = -qJD(1) * t268 - t238 * t278 + (t235 * t261 + t260) * t236;
t170 = t175 * t234 - t195 * t232;
t286 = t247 * t232;
t188 = t209 * t234 - t286;
t180 = 0.1e1 / t188;
t181 = 0.1e1 / t188 ^ 2;
t297 = t170 * t180 * t181;
t259 = 0.2e1 * t187 * t297;
t203 = -t229 * t249 + t230 * t282;
t284 = t233 * t236;
t213 = t229 * t284 - t230 * t235;
t192 = atan2(-t203, t213);
t183 = sin(t192);
t184 = cos(t192);
t168 = -t183 * t203 + t184 * t213;
t166 = 0.1e1 / t168 ^ 2;
t202 = t248 ^ 2;
t164 = t166 * t202 + 0.1e1;
t174 = t209 * qJD(4) + t196 * t229 - t230 * t264;
t296 = t174 * t166;
t201 = t203 ^ 2;
t211 = 0.1e1 / t213 ^ 2;
t191 = t201 * t211 + 0.1e1;
t189 = 0.1e1 / t191;
t255 = t233 * t260;
t277 = qJD(4) * t230;
t176 = t307 * t229 - t230 * t255 - t249 * t277;
t214 = t229 * t235 + t230 * t284;
t263 = t233 * t278;
t199 = t214 * qJD(4) + t229 * t263;
t210 = 0.1e1 / t213;
t288 = t203 * t211;
t253 = -t176 * t210 + t199 * t288;
t158 = t253 * t189;
t254 = -t183 * t213 - t184 * t203;
t153 = t254 * t158 - t176 * t183 + t184 * t199;
t165 = 0.1e1 / t168;
t167 = t165 * t166;
t300 = t153 * t167;
t276 = 0.2e1 * (-t202 * t300 - t248 * t296) / t164 ^ 2;
t306 = t199 * t211;
t215 = -t266 + t268;
t283 = t233 * t237;
t250 = -t210 * t215 + t283 * t288;
t305 = t229 * t250;
t177 = (qJD(4) * t249 + t255) * t229 + t307 * t230;
t304 = -0.2e1 * t203;
t303 = -0.2e1 * t248;
t290 = t210 * t306;
t299 = (t176 * t288 - t201 * t290) / t191 ^ 2;
t298 = t166 * t248;
t295 = t180 * t232;
t294 = t181 * t187;
t293 = t183 * t248;
t292 = t184 * t248;
t291 = t187 * t234;
t289 = t203 * t210;
t287 = t247 * t229;
t279 = qJD(2) * t236;
t169 = t175 * t232 + t195 * t234;
t179 = t187 ^ 2;
t173 = t179 * t181 + 0.1e1;
t275 = 0.2e1 * (t169 * t294 - t179 * t297) / t173 ^ 2;
t274 = -0.2e1 * t299;
t273 = t167 * t303;
t272 = t210 * t299;
t271 = t166 * t293;
t270 = t166 * t292;
t258 = t290 * t304;
t205 = -t229 * t282 - t230 * t249;
t252 = -t205 * t210 + t214 * t288;
t251 = -qJD(4) * t287 + t195 * t230;
t245 = -t183 + (t184 * t289 + t183) * t189;
t200 = -t213 * qJD(4) + t230 * t263;
t197 = t247 * qJD(1) + t249 * qJD(2);
t194 = t218 * t232 + t230 * t285;
t193 = -t218 * t234 + t230 * t286;
t186 = -t205 * t234 + t215 * t232;
t185 = -t205 * t232 - t215 * t234;
t171 = 0.1e1 / t173;
t162 = 0.1e1 / t164;
t161 = t189 * t305;
t159 = t252 * t189;
t157 = t245 * t248;
t155 = (-t183 * t215 + t184 * t283) * t229 + t254 * t161;
t154 = t254 * t159 - t183 * t205 + t184 * t214;
t152 = t252 * t274 + (t214 * t258 - t177 * t210 + (t176 * t214 + t199 * t205 + t200 * t203) * t211) * t189;
t150 = t274 * t305 + (t250 * t277 + (t258 * t283 - t197 * t210 + (t199 * t215 + (t176 * t237 - t203 * t279) * t233) * t211) * t229) * t189;
t1 = [t272 * t303 + (-t174 * t210 - t248 * t306) * t189, t150, 0, t152, 0, 0; t203 * t165 * t276 + (-t176 * t165 + (t153 * t203 + t157 * t174) * t166) * t162 - (-t157 * t166 * t276 + (-0.2e1 * t157 * t300 + (-t158 * t189 * t289 + t274) * t271 + (t272 * t304 - t158 + (t158 - t253) * t189) * t270 - t245 * t296) * t162) * t248 (-t155 * t298 - t165 * t287) * t276 + (-t155 * t296 + (t195 * t229 + t247 * t277) * t165 + (t155 * t273 - t166 * t287) * t153 + (-t150 * t203 - t161 * t176 + (-t229 * t279 + t237 * t277) * t233 + (-t161 * t213 - t215 * t229) * t158) * t270 + (-t215 * t277 - t150 * t213 - t161 * t199 - t197 * t229 + (t161 * t203 - t229 * t283) * t158) * t271) * t162, 0 (-t154 * t298 - t165 * t209) * t276 + (t154 * t153 * t273 + t175 * t165 + (-t209 * t153 - t154 * t174 + (-t152 * t203 - t159 * t176 + t200 + (-t159 * t213 - t205) * t158) * t292 + (-t152 * t213 - t159 * t199 - t177 + (t159 * t203 - t214) * t158) * t293) * t166) * t162, 0, 0; (-t180 * t185 + t186 * t294) * t275 + ((-t177 * t232 - t197 * t234) * t180 + t186 * t259 + (-t185 * t170 - (-t177 * t234 + t197 * t232) * t187 - t186 * t169) * t181) * t171 (-t180 * t193 + t194 * t294) * t275 + ((-t196 * t234 + t251 * t232) * t180 + t194 * t259 + (-t193 * t170 - (t196 * t232 + t251 * t234) * t187 - t194 * t169) * t181) * t171, 0 -(-t181 * t291 + t295) * t248 * t275 + (t248 * t234 * t259 - t174 * t295 + (t174 * t291 - (t169 * t234 + t170 * t232) * t248) * t181) * t171, 0, 0;];
JaD_rot  = t1;
