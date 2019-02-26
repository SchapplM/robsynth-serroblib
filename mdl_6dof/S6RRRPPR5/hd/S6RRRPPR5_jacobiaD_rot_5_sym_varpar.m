% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPPR5_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR5_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:05:46
% EndTime: 2019-02-26 22:05:47
% DurationCPUTime: 1.31s
% Computational Cost: add. (7870->139), mult. (12381->281), div. (702->12), fcn. (15714->13), ass. (0->122)
t236 = cos(pkin(6));
t237 = sin(qJ(2));
t303 = sin(qJ(1));
t267 = t303 * t237;
t257 = t236 * t267;
t262 = qJD(2) * t303;
t238 = cos(qJ(2));
t239 = cos(qJ(1));
t281 = t239 * t238;
t234 = sin(pkin(6));
t283 = t234 * t239;
t308 = -qJD(1) * t257 - t237 * t262 + (qJD(2) * t236 + qJD(1)) * t281 - qJD(3) * t283;
t219 = -t257 + t281;
t232 = qJ(3) + pkin(11);
t230 = sin(t232);
t231 = cos(t232);
t268 = t234 * t303;
t210 = t219 * t231 + t230 * t268;
t233 = sin(pkin(12));
t235 = cos(pkin(12));
t266 = t303 * t238;
t282 = t239 * t237;
t248 = -t236 * t266 - t282;
t286 = t248 * t235;
t188 = t210 * t233 + t286;
t250 = -t236 * t282 - t266;
t197 = t250 * qJD(1) + t248 * qJD(2);
t249 = -t219 * t230 + t231 * t268;
t265 = qJD(1) * t283;
t176 = t249 * qJD(3) + t197 * t231 + t230 * t265;
t261 = t303 * qJD(1);
t269 = t236 * t281;
t279 = qJD(2) * t238;
t196 = -qJD(1) * t269 - t239 * t279 + (t236 * t262 + t261) * t237;
t171 = t176 * t235 - t196 * t233;
t287 = t248 * t233;
t189 = t210 * t235 - t287;
t181 = 0.1e1 / t189;
t182 = 0.1e1 / t189 ^ 2;
t298 = t171 * t181 * t182;
t260 = 0.2e1 * t188 * t298;
t204 = -t230 * t250 + t231 * t283;
t285 = t234 * t237;
t214 = t230 * t285 - t231 * t236;
t193 = atan2(-t204, t214);
t184 = sin(t193);
t185 = cos(t193);
t169 = -t184 * t204 + t185 * t214;
t167 = 0.1e1 / t169 ^ 2;
t203 = t249 ^ 2;
t165 = t167 * t203 + 0.1e1;
t175 = t210 * qJD(3) + t197 * t230 - t231 * t265;
t297 = t175 * t167;
t202 = t204 ^ 2;
t212 = 0.1e1 / t214 ^ 2;
t192 = t202 * t212 + 0.1e1;
t190 = 0.1e1 / t192;
t256 = t234 * t261;
t278 = qJD(3) * t231;
t177 = t308 * t230 - t231 * t256 - t250 * t278;
t215 = t230 * t236 + t231 * t285;
t264 = t234 * t279;
t200 = t215 * qJD(3) + t230 * t264;
t211 = 0.1e1 / t214;
t289 = t204 * t212;
t254 = -t177 * t211 + t200 * t289;
t159 = t254 * t190;
t255 = -t184 * t214 - t185 * t204;
t154 = t255 * t159 - t177 * t184 + t185 * t200;
t166 = 0.1e1 / t169;
t168 = t166 * t167;
t301 = t154 * t168;
t277 = 0.2e1 * (-t203 * t301 - t249 * t297) / t165 ^ 2;
t307 = t200 * t212;
t216 = -t267 + t269;
t284 = t234 * t238;
t251 = -t211 * t216 + t284 * t289;
t306 = t230 * t251;
t178 = (qJD(3) * t250 + t256) * t230 + t308 * t231;
t305 = -0.2e1 * t204;
t304 = -0.2e1 * t249;
t291 = t211 * t307;
t300 = (t177 * t289 - t202 * t291) / t192 ^ 2;
t299 = t167 * t249;
t296 = t181 * t233;
t295 = t182 * t188;
t294 = t184 * t249;
t293 = t185 * t249;
t292 = t188 * t235;
t290 = t204 * t211;
t288 = t248 * t230;
t280 = qJD(2) * t237;
t170 = t176 * t233 + t196 * t235;
t180 = t188 ^ 2;
t174 = t180 * t182 + 0.1e1;
t276 = 0.2e1 * (t170 * t295 - t180 * t298) / t174 ^ 2;
t275 = -0.2e1 * t300;
t274 = t168 * t304;
t273 = t211 * t300;
t272 = t167 * t294;
t271 = t167 * t293;
t259 = t291 * t305;
t206 = -t230 * t283 - t231 * t250;
t253 = -t206 * t211 + t215 * t289;
t252 = -qJD(3) * t288 + t196 * t231;
t246 = -t184 + (t185 * t290 + t184) * t190;
t201 = -t214 * qJD(3) + t231 * t264;
t198 = t248 * qJD(1) + t250 * qJD(2);
t195 = t219 * t233 + t231 * t286;
t194 = -t219 * t235 + t231 * t287;
t187 = -t206 * t235 + t216 * t233;
t186 = -t206 * t233 - t216 * t235;
t172 = 0.1e1 / t174;
t163 = 0.1e1 / t165;
t162 = t190 * t306;
t160 = t253 * t190;
t158 = t246 * t249;
t156 = (-t184 * t216 + t185 * t284) * t230 + t255 * t162;
t155 = t255 * t160 - t184 * t206 + t185 * t215;
t153 = t253 * t275 + (t215 * t259 - t178 * t211 + (t177 * t215 + t200 * t206 + t201 * t204) * t212) * t190;
t151 = t275 * t306 + (t251 * t278 + (t259 * t284 - t198 * t211 + (t200 * t216 + (t177 * t238 - t204 * t280) * t234) * t212) * t230) * t190;
t1 = [t273 * t304 + (-t175 * t211 - t249 * t307) * t190, t151, t153, 0, 0, 0; t204 * t166 * t277 + (-t177 * t166 + (t154 * t204 + t158 * t175) * t167) * t163 - (-t158 * t167 * t277 + (-0.2e1 * t158 * t301 + (-t159 * t190 * t290 + t275) * t272 + (t273 * t305 - t159 + (t159 - t254) * t190) * t271 - t246 * t297) * t163) * t249 (-t156 * t299 - t166 * t288) * t277 + (-t156 * t297 + (t196 * t230 + t248 * t278) * t166 + (t156 * t274 - t167 * t288) * t154 + (-t151 * t204 - t162 * t177 + (-t230 * t280 + t238 * t278) * t234 + (-t162 * t214 - t216 * t230) * t159) * t271 + (-t216 * t278 - t151 * t214 - t162 * t200 - t198 * t230 + (t162 * t204 - t230 * t284) * t159) * t272) * t163 (-t155 * t299 - t166 * t210) * t277 + (t155 * t154 * t274 + t176 * t166 + (-t210 * t154 - t155 * t175 + (-t153 * t204 - t160 * t177 + t201 + (-t160 * t214 - t206) * t159) * t293 + (-t153 * t214 - t160 * t200 - t178 + (t160 * t204 - t215) * t159) * t294) * t167) * t163, 0, 0, 0; (-t181 * t186 + t187 * t295) * t276 + ((-t178 * t233 - t198 * t235) * t181 + t187 * t260 + (-t186 * t171 - (-t178 * t235 + t198 * t233) * t188 - t187 * t170) * t182) * t172 (-t181 * t194 + t195 * t295) * t276 + ((-t197 * t235 + t252 * t233) * t181 + t195 * t260 + (-t194 * t171 - (t197 * t233 + t252 * t235) * t188 - t195 * t170) * t182) * t172 -(-t182 * t292 + t296) * t249 * t276 + (t249 * t235 * t260 - t175 * t296 + (t175 * t292 - (t170 * t235 + t171 * t233) * t249) * t182) * t172, 0, 0, 0;];
JaD_rot  = t1;
