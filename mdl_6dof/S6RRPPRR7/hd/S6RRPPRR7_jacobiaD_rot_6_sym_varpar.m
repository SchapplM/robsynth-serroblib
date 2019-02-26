% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR7_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:31:57
% EndTime: 2019-02-26 21:31:58
% DurationCPUTime: 1.39s
% Computational Cost: add. (4522->147), mult. (13478->295), div. (726->12), fcn. (17045->13), ass. (0->127)
t250 = cos(pkin(6));
t252 = sin(qJ(5));
t256 = cos(qJ(5));
t249 = sin(pkin(6));
t257 = cos(qJ(2));
t306 = t249 * t257;
t269 = -t250 * t256 + t252 * t306;
t258 = cos(qJ(1));
t300 = t258 * t257;
t283 = t250 * t300;
t253 = sin(qJ(2));
t254 = sin(qJ(1));
t303 = t254 * t253;
t237 = -t283 + t303;
t305 = t249 * t258;
t270 = -t237 * t252 + t256 * t305;
t213 = atan2(t270, -t269);
t208 = sin(t213);
t209 = cos(t213);
t191 = t208 * t270 - t209 * t269;
t189 = 0.1e1 / t191 ^ 2;
t301 = t258 * t253;
t302 = t254 * t257;
t239 = t250 * t302 + t301;
t307 = t249 * t254;
t229 = t239 * t252 + t256 * t307;
t221 = t229 ^ 2;
t187 = t221 * t189 + 0.1e1;
t276 = qJD(2) * t250 + qJD(1);
t297 = qJD(2) * t257;
t216 = -qJD(1) * t283 - t258 * t297 + t276 * t303;
t299 = qJD(1) * t249;
t281 = t258 * t299;
t285 = t252 * t307;
t195 = -t216 * t252 - qJD(5) * t285 + (qJD(5) * t239 + t281) * t256;
t319 = t195 * t189;
t220 = t270 ^ 2;
t233 = 0.1e1 / t269 ^ 2;
t212 = t220 * t233 + 0.1e1;
t210 = 0.1e1 / t212;
t238 = t250 * t301 + t302;
t218 = qJD(1) * t239 + qJD(2) * t238;
t226 = t237 * t256 + t252 * t305;
t282 = t254 * t299;
t197 = qJD(5) * t226 + t218 * t252 + t256 * t282;
t236 = -t250 * t252 - t256 * t306;
t298 = qJD(2) * t253;
t280 = t249 * t298;
t222 = qJD(5) * t236 + t252 * t280;
t232 = 0.1e1 / t269;
t312 = t270 * t233;
t273 = t197 * t232 - t222 * t312;
t179 = t273 * t210;
t274 = t208 * t269 + t209 * t270;
t174 = t179 * t274 - t208 * t197 + t209 * t222;
t188 = 0.1e1 / t191;
t190 = t188 * t189;
t324 = t174 * t190;
t294 = 0.2e1 * (-t221 * t324 + t229 * t319) / t187 ^ 2;
t329 = t222 * t233;
t308 = t249 * t253;
t268 = t232 * t238 - t308 * t312;
t328 = t252 * t268;
t198 = qJD(5) * t270 + t218 * t256 - t252 * t282;
t230 = t239 * t256 - t285;
t284 = t250 * t303;
t240 = -t284 + t300;
t251 = sin(qJ(6));
t255 = cos(qJ(6));
t207 = t230 * t255 + t240 * t251;
t201 = 0.1e1 / t207;
t202 = 0.1e1 / t207 ^ 2;
t327 = 0.2e1 * t270;
t326 = 0.2e1 * t229;
t196 = -qJD(5) * t229 - t216 * t256 - t252 * t281;
t217 = -qJD(1) * t238 - qJD(2) * t239;
t183 = qJD(6) * t207 + t196 * t251 - t217 * t255;
t206 = t230 * t251 - t240 * t255;
t200 = t206 ^ 2;
t194 = t200 * t202 + 0.1e1;
t318 = t202 * t206;
t295 = qJD(6) * t206;
t184 = t196 * t255 + t217 * t251 - t295;
t321 = t184 * t201 * t202;
t323 = (t183 * t318 - t200 * t321) / t194 ^ 2;
t314 = t232 * t329;
t322 = (-t197 * t312 + t220 * t314) / t212 ^ 2;
t320 = t189 * t229;
t317 = t206 * t255;
t316 = t208 * t229;
t315 = t209 * t229;
t313 = t270 * t232;
t310 = t240 * t252;
t309 = t240 * t256;
t304 = t251 * t201;
t296 = qJD(5) * t256;
t293 = -0.2e1 * t323;
t292 = 0.2e1 * t323;
t291 = -0.2e1 * t322;
t290 = t190 * t326;
t289 = t232 * t322;
t288 = t189 * t316;
t287 = t189 * t315;
t286 = t206 * t321;
t279 = 0.2e1 * t286;
t278 = t314 * t327;
t275 = qJD(6) * t309 - t216;
t205 = -t226 * t255 - t238 * t251;
t204 = -t226 * t251 + t238 * t255;
t272 = t202 * t317 - t304;
t271 = t226 * t232 - t236 * t312;
t266 = -t208 + (t209 * t313 + t208) * t210;
t265 = -qJD(5) * t310 - qJD(6) * t239 + t217 * t256;
t223 = qJD(5) * t269 + t256 * t280;
t219 = -qJD(1) * t284 - t254 * t298 + t276 * t300;
t215 = -t239 * t251 + t255 * t309;
t214 = t239 * t255 + t251 * t309;
t192 = 0.1e1 / t194;
t185 = 0.1e1 / t187;
t182 = t210 * t328;
t181 = t271 * t210;
t178 = t266 * t229;
t176 = (-t208 * t238 + t209 * t308) * t252 + t274 * t182;
t175 = t181 * t274 - t208 * t226 + t209 * t236;
t173 = t271 * t291 + (-t236 * t278 + t198 * t232 + (t197 * t236 + t222 * t226 - t223 * t270) * t233) * t210;
t171 = t291 * t328 + (t268 * t296 + (-t278 * t308 + t219 * t232 + (t222 * t238 + (t197 * t253 - t270 * t297) * t249) * t233) * t252) * t210;
t1 = [-t289 * t326 + (t195 * t232 + t229 * t329) * t210, t171, 0, 0, t173, 0; -t270 * t188 * t294 + (-t197 * t188 + (-t174 * t270 - t178 * t195) * t189) * t185 + (t178 * t189 * t294 + (0.2e1 * t178 * t324 - (-t179 * t210 * t313 + t291) * t288 - (-t289 * t327 - t179 + (t179 - t273) * t210) * t287 - t266 * t319) * t185) * t229 (t176 * t320 - t188 * t310) * t294 + (-t176 * t319 + (t217 * t252 + t240 * t296) * t188 + (t176 * t290 - t189 * t310) * t174 - (t171 * t270 - t182 * t197 + (t252 * t297 + t253 * t296) * t249 + (t182 * t269 - t238 * t252) * t179) * t287 - (-t238 * t296 + t171 * t269 - t182 * t222 - t219 * t252 + (-t182 * t270 - t252 * t308) * t179) * t288) * t185, 0, 0 (t175 * t320 - t188 * t230) * t294 + (t175 * t174 * t290 + t196 * t188 + (-t230 * t174 - t175 * t195 - (t173 * t270 - t181 * t197 + t223 + (t181 * t269 - t226) * t179) * t315 - (t173 * t269 - t181 * t222 - t198 + (-t181 * t270 - t236) * t179) * t316) * t189) * t185, 0; (-t201 * t204 + t205 * t318) * t292 + ((qJD(6) * t205 - t198 * t251 + t219 * t255) * t201 + t205 * t279 + (-t204 * t184 - (-qJD(6) * t204 - t198 * t255 - t219 * t251) * t206 - t205 * t183) * t202) * t192 (-t201 * t214 + t215 * t318) * t292 + (t215 * t279 + t275 * t201 * t255 + t265 * t304 + (t206 * t251 * t275 - t215 * t183 - t214 * t184 - t265 * t317) * t202) * t192, 0, 0, t272 * t229 * t293 + (t272 * t195 + ((-qJD(6) * t201 - 0.2e1 * t286) * t255 + (t183 * t255 + (t184 - t295) * t251) * t202) * t229) * t192, t293 + 0.2e1 * (t183 * t202 * t192 + (-t192 * t321 - t202 * t323) * t206) * t206;];
JaD_rot  = t1;
