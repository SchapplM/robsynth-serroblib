% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR2
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
% Datum: 2019-02-26 20:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRPR2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:11:07
% EndTime: 2019-02-26 20:11:09
% DurationCPUTime: 1.22s
% Computational Cost: add. (9451->115), mult. (13312->231), div. (822->12), fcn. (17117->13), ass. (0->113)
t270 = sin(pkin(11));
t272 = cos(pkin(11));
t274 = sin(qJ(2));
t273 = cos(pkin(6));
t275 = cos(qJ(2));
t302 = t273 * t275;
t255 = -t270 * t274 + t272 * t302;
t251 = t255 * qJD(2);
t303 = t273 * t274;
t256 = t270 * t275 + t272 * t303;
t269 = qJ(3) + qJ(4);
t265 = sin(t269);
t268 = qJD(3) + qJD(4);
t271 = sin(pkin(6));
t306 = t271 * t272;
t291 = t265 * t306;
t266 = cos(t269);
t308 = t266 * t268;
t218 = t251 * t265 + t256 * t308 - t268 * t291;
t240 = t256 * t265 + t266 * t306;
t238 = t240 ^ 2;
t305 = t271 * t274;
t293 = t265 * t305;
t249 = -t273 * t266 + t293;
t247 = 0.1e1 / t249 ^ 2;
t232 = t238 * t247 + 0.1e1;
t230 = 0.1e1 / t232;
t300 = qJD(2) * t275;
t284 = t268 * t273 + t271 * t300;
t292 = t266 * t305;
t236 = t284 * t265 + t268 * t292;
t246 = 0.1e1 / t249;
t312 = t240 * t247;
t202 = (-t218 * t246 + t236 * t312) * t230;
t233 = atan2(-t240, t249);
t228 = sin(t233);
t229 = cos(t233);
t287 = -t228 * t249 - t229 * t240;
t198 = t287 * t202 - t228 * t218 + t229 * t236;
t212 = -t228 * t240 + t229 * t249;
t209 = 0.1e1 / t212;
t210 = 0.1e1 / t212 ^ 2;
t326 = t198 * t209 * t210;
t294 = t270 * t303;
t258 = t272 * t275 - t294;
t307 = t270 * t271;
t243 = t258 * t265 - t266 * t307;
t325 = 0.2e1 * t243 * t326;
t304 = t271 * t275;
t283 = -t246 * t255 + t304 * t312;
t324 = t265 * t283;
t313 = t236 * t246 * t247;
t323 = -0.2e1 * (t218 * t312 - t238 * t313) / t232 ^ 2;
t244 = t258 * t266 + t265 * t307;
t257 = t270 * t302 + t272 * t274;
t267 = pkin(12) + qJ(6);
t263 = sin(t267);
t264 = cos(t267);
t227 = t244 * t264 + t257 * t263;
t223 = 0.1e1 / t227;
t224 = 0.1e1 / t227 ^ 2;
t253 = t257 * qJD(2);
t289 = t268 * t307 - t253;
t309 = t265 * t268;
t221 = -t258 * t309 + t289 * t266;
t254 = -qJD(2) * t294 + t272 * t300;
t213 = t227 * qJD(6) + t221 * t263 - t254 * t264;
t226 = t244 * t263 - t257 * t264;
t222 = t226 ^ 2;
t217 = t222 * t224 + 0.1e1;
t317 = t224 * t226;
t299 = qJD(6) * t226;
t214 = t221 * t264 + t254 * t263 - t299;
t320 = t214 * t223 * t224;
t322 = (t213 * t317 - t222 * t320) / t217 ^ 2;
t321 = t210 * t243;
t220 = t258 * t308 + t289 * t265;
t319 = t220 * t210;
t318 = t223 * t263;
t316 = t226 * t264;
t315 = t228 * t243;
t314 = t229 * t243;
t311 = t257 * t265;
t310 = t257 * t266;
t301 = qJD(2) * t274;
t239 = t243 ^ 2;
t208 = t239 * t210 + 0.1e1;
t298 = 0.2e1 * (-t239 * t326 + t243 * t319) / t208 ^ 2;
t297 = -0.2e1 * t322;
t295 = t226 * t320;
t290 = -0.2e1 * t240 * t313;
t288 = qJD(6) * t310 - t253;
t286 = t224 * t316 - t318;
t242 = t256 * t266 - t291;
t250 = t273 * t265 + t292;
t285 = -t242 * t246 + t250 * t312;
t282 = qJD(6) * t258 - t254 * t266 + t257 * t309;
t252 = t256 * qJD(2);
t237 = t284 * t266 - t268 * t293;
t235 = t258 * t263 - t264 * t310;
t234 = -t258 * t264 - t263 * t310;
t219 = -t256 * t309 + (-t268 * t306 + t251) * t266;
t215 = 0.1e1 / t217;
t206 = 0.1e1 / t208;
t204 = t230 * t324;
t203 = t285 * t230;
t200 = (-t228 * t255 + t229 * t304) * t265 + t287 * t204;
t199 = t287 * t203 - t228 * t242 + t229 * t250;
t197 = t285 * t323 + (t250 * t290 - t219 * t246 + (t218 * t250 + t236 * t242 + t237 * t240) * t247) * t230;
t195 = t323 * t324 + (t283 * t308 + (t290 * t304 + t246 * t252 + (t236 * t255 + (t218 * t275 - t240 * t301) * t271) * t247) * t265) * t230;
t194 = t286 * t243 * t297 + (t286 * t220 + ((-qJD(6) * t223 - 0.2e1 * t295) * t264 + (t213 * t264 + (t214 - t299) * t263) * t224) * t243) * t215;
t193 = (t199 * t321 - t209 * t244) * t298 + (t199 * t325 + t221 * t209 + (-t244 * t198 - t199 * t220 - (-t197 * t240 - t203 * t218 + t237 + (-t203 * t249 - t242) * t202) * t314 - (-t197 * t249 - t203 * t236 - t219 + (t203 * t240 - t250) * t202) * t315) * t210) * t206;
t1 = [0, t195, t197, t197, 0, 0; 0 (t200 * t321 + t209 * t311) * t298 + ((-t254 * t265 - t257 * t308) * t209 + (-t319 + t325) * t200 + (t311 * t198 - (-t195 * t240 - t204 * t218 + (-t265 * t301 + t275 * t308) * t271 + (-t204 * t249 - t255 * t265) * t202) * t314 - (-t255 * t308 - t195 * t249 - t204 * t236 + t252 * t265 + (t204 * t240 - t265 * t304) * t202) * t315) * t210) * t206, t193, t193, 0, 0; 0, 0.2e1 * (-t223 * t234 + t235 * t317) * t322 + (0.2e1 * t235 * t295 - t288 * t223 * t264 + t282 * t318 + (-t288 * t226 * t263 - t235 * t213 - t234 * t214 - t282 * t316) * t224) * t215, t194, t194, 0, t297 + 0.2e1 * (t213 * t224 * t215 + (-t215 * t320 - t224 * t322) * t226) * t226;];
JaD_rot  = t1;
