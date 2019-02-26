% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRP5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:03:31
% EndTime: 2019-02-26 20:03:33
% DurationCPUTime: 1.64s
% Computational Cost: add. (7266->159), mult. (21195->316), div. (787->12), fcn. (27318->13), ass. (0->127)
t252 = sin(qJ(5));
t253 = sin(qJ(3));
t256 = cos(qJ(3));
t254 = sin(qJ(2));
t257 = cos(qJ(2));
t317 = cos(pkin(10));
t318 = cos(pkin(6));
t280 = t318 * t317;
t316 = sin(pkin(10));
t271 = t254 * t280 + t316 * t257;
t251 = sin(pkin(6));
t287 = t251 * t317;
t266 = t271 * t253 + t256 * t287;
t255 = cos(qJ(5));
t270 = -t316 * t254 + t257 * t280;
t306 = t270 * t255;
t219 = t266 * t252 - t306;
t242 = t271 * qJD(2);
t267 = qJD(3) * t271;
t268 = qJD(2) * t270;
t322 = -qJD(3) * t287 + t268;
t264 = t322 * t253 + t256 * t267;
t191 = t219 * qJD(5) + t242 * t252 - t264 * t255;
t265 = t266 * t255;
t217 = -t252 * t270 - t265;
t212 = t217 ^ 2;
t305 = t251 * t254;
t248 = t253 * t305 - t318 * t256;
t303 = t252 * t257;
t272 = t248 * t255 + t251 * t303;
t233 = 0.1e1 / t272 ^ 2;
t203 = t212 * t233 + 0.1e1;
t201 = 0.1e1 / t203;
t249 = t318 * t253 + t256 * t305;
t300 = qJD(2) * t251;
t289 = t257 * t300;
t235 = t249 * qJD(3) + t253 * t289;
t301 = t255 * t257;
t238 = t248 * t252 - t251 * t301;
t290 = t254 * t300;
t205 = t238 * qJD(5) - t235 * t255 + t252 * t290;
t232 = 0.1e1 / t272;
t308 = t217 * t233;
t178 = (t191 * t232 + t205 * t308) * t201;
t204 = atan2(-t217, -t272);
t198 = sin(t204);
t199 = cos(t204);
t278 = t198 * t272 - t199 * t217;
t174 = t278 * t178 - t198 * t191 + t199 * t205;
t190 = -t198 * t217 - t199 * t272;
t187 = 0.1e1 / t190;
t188 = 0.1e1 / t190 ^ 2;
t326 = t174 * t187 * t188;
t279 = t318 * t316;
t247 = -t254 * t279 + t317 * t257;
t286 = t251 * t316;
t229 = t247 * t253 - t256 * t286;
t269 = -t317 * t254 - t257 * t279;
t277 = t229 * t255 + t252 * t269;
t320 = -0.2e1 * t277;
t284 = t320 * t326;
t310 = t205 * t232 * t233;
t325 = (t191 * t308 + t212 * t310) / t203 ^ 2;
t230 = t247 * t256 + t253 * t286;
t243 = t269 * qJD(2);
t210 = t230 * qJD(3) + t243 * t253;
t244 = t247 * qJD(2);
t194 = t277 * qJD(5) + t210 * t252 + t244 * t255;
t221 = t229 * t252 - t255 * t269;
t215 = 0.1e1 / t221 ^ 2;
t324 = t194 * t215;
t228 = -t253 * t287 + t271 * t256;
t273 = t228 * t232 + t249 * t308;
t323 = t255 * t273;
t321 = -0.2e1 * t325;
t214 = 0.1e1 / t221;
t319 = -0.2e1 * t252;
t227 = t230 ^ 2;
t200 = t227 * t215 + 0.1e1;
t299 = qJD(3) * t253;
t211 = -t247 * t299 + (qJD(3) * t286 + t243) * t256;
t309 = t215 * t230;
t313 = t214 * t324;
t315 = (t211 * t309 - t227 * t313) / t200 ^ 2;
t314 = t188 * t277;
t312 = t198 * t277;
t311 = t199 * t277;
t307 = t230 * t255;
t304 = t252 * t253;
t302 = t253 * t255;
t298 = qJD(3) * t256;
t297 = qJD(5) * t252;
t296 = qJD(5) * t253;
t295 = qJD(5) * t255;
t213 = t277 ^ 2;
t186 = t213 * t188 + 0.1e1;
t193 = t221 * qJD(5) - t210 * t255 + t244 * t252;
t294 = 0.2e1 * (-t193 * t314 - t213 * t326) / t186 ^ 2;
t292 = t230 * t313;
t291 = t217 * t310;
t288 = t270 * t297;
t283 = 0.2e1 * t291;
t282 = t309 * t315;
t275 = t219 * t232 + t238 * t308;
t222 = t271 * t252 - t270 * t302;
t240 = (t252 * t254 - t253 * t301) * t251;
t274 = t222 * t232 + t240 * t308;
t236 = -t248 * qJD(3) + t256 * t289;
t224 = t247 * t255 + t269 * t304;
t223 = t247 * t252 - t269 * t302;
t209 = -t253 * t267 + t322 * t256;
t207 = ((qJD(2) + t296) * t303 + (-t257 * t298 + (qJD(2) * t253 + qJD(5)) * t254) * t255) * t251;
t206 = t272 * qJD(5) + t235 * t252 + t255 * t290;
t196 = 0.1e1 / t200;
t195 = t242 * t302 + t252 * t268 + t253 * t288 + t271 * t295 - t298 * t306;
t192 = qJD(5) * t265 + t242 * t255 + t264 * t252 + t288;
t184 = 0.1e1 / t186;
t183 = t201 * t323;
t182 = t274 * t201;
t180 = t275 * t201;
t177 = (t198 * t228 - t199 * t249) * t255 - t278 * t183;
t176 = t278 * t182 - t198 * t222 + t199 * t240;
t175 = t278 * t180 - t198 * t219 + t199 * t238;
t173 = t274 * t321 + (t240 * t283 + t195 * t232 + (t191 * t240 + t205 * t222 + t207 * t217) * t233) * t201;
t171 = t275 * t321 + (t238 * t283 + t192 * t232 + (t191 * t238 + t205 * t219 + t206 * t217) * t233) * t201;
t170 = 0.2e1 * t323 * t325 + (t273 * t297 + (-0.2e1 * t249 * t291 - t209 * t232 + (-t191 * t249 - t205 * t228 - t217 * t236) * t233) * t255) * t201;
t1 = [0, t173, t170, 0, t171, 0; 0 (-t176 * t314 - t187 * t223) * t294 + (t176 * t284 + (-t223 * t174 - t176 * t193 + (-t173 * t217 - t182 * t191 + t207 + (t182 * t272 - t222) * t178) * t311 + (t173 * t272 - t182 * t205 - t195 + (t182 * t217 - t240) * t178) * t312) * t188 + ((t269 * t296 + t243) * t252 + (qJD(5) * t247 + t244 * t253 - t269 * t298) * t255) * t187) * t184 (-t177 * t314 + t187 * t307) * t294 + ((-t211 * t255 + t230 * t297) * t187 + t177 * t284 + (-t177 * t193 + t307 * t174 + (t249 * t297 - t170 * t217 + t183 * t191 - t236 * t255 + (-t183 * t272 + t228 * t255) * t178) * t311 + (-t228 * t297 + t170 * t272 + t183 * t205 + t209 * t255 + (-t183 * t217 + t249 * t255) * t178) * t312) * t188) * t184, 0 (-t175 * t314 - t187 * t221) * t294 + (t175 * t284 + t194 * t187 + (-t221 * t174 - t175 * t193 + (-t171 * t217 - t180 * t191 + t206 + (t180 * t272 - t219) * t178) * t311 + (t171 * t272 - t180 * t205 - t192 + (t180 * t217 - t238) * t178) * t312) * t188) * t184, 0; 0, 0.2e1 * (t214 * t256 * t269 - t224 * t309) * t315 + (-0.2e1 * t224 * t292 + (t244 * t256 + t269 * t299) * t214 + ((t243 * t255 - t244 * t304 - t247 * t297) * t230 + t224 * t211 - (-t256 * t194 + (-t252 * t298 - t253 * t295) * t230) * t269) * t215) * t196, -0.2e1 * t214 * t229 * t315 + (t210 * t214 - t229 * t324) * t196 + (t282 * t319 + (0.2e1 * t211 * t215 * t252 + (t215 * t295 + t313 * t319) * t230) * t196) * t230, 0, t282 * t320 + (t292 * t320 + (-t193 * t230 + t211 * t277) * t215) * t196, 0;];
JaD_rot  = t1;
