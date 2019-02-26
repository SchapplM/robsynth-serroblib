% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRRR1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:18:37
% EndTime: 2019-02-26 20:18:38
% DurationCPUTime: 1.35s
% Computational Cost: add. (18064->114), mult. (17539->232), div. (1085->12), fcn. (22566->13), ass. (0->112)
t266 = sin(pkin(12));
t268 = cos(pkin(12));
t271 = sin(qJ(2));
t269 = cos(pkin(6));
t273 = cos(qJ(2));
t299 = t269 * t273;
t255 = -t266 * t271 + t268 * t299;
t250 = t255 * qJD(2);
t300 = t269 * t271;
t256 = t266 * t273 + t268 * t300;
t265 = qJ(3) + qJ(4) + qJ(5);
t262 = sin(t265);
t264 = qJD(3) + qJD(4) + qJD(5);
t267 = sin(pkin(6));
t303 = t267 * t268;
t289 = t262 * t303;
t263 = cos(t265);
t305 = t263 * t264;
t217 = t250 * t262 + t256 * t305 - t264 * t289;
t239 = t256 * t262 + t263 * t303;
t237 = t239 ^ 2;
t302 = t267 * t271;
t247 = t262 * t302 - t263 * t269;
t245 = 0.1e1 / t247 ^ 2;
t225 = t237 * t245 + 0.1e1;
t223 = 0.1e1 / t225;
t297 = qJD(2) * t273;
t282 = t264 * t269 + t267 * t297;
t291 = t264 * t302;
t235 = t282 * t262 + t263 * t291;
t244 = 0.1e1 / t247;
t310 = t239 * t245;
t201 = (-t217 * t244 + t235 * t310) * t223;
t230 = atan2(-t239, t247);
t221 = sin(t230);
t222 = cos(t230);
t285 = -t221 * t247 - t222 * t239;
t197 = t285 * t201 - t217 * t221 + t222 * t235;
t211 = -t221 * t239 + t222 * t247;
t208 = 0.1e1 / t211;
t209 = 0.1e1 / t211 ^ 2;
t324 = t197 * t208 * t209;
t290 = t266 * t300;
t258 = t268 * t273 - t290;
t304 = t266 * t267;
t242 = t258 * t262 - t263 * t304;
t323 = 0.2e1 * t242 * t324;
t301 = t267 * t273;
t281 = -t244 * t255 + t301 * t310;
t322 = t262 * t281;
t311 = t235 * t244 * t245;
t321 = -0.2e1 * (t217 * t310 - t237 * t311) / t225 ^ 2;
t243 = t258 * t263 + t262 * t304;
t272 = cos(qJ(6));
t257 = t266 * t299 + t268 * t271;
t270 = sin(qJ(6));
t308 = t257 * t270;
t232 = t243 * t272 + t308;
t227 = 0.1e1 / t232;
t228 = 0.1e1 / t232 ^ 2;
t252 = t257 * qJD(2);
t287 = t264 * t304 - t252;
t306 = t262 * t264;
t220 = -t258 * t306 + t287 * t263;
t253 = -qJD(2) * t290 + t268 * t297;
t212 = t232 * qJD(6) + t220 * t270 - t253 * t272;
t307 = t257 * t272;
t231 = t243 * t270 - t307;
t226 = t231 ^ 2;
t216 = t226 * t228 + 0.1e1;
t313 = t228 * t231;
t296 = qJD(6) * t231;
t213 = t220 * t272 + t253 * t270 - t296;
t318 = t213 * t227 * t228;
t320 = (t212 * t313 - t226 * t318) / t216 ^ 2;
t319 = t209 * t242;
t219 = t258 * t305 + t287 * t262;
t317 = t219 * t209;
t316 = t221 * t242;
t315 = t222 * t242;
t314 = t227 * t270;
t312 = t231 * t272;
t309 = t257 * t262;
t298 = qJD(2) * t271;
t238 = t242 ^ 2;
t207 = t209 * t238 + 0.1e1;
t295 = 0.2e1 * (-t238 * t324 + t242 * t317) / t207 ^ 2;
t294 = -0.2e1 * t320;
t292 = t231 * t318;
t288 = -0.2e1 * t239 * t311;
t286 = qJD(6) * t257 * t263 - t252;
t284 = t228 * t312 - t314;
t241 = t256 * t263 - t289;
t248 = t262 * t269 + t263 * t302;
t283 = -t241 * t244 + t248 * t310;
t280 = qJD(6) * t258 - t253 * t263 + t257 * t306;
t251 = t256 * qJD(2);
t236 = -t262 * t291 + t282 * t263;
t234 = t258 * t270 - t263 * t307;
t233 = -t258 * t272 - t263 * t308;
t218 = -t256 * t306 + (-t264 * t303 + t250) * t263;
t214 = 0.1e1 / t216;
t205 = 0.1e1 / t207;
t203 = t223 * t322;
t202 = t283 * t223;
t199 = (-t221 * t255 + t222 * t301) * t262 + t285 * t203;
t198 = t285 * t202 - t221 * t241 + t222 * t248;
t195 = t283 * t321 + (t248 * t288 - t218 * t244 + (t217 * t248 + t235 * t241 + t236 * t239) * t245) * t223;
t194 = t321 * t322 + (t281 * t305 + (t288 * t301 + t244 * t251 + (t235 * t255 + (t217 * t273 - t239 * t298) * t267) * t245) * t262) * t223;
t193 = t284 * t242 * t294 + (t284 * t219 + ((-qJD(6) * t227 - 0.2e1 * t292) * t272 + (t212 * t272 + (t213 - t296) * t270) * t228) * t242) * t214;
t192 = (t198 * t319 - t208 * t243) * t295 + (t198 * t323 + t220 * t208 + (-t243 * t197 - t198 * t219 - (-t195 * t239 - t202 * t217 + t236 + (-t202 * t247 - t241) * t201) * t315 - (-t195 * t247 - t202 * t235 - t218 + (t202 * t239 - t248) * t201) * t316) * t209) * t205;
t1 = [0, t194, t195, t195, t195, 0; 0 (t199 * t319 + t208 * t309) * t295 + ((-t253 * t262 - t257 * t305) * t208 + (-t317 + t323) * t199 + (t309 * t197 - (-t194 * t239 - t203 * t217 + (-t262 * t298 + t273 * t305) * t267 + (-t203 * t247 - t255 * t262) * t201) * t315 - (-t255 * t305 - t194 * t247 - t203 * t235 + t251 * t262 + (t203 * t239 - t262 * t301) * t201) * t316) * t209) * t205, t192, t192, t192, 0; 0, 0.2e1 * (-t227 * t233 + t234 * t313) * t320 + (0.2e1 * t234 * t292 - t286 * t227 * t272 + t280 * t314 + (-t286 * t231 * t270 - t234 * t212 - t233 * t213 - t280 * t312) * t228) * t214, t193, t193, t193, t294 + 0.2e1 * (t212 * t228 * t214 + (-t214 * t318 - t228 * t320) * t231) * t231;];
JaD_rot  = t1;
