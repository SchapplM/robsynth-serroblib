% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPPRR1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:44:35
% EndTime: 2019-02-26 19:44:37
% DurationCPUTime: 1.25s
% Computational Cost: add. (8444->115), mult. (16325->235), div. (559->12), fcn. (21250->15), ass. (0->111)
t281 = sin(pkin(11));
t284 = cos(pkin(11));
t288 = sin(qJ(2));
t290 = cos(qJ(2));
t271 = t288 * t281 - t290 * t284;
t286 = cos(pkin(6));
t299 = t271 * t286;
t265 = qJD(2) * t299;
t304 = t290 * t281 + t288 * t284;
t270 = t304 * qJD(2);
t282 = sin(pkin(10));
t285 = cos(pkin(10));
t246 = -t285 * t265 - t282 * t270;
t268 = t304 * t286;
t252 = t285 * t268 - t282 * t271;
t280 = pkin(12) + qJ(5);
t278 = sin(t280);
t283 = sin(pkin(6));
t320 = t283 * t285;
t310 = t278 * t320;
t279 = cos(t280);
t316 = qJD(5) * t279;
t218 = -qJD(5) * t310 + t246 * t278 + t252 * t316;
t240 = t252 * t278 + t279 * t320;
t238 = t240 ^ 2;
t267 = t304 * t283;
t259 = t267 * t278 - t286 * t279;
t257 = 0.1e1 / t259 ^ 2;
t232 = t238 * t257 + 0.1e1;
t230 = 0.1e1 / t232;
t260 = t267 * t279 + t286 * t278;
t266 = t271 * t283;
t264 = qJD(2) * t266;
t236 = t260 * qJD(5) - t264 * t278;
t256 = 0.1e1 / t259;
t325 = t240 * t257;
t202 = (-t218 * t256 + t236 * t325) * t230;
t233 = atan2(-t240, t259);
t228 = sin(t233);
t229 = cos(t233);
t307 = -t228 * t259 - t229 * t240;
t198 = t307 * t202 - t228 * t218 + t229 * t236;
t212 = -t228 * t240 + t229 * t259;
t209 = 0.1e1 / t212;
t210 = 0.1e1 / t212 ^ 2;
t338 = t198 * t209 * t210;
t305 = -t282 * t268 - t285 * t271;
t321 = t282 * t283;
t300 = -t278 * t305 + t279 * t321;
t337 = -0.2e1 * t300 * t338;
t251 = -t282 * t304 - t285 * t299;
t301 = -t251 * t256 - t266 * t325;
t336 = t278 * t301;
t326 = t236 * t256 * t257;
t335 = -0.2e1 * (t218 * t325 - t238 * t326) / t232 ^ 2;
t244 = t278 * t321 + t279 * t305;
t289 = cos(qJ(6));
t254 = t282 * t299 - t285 * t304;
t287 = sin(qJ(6));
t323 = t254 * t287;
t227 = t244 * t289 - t323;
t223 = 0.1e1 / t227;
t224 = 0.1e1 / t227 ^ 2;
t306 = t282 * t265 - t285 * t270;
t221 = t300 * qJD(5) + t279 * t306;
t269 = t271 * qJD(2);
t298 = t286 * t270;
t247 = t285 * t269 + t282 * t298;
t213 = t227 * qJD(6) + t221 * t287 + t247 * t289;
t322 = t254 * t289;
t226 = t244 * t287 + t322;
t222 = t226 ^ 2;
t217 = t222 * t224 + 0.1e1;
t330 = t224 * t226;
t315 = qJD(6) * t226;
t214 = t221 * t289 - t247 * t287 - t315;
t332 = t214 * t223 * t224;
t334 = (t213 * t330 - t222 * t332) / t217 ^ 2;
t333 = t210 * t300;
t331 = t223 * t287;
t329 = t226 * t289;
t328 = t228 * t300;
t327 = t229 * t300;
t324 = t254 * t278;
t239 = t300 ^ 2;
t208 = t239 * t210 + 0.1e1;
t220 = t244 * qJD(5) + t278 * t306;
t314 = 0.2e1 * (-t220 * t333 - t239 * t338) / t208 ^ 2;
t313 = -0.2e1 * t334;
t311 = t226 * t332;
t309 = -0.2e1 * t240 * t326;
t308 = qJD(6) * t254 * t279 - t306;
t303 = t224 * t329 - t331;
t242 = t252 * t279 - t310;
t302 = -t242 * t256 + t260 * t325;
t297 = -qJD(5) * t324 + qJD(6) * t305 + t247 * t279;
t263 = t283 * t270;
t245 = t282 * t269 - t285 * t298;
t237 = -t259 * qJD(5) - t264 * t279;
t235 = t279 * t322 + t287 * t305;
t234 = t279 * t323 - t289 * t305;
t219 = -t240 * qJD(5) + t246 * t279;
t215 = 0.1e1 / t217;
t206 = 0.1e1 / t208;
t204 = t230 * t336;
t203 = t302 * t230;
t200 = (-t228 * t251 - t229 * t266) * t278 + t307 * t204;
t199 = t307 * t203 - t228 * t242 + t229 * t260;
t196 = t302 * t335 + (t260 * t309 - t219 * t256 + (t218 * t260 + t236 * t242 + t237 * t240) * t257) * t230;
t195 = t335 * t336 + (t301 * t316 + (-t266 * t309 - t245 * t256 + (-t218 * t266 + t236 * t251 - t240 * t263) * t257) * t278) * t230;
t1 = [0, t195, 0, 0, t196, 0; 0 (-t200 * t333 - t209 * t324) * t314 + ((t247 * t278 + t254 * t316) * t209 + t200 * t337 + (-t200 * t220 - t324 * t198 + (-t266 * t316 - t195 * t240 - t204 * t218 - t263 * t278 + (-t204 * t259 - t251 * t278) * t202) * t327 + (-t251 * t316 - t195 * t259 - t204 * t236 - t245 * t278 + (t204 * t240 + t266 * t278) * t202) * t328) * t210) * t206, 0, 0 (-t199 * t333 - t209 * t244) * t314 + (t199 * t337 + t221 * t209 + (-t244 * t198 - t199 * t220 + (-t196 * t240 - t203 * t218 + t237 + (-t203 * t259 - t242) * t202) * t327 + (-t196 * t259 - t203 * t236 - t219 + (t203 * t240 - t260) * t202) * t328) * t210) * t206, 0; 0, 0.2e1 * (-t223 * t234 + t235 * t330) * t334 + (0.2e1 * t235 * t311 + t308 * t223 * t289 + t297 * t331 + (t308 * t226 * t287 - t235 * t213 - t234 * t214 - t297 * t329) * t224) * t215, 0, 0, -t303 * t300 * t313 + (t303 * t220 - ((-qJD(6) * t223 - 0.2e1 * t311) * t289 + (t213 * t289 + (t214 - t315) * t287) * t224) * t300) * t215, t313 + 0.2e1 * (t213 * t224 * t215 + (-t215 * t332 - t224 * t334) * t226) * t226;];
JaD_rot  = t1;
