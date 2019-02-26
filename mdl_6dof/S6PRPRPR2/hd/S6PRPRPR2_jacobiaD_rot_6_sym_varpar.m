% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRPR2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:47:04
% EndTime: 2019-02-26 19:47:05
% DurationCPUTime: 1.23s
% Computational Cost: add. (5953->115), mult. (16325->235), div. (559->12), fcn. (21250->15), ass. (0->111)
t282 = sin(pkin(11));
t285 = cos(pkin(11));
t289 = sin(qJ(2));
t291 = cos(qJ(2));
t272 = t289 * t282 - t291 * t285;
t287 = cos(pkin(6));
t300 = t272 * t287;
t266 = qJD(2) * t300;
t305 = t291 * t282 + t289 * t285;
t271 = t305 * qJD(2);
t283 = sin(pkin(10));
t286 = cos(pkin(10));
t247 = -t286 * t266 - t283 * t271;
t269 = t305 * t287;
t253 = t286 * t269 - t283 * t272;
t288 = sin(qJ(4));
t284 = sin(pkin(6));
t322 = t284 * t288;
t311 = t286 * t322;
t290 = cos(qJ(4));
t317 = qJD(4) * t290;
t225 = -qJD(4) * t311 + t247 * t288 + t253 * t317;
t321 = t284 * t290;
t241 = t253 * t288 + t286 * t321;
t239 = t241 ^ 2;
t268 = t305 * t284;
t260 = t268 * t288 - t287 * t290;
t258 = 0.1e1 / t260 ^ 2;
t235 = t239 * t258 + 0.1e1;
t233 = 0.1e1 / t235;
t261 = t268 * t290 + t287 * t288;
t267 = t272 * t284;
t265 = qJD(2) * t267;
t237 = t261 * qJD(4) - t265 * t288;
t257 = 0.1e1 / t260;
t325 = t241 * t258;
t203 = (-t225 * t257 + t237 * t325) * t233;
t236 = atan2(-t241, t260);
t231 = sin(t236);
t232 = cos(t236);
t308 = -t231 * t260 - t232 * t241;
t199 = t203 * t308 - t231 * t225 + t232 * t237;
t215 = -t231 * t241 + t232 * t260;
t212 = 0.1e1 / t215;
t213 = 0.1e1 / t215 ^ 2;
t339 = t199 * t212 * t213;
t306 = -t283 * t269 - t286 * t272;
t301 = t283 * t321 - t288 * t306;
t338 = -0.2e1 * t301 * t339;
t252 = -t283 * t305 - t286 * t300;
t302 = -t252 * t257 - t267 * t325;
t337 = t288 * t302;
t326 = t237 * t257 * t258;
t336 = -0.2e1 * (t225 * t325 - t239 * t326) / t235 ^ 2;
t245 = t283 * t322 + t290 * t306;
t255 = t283 * t300 - t286 * t305;
t281 = pkin(12) + qJ(6);
t279 = sin(t281);
t280 = cos(t281);
t224 = t245 * t280 - t255 * t279;
t220 = 0.1e1 / t224;
t221 = 0.1e1 / t224 ^ 2;
t307 = t283 * t266 - t286 * t271;
t228 = qJD(4) * t301 + t290 * t307;
t270 = t272 * qJD(2);
t299 = t287 * t271;
t248 = t286 * t270 + t283 * t299;
t210 = qJD(6) * t224 + t228 * t279 + t248 * t280;
t223 = t245 * t279 + t255 * t280;
t219 = t223 ^ 2;
t218 = t219 * t221 + 0.1e1;
t331 = t221 * t223;
t316 = qJD(6) * t223;
t211 = t228 * t280 - t248 * t279 - t316;
t334 = t211 * t220 * t221;
t335 = (t210 * t331 - t219 * t334) / t218 ^ 2;
t333 = t213 * t301;
t332 = t220 * t279;
t330 = t223 * t280;
t227 = qJD(4) * t245 + t288 * t307;
t329 = t227 * t213;
t328 = t231 * t301;
t327 = t232 * t301;
t324 = t255 * t288;
t323 = t255 * t290;
t240 = t301 ^ 2;
t209 = t240 * t213 + 0.1e1;
t315 = 0.2e1 * (-t240 * t339 - t301 * t329) / t209 ^ 2;
t314 = -0.2e1 * t335;
t312 = t223 * t334;
t310 = -0.2e1 * t241 * t326;
t309 = qJD(6) * t323 - t307;
t304 = t221 * t330 - t332;
t243 = t253 * t290 - t311;
t303 = -t243 * t257 + t261 * t325;
t298 = -qJD(4) * t324 + qJD(6) * t306 + t248 * t290;
t264 = t284 * t271;
t246 = t283 * t270 - t286 * t299;
t238 = -t260 * qJD(4) - t265 * t290;
t230 = t279 * t306 + t280 * t323;
t229 = t279 * t323 - t280 * t306;
t226 = -qJD(4) * t241 + t247 * t290;
t216 = 0.1e1 / t218;
t207 = 0.1e1 / t209;
t205 = t233 * t337;
t204 = t303 * t233;
t201 = (-t231 * t252 - t232 * t267) * t288 + t308 * t205;
t200 = t204 * t308 - t231 * t243 + t232 * t261;
t198 = t303 * t336 + (t261 * t310 - t226 * t257 + (t225 * t261 + t237 * t243 + t238 * t241) * t258) * t233;
t196 = t336 * t337 + (t302 * t317 + (-t267 * t310 - t246 * t257 + (-t225 * t267 + t237 * t252 - t241 * t264) * t258) * t288) * t233;
t1 = [0, t196, 0, t198, 0, 0; 0 (-t201 * t333 - t212 * t324) * t315 + ((t248 * t288 + t255 * t317) * t212 + (-t329 + t338) * t201 + (-t324 * t199 + (-t267 * t317 - t196 * t241 - t205 * t225 - t264 * t288 + (-t205 * t260 - t252 * t288) * t203) * t327 + (-t252 * t317 - t196 * t260 - t205 * t237 - t246 * t288 + (t205 * t241 + t267 * t288) * t203) * t328) * t213) * t207, 0 (-t200 * t333 - t212 * t245) * t315 + (t200 * t338 + t228 * t212 + (-t245 * t199 - t200 * t227 + (-t198 * t241 - t204 * t225 + t238 + (-t204 * t260 - t243) * t203) * t327 + (-t198 * t260 - t204 * t237 - t226 + (t204 * t241 - t261) * t203) * t328) * t213) * t207, 0, 0; 0, 0.2e1 * (-t220 * t229 + t230 * t331) * t335 + (0.2e1 * t230 * t312 + t309 * t220 * t280 + t298 * t332 + (t223 * t279 * t309 - t230 * t210 - t229 * t211 - t298 * t330) * t221) * t216, 0, -t304 * t301 * t314 + (t304 * t227 - ((-qJD(6) * t220 - 0.2e1 * t312) * t280 + (t210 * t280 + (t211 - t316) * t279) * t221) * t301) * t216, 0, t314 + 0.2e1 * (t210 * t221 * t216 + (-t216 * t334 - t221 * t335) * t223) * t223;];
JaD_rot  = t1;
