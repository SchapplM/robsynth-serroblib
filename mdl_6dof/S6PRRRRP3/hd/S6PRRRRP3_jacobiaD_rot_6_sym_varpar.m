% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRRP3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:16:21
% EndTime: 2019-02-26 20:16:22
% DurationCPUTime: 1.02s
% Computational Cost: add. (3659->111), mult. (9671->227), div. (577->12), fcn. (12382->13), ass. (0->107)
t248 = sin(pkin(11));
t250 = cos(pkin(11));
t253 = sin(qJ(2));
t251 = cos(pkin(6));
t255 = cos(qJ(2));
t282 = t251 * t255;
t234 = -t248 * t253 + t250 * t282;
t227 = t234 * qJD(2);
t283 = t251 * t253;
t235 = t248 * t255 + t250 * t283;
t252 = sin(qJ(3));
t249 = sin(pkin(6));
t286 = t249 * t252;
t273 = t250 * t286;
t254 = cos(qJ(3));
t279 = qJD(3) * t254;
t205 = -qJD(3) * t273 + t227 * t252 + t235 * t279;
t285 = t249 * t254;
t219 = t235 * t252 + t250 * t285;
t217 = t219 ^ 2;
t238 = -t251 * t254 + t253 * t286;
t232 = 0.1e1 / t238 ^ 2;
t213 = t217 * t232 + 0.1e1;
t211 = 0.1e1 / t213;
t239 = t251 * t252 + t253 * t285;
t280 = qJD(2) * t255;
t272 = t249 * t280;
t224 = t239 * qJD(3) + t252 * t272;
t231 = 0.1e1 / t238;
t290 = t219 * t232;
t183 = (-t205 * t231 + t224 * t290) * t211;
t214 = atan2(-t219, t238);
t209 = sin(t214);
t210 = cos(t214);
t267 = -t209 * t238 - t210 * t219;
t179 = t267 * t183 - t205 * t209 + t210 * t224;
t195 = -t209 * t219 + t210 * t238;
t192 = 0.1e1 / t195;
t193 = 0.1e1 / t195 ^ 2;
t303 = t179 * t192 * t193;
t274 = t248 * t283;
t237 = t250 * t255 - t274;
t264 = -t237 * t252 + t248 * t285;
t302 = -0.2e1 * t264 * t303;
t284 = t249 * t255;
t263 = -t231 * t234 + t284 * t290;
t301 = t252 * t263;
t289 = t224 * t231 * t232;
t300 = -0.2e1 * (t205 * t290 - t217 * t289) / t213 ^ 2;
t223 = t237 * t254 + t248 * t286;
t236 = t248 * t282 + t250 * t253;
t247 = qJ(4) + qJ(5);
t244 = sin(t247);
t245 = cos(t247);
t204 = t223 * t245 + t236 * t244;
t200 = 0.1e1 / t204;
t201 = 0.1e1 / t204 ^ 2;
t230 = -qJD(2) * t274 + t250 * t280;
t246 = qJD(4) + qJD(5);
t269 = t223 * t246 - t230;
t229 = t236 * qJD(2);
t208 = t264 * qJD(3) - t229 * t254;
t270 = t236 * t246 + t208;
t190 = t270 * t244 + t269 * t245;
t203 = t223 * t244 - t236 * t245;
t199 = t203 ^ 2;
t198 = t199 * t201 + 0.1e1;
t295 = t201 * t203;
t191 = -t269 * t244 + t270 * t245;
t298 = t191 * t200 * t201;
t299 = (t190 * t295 - t199 * t298) / t198 ^ 2;
t297 = t193 * t264;
t296 = t200 * t244;
t294 = t203 * t245;
t207 = t223 * qJD(3) - t229 * t252;
t293 = t207 * t193;
t292 = t209 * t264;
t291 = t210 * t264;
t288 = t236 * t252;
t287 = t236 * t254;
t281 = qJD(2) * t253;
t218 = t264 ^ 2;
t189 = t193 * t218 + 0.1e1;
t278 = 0.2e1 * (-t218 * t303 - t264 * t293) / t189 ^ 2;
t277 = -0.2e1 * t299;
t275 = t203 * t298;
t271 = -0.2e1 * t219 * t289;
t268 = t246 * t287 - t229;
t266 = t201 * t294 - t296;
t221 = t235 * t254 - t273;
t265 = -t221 * t231 + t239 * t290;
t262 = qJD(3) * t288 - t230 * t254 + t237 * t246;
t228 = t235 * qJD(2);
t225 = -t238 * qJD(3) + t254 * t272;
t216 = t237 * t244 - t245 * t287;
t215 = -t237 * t245 - t244 * t287;
t206 = -t219 * qJD(3) + t227 * t254;
t196 = 0.1e1 / t198;
t186 = 0.1e1 / t189;
t185 = t211 * t301;
t184 = t265 * t211;
t181 = (-t209 * t234 + t210 * t284) * t252 + t267 * t185;
t180 = t267 * t184 - t209 * t221 + t210 * t239;
t178 = t265 * t300 + (t239 * t271 - t206 * t231 + (t205 * t239 + t219 * t225 + t221 * t224) * t232) * t211;
t176 = t300 * t301 + (t263 * t279 + (t271 * t284 + t228 * t231 + (t224 * t234 + (t205 * t255 - t219 * t281) * t249) * t232) * t252) * t211;
t175 = t277 + 0.2e1 * (t190 * t201 * t196 + (-t196 * t298 - t201 * t299) * t203) * t203;
t1 = [0, t176, t178, 0, 0, 0; 0 (-t181 * t297 + t192 * t288) * t278 + ((-t230 * t252 - t236 * t279) * t192 + (-t293 + t302) * t181 + (t288 * t179 + (-t176 * t219 - t185 * t205 + (-t252 * t281 + t255 * t279) * t249 + (-t185 * t238 - t234 * t252) * t183) * t291 + (-t234 * t279 - t176 * t238 - t185 * t224 + t228 * t252 + (t185 * t219 - t252 * t284) * t183) * t292) * t193) * t186 (-t180 * t297 - t192 * t223) * t278 + (t180 * t302 + t208 * t192 + (-t223 * t179 - t180 * t207 + (-t178 * t219 - t184 * t205 + t225 + (-t184 * t238 - t221) * t183) * t291 + (-t178 * t238 - t184 * t224 - t206 + (t184 * t219 - t239) * t183) * t292) * t193) * t186, 0, 0, 0; 0, 0.2e1 * (-t200 * t215 + t216 * t295) * t299 + (0.2e1 * t216 * t275 - t268 * t200 * t245 + t262 * t296 + (-t268 * t203 * t244 - t216 * t190 - t215 * t191 - t262 * t294) * t201) * t196, -t266 * t264 * t277 + (t266 * t207 - ((-t200 * t246 - 0.2e1 * t275) * t245 + (t190 * t245 + (-t203 * t246 + t191) * t244) * t201) * t264) * t196, t175, t175, 0;];
JaD_rot  = t1;
