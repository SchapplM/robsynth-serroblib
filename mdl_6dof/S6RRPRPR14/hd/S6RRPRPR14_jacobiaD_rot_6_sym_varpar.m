% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR14_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR14_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:45:26
% EndTime: 2019-02-26 21:45:28
% DurationCPUTime: 1.59s
% Computational Cost: add. (4522->147), mult. (13478->287), div. (726->12), fcn. (17045->13), ass. (0->124)
t256 = sin(qJ(4));
t254 = cos(pkin(6));
t261 = cos(qJ(2));
t262 = cos(qJ(1));
t303 = t262 * t261;
t250 = t254 * t303;
t257 = sin(qJ(2));
t258 = sin(qJ(1));
t307 = t258 * t257;
t282 = -t250 + t307;
t260 = cos(qJ(4));
t253 = sin(pkin(6));
t309 = t253 * t262;
t291 = t260 * t309;
t330 = -t282 * t256 + t291;
t227 = t330 ^ 2;
t310 = t253 * t261;
t273 = -t254 * t260 + t256 * t310;
t239 = 0.1e1 / t273 ^ 2;
t217 = t227 * t239 + 0.1e1;
t215 = 0.1e1 / t217;
t304 = t262 * t257;
t306 = t258 * t261;
t243 = t254 * t304 + t306;
t244 = t254 * t306 + t304;
t223 = t244 * qJD(1) + t243 * qJD(2);
t302 = qJD(1) * t253;
t289 = t258 * t302;
t329 = t256 * t309 + t282 * t260;
t202 = qJD(4) * t329 + t223 * t256 + t260 * t289;
t241 = -t254 * t256 - t260 * t310;
t301 = qJD(2) * t257;
t288 = t253 * t301;
t228 = t241 * qJD(4) + t256 * t288;
t238 = 0.1e1 / t273;
t315 = t330 * t239;
t276 = t202 * t238 - t228 * t315;
t184 = t276 * t215;
t218 = atan2(t330, -t273);
t213 = sin(t218);
t214 = cos(t218);
t279 = t213 * t273 + t214 * t330;
t179 = t279 * t184 - t213 * t202 + t214 * t228;
t196 = t213 * t330 - t214 * t273;
t194 = 0.1e1 / t196 ^ 2;
t337 = t194 * t179;
t336 = -0.2e1 * t330;
t193 = 0.1e1 / t196;
t335 = t193 * t337;
t311 = t253 * t258;
t231 = t244 * t256 + t260 * t311;
t287 = 0.2e1 * t231 * t335;
t281 = qJD(2) * t254 + qJD(1);
t300 = qJD(2) * t261;
t221 = -qJD(1) * t250 - t262 * t300 + t281 * t307;
t290 = t262 * t302;
t294 = t256 * t311;
t204 = -t221 * t256 - qJD(4) * t294 + (qJD(4) * t244 + t290) * t260;
t321 = t204 * t194;
t334 = -t321 + t287;
t230 = -t244 * t260 + t294;
t293 = t254 * t307;
t245 = -t293 + t303;
t255 = sin(qJ(6));
t259 = cos(qJ(6));
t278 = t230 * t259 - t245 * t255;
t333 = t278 * qJD(6);
t200 = -qJD(4) * t291 - t223 * t260 + (t282 * qJD(4) + t289) * t256;
t332 = t228 * t239;
t312 = t253 * t257;
t272 = -t238 * t243 + t312 * t315;
t331 = t256 * t272;
t212 = t230 * t255 + t245 * t259;
t206 = 0.1e1 / t212;
t207 = 0.1e1 / t212 ^ 2;
t226 = t231 ^ 2;
t192 = t226 * t194 + 0.1e1;
t327 = (-t226 * t335 + t231 * t321) / t192 ^ 2;
t203 = t231 * qJD(4) + t221 * t260 + t256 * t290;
t222 = -t243 * qJD(1) - t244 * qJD(2);
t188 = t212 * qJD(6) - t203 * t259 + t222 * t255;
t205 = t278 ^ 2;
t199 = t205 * t207 + 0.1e1;
t320 = t207 * t278;
t189 = t203 * t255 + t222 * t259 + t333;
t323 = t189 * t206 * t207;
t326 = (-t188 * t320 - t205 * t323) / t199 ^ 2;
t317 = t238 * t332;
t324 = (-t202 * t315 + t227 * t317) / t217 ^ 2;
t322 = t194 * t231;
t319 = t213 * t231;
t318 = t214 * t231;
t316 = t330 * t238;
t314 = t245 * t256;
t313 = t245 * t260;
t308 = t255 * t278;
t305 = t259 * t206;
t299 = qJD(4) * t260;
t298 = 0.2e1 * t327;
t297 = 0.2e1 * t326;
t296 = 0.2e1 * t324;
t286 = t238 * t296;
t285 = -0.2e1 * t278 * t323;
t284 = t317 * t336;
t280 = -qJD(6) * t313 + t221;
t277 = t243 * t255 + t259 * t329;
t210 = -t243 * t259 + t255 * t329;
t275 = -t207 * t308 + t305;
t274 = -t238 * t329 + t241 * t315;
t270 = -t213 + (t214 * t316 + t213) * t215;
t269 = qJD(4) * t314 + qJD(6) * t244 - t222 * t260;
t229 = t273 * qJD(4) + t260 * t288;
t224 = -qJD(1) * t293 - t258 * t301 + t281 * t303;
t220 = -t244 * t259 - t255 * t313;
t219 = -t244 * t255 + t259 * t313;
t197 = 0.1e1 / t199;
t190 = 0.1e1 / t192;
t187 = t215 * t331;
t186 = t274 * t215;
t181 = (-t213 * t243 + t214 * t312) * t256 - t279 * t187;
t180 = -t279 * t186 - t213 * t329 + t214 * t241;
t178 = t274 * t296 + (t241 * t284 - t200 * t238 + (t202 * t241 + t228 * t329 - t229 * t330) * t239) * t215;
t176 = t296 * t331 + (-t272 * t299 + (t284 * t312 + t224 * t238 + (t228 * t243 + (t202 * t257 - t300 * t330) * t253) * t239) * t256) * t215;
t1 = [-t231 * t286 + (t204 * t238 + t231 * t332) * t215, t176, 0, t178, 0, 0; t193 * t327 * t336 + (-t202 * t193 - t330 * t337 - (t270 * t204 + ((-t184 * t215 * t316 - 0.2e1 * t324) * t213 + (-t330 * t286 - t184 + (t184 - t276) * t215) * t214) * t231) * t322) * t190 + (t190 * t334 + t322 * t298) * t270 * t231 (t181 * t322 - t193 * t314) * t298 + ((t222 * t256 + t245 * t299) * t193 + t334 * t181 + (-t314 * t179 - (t176 * t330 + t187 * t202 + (t256 * t300 + t257 * t299) * t253 + (-t187 * t273 - t243 * t256) * t184) * t318 - (-t243 * t299 + t176 * t273 + t187 * t228 - t224 * t256 + (t187 * t330 - t256 * t312) * t184) * t319) * t194) * t190, 0 (t180 * t322 + t193 * t230) * t298 + (t180 * t287 - t203 * t193 + (t230 * t179 - t180 * t204 - (t178 * t330 + t186 * t202 + t229 + (-t186 * t273 - t329) * t184) * t318 - (t178 * t273 + t186 * t228 + t200 + (t186 * t330 - t241) * t184) * t319) * t194) * t190, 0, 0; (t206 * t277 - t210 * t320) * t297 + ((t210 * qJD(6) + t200 * t259 - t224 * t255) * t206 + t210 * t285 + (t277 * t189 + (t277 * qJD(6) - t200 * t255 - t224 * t259) * t278 - t210 * t188) * t207) * t197 (-t206 * t219 - t220 * t320) * t297 + (t220 * t285 + t280 * t206 * t255 - t269 * t305 + (t259 * t278 * t280 - t220 * t188 - t219 * t189 + t269 * t308) * t207) * t197, 0, t275 * t231 * t297 + (-t275 * t204 + ((qJD(6) * t206 + t285) * t255 + (-t188 * t255 + (t189 + t333) * t259) * t207) * t231) * t197, 0, -0.2e1 * t326 - 0.2e1 * (t188 * t207 * t197 - (-t197 * t323 - t207 * t326) * t278) * t278;];
JaD_rot  = t1;
