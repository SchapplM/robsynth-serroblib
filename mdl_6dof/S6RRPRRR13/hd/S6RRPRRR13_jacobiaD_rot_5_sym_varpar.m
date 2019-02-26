% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR13_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR13_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:01:00
% EndTime: 2019-02-26 22:01:01
% DurationCPUTime: 1.73s
% Computational Cost: add. (4522->152), mult. (13478->303), div. (726->12), fcn. (17045->13), ass. (0->127)
t255 = sin(qJ(2));
t256 = sin(qJ(1));
t259 = cos(qJ(2));
t260 = cos(qJ(1));
t338 = cos(pkin(6));
t290 = t260 * t338;
t243 = t255 * t256 - t259 * t290;
t254 = sin(qJ(4));
t258 = cos(qJ(4));
t252 = sin(pkin(6));
t318 = t252 * t260;
t279 = t243 * t258 + t254 * t318;
t226 = t279 ^ 2;
t319 = t252 * t259;
t276 = -t338 * t254 - t258 * t319;
t239 = 0.1e1 / t276 ^ 2;
t217 = t226 * t239 + 0.1e1;
t215 = 0.1e1 / t217;
t313 = qJD(1) * t256;
t297 = t252 * t313;
t275 = -t255 * t290 - t256 * t259;
t291 = t256 * t338;
t277 = t260 * t255 + t259 * t291;
t268 = t277 * qJD(1) - t275 * qJD(2);
t298 = t258 * t318;
t346 = -qJD(4) * t298 - t268 * t258;
t200 = (qJD(4) * t243 + t297) * t254 + t346;
t242 = -t254 * t319 + t338 * t258;
t312 = qJD(2) * t255;
t295 = t252 * t312;
t228 = t242 * qJD(4) - t258 * t295;
t238 = 0.1e1 / t276;
t324 = t279 * t239;
t283 = t200 * t238 - t228 * t324;
t184 = t283 * t215;
t218 = atan2(t279, -t276);
t213 = sin(t218);
t214 = cos(t218);
t285 = t213 * t276 + t214 * t279;
t179 = t285 * t184 - t213 * t200 + t214 * t228;
t196 = t213 * t279 - t214 * t276;
t194 = 0.1e1 / t196 ^ 2;
t347 = t179 * t194;
t310 = qJD(4) * t254;
t292 = t252 * t310;
t309 = qJD(4) * t258;
t201 = t243 * t309 + t268 * t254 + t258 * t297 + t260 * t292;
t320 = t252 * t258;
t230 = t277 * t254 + t256 * t320;
t286 = t255 * t291;
t315 = t260 * t259;
t245 = -t286 + t315;
t253 = sin(qJ(5));
t257 = cos(qJ(5));
t209 = t230 * t253 - t245 * t257;
t345 = 0.2e1 * t209;
t193 = 0.1e1 / t196;
t344 = t193 * t347;
t271 = (t338 * qJD(1) + qJD(2)) * t315 - qJD(2) * t286 - t255 * t313;
t296 = qJD(1) * t318;
t203 = t230 * qJD(4) + t254 * t296 - t271 * t258;
t273 = t277 * t258;
t229 = t256 * t252 * t254 - t273;
t289 = 0.2e1 * t229 * t344;
t343 = -t194 * t203 + t289;
t342 = t228 * t239;
t321 = t252 * t255;
t299 = t279 * t321;
t278 = t238 * t275 + t239 * t299;
t341 = t258 * t278;
t340 = -t245 * t254 * qJD(5) - t271;
t339 = -qJD(5) * t277 + t245 * t309;
t210 = t230 * t257 + t245 * t253;
t206 = 0.1e1 / t210;
t207 = 0.1e1 / t210 ^ 2;
t225 = t229 ^ 2;
t192 = t194 * t225 + 0.1e1;
t331 = t194 * t229;
t337 = (t203 * t331 - t225 * t344) / t192 ^ 2;
t204 = qJD(4) * t273 + t271 * t254 - t256 * t292 + t258 * t296;
t223 = t275 * qJD(1) - t277 * qJD(2);
t188 = t210 * qJD(5) + t204 * t253 - t223 * t257;
t205 = t209 ^ 2;
t199 = t205 * t207 + 0.1e1;
t329 = t207 * t209;
t308 = qJD(5) * t209;
t189 = t204 * t257 + t223 * t253 - t308;
t333 = t189 * t206 * t207;
t336 = (t188 * t329 - t205 * t333) / t199 ^ 2;
t326 = t238 * t342;
t334 = (-t200 * t324 + t226 * t326) / t217 ^ 2;
t197 = 0.1e1 / t199;
t330 = t197 * t207;
t328 = t213 * t229;
t327 = t214 * t229;
t325 = t279 * t238;
t323 = t279 * t242;
t322 = t245 * t258;
t317 = t253 * t254;
t316 = t254 * t257;
t311 = qJD(2) * t259;
t307 = 0.2e1 * t337;
t306 = -0.2e1 * t336;
t305 = -0.2e1 * t334;
t304 = 0.2e1 * t334;
t302 = t207 * t336;
t301 = t188 * t330;
t300 = t209 * t333;
t288 = t238 * t304;
t287 = 0.2e1 * t300;
t280 = -t243 * t254 + t298;
t212 = t253 * t275 + t257 * t280;
t211 = t253 * t280 - t257 * t275;
t282 = -t206 * t253 + t257 * t329;
t281 = t238 * t280 + t239 * t323;
t274 = -t213 + (t214 * t325 + t213) * t215;
t227 = t276 * qJD(4) + t254 * t295;
t224 = -qJD(1) * t286 - t256 * t312 + (qJD(2) * t338 + qJD(1)) * t315;
t220 = t245 * t316 - t277 * t253;
t190 = 0.1e1 / t192;
t187 = t215 * t341;
t186 = t281 * t215;
t181 = (-t213 * t275 - t214 * t321) * t258 + t285 * t187;
t180 = -t285 * t186 + t213 * t280 + t214 * t242;
t178 = t281 * t304 + (-0.2e1 * t323 * t326 + t201 * t238 + (t200 * t242 - t227 * t279 - t228 * t280) * t239) * t215;
t176 = t305 * t341 + (-t278 * t310 + (0.2e1 * t299 * t326 - t224 * t238 + (t228 * t275 + (-t200 * t255 + t279 * t311) * t252) * t239) * t258) * t215;
t1 = [-t229 * t288 + (t203 * t238 + t229 * t342) * t215, t176, 0, t178, 0, 0; -0.2e1 * t279 * t193 * t337 + ((-t243 * t310 - t254 * t297 - t346) * t193 - t279 * t347 - (t274 * t203 + ((-t184 * t215 * t325 + t305) * t213 + (-t279 * t288 - t184 + (t184 - t283) * t215) * t214) * t229) * t331) * t190 + (t343 * t190 + t331 * t307) * t274 * t229 (t181 * t331 + t193 * t322) * t307 + ((-t223 * t258 + t245 * t310) * t193 + t343 * t181 + (t322 * t179 - (t176 * t279 - t187 * t200 + (t255 * t310 - t258 * t311) * t252 + (t187 * t276 - t258 * t275) * t184) * t327 - (t275 * t310 + t176 * t276 - t187 * t228 + t224 * t258 + (-t187 * t279 + t255 * t320) * t184) * t328) * t194) * t190, 0 (t180 * t331 - t193 * t230) * t307 + (t180 * t289 + t204 * t193 + (-t230 * t179 - t180 * t203 - (t178 * t279 + t186 * t200 + t227 + (-t186 * t276 + t280) * t184) * t327 - (t178 * t276 + t186 * t228 - t201 + (t186 * t279 - t242) * t184) * t328) * t194) * t190, 0, 0; 0.2e1 * (-t206 * t211 + t212 * t329) * t336 + ((t212 * qJD(5) - t201 * t253 + t224 * t257) * t206 + t212 * t287 + (-t211 * t189 - (-t211 * qJD(5) - t201 * t257 - t224 * t253) * t209 - t212 * t188) * t207) * t197 (t302 * t345 - t301) * t220 + (-t189 * t330 + t206 * t306) * (t245 * t317 + t277 * t257) + (t220 * t287 + (t317 * t206 - t316 * t329) * t223 + (-t206 * t340 - t329 * t339) * t257 + (t206 * t339 - t329 * t340) * t253) * t197, 0, t282 * t229 * t306 + (t282 * t203 + ((-qJD(5) * t206 - 0.2e1 * t300) * t257 + (t188 * t257 + (t189 - t308) * t253) * t207) * t229) * t197, t306 + (t301 + (-t197 * t333 - t302) * t209) * t345, 0;];
JaD_rot  = t1;
