% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR11_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR11_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_jacobiaD_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:36:22
% EndTime: 2019-02-26 22:36:23
% DurationCPUTime: 1.44s
% Computational Cost: add. (4522->148), mult. (13478->297), div. (726->12), fcn. (17045->13), ass. (0->126)
t256 = cos(pkin(6));
t259 = sin(qJ(2));
t331 = sin(qJ(1));
t293 = t331 * t259;
t283 = t256 * t293;
t287 = t331 * qJD(2);
t262 = cos(qJ(2));
t263 = cos(qJ(1));
t309 = t263 * t262;
t255 = sin(pkin(6));
t313 = t255 * t263;
t336 = -qJD(1) * t283 - t259 * t287 + (qJD(2) * t256 + qJD(1)) * t309 - qJD(3) * t313;
t258 = sin(qJ(3));
t261 = cos(qJ(3));
t292 = t331 * t262;
t310 = t263 * t259;
t275 = -t256 * t310 - t292;
t228 = -t258 * t275 + t261 * t313;
t315 = t255 * t259;
t239 = -t256 * t261 + t258 * t315;
t217 = atan2(-t228, t239);
t212 = sin(t217);
t213 = cos(t217);
t195 = -t212 * t228 + t213 * t239;
t193 = 0.1e1 / t195 ^ 2;
t244 = -t283 + t309;
t294 = t255 * t331;
t274 = -t244 * t258 + t261 * t294;
t225 = t274 ^ 2;
t191 = t225 * t193 + 0.1e1;
t273 = -t256 * t292 - t310;
t221 = qJD(1) * t275 + qJD(2) * t273;
t234 = t244 * t261 + t258 * t294;
t291 = qJD(1) * t313;
t199 = qJD(3) * t234 + t221 * t258 - t261 * t291;
t324 = t199 * t193;
t224 = t228 ^ 2;
t237 = 0.1e1 / t239 ^ 2;
t216 = t224 * t237 + 0.1e1;
t214 = 0.1e1 / t216;
t288 = t331 * qJD(1);
t282 = t255 * t288;
t306 = qJD(3) * t261;
t201 = t336 * t258 - t261 * t282 - t275 * t306;
t240 = t256 * t258 + t261 * t315;
t307 = qJD(2) * t262;
t290 = t255 * t307;
t226 = qJD(3) * t240 + t258 * t290;
t236 = 0.1e1 / t239;
t318 = t228 * t237;
t279 = -t201 * t236 + t226 * t318;
t183 = t279 * t214;
t280 = -t212 * t239 - t213 * t228;
t178 = t183 * t280 - t212 * t201 + t213 * t226;
t192 = 0.1e1 / t195;
t194 = t192 * t193;
t329 = t178 * t194;
t304 = 0.2e1 * (-t225 * t329 - t274 * t324) / t191 ^ 2;
t335 = t226 * t237;
t295 = t256 * t309;
t241 = -t293 + t295;
t314 = t255 * t262;
t276 = -t236 * t241 + t314 * t318;
t334 = t258 * t276;
t202 = t258 * (qJD(3) * t275 + t282) + t336 * t261;
t257 = sin(qJ(4));
t260 = cos(qJ(4));
t211 = t234 * t260 - t257 * t273;
t205 = 0.1e1 / t211;
t206 = 0.1e1 / t211 ^ 2;
t333 = -0.2e1 * t228;
t332 = -0.2e1 * t274;
t200 = qJD(3) * t274 + t221 * t261 + t258 * t291;
t220 = -qJD(1) * t295 - t263 * t307 + (t256 * t287 + t288) * t259;
t187 = qJD(4) * t211 + t200 * t257 + t220 * t260;
t210 = t234 * t257 + t260 * t273;
t204 = t210 ^ 2;
t198 = t204 * t206 + 0.1e1;
t323 = t206 * t210;
t305 = qJD(4) * t210;
t188 = t200 * t260 - t220 * t257 - t305;
t326 = t188 * t205 * t206;
t328 = (t187 * t323 - t204 * t326) / t198 ^ 2;
t320 = t236 * t335;
t327 = (t201 * t318 - t224 * t320) / t216 ^ 2;
t325 = t193 * t274;
t322 = t212 * t274;
t321 = t213 * t274;
t319 = t228 * t236;
t317 = t273 * t258;
t316 = t273 * t261;
t312 = t257 * t205;
t311 = t260 * t210;
t308 = qJD(2) * t259;
t303 = -0.2e1 * t328;
t302 = 0.2e1 * t328;
t301 = -0.2e1 * t327;
t300 = t194 * t332;
t299 = t236 * t327;
t298 = t193 * t322;
t297 = t193 * t321;
t296 = t210 * t326;
t286 = 0.2e1 * t296;
t285 = t320 * t333;
t230 = -t258 * t313 - t261 * t275;
t281 = -qJD(4) * t316 + t221;
t209 = -t230 * t260 + t241 * t257;
t208 = -t230 * t257 - t241 * t260;
t278 = t206 * t311 - t312;
t277 = -t230 * t236 + t240 * t318;
t271 = -t212 + (t213 * t319 + t212) * t214;
t270 = -qJD(3) * t317 + qJD(4) * t244 + t220 * t261;
t227 = -qJD(3) * t239 + t261 * t290;
t222 = qJD(1) * t273 + qJD(2) * t275;
t219 = t244 * t257 + t260 * t316;
t218 = -t244 * t260 + t257 * t316;
t196 = 0.1e1 / t198;
t189 = 0.1e1 / t191;
t186 = t214 * t334;
t185 = t277 * t214;
t182 = t271 * t274;
t180 = (-t212 * t241 + t213 * t314) * t258 + t280 * t186;
t179 = t185 * t280 - t212 * t230 + t213 * t240;
t177 = t277 * t301 + (t240 * t285 - t202 * t236 + (t201 * t240 + t226 * t230 + t227 * t228) * t237) * t214;
t175 = t301 * t334 + (t276 * t306 + (t285 * t314 - t222 * t236 + (t226 * t241 + (t201 * t262 - t228 * t308) * t255) * t237) * t258) * t214;
t1 = [t299 * t332 + (-t199 * t236 - t274 * t335) * t214, t175, t177, 0, 0, 0; t228 * t192 * t304 + (-t201 * t192 + (t178 * t228 + t182 * t199) * t193) * t189 - (-t182 * t193 * t304 + (-0.2e1 * t182 * t329 + (-t183 * t214 * t319 + t301) * t298 + (t299 * t333 - t183 + (t183 - t279) * t214) * t297 - t271 * t324) * t189) * t274 (-t180 * t325 - t192 * t317) * t304 + (-t180 * t324 + (t220 * t258 + t273 * t306) * t192 + (t180 * t300 - t193 * t317) * t178 + (-t175 * t228 - t186 * t201 + (-t258 * t308 + t262 * t306) * t255 + (-t186 * t239 - t241 * t258) * t183) * t297 + (-t241 * t306 - t175 * t239 - t186 * t226 - t222 * t258 + (t186 * t228 - t258 * t314) * t183) * t298) * t189 (-t179 * t325 - t192 * t234) * t304 + (t179 * t178 * t300 + t200 * t192 + (-t234 * t178 - t179 * t199 + (-t177 * t228 - t185 * t201 + t227 + (-t185 * t239 - t230) * t183) * t321 + (-t177 * t239 - t185 * t226 - t202 + (t185 * t228 - t240) * t183) * t322) * t193) * t189, 0, 0, 0; (-t205 * t208 + t209 * t323) * t302 + ((qJD(4) * t209 - t202 * t257 - t222 * t260) * t205 + t209 * t286 + (-t208 * t188 - (-qJD(4) * t208 - t202 * t260 + t222 * t257) * t210 - t209 * t187) * t206) * t196 (-t205 * t218 + t219 * t323) * t302 + (t219 * t286 - t281 * t205 * t260 + t270 * t312 + (-t210 * t257 * t281 - t219 * t187 - t218 * t188 - t270 * t311) * t206) * t196, -t278 * t274 * t303 + (t278 * t199 - ((-qJD(4) * t205 - 0.2e1 * t296) * t260 + (t187 * t260 + (t188 - t305) * t257) * t206) * t274) * t196, t303 + 0.2e1 * (t187 * t206 * t196 + (-t196 * t326 - t206 * t328) * t210) * t210, 0, 0;];
JaD_rot  = t1;
