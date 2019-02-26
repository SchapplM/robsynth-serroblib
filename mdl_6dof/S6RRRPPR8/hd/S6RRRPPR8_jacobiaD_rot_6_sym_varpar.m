% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPPR8_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR8_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:07:40
% EndTime: 2019-02-26 22:07:42
% DurationCPUTime: 1.37s
% Computational Cost: add. (4522->149), mult. (13478->297), div. (726->12), fcn. (17045->13), ass. (0->124)
t258 = cos(pkin(6));
t262 = sin(qJ(1));
t265 = cos(qJ(2));
t328 = cos(qJ(1));
t287 = t328 * qJD(1);
t257 = sin(pkin(6));
t292 = t257 * t328;
t261 = sin(qJ(2));
t307 = t262 * t261;
t293 = t258 * t307;
t304 = qJD(2) * t261;
t333 = -qJD(1) * t293 - t262 * t304 + (t328 * qJD(2) * t258 + t287) * t265 - qJD(3) * t292;
t291 = t328 * t261;
t306 = t262 * t265;
t245 = t258 * t291 + t306;
t260 = sin(qJ(3));
t264 = cos(qJ(3));
t230 = t245 * t264 - t260 * t292;
t310 = t257 * t264;
t243 = t258 * t260 + t261 * t310;
t218 = atan2(-t230, t243);
t213 = sin(t218);
t214 = cos(t218);
t196 = -t213 * t230 + t214 * t243;
t194 = 0.1e1 / t196 ^ 2;
t290 = t328 * t265;
t275 = -t290 + t293;
t311 = t257 * t262;
t235 = t260 * t311 - t264 * t275;
t226 = t235 ^ 2;
t192 = t226 * t194 + 0.1e1;
t274 = t258 * t306 + t291;
t222 = t245 * qJD(1) + t274 * qJD(2);
t302 = qJD(3) * t264;
t303 = qJD(3) * t260;
t201 = -t222 * t264 + t275 * t303 + (t260 * t287 + t262 * t302) * t257;
t321 = t201 * t194;
t225 = t230 ^ 2;
t240 = 0.1e1 / t243 ^ 2;
t217 = t225 * t240 + 0.1e1;
t215 = 0.1e1 / t217;
t284 = t333 * t264;
t289 = qJD(1) * t311;
t203 = -t245 * t303 + t260 * t289 + t284;
t242 = -t257 * t261 * t260 + t258 * t264;
t309 = t257 * t265;
t288 = qJD(2) * t309;
t228 = t242 * qJD(3) + t264 * t288;
t239 = 0.1e1 / t243;
t315 = t230 * t240;
t279 = -t203 * t239 + t228 * t315;
t184 = t279 * t215;
t280 = -t213 * t243 - t214 * t230;
t179 = t280 * t184 - t213 * t203 + t214 * t228;
t193 = 0.1e1 / t196;
t195 = t193 * t194;
t326 = t179 * t195;
t300 = 0.2e1 * (-t226 * t326 + t235 * t321) / t192 ^ 2;
t332 = t228 * t240;
t244 = -t258 * t290 + t307;
t276 = t239 * t244 + t309 * t315;
t331 = t264 * t276;
t234 = -t260 * t275 - t262 * t310;
t259 = sin(qJ(6));
t263 = cos(qJ(6));
t212 = t234 * t263 - t259 * t274;
t206 = 0.1e1 / t212;
t207 = 0.1e1 / t212 ^ 2;
t330 = -0.2e1 * t230;
t329 = 0.2e1 * t235;
t283 = t264 * t292;
t200 = -qJD(1) * t283 + t235 * qJD(3) - t222 * t260;
t221 = t244 * qJD(1) + t275 * qJD(2);
t188 = t212 * qJD(6) + t200 * t259 - t221 * t263;
t313 = t274 * t263;
t211 = t234 * t259 + t313;
t205 = t211 ^ 2;
t199 = t205 * t207 + 0.1e1;
t320 = t207 * t211;
t301 = qJD(6) * t211;
t189 = t200 * t263 + t221 * t259 - t301;
t323 = t189 * t206 * t207;
t325 = (t188 * t320 - t205 * t323) / t199 ^ 2;
t317 = t239 * t332;
t324 = (t203 * t315 - t225 * t317) / t217 ^ 2;
t322 = t194 * t235;
t319 = t213 * t235;
t318 = t214 * t235;
t316 = t230 * t239;
t314 = t274 * t260;
t312 = t274 * t264;
t308 = t259 * t206;
t305 = t263 * t211;
t299 = 0.2e1 * t325;
t298 = -0.2e1 * t324;
t297 = t195 * t329;
t296 = t239 * t324;
t295 = t194 * t319;
t294 = t194 * t318;
t286 = 0.2e1 * t211 * t323;
t285 = t317 * t330;
t281 = -qJD(6) * t314 - t222;
t229 = t245 * t260 + t283;
t210 = -t229 * t263 + t244 * t259;
t209 = -t229 * t259 - t244 * t263;
t278 = t207 * t305 - t308;
t277 = t229 * t239 + t242 * t315;
t273 = -t213 + (t214 * t316 + t213) * t215;
t202 = t245 * t302 + t333 * t260 - t264 * t289;
t272 = qJD(6) * t275 + t221 * t260 - t274 * t302;
t227 = -t243 * qJD(3) - t260 * t288;
t223 = t274 * qJD(1) + t245 * qJD(2);
t220 = t259 * t275 - t260 * t313;
t219 = -t259 * t314 - t263 * t275;
t197 = 0.1e1 / t199;
t190 = 0.1e1 / t192;
t187 = t215 * t331;
t186 = t277 * t215;
t183 = t273 * t235;
t181 = (t213 * t244 + t214 * t309) * t264 + t280 * t187;
t180 = t280 * t186 + t213 * t229 + t214 * t242;
t178 = t277 * t298 + (t242 * t285 + t202 * t239 + (t203 * t242 + t227 * t230 - t228 * t229) * t240) * t215;
t176 = t298 * t331 + (-t276 * t303 + (t285 * t309 + t223 * t239 + (-t228 * t244 + (t203 * t265 - t230 * t304) * t257) * t240) * t264) * t215;
t1 = [t296 * t329 + (-t201 * t239 + t235 * t332) * t215, t176, t178, 0, 0, 0; t230 * t193 * t300 + (((qJD(3) * t245 - t289) * t260 - t284) * t193 + (t179 * t230 - t183 * t201) * t194) * t190 + (t183 * t194 * t300 + (0.2e1 * t183 * t326 - (-t184 * t215 * t316 + t298) * t295 - (t296 * t330 - t184 + (t184 - t279) * t215) * t294 - t273 * t321) * t190) * t235 (t181 * t322 + t193 * t312) * t300 + (-t181 * t321 + (t221 * t264 + t274 * t303) * t193 + (t181 * t297 + t194 * t312) * t179 - (-t176 * t230 - t187 * t203 + (-t264 * t304 - t265 * t303) * t257 + (-t187 * t243 + t244 * t264) * t184) * t294 - (-t244 * t303 - t176 * t243 - t187 * t228 + t223 * t264 + (t187 * t230 - t264 * t309) * t184) * t295) * t190 (t180 * t322 + t193 * t234) * t300 + (t180 * t179 * t297 - t200 * t193 + (t234 * t179 - t180 * t201 - (-t178 * t230 - t186 * t203 + t227 + (-t186 * t243 + t229) * t184) * t318 - (-t178 * t243 - t186 * t228 + t202 + (t186 * t230 - t242) * t184) * t319) * t194) * t190, 0, 0, 0; (-t206 * t209 + t210 * t320) * t299 + ((t210 * qJD(6) - t202 * t259 - t223 * t263) * t206 + t210 * t286 + (-t209 * t189 - (-t209 * qJD(6) - t202 * t263 + t223 * t259) * t211 - t210 * t188) * t207) * t197 (-t206 * t219 + t220 * t320) * t299 + (t220 * t286 + t281 * t206 * t263 + t272 * t308 + (t281 * t211 * t259 - t220 * t188 - t219 * t189 - t272 * t305) * t207) * t197, t278 * t235 * t299 + (-t278 * t201 + ((qJD(6) * t206 + t286) * t263 + (-t188 * t263 + (-t189 + t301) * t259) * t207) * t235) * t197, 0, 0, -0.2e1 * t325 + 0.2e1 * (t188 * t207 * t197 + (-t197 * t323 - t207 * t325) * t211) * t211;];
JaD_rot  = t1;
