% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:48
% EndTime: 2019-02-26 21:30:50
% DurationCPUTime: 1.40s
% Computational Cost: add. (4522->145), mult. (13478->293), div. (726->12), fcn. (17045->13), ass. (0->124)
t251 = cos(pkin(6));
t254 = sin(qJ(2));
t257 = cos(qJ(5));
t250 = sin(pkin(6));
t253 = sin(qJ(5));
t307 = t250 * t253;
t237 = t251 * t257 + t254 * t307;
t255 = sin(qJ(1));
t258 = cos(qJ(2));
t300 = t255 * t258;
t259 = cos(qJ(1));
t301 = t254 * t259;
t240 = t251 * t301 + t300;
t304 = t250 * t259;
t272 = -t240 * t253 + t257 * t304;
t215 = atan2(t272, t237);
t210 = sin(t215);
t211 = cos(t215);
t193 = t210 * t272 + t211 * t237;
t191 = 0.1e1 / t193 ^ 2;
t302 = t254 * t255;
t284 = t251 * t302;
t299 = t258 * t259;
t270 = t284 - t299;
t306 = t250 * t257;
t231 = -t253 * t270 + t255 * t306;
t223 = t231 ^ 2;
t189 = t191 * t223 + 0.1e1;
t271 = t251 * t300 + t301;
t219 = qJD(1) * t240 + qJD(2) * t271;
t298 = qJD(1) * t250;
t282 = t259 * t298;
t285 = t255 * t307;
t197 = -t219 * t253 - qJD(5) * t285 + (-qJD(5) * t270 + t282) * t257;
t318 = t191 * t231;
t222 = t272 ^ 2;
t235 = 0.1e1 / t237 ^ 2;
t214 = t222 * t235 + 0.1e1;
t212 = 0.1e1 / t214;
t297 = qJD(2) * t254;
t221 = -qJD(1) * t284 - t255 * t297 + (qJD(2) * t251 + qJD(1)) * t299;
t228 = t240 * t257 + t253 * t304;
t283 = t255 * t298;
t199 = qJD(5) * t228 + t221 * t253 + t257 * t283;
t238 = -t251 * t253 + t254 * t306;
t305 = t250 * t258;
t281 = qJD(2) * t305;
t224 = qJD(5) * t238 + t253 * t281;
t234 = 0.1e1 / t237;
t311 = t272 * t235;
t275 = -t199 * t234 - t224 * t311;
t181 = t275 * t212;
t276 = -t210 * t237 + t211 * t272;
t176 = t181 * t276 - t199 * t210 + t211 * t224;
t190 = 0.1e1 / t193;
t192 = t190 * t191;
t323 = t176 * t192;
t294 = 0.2e1 * (t197 * t318 - t223 * t323) / t189 ^ 2;
t328 = t224 * t235;
t239 = -t251 * t299 + t302;
t269 = t234 * t239 - t305 * t311;
t327 = t253 * t269;
t200 = qJD(5) * t272 + t221 * t257 - t253 * t283;
t232 = -t257 * t270 - t285;
t252 = sin(qJ(6));
t256 = cos(qJ(6));
t209 = t232 * t256 - t252 * t271;
t203 = 0.1e1 / t209;
t204 = 0.1e1 / t209 ^ 2;
t326 = 0.2e1 * t272;
t325 = 0.2e1 * t231;
t198 = -qJD(5) * t231 - t219 * t257 - t253 * t282;
t218 = qJD(1) * t239 + qJD(2) * t270;
t185 = qJD(6) * t209 + t198 * t252 - t218 * t256;
t208 = t232 * t252 + t256 * t271;
t202 = t208 ^ 2;
t196 = t202 * t204 + 0.1e1;
t317 = t204 * t208;
t295 = qJD(6) * t208;
t186 = t198 * t256 + t218 * t252 - t295;
t320 = t186 * t203 * t204;
t322 = (t185 * t317 - t202 * t320) / t196 ^ 2;
t313 = t234 * t328;
t321 = (-t199 * t311 - t222 * t313) / t214 ^ 2;
t319 = t191 * t197;
t316 = t208 * t256;
t315 = t210 * t231;
t314 = t211 * t231;
t312 = t272 * t234;
t309 = t271 * t253;
t308 = t271 * t257;
t303 = t252 * t203;
t296 = qJD(5) * t257;
t293 = -0.2e1 * t322;
t292 = 0.2e1 * t322;
t291 = -0.2e1 * t321;
t290 = t192 * t325;
t289 = t234 * t321;
t288 = t191 * t315;
t287 = t191 * t314;
t286 = t208 * t320;
t280 = 0.2e1 * t286;
t279 = t313 * t326;
t277 = -qJD(6) * t308 - t219;
t207 = -t228 * t256 + t239 * t252;
t206 = -t228 * t252 - t239 * t256;
t274 = t204 * t316 - t303;
t273 = -t228 * t234 - t238 * t311;
t267 = -t210 + (-t211 * t312 + t210) * t212;
t266 = qJD(5) * t309 + qJD(6) * t270 + t218 * t257;
t225 = -qJD(5) * t237 + t257 * t281;
t220 = qJD(1) * t271 + qJD(2) * t240;
t217 = t252 * t270 - t256 * t308;
t216 = -t252 * t308 - t256 * t270;
t194 = 0.1e1 / t196;
t187 = 0.1e1 / t189;
t184 = t212 * t327;
t183 = t273 * t212;
t180 = t267 * t231;
t178 = (t210 * t239 + t211 * t305) * t253 + t276 * t184;
t177 = t183 * t276 - t210 * t228 + t211 * t238;
t175 = t273 * t291 + (t238 * t279 - t200 * t234 + (t199 * t238 + t224 * t228 - t225 * t272) * t235) * t212;
t173 = t291 * t327 + (t269 * t296 + (t279 * t305 + t220 * t234 + (-t224 * t239 + (t199 * t258 + t272 * t297) * t250) * t235) * t253) * t212;
t1 = [t289 * t325 + (-t197 * t234 + t231 * t328) * t212, t173, 0, 0, t175, 0; -t272 * t190 * t294 + (-t199 * t190 + (-t176 * t272 - t180 * t197) * t191) * t187 + (t180 * t191 * t294 + (0.2e1 * t180 * t323 - (t181 * t212 * t312 + t291) * t288 - (t289 * t326 - t181 + (t181 - t275) * t212) * t287 - t267 * t319) * t187) * t231 (t178 * t318 + t190 * t309) * t294 + (-t178 * t319 + (t218 * t253 - t271 * t296) * t190 + (t178 * t290 + t191 * t309) * t176 - (t173 * t272 - t184 * t199 + (-t253 * t297 + t258 * t296) * t250 + (-t184 * t237 + t239 * t253) * t181) * t287 - (t239 * t296 - t173 * t237 - t184 * t224 + t220 * t253 + (-t184 * t272 - t253 * t305) * t181) * t288) * t187, 0, 0 (t177 * t318 - t190 * t232) * t294 + (t177 * t176 * t290 + t198 * t190 + (-t232 * t176 - t177 * t197 - (t175 * t272 - t183 * t199 + t225 + (-t183 * t237 - t228) * t181) * t314 - (-t175 * t237 - t183 * t224 - t200 + (-t183 * t272 - t238) * t181) * t315) * t191) * t187, 0; (-t203 * t206 + t207 * t317) * t292 + ((qJD(6) * t207 - t200 * t252 - t220 * t256) * t203 + t207 * t280 + (-t206 * t186 - (-qJD(6) * t206 - t200 * t256 + t220 * t252) * t208 - t207 * t185) * t204) * t194 (-t203 * t216 + t217 * t317) * t292 + (t217 * t280 + t277 * t203 * t256 + t266 * t303 + (t208 * t252 * t277 - t217 * t185 - t216 * t186 - t266 * t316) * t204) * t194, 0, 0, t274 * t231 * t293 + (t274 * t197 + ((-qJD(6) * t203 - 0.2e1 * t286) * t256 + (t185 * t256 + (t186 - t295) * t252) * t204) * t231) * t194, t293 + 0.2e1 * (t185 * t204 * t194 + (-t194 * t320 - t204 * t322) * t208) * t208;];
JaD_rot  = t1;
