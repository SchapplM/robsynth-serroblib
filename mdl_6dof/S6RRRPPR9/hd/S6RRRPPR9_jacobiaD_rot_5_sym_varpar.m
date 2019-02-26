% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPPR9_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR9_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:08:10
% EndTime: 2019-02-26 22:08:12
% DurationCPUTime: 1.83s
% Computational Cost: add. (6351->137), mult. (18580->286), div. (681->12), fcn. (23728->13), ass. (0->126)
t247 = cos(pkin(6));
t249 = sin(qJ(2));
t250 = sin(qJ(1));
t301 = t250 * t249;
t288 = t247 * t301;
t252 = cos(qJ(2));
t253 = cos(qJ(1));
t298 = t253 * t252;
t223 = -qJD(1) * t288 - qJD(2) * t301 + (qJD(2) * t247 + qJD(1)) * t298;
t299 = t253 * t249;
t300 = t250 * t252;
t240 = t247 * t299 + t300;
t248 = sin(qJ(3));
t251 = cos(qJ(3));
t246 = sin(pkin(6));
t304 = t246 * t253;
t231 = t240 * t248 + t251 * t304;
t297 = qJD(1) * t246;
t287 = t250 * t297;
t202 = t231 * qJD(3) - t223 * t251 - t248 * t287;
t245 = sin(pkin(11));
t321 = cos(pkin(11));
t268 = t247 * t300 + t299;
t323 = t268 * qJD(1) + t240 * qJD(2);
t328 = t202 * t245 + t323 * t321;
t232 = -t240 * t251 + t248 * t304;
t276 = -t247 * t298 + t301;
t210 = t232 * t245 + t276 * t321;
t302 = t249 * t251;
t269 = t246 * t302 + t247 * t248;
t282 = t252 * t321;
t227 = t269 * t245 + t246 * t282;
t195 = atan2(t210, t227);
t190 = sin(t195);
t191 = cos(t195);
t273 = -t190 * t210 - t191 * t227;
t183 = 0.1e1 / t273 ^ 2;
t267 = t288 - t298;
t305 = t246 * t250;
t234 = t248 * t305 - t251 * t267;
t261 = t268 * t321;
t212 = t234 * t245 - t261;
t204 = t212 ^ 2;
t181 = t204 * t183 + 0.1e1;
t222 = -t240 * qJD(1) - t268 * qJD(2);
t233 = t248 * t267 + t251 * t305;
t286 = t253 * t297;
t200 = t233 * qJD(3) + t222 * t251 + t248 * t286;
t221 = t276 * qJD(1) + t267 * qJD(2);
t285 = t221 * t321;
t186 = t200 * t245 + t285;
t315 = t186 * t183;
t182 = 0.1e1 / t273;
t203 = t210 ^ 2;
t225 = 0.1e1 / t227 ^ 2;
t194 = t203 * t225 + 0.1e1;
t192 = 0.1e1 / t194;
t224 = 0.1e1 / t227;
t283 = t249 * t321;
t295 = qJD(3) * t248;
t296 = qJD(2) * t252;
t303 = t247 * t251;
t214 = qJD(3) * t245 * t303 + ((-t249 * t295 + t251 * t296) * t245 - qJD(2) * t283) * t246;
t310 = t214 * t225;
t272 = -t210 * t310 + t224 * t328;
t174 = t272 * t192;
t274 = -t190 * t227 + t191 * t210;
t170 = t274 * t174 + t190 * t328 + t191 * t214;
t327 = t170 * t183;
t319 = t182 * t327;
t294 = 0.2e1 * (t204 * t319 + t212 * t315) / t181 ^ 2;
t325 = t268 * t245;
t293 = -0.2e1 * t319;
t324 = t212 * t293 - t315;
t201 = t232 * qJD(3) - t223 * t248 + t251 * t287;
t213 = t234 * t321 + t325;
t205 = 0.1e1 / t213;
t206 = 0.1e1 / t213 ^ 2;
t322 = 0.2e1 * t210;
t309 = t224 * t310;
t311 = t210 * t225;
t317 = (-t203 * t309 + t311 * t328) / t194 ^ 2;
t316 = t183 * t212;
t187 = t200 * t321 - t221 * t245;
t207 = t205 * t206;
t314 = t187 * t207;
t313 = t206 * t233;
t312 = t210 * t224;
t307 = t245 * t252;
t306 = t246 * t248;
t228 = t233 ^ 2;
t198 = t228 * t206 + 0.1e1;
t199 = -t234 * qJD(3) - t222 * t248 + t251 * t286;
t289 = t199 * t313;
t292 = 0.2e1 * (-t228 * t314 + t289) / t198 ^ 2;
t291 = -0.2e1 * t317;
t290 = t224 * t317;
t284 = t228 * t321;
t281 = t182 * t294;
t280 = t183 * t294;
t279 = t309 * t322;
t277 = 0.2e1 * t233 * t314;
t275 = t276 * t245;
t215 = -t240 * t321 - t251 * t275;
t235 = (t251 * t307 - t283) * t246;
t271 = -t215 * t224 - t235 * t311;
t239 = -t249 * t306 + t303;
t270 = t224 * t231 - t239 * t311;
t264 = t248 * t268;
t263 = -t190 + (-t191 * t312 + t190) * t192;
t262 = qJD(3) * t268;
t229 = -t269 * qJD(3) - t296 * t306;
t220 = (-t295 * t307 + (-t245 * t302 - t282) * qJD(2)) * t246;
t217 = -t245 * t267 - t251 * t261;
t216 = -t251 * t325 + t267 * t321;
t211 = t232 * t321 - t275;
t196 = 0.1e1 / t198;
t189 = -t223 * t321 + (-t251 * t323 + t276 * t295) * t245;
t179 = 0.1e1 / t181;
t178 = t270 * t245 * t192;
t176 = t271 * t192;
t173 = t263 * t212;
t172 = (t190 * t231 + t191 * t239) * t245 + t274 * t178;
t169 = (t270 * t291 + (t239 * t279 - t201 * t224 + (-t210 * t229 - t214 * t231 - t239 * t328) * t225) * t192) * t245;
t168 = t271 * t291 + (t235 * t279 - t189 * t224 + (-t210 * t220 + t214 * t215 - t235 * t328) * t225) * t192;
t1 = [0.2e1 * t212 * t290 + (-t186 * t224 + t212 * t310) * t192, t168, t169, 0, 0, 0; t210 * t281 + (-t328 * t182 + (-t210 * t170 - t173 * t186) * t183) * t179 + (t173 * t280 + (t173 * t293 - t263 * t315 + (-(t174 * t192 * t312 + t291) * t190 - (t290 * t322 - t174 + (t174 - t272) * t192) * t191) * t316) * t179) * t212, t216 * t281 + (-(-t222 * t321 + (t221 * t251 + t248 * t262) * t245) * t182 - t216 * t327 - ((t168 * t210 + t176 * t328 + t220 + (-t176 * t227 - t215) * t174) * t191 + (-t168 * t227 - t176 * t214 - t189 + (-t176 * t210 - t235) * t174) * t190) * t316) * t179 + (t324 * t179 + t212 * t280) * (t274 * t176 - t190 * t215 + t191 * t235) (t182 * t233 * t245 + t172 * t316) * t294 + (-(t274 * t169 + (t273 * t174 - t190 * t214 + t191 * t328) * t178) * t316 + t324 * t172 + (-t199 * t182 + (-t233 * t170 - (-t190 * t201 + t191 * t229 + (-t190 * t239 + t191 * t231) * t174) * t212) * t183) * t245) * t179, 0, 0, 0; (-t205 * t231 + t211 * t313) * t292 + (-t201 * t205 + t211 * t277 + (-t231 * t187 - (t202 * t321 - t245 * t323) * t233 - t211 * t199) * t206) * t196 (-t205 * t264 + t217 * t313) * t292 + (t217 * t277 + (-t221 * t248 + t251 * t262) * t205 + (-t187 * t264 - (t222 * t245 + t251 * t285 + t261 * t295) * t233 - t217 * t199) * t206) * t196 (t205 * t234 + t206 * t284) * t292 + (-0.2e1 * t321 * t289 - t200 * t205 + (t206 * t234 + 0.2e1 * t207 * t284) * t187) * t196, 0, 0, 0;];
JaD_rot  = t1;
