% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR10
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRR10_jacobia_rot_6_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobia_rot_6_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobia_rot_6_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:27:18
% EndTime: 2018-11-23 11:27:20
% DurationCPUTime: 1.97s
% Computational Cost: add. (25558->128), mult. (25589->247), div. (145->9), fcn. (26142->33), ass. (0->125)
t272 = sin(qJ(2));
t273 = sin(qJ(1));
t279 = cos(qJ(1));
t306 = pkin(6) + qJ(2);
t292 = cos(t306) / 0.2e1;
t307 = pkin(6) - qJ(2);
t298 = cos(t307);
t286 = t298 / 0.2e1 + t292;
t234 = t273 * t272 - t279 * t286;
t289 = sin(t306) / 0.2e1;
t295 = sin(t307);
t246 = t289 - t295 / 0.2e1;
t278 = cos(qJ(2));
t235 = t279 * t246 + t273 * t278;
t304 = pkin(7) + qJ(3);
t288 = sin(t304) / 0.2e1;
t305 = pkin(7) - qJ(3);
t294 = sin(t305);
t244 = t288 - t294 / 0.2e1;
t291 = cos(t304) / 0.2e1;
t297 = cos(t305);
t249 = t291 - t297 / 0.2e1;
t277 = cos(qJ(3));
t266 = sin(pkin(6));
t308 = t266 * t279;
t218 = -t234 * t244 + t235 * t277 + t249 * t308;
t243 = t288 + t294 / 0.2e1;
t250 = t297 / 0.2e1 + t291;
t271 = sin(qJ(3));
t220 = t234 * t250 + t235 * t271 + t243 * t308;
t303 = pkin(8) - qJ(4);
t293 = sin(t303);
t302 = pkin(8) + qJ(4);
t323 = sin(t302) / 0.2e1;
t242 = t323 - t293 / 0.2e1;
t290 = cos(t302) / 0.2e1;
t296 = cos(t303);
t247 = t290 - t296 / 0.2e1;
t276 = cos(qJ(4));
t265 = sin(pkin(7));
t322 = cos(pkin(7));
t299 = t266 * t322;
t285 = t234 * t265 - t279 * t299;
t197 = t218 * t276 - t220 * t242 - t247 * t285;
t269 = sin(qJ(5));
t275 = cos(qJ(5));
t264 = sin(pkin(8));
t321 = cos(pkin(8));
t281 = t220 * t264 + t285 * t321;
t184 = t197 * t275 + t269 * t281;
t182 = t197 * t269 - t275 * t281;
t241 = t323 + t293 / 0.2e1;
t248 = t296 / 0.2e1 + t290;
t270 = sin(qJ(4));
t324 = -t218 * t270 - t220 * t248 + t241 * t285;
t245 = t289 + t295 / 0.2e1;
t251 = t292 - t298 / 0.2e1;
t267 = cos(pkin(6));
t227 = t267 * t243 + t245 * t250 + t251 * t271;
t228 = t245 * t244 - t267 * t249 - t251 * t277;
t233 = -t245 * t265 + t267 * t322;
t210 = t227 * t242 + t228 * t276 - t233 * t247;
t216 = -t227 * t264 + t233 * t321;
t194 = t210 * t269 - t216 * t275;
t177 = atan2(-t182, t194);
t174 = sin(t177);
t175 = cos(t177);
t172 = -t174 * t182 + t175 * t194;
t171 = 0.1e1 / t172 ^ 2;
t237 = -t279 * t272 - t273 * t286;
t239 = t273 * t246 - t279 * t278;
t309 = t266 * t273;
t283 = t237 * t250 + t239 * t271 + t243 * t309;
t284 = t237 * t265 - t273 * t299;
t280 = -t264 * t283 - t284 * t321;
t223 = -t237 * t244 + t239 * t277 + t249 * t309;
t282 = -t223 * t276 + t242 * t283 + t247 * t284;
t186 = t269 * t282 - t275 * t280;
t320 = t171 * t186;
t319 = t175 * t182;
t187 = t269 * t280 + t275 * t282;
t200 = -t223 * t270 + t241 * t284 - t248 * t283;
t268 = sin(qJ(6));
t274 = cos(qJ(6));
t181 = t187 * t274 + t200 * t268;
t179 = 0.1e1 / t181 ^ 2;
t180 = t187 * t268 - t200 * t274;
t318 = t179 * t180;
t193 = 0.1e1 / t194 ^ 2;
t317 = t182 * t193;
t316 = t186 ^ 2 * t171;
t315 = t200 * t275;
t311 = t264 * t275;
t310 = t265 * t247;
t301 = t180 ^ 2 * t179 + 0.1e1;
t300 = t265 * t321;
t287 = -t174 * t194 - t319;
t231 = -t245 * t271 + t251 * t250;
t226 = t237 * t277 + t239 * t244;
t225 = -t237 * t271 + t239 * t250;
t224 = t234 * t271 - t235 * t250;
t214 = -t225 * t264 - t239 * t300;
t209 = t227 * t248 - t228 * t270 + t233 * t241;
t208 = t223 * t242 + t276 * t283;
t207 = -t223 * t248 + t270 * t283;
t206 = t225 * t242 + t226 * t276 + t239 * t310;
t205 = t239 * t265 * t241 - t225 * t248 + t226 * t270;
t204 = (t227 * t276 - t228 * t242) * t269 - t228 * t311;
t203 = ((t251 * t244 + t245 * t277) * t276 + t231 * t242 + t251 * t310) * t269 - (-t231 * t264 - t251 * t300) * t275;
t195 = t210 * t275 + t216 * t269;
t192 = 0.1e1 / t194;
t191 = -t223 * t264 * t269 + t208 * t275;
t190 = (-t218 * t242 - t220 * t276) * t269 - t218 * t311;
t189 = t206 * t275 + t214 * t269;
t188 = ((-t234 * t277 - t235 * t244) * t276 + t224 * t242 - t235 * t310) * t269 - (-t224 * t264 + t235 * t300) * t275;
t178 = 0.1e1 / t181;
t176 = 0.1e1 / (t182 ^ 2 * t193 + 0.1e1);
t173 = 0.1e1 / t301;
t170 = 0.1e1 / t172;
t169 = 0.1e1 / (0.1e1 + t316);
t168 = (-t192 * t324 + t209 * t317) * t269 * t176;
t167 = (-t190 * t192 + t204 * t317) * t176;
t166 = (-t188 * t192 + t203 * t317) * t176;
t165 = (-t184 * t192 + t195 * t317) * t176;
t1 = [-t186 * t192 * t176, t166, t167, t168, t165, 0; (-t182 * t170 - (-t174 + (t192 * t319 + t174) * t176) * t316) * t169 ((t206 * t269 - t214 * t275) * t170 - (t166 * t287 - t174 * t188 + t175 * t203) * t320) * t169 ((t208 * t269 + t223 * t311) * t170 - (t167 * t287 - t174 * t190 + t175 * t204) * t320) * t169 (-t200 * t269 * t170 - ((-t174 * t324 + t175 * t209) * t269 + t287 * t168) * t320) * t169 (t187 * t170 - (t165 * t287 - t174 * t184 + t175 * t195) * t320) * t169, 0; ((-t184 * t268 - t274 * t324) * t178 - (-t184 * t274 + t268 * t324) * t318) * t173 ((t189 * t268 - t205 * t274) * t178 - (t189 * t274 + t205 * t268) * t318) * t173 ((t191 * t268 - t207 * t274) * t178 - (t191 * t274 + t207 * t268) * t318) * t173 ((-t268 * t315 - t274 * t282) * t178 - (t268 * t282 - t274 * t315) * t318) * t173 (-t268 * t178 + t274 * t318) * t186 * t173, t301 * t173;];
Ja_rot  = t1;
