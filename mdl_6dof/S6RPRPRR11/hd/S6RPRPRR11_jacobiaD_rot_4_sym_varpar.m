% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRR11
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR11_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR11_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_jacobiaD_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:54:31
% EndTime: 2019-02-26 20:54:33
% DurationCPUTime: 1.27s
% Computational Cost: add. (4301->104), mult. (13743->216), div. (424->12), fcn. (17471->15), ass. (0->107)
t307 = sin(pkin(12));
t311 = cos(pkin(6));
t279 = t311 * t307;
t309 = cos(pkin(12));
t312 = sin(qJ(1));
t314 = cos(qJ(1));
t237 = t314 * t279 + t312 * t309;
t250 = sin(qJ(3));
t248 = sin(pkin(6));
t308 = sin(pkin(7));
t285 = t248 * t308;
t275 = t314 * t285;
t244 = t250 * t275;
t281 = t311 * t309;
t236 = -t314 * t281 + t312 * t307;
t310 = cos(pkin(7));
t287 = t236 * t310;
t313 = cos(qJ(3));
t321 = -t237 * t313 + t250 * t287 + t244;
t238 = -t312 * t279 + t314 * t309;
t261 = t312 * t281 + t314 * t307;
t259 = t261 * t310;
t273 = t312 * t285;
t217 = t238 * t313 + (-t259 + t273) * t250;
t286 = t248 * t310;
t274 = t312 * t286;
t228 = t261 * t308 + t274;
t247 = sin(pkin(13));
t249 = cos(pkin(13));
t201 = t217 * t247 - t228 * t249;
t233 = t237 * qJD(1);
t258 = qJD(1) * t287;
t267 = t313 * t273;
t318 = -t238 * t250 - t313 * t259 + t267;
t191 = qJD(1) * t244 + t318 * qJD(3) - t233 * t313 + t250 * t258;
t227 = -t236 * t308 + t314 * t286;
t220 = t227 * qJD(1);
t189 = t191 * t249 + t220 * t247;
t202 = t217 * t249 + t228 * t247;
t196 = 0.1e1 / t202;
t197 = 0.1e1 / t202 ^ 2;
t301 = t189 * t196 * t197;
t320 = 0.2e1 * t201 * t301;
t268 = t313 * t275;
t283 = t310 * t313;
t319 = -t236 * t283 - t268;
t278 = t310 * t309;
t280 = t311 * t308;
t226 = t250 * t280 + (t250 * t278 + t313 * t307) * t248;
t219 = t226 * qJD(3);
t317 = (-t307 * t250 + t313 * t278) * t248 + t313 * t280;
t223 = 0.1e1 / t317 ^ 2;
t296 = t219 * t223;
t234 = t261 * qJD(1);
t235 = t238 * qJD(1);
t193 = (qJD(1) * t273 - qJD(3) * t237 - t310 * t234) * t250 + t235 * t313 + t319 * qJD(3);
t316 = qJD(1) * t267 + t321 * qJD(3) - t234 * t283 - t235 * t250;
t211 = t237 * t250 - t319;
t208 = atan2(-t211, -t317);
t203 = sin(t208);
t204 = cos(t208);
t184 = -t203 * t211 - t204 * t317;
t181 = 0.1e1 / t184;
t222 = 0.1e1 / t317;
t182 = 0.1e1 / t184 ^ 2;
t315 = -0.2e1 * t318;
t210 = t318 ^ 2;
t180 = t210 * t182 + 0.1e1;
t190 = -qJD(1) * t268 + t217 * qJD(3) - t233 * t250 - t313 * t258;
t302 = t182 * t318;
t209 = t211 ^ 2;
t207 = t209 * t223 + 0.1e1;
t205 = 0.1e1 / t207;
t272 = t211 * t296 - t222 * t316;
t175 = t272 * t205;
t277 = t203 * t317 - t204 * t211;
t171 = t277 * t175 + t203 * t316 + t204 * t219;
t305 = t171 * t181 * t182;
t306 = (-t190 * t302 - t210 * t305) / t180 ^ 2;
t295 = t222 * t296;
t304 = (-t211 * t223 * t316 + t209 * t295) / t207 ^ 2;
t178 = 0.1e1 / t180;
t303 = t178 * t182;
t300 = t197 * t201;
t299 = t201 * t249;
t298 = t211 * t222;
t297 = t211 * t226;
t293 = t247 * t196;
t292 = 0.2e1 * t306;
t195 = t201 ^ 2;
t187 = t195 * t197 + 0.1e1;
t188 = t191 * t247 - t220 * t249;
t291 = 0.2e1 * (t188 * t300 - t195 * t301) / t187 ^ 2;
t290 = -0.2e1 * t304;
t289 = t222 * t304;
t284 = t305 * t315;
t271 = -t222 * t321 + t223 * t297;
t266 = -t203 + (-t204 * t298 + t203) * t205;
t221 = -qJD(1) * t274 - t234 * t308;
t218 = t317 * qJD(3);
t200 = t227 * t247 + t249 * t321;
t199 = -t227 * t249 + t247 * t321;
t185 = 0.1e1 / t187;
t176 = t271 * t205;
t172 = t277 * t176 + t203 * t321 + t204 * t226;
t170 = t271 * t290 + (0.2e1 * t295 * t297 + t193 * t222 + (t211 * t218 - t219 * t321 - t226 * t316) * t223) * t205;
t1 = [-t289 * t315 + (t190 * t222 - t296 * t318) * t205, 0, t170, 0, 0, 0; -(-t171 * t303 - 0.2e1 * t181 * t306) * t211 + (t316 * t181 + (t266 * t190 - ((t175 * t205 * t298 + t290) * t203 + (0.2e1 * t211 * t289 - t175 + (t175 - t272) * t205) * t204) * t318) * t302) * t178 - (t178 * t284 - t190 * t303 - t302 * t292) * t266 * t318, 0 (-t172 * t302 - t181 * t217) * t292 + (t172 * t284 + t191 * t181 + (-t217 * t171 - t172 * t190 - (-(-t170 * t211 + t176 * t316 + t218 + (t176 * t317 + t321) * t175) * t204 - (t170 * t317 - t176 * t219 - t193 + (t176 * t211 - t226) * t175) * t203) * t318) * t182) * t178, 0, 0, 0; (-t196 * t199 + t200 * t300) * t291 + ((-t193 * t247 - t221 * t249) * t196 + t200 * t320 + (-t199 * t189 - (-t193 * t249 + t221 * t247) * t201 - t200 * t188) * t197) * t185, 0 -(-t197 * t299 + t293) * t318 * t291 + (t318 * t249 * t320 - t190 * t293 + (t190 * t299 - (t188 * t249 + t189 * t247) * t318) * t197) * t185, 0, 0, 0;];
JaD_rot  = t1;
