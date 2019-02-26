% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function JaD_rot = S6PRPRPR2_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:47:03
% EndTime: 2019-02-26 19:47:05
% DurationCPUTime: 1.14s
% Computational Cost: add. (5050->106), mult. (14713->222), div. (535->12), fcn. (19181->15), ass. (0->105)
t244 = sin(pkin(11));
t248 = cos(pkin(11));
t252 = sin(qJ(2));
t254 = cos(qJ(2));
t236 = t244 * t252 - t254 * t248;
t250 = cos(pkin(6));
t262 = t236 * t250;
t230 = qJD(2) * t262;
t267 = t244 * t254 + t252 * t248;
t235 = t267 * qJD(2);
t245 = sin(pkin(10));
t249 = cos(pkin(10));
t212 = -t230 * t249 - t235 * t245;
t233 = t267 * t250;
t217 = t233 * t249 - t236 * t245;
t251 = sin(qJ(4));
t246 = sin(pkin(6));
t281 = t246 * t251;
t272 = t249 * t281;
t253 = cos(qJ(4));
t277 = qJD(4) * t253;
t184 = -qJD(4) * t272 + t212 * t251 + t217 * t277;
t280 = t246 * t253;
t206 = t217 * t251 + t249 * t280;
t204 = t206 ^ 2;
t232 = t267 * t246;
t224 = t232 * t251 - t250 * t253;
t222 = 0.1e1 / t224 ^ 2;
t200 = t204 * t222 + 0.1e1;
t198 = 0.1e1 / t200;
t225 = t232 * t253 + t250 * t251;
t231 = t236 * t246;
t229 = qJD(2) * t231;
t202 = t225 * qJD(4) - t229 * t251;
t221 = 0.1e1 / t224;
t285 = t206 * t222;
t168 = (-t184 * t221 + t202 * t285) * t198;
t201 = atan2(-t206, t224);
t196 = sin(t201);
t197 = cos(t201);
t270 = -t196 * t224 - t197 * t206;
t164 = t270 * t168 - t196 * t184 + t197 * t202;
t178 = -t196 * t206 + t197 * t224;
t175 = 0.1e1 / t178;
t176 = 0.1e1 / t178 ^ 2;
t298 = t164 * t175 * t176;
t268 = -t233 * t245 - t236 * t249;
t210 = t245 * t281 + t253 * t268;
t219 = t245 * t262 - t249 * t267;
t243 = sin(pkin(12));
t247 = cos(pkin(12));
t192 = t210 * t243 + t219 * t247;
t264 = t245 * t280 - t251 * t268;
t269 = t230 * t245 - t235 * t249;
t187 = t264 * qJD(4) + t253 * t269;
t234 = t236 * qJD(2);
t261 = t250 * t235;
t213 = t234 * t249 + t245 * t261;
t183 = t187 * t247 - t213 * t243;
t193 = t210 * t247 - t219 * t243;
t189 = 0.1e1 / t193;
t190 = 0.1e1 / t193 ^ 2;
t292 = t183 * t189 * t190;
t297 = 0.2e1 * t192 * t292;
t296 = -0.2e1 * t264 * t298;
t216 = -t245 * t267 - t249 * t262;
t265 = -t216 * t221 - t231 * t285;
t295 = t251 * t265;
t286 = t202 * t221 * t222;
t294 = -0.2e1 * (t184 * t285 - t204 * t286) / t200 ^ 2;
t293 = t176 * t264;
t291 = t189 * t243;
t290 = t190 * t192;
t289 = t192 * t247;
t288 = t196 * t264;
t287 = t197 * t264;
t284 = t219 * t251;
t283 = t219 * t253;
t205 = t264 ^ 2;
t174 = t176 * t205 + 0.1e1;
t186 = t210 * qJD(4) + t251 * t269;
t276 = 0.2e1 * (-t186 * t293 - t205 * t298) / t174 ^ 2;
t188 = t192 ^ 2;
t181 = t188 * t190 + 0.1e1;
t182 = t187 * t243 + t213 * t247;
t275 = 0.2e1 * (t182 * t290 - t188 * t292) / t181 ^ 2;
t271 = -0.2e1 * t206 * t286;
t208 = t217 * t253 - t272;
t266 = -t208 * t221 + t225 * t285;
t263 = -qJD(4) * t284 + t213 * t253;
t228 = t246 * t235;
t211 = t245 * t234 - t249 * t261;
t203 = -t224 * qJD(4) - t229 * t253;
t195 = t243 * t268 + t247 * t283;
t194 = t243 * t283 - t247 * t268;
t185 = -t206 * qJD(4) + t212 * t253;
t179 = 0.1e1 / t181;
t172 = 0.1e1 / t174;
t170 = t198 * t295;
t169 = t266 * t198;
t166 = (-t196 * t216 - t197 * t231) * t251 + t270 * t170;
t165 = t270 * t169 - t196 * t208 + t197 * t225;
t163 = t266 * t294 + (t225 * t271 - t185 * t221 + (t184 * t225 + t202 * t208 + t203 * t206) * t222) * t198;
t161 = t294 * t295 + (t265 * t277 + (-t231 * t271 - t211 * t221 + (-t184 * t231 + t202 * t216 - t206 * t228) * t222) * t251) * t198;
t1 = [0, t161, 0, t163, 0, 0; 0 (-t166 * t293 - t175 * t284) * t276 + ((t213 * t251 + t219 * t277) * t175 + t166 * t296 + (-t166 * t186 - t284 * t164 + (-t231 * t277 - t161 * t206 - t170 * t184 - t228 * t251 + (-t170 * t224 - t216 * t251) * t168) * t287 + (-t216 * t277 - t161 * t224 - t170 * t202 - t211 * t251 + (t170 * t206 + t231 * t251) * t168) * t288) * t176) * t172, 0 (-t165 * t293 - t175 * t210) * t276 + (t165 * t296 + t187 * t175 + (-t210 * t164 - t165 * t186 + (-t163 * t206 - t169 * t184 + t203 + (-t169 * t224 - t208) * t168) * t287 + (-t163 * t224 - t169 * t202 - t185 + (t169 * t206 - t225) * t168) * t288) * t176) * t172, 0, 0; 0 (-t189 * t194 + t195 * t290) * t275 + ((t263 * t243 - t247 * t269) * t189 + t195 * t297 + (-t194 * t183 - (t243 * t269 + t263 * t247) * t192 - t195 * t182) * t190) * t179, 0 -(-t190 * t289 + t291) * t264 * t275 + (t264 * t247 * t297 - t186 * t291 + (t186 * t289 - (t182 * t247 + t183 * t243) * t264) * t190) * t179, 0, 0;];
JaD_rot  = t1;
