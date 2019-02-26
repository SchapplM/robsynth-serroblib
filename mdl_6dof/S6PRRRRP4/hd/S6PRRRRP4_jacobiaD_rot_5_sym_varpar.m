% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRP4
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
% Datum: 2019-02-26 20:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRRP4_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:16:59
% EndTime: 2019-02-26 20:17:00
% DurationCPUTime: 0.92s
% Computational Cost: add. (3659->111), mult. (9671->227), div. (577->12), fcn. (12382->13), ass. (0->107)
t241 = sin(pkin(11));
t243 = cos(pkin(11));
t246 = sin(qJ(2));
t244 = cos(pkin(6));
t248 = cos(qJ(2));
t275 = t244 * t248;
t227 = -t241 * t246 + t243 * t275;
t220 = t227 * qJD(2);
t276 = t244 * t246;
t228 = t241 * t248 + t243 * t276;
t245 = sin(qJ(3));
t242 = sin(pkin(6));
t279 = t242 * t245;
t266 = t243 * t279;
t247 = cos(qJ(3));
t272 = qJD(3) * t247;
t198 = -qJD(3) * t266 + t220 * t245 + t228 * t272;
t278 = t242 * t247;
t212 = t228 * t245 + t243 * t278;
t210 = t212 ^ 2;
t231 = -t244 * t247 + t246 * t279;
t225 = 0.1e1 / t231 ^ 2;
t206 = t210 * t225 + 0.1e1;
t204 = 0.1e1 / t206;
t232 = t244 * t245 + t246 * t278;
t273 = qJD(2) * t248;
t265 = t242 * t273;
t217 = t232 * qJD(3) + t245 * t265;
t224 = 0.1e1 / t231;
t283 = t212 * t225;
t176 = (-t198 * t224 + t217 * t283) * t204;
t207 = atan2(-t212, t231);
t202 = sin(t207);
t203 = cos(t207);
t260 = -t202 * t231 - t203 * t212;
t172 = t260 * t176 - t202 * t198 + t203 * t217;
t188 = -t202 * t212 + t203 * t231;
t185 = 0.1e1 / t188;
t186 = 0.1e1 / t188 ^ 2;
t296 = t172 * t185 * t186;
t267 = t241 * t276;
t230 = t243 * t248 - t267;
t257 = -t230 * t245 + t241 * t278;
t295 = -0.2e1 * t257 * t296;
t277 = t242 * t248;
t256 = -t224 * t227 + t277 * t283;
t294 = t245 * t256;
t282 = t217 * t224 * t225;
t293 = -0.2e1 * (t198 * t283 - t210 * t282) / t206 ^ 2;
t216 = t230 * t247 + t241 * t279;
t229 = t241 * t275 + t243 * t246;
t240 = qJ(4) + qJ(5);
t237 = sin(t240);
t238 = cos(t240);
t197 = t216 * t238 + t229 * t237;
t193 = 0.1e1 / t197;
t194 = 0.1e1 / t197 ^ 2;
t223 = -qJD(2) * t267 + t243 * t273;
t239 = qJD(4) + qJD(5);
t262 = t216 * t239 - t223;
t222 = t229 * qJD(2);
t201 = t257 * qJD(3) - t222 * t247;
t263 = t229 * t239 + t201;
t183 = t263 * t237 + t262 * t238;
t196 = t216 * t237 - t229 * t238;
t192 = t196 ^ 2;
t191 = t192 * t194 + 0.1e1;
t288 = t194 * t196;
t184 = -t262 * t237 + t263 * t238;
t291 = t184 * t193 * t194;
t292 = (t183 * t288 - t192 * t291) / t191 ^ 2;
t290 = t186 * t257;
t289 = t193 * t237;
t287 = t196 * t238;
t200 = t216 * qJD(3) - t222 * t245;
t286 = t200 * t186;
t285 = t202 * t257;
t284 = t203 * t257;
t281 = t229 * t245;
t280 = t229 * t247;
t274 = qJD(2) * t246;
t211 = t257 ^ 2;
t182 = t211 * t186 + 0.1e1;
t271 = 0.2e1 * (-t211 * t296 - t257 * t286) / t182 ^ 2;
t270 = -0.2e1 * t292;
t268 = t196 * t291;
t264 = -0.2e1 * t212 * t282;
t261 = t239 * t280 - t222;
t259 = t194 * t287 - t289;
t214 = t228 * t247 - t266;
t258 = -t214 * t224 + t232 * t283;
t255 = qJD(3) * t281 - t223 * t247 + t230 * t239;
t221 = t228 * qJD(2);
t218 = -t231 * qJD(3) + t247 * t265;
t209 = t230 * t237 - t238 * t280;
t208 = -t230 * t238 - t237 * t280;
t199 = -t212 * qJD(3) + t220 * t247;
t189 = 0.1e1 / t191;
t179 = 0.1e1 / t182;
t178 = t204 * t294;
t177 = t258 * t204;
t174 = (-t202 * t227 + t203 * t277) * t245 + t260 * t178;
t173 = t260 * t177 - t202 * t214 + t203 * t232;
t171 = t258 * t293 + (t232 * t264 - t199 * t224 + (t198 * t232 + t212 * t218 + t214 * t217) * t225) * t204;
t169 = t293 * t294 + (t256 * t272 + (t264 * t277 + t221 * t224 + (t217 * t227 + (t198 * t248 - t212 * t274) * t242) * t225) * t245) * t204;
t168 = t270 + 0.2e1 * (t183 * t194 * t189 + (-t189 * t291 - t194 * t292) * t196) * t196;
t1 = [0, t169, t171, 0, 0, 0; 0 (-t174 * t290 + t185 * t281) * t271 + ((-t223 * t245 - t229 * t272) * t185 + (-t286 + t295) * t174 + (t281 * t172 + (-t169 * t212 - t178 * t198 + (-t245 * t274 + t248 * t272) * t242 + (-t178 * t231 - t227 * t245) * t176) * t284 + (-t227 * t272 - t169 * t231 - t178 * t217 + t221 * t245 + (t178 * t212 - t245 * t277) * t176) * t285) * t186) * t179 (-t173 * t290 - t185 * t216) * t271 + (t173 * t295 + t201 * t185 + (-t216 * t172 - t173 * t200 + (-t171 * t212 - t177 * t198 + t218 + (-t177 * t231 - t214) * t176) * t284 + (-t171 * t231 - t177 * t217 - t199 + (t177 * t212 - t232) * t176) * t285) * t186) * t179, 0, 0, 0; 0, 0.2e1 * (-t193 * t208 + t209 * t288) * t292 + (0.2e1 * t209 * t268 - t261 * t193 * t238 + t255 * t289 + (-t261 * t196 * t237 - t209 * t183 - t208 * t184 - t255 * t287) * t194) * t189, -t259 * t257 * t270 + (t259 * t200 - ((-t193 * t239 - 0.2e1 * t268) * t238 + (t183 * t238 + (-t196 * t239 + t184) * t237) * t194) * t257) * t189, t168, t168, 0;];
JaD_rot  = t1;
