% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRR6_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:56:55
% EndTime: 2019-02-26 19:56:56
% DurationCPUTime: 0.94s
% Computational Cost: add. (3659->110), mult. (9671->232), div. (577->12), fcn. (12382->13), ass. (0->108)
t241 = sin(pkin(11));
t243 = cos(pkin(11));
t248 = cos(qJ(2));
t244 = cos(pkin(6));
t246 = sin(qJ(2));
t275 = t244 * t246;
t231 = t241 * t248 + t243 * t275;
t225 = t231 * qJD(2);
t242 = sin(pkin(6));
t247 = cos(qJ(4));
t245 = sin(qJ(4));
t274 = t244 * t248;
t258 = -t241 * t246 + t243 * t274;
t256 = t258 * t245;
t202 = qJD(4) * t256 + (t243 * t242 * qJD(4) + t225) * t247;
t279 = t242 * t245;
t215 = t243 * t279 - t247 * t258;
t212 = t215 ^ 2;
t276 = t242 * t248;
t234 = t244 * t245 + t247 * t276;
t229 = 0.1e1 / t234 ^ 2;
t207 = t212 * t229 + 0.1e1;
t205 = 0.1e1 / t207;
t235 = t244 * t247 - t245 * t276;
t278 = t242 * t246;
t265 = qJD(2) * t278;
t218 = qJD(4) * t235 - t247 * t265;
t228 = 0.1e1 / t234;
t284 = t215 * t229;
t177 = (t202 * t228 - t218 * t284) * t205;
t208 = atan2(t215, t234);
t203 = sin(t208);
t204 = cos(t208);
t261 = -t203 * t234 + t204 * t215;
t173 = t177 * t261 + t203 * t202 + t204 * t218;
t189 = t203 * t215 + t204 * t234;
t186 = 0.1e1 / t189;
t187 = 0.1e1 / t189 ^ 2;
t297 = t173 * t186 * t187;
t232 = t241 * t274 + t243 * t246;
t213 = -t232 * t247 + t241 * t279;
t296 = 0.2e1 * t213 * t297;
t282 = t218 * t228 * t229;
t295 = (t202 * t284 - t212 * t282) / t207 ^ 2;
t267 = t215 * t278;
t257 = t228 * t231 + t229 * t267;
t294 = t247 * t257;
t277 = t242 * t247;
t214 = t232 * t245 + t241 * t277;
t266 = t241 * t275;
t233 = t243 * t248 - t266;
t240 = qJ(5) + qJ(6);
t237 = sin(t240);
t238 = cos(t240);
t198 = t214 * t238 + t233 * t237;
t194 = 0.1e1 / t198;
t195 = 0.1e1 / t198 ^ 2;
t226 = t232 * qJD(2);
t239 = qJD(5) + qJD(6);
t263 = t214 * t239 + t226;
t273 = qJD(2) * t248;
t227 = -qJD(2) * t266 + t243 * t273;
t199 = -qJD(4) * t213 + t227 * t245;
t264 = t233 * t239 + t199;
t184 = t237 * t264 + t238 * t263;
t197 = t214 * t237 - t233 * t238;
t193 = t197 ^ 2;
t192 = t193 * t195 + 0.1e1;
t289 = t195 * t197;
t185 = -t237 * t263 + t238 * t264;
t292 = t185 * t194 * t195;
t293 = (t184 * t289 - t193 * t292) / t192 ^ 2;
t291 = t187 * t213;
t290 = t194 * t237;
t288 = t197 * t238;
t200 = qJD(4) * t214 - t227 * t247;
t287 = t200 * t187;
t286 = t203 * t213;
t285 = t204 * t213;
t283 = t215 * t235;
t281 = t233 * t245;
t280 = t233 * t247;
t272 = qJD(4) * t245;
t211 = t213 ^ 2;
t183 = t211 * t187 + 0.1e1;
t271 = 0.2e1 * (-t211 * t297 + t213 * t287) / t183 ^ 2;
t270 = -0.2e1 * t293;
t268 = t197 * t292;
t262 = t239 * t281 + t227;
t260 = t195 * t288 - t290;
t216 = t243 * t277 + t256;
t259 = -t216 * t228 + t229 * t283;
t255 = qJD(4) * t280 - t226 * t245 - t232 * t239;
t224 = t258 * qJD(2);
t217 = -qJD(4) * t234 + t245 * t265;
t210 = -t232 * t237 + t238 * t281;
t209 = t232 * t238 + t237 * t281;
t201 = qJD(4) * t215 + t225 * t245;
t190 = 0.1e1 / t192;
t180 = 0.1e1 / t183;
t179 = t205 * t294;
t178 = t259 * t205;
t175 = (t203 * t231 - t204 * t278) * t247 + t261 * t179;
t174 = -t178 * t261 + t203 * t216 + t204 * t235;
t172 = 0.2e1 * t259 * t295 + (0.2e1 * t282 * t283 - t201 * t228 + (-t202 * t235 - t215 * t217 - t216 * t218) * t229) * t205;
t170 = -0.2e1 * t294 * t295 + (-t257 * t272 + (-0.2e1 * t267 * t282 + t224 * t228 + (-t218 * t231 + (t202 * t246 + t215 * t273) * t242) * t229) * t247) * t205;
t169 = t270 + 0.2e1 * (t184 * t195 * t190 + (-t190 * t292 - t195 * t293) * t197) * t197;
t1 = [0, t170, 0, t172, 0, 0; 0 (t175 * t291 + t186 * t280) * t271 + ((t226 * t247 + t233 * t272) * t186 + (-t287 + t296) * t175 + (t280 * t173 - (t170 * t215 + t179 * t202 + (t246 * t272 - t247 * t273) * t242 + (-t179 * t234 + t231 * t247) * t177) * t285 - (-t231 * t272 - t170 * t234 - t179 * t218 + t224 * t247 + (-t179 * t215 + t246 * t277) * t177) * t286) * t187) * t180, 0 (t174 * t291 - t186 * t214) * t271 + (t174 * t296 + t199 * t186 + (-t214 * t173 - t174 * t200 - (t172 * t215 - t178 * t202 + t217 + (t178 * t234 + t216) * t177) * t285 - (-t172 * t234 + t178 * t218 - t201 + (t178 * t215 - t235) * t177) * t286) * t187) * t180, 0, 0; 0, 0.2e1 * (-t194 * t209 + t210 * t289) * t293 + (0.2e1 * t210 * t268 + t262 * t194 * t238 + t255 * t290 + (t197 * t237 * t262 - t210 * t184 - t209 * t185 - t255 * t288) * t195) * t190, 0, t260 * t213 * t270 + (t260 * t200 + ((-t194 * t239 - 0.2e1 * t268) * t238 + (t184 * t238 + (-t197 * t239 + t185) * t237) * t195) * t213) * t190, t169, t169;];
JaD_rot  = t1;
