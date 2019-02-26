% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRP6_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:53:16
% EndTime: 2019-02-26 19:53:17
% DurationCPUTime: 1.57s
% Computational Cost: add. (7242->157), mult. (21125->305), div. (784->12), fcn. (27225->13), ass. (0->123)
t256 = sin(qJ(4));
t259 = cos(qJ(4));
t257 = sin(qJ(2));
t260 = cos(qJ(2));
t319 = cos(pkin(10));
t320 = cos(pkin(6));
t284 = t320 * t319;
t318 = sin(pkin(10));
t273 = -t318 * t257 + t260 * t284;
t254 = sin(pkin(6));
t290 = t254 * t319;
t229 = t256 * t290 - t259 * t273;
t274 = t257 * t284 + t318 * t260;
t245 = t274 * qJD(2);
t211 = t229 * qJD(4) + t245 * t256;
t258 = cos(qJ(5));
t255 = sin(qJ(5));
t271 = t274 * t255;
t275 = t256 * t273 + t259 * t290;
t222 = -t258 * t275 + t271;
t269 = qJD(2) * t273;
t194 = t222 * qJD(5) + t211 * t255 - t258 * t269;
t270 = t274 * t258;
t220 = -t255 * t275 - t270;
t217 = t220 ^ 2;
t306 = t254 * t260;
t276 = t256 * t306 - t320 * t259;
t303 = t257 * t258;
t277 = t254 * t303 + t255 * t276;
t232 = 0.1e1 / t277 ^ 2;
t204 = t217 * t232 + 0.1e1;
t202 = 0.1e1 / t204;
t251 = -t320 * t256 - t259 * t306;
t302 = qJD(2) * t254;
t292 = t257 * t302;
t234 = t251 * qJD(4) + t256 * t292;
t305 = t255 * t257;
t237 = t254 * t305 - t258 * t276;
t291 = t260 * t302;
t206 = t237 * qJD(5) + t234 * t255 - t258 * t291;
t231 = 0.1e1 / t277;
t309 = t220 * t232;
t179 = (t194 * t231 + t206 * t309) * t202;
t205 = atan2(-t220, -t277);
t199 = sin(t205);
t200 = cos(t205);
t282 = t199 * t277 - t200 * t220;
t175 = t282 * t179 - t199 * t194 + t200 * t206;
t191 = -t199 * t220 - t200 * t277;
t188 = 0.1e1 / t191;
t189 = 0.1e1 / t191 ^ 2;
t325 = t175 * t188 * t189;
t283 = t320 * t318;
t272 = -t319 * t257 - t260 * t283;
t289 = t254 * t318;
t228 = -t256 * t272 + t259 * t289;
t250 = -t257 * t283 + t319 * t260;
t218 = t228 * t255 - t250 * t258;
t324 = -0.2e1 * t218;
t287 = 0.2e1 * t218 * t325;
t278 = t229 * t231 + t251 * t309;
t323 = t255 * t278;
t322 = qJD(5) * t270 + t255 * t269;
t312 = t206 * t231 * t232;
t321 = -0.2e1 * (t194 * t309 + t217 * t312) / t204 ^ 2;
t219 = t228 * t258 + t250 * t255;
t214 = 0.1e1 / t219;
t215 = 0.1e1 / t219 ^ 2;
t285 = t256 * t289;
t227 = -t259 * t272 - t285;
t226 = t227 ^ 2;
t311 = t215 * t226;
t201 = 0.1e1 + t311;
t247 = t250 * qJD(2);
t210 = -t228 * qJD(4) + t247 * t259;
t301 = qJD(4) * t259;
t209 = qJD(4) * t285 - t247 * t256 + t272 * t301;
t246 = t272 * qJD(2);
t193 = -t218 * qJD(5) - t209 * t258 + t246 * t255;
t315 = t193 * t214 * t215;
t294 = t226 * t315;
t310 = t215 * t227;
t317 = (t210 * t310 - t294) / t201 ^ 2;
t316 = t189 * t218;
t314 = t199 * t218;
t313 = t200 * t218;
t308 = t227 * t255;
t307 = t250 * t256;
t304 = t256 * t258;
t300 = qJD(5) * t255;
t299 = qJD(5) * t256;
t298 = qJD(5) * t258;
t213 = t218 ^ 2;
t187 = t189 * t213 + 0.1e1;
t192 = t219 * qJD(5) - t209 * t255 - t246 * t258;
t297 = 0.2e1 * (t192 * t316 - t213 * t325) / t187 ^ 2;
t296 = 0.2e1 * t317;
t293 = t227 * t315;
t286 = 0.2e1 * t220 * t312;
t280 = t222 * t231 + t237 * t309;
t223 = t256 * t271 - t258 * t273;
t240 = (t256 * t305 - t258 * t260) * t254;
t279 = t223 * t231 + t240 * t309;
t235 = t276 * qJD(4) + t259 * t292;
t225 = t250 * t304 + t255 * t272;
t224 = t255 * t307 - t258 * t272;
t212 = t275 * qJD(4) + t245 * t259;
t208 = ((qJD(2) + t299) * t303 + (t257 * t301 + (qJD(2) * t256 + qJD(5)) * t260) * t255) * t254;
t207 = t277 * qJD(5) + t234 * t258 + t255 * t291;
t197 = 0.1e1 / t201;
t196 = t245 * t258 + t256 * t322 + t271 * t301 + t273 * t300;
t195 = t211 * t258 + t275 * t300 + t322;
t185 = 0.1e1 / t187;
t184 = t202 * t323;
t183 = t279 * t202;
t181 = t280 * t202;
t178 = (-t199 * t229 + t200 * t251) * t255 + t282 * t184;
t177 = t282 * t183 - t199 * t223 + t200 * t240;
t176 = t282 * t181 - t199 * t222 + t200 * t237;
t174 = t279 * t321 + (t240 * t286 + t196 * t231 + (t194 * t240 + t206 * t223 + t208 * t220) * t232) * t202;
t172 = t280 * t321 + (t237 * t286 + t195 * t231 + (t194 * t237 + t206 * t222 + t207 * t220) * t232) * t202;
t171 = t321 * t323 + (t278 * t298 + (t251 * t286 + t212 * t231 + (t194 * t251 + t206 * t229 + t220 * t235) * t232) * t255) * t202;
t1 = [0, t174, 0, t171, t172, 0; 0 (t177 * t316 - t188 * t224) * t297 + (t177 * t287 + (-t224 * t175 - t177 * t192 - (-t174 * t220 - t183 * t194 + t208 + (t183 * t277 - t223) * t179) * t313 - (t174 * t277 - t183 * t206 - t196 + (t183 * t220 - t240) * t179) * t314) * t189 + ((t250 * t299 + t247) * t258 + (qJD(5) * t272 + t246 * t256 + t250 * t301) * t255) * t188) * t185, 0 (t178 * t316 - t188 * t308) * t297 + ((t210 * t255 + t227 * t298) * t188 + t178 * t287 + (-t178 * t192 - t308 * t175 - (t251 * t298 - t171 * t220 - t184 * t194 + t235 * t255 + (t184 * t277 - t229 * t255) * t179) * t313 - (-t229 * t298 + t171 * t277 - t184 * t206 - t212 * t255 + (t184 * t220 - t251 * t255) * t179) * t314) * t189) * t185 (t176 * t316 - t188 * t219) * t297 + (t176 * t287 + t193 * t188 + (-t219 * t175 - t176 * t192 - (-t172 * t220 - t181 * t194 + t207 + (t181 * t277 - t222) * t179) * t313 - (t172 * t277 - t181 * t206 - t195 + (t181 * t220 - t237) * t179) * t314) * t189) * t185, 0; 0 (-t214 * t250 * t259 + t225 * t310) * t296 + (0.2e1 * t225 * t293 + (-qJD(4) * t307 + t246 * t259) * t214 + (-(t246 * t304 - t247 * t255 + t272 * t298) * t227 - t225 * t210 + (-t259 * t193 - (-t255 * t299 + t258 * t301) * t227) * t250) * t215) * t197, 0 (t214 * t228 + t258 * t311) * t296 + (0.2e1 * t258 * t294 + t209 * t214 + (-0.2e1 * t210 * t227 * t258 + t193 * t228 + t226 * t300) * t215) * t197, t310 * t317 * t324 + (t293 * t324 + (t192 * t227 + t210 * t218) * t215) * t197, 0;];
JaD_rot  = t1;
