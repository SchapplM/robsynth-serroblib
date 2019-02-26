% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRP2_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:50:55
% EndTime: 2019-02-26 19:50:56
% DurationCPUTime: 1.18s
% Computational Cost: add. (5642->114), mult. (16325->235), div. (559->12), fcn. (21250->15), ass. (0->110)
t270 = sin(pkin(11));
t273 = cos(pkin(11));
t278 = sin(qJ(2));
t281 = cos(qJ(2));
t263 = t278 * t270 - t281 * t273;
t275 = cos(pkin(6));
t290 = t263 * t275;
t257 = qJD(2) * t290;
t295 = t281 * t270 + t278 * t273;
t262 = t295 * qJD(2);
t271 = sin(pkin(10));
t274 = cos(pkin(10));
t238 = -t274 * t257 - t271 * t262;
t260 = t295 * t275;
t244 = t274 * t260 - t271 * t263;
t277 = sin(qJ(4));
t272 = sin(pkin(6));
t312 = t272 * t277;
t301 = t274 * t312;
t280 = cos(qJ(4));
t307 = qJD(4) * t280;
t210 = -qJD(4) * t301 + t238 * t277 + t244 * t307;
t311 = t272 * t280;
t232 = t244 * t277 + t274 * t311;
t230 = t232 ^ 2;
t259 = t295 * t272;
t251 = t259 * t277 - t275 * t280;
t249 = 0.1e1 / t251 ^ 2;
t226 = t230 * t249 + 0.1e1;
t224 = 0.1e1 / t226;
t252 = t259 * t280 + t275 * t277;
t258 = t263 * t272;
t256 = qJD(2) * t258;
t228 = t252 * qJD(4) - t256 * t277;
t248 = 0.1e1 / t251;
t315 = t232 * t249;
t194 = (-t210 * t248 + t228 * t315) * t224;
t227 = atan2(-t232, t251);
t222 = sin(t227);
t223 = cos(t227);
t298 = -t222 * t251 - t223 * t232;
t190 = t298 * t194 - t222 * t210 + t223 * t228;
t206 = -t222 * t232 + t223 * t251;
t203 = 0.1e1 / t206;
t204 = 0.1e1 / t206 ^ 2;
t329 = t190 * t203 * t204;
t296 = -t271 * t260 - t274 * t263;
t291 = t271 * t311 - t277 * t296;
t328 = -0.2e1 * t291 * t329;
t243 = -t271 * t295 - t274 * t290;
t292 = -t243 * t248 - t258 * t315;
t327 = t277 * t292;
t316 = t228 * t248 * t249;
t326 = -0.2e1 * (t210 * t315 - t230 * t316) / t226 ^ 2;
t236 = t271 * t312 + t280 * t296;
t246 = t271 * t290 - t274 * t295;
t276 = sin(qJ(5));
t279 = cos(qJ(5));
t219 = t236 * t279 - t246 * t276;
t215 = 0.1e1 / t219;
t216 = 0.1e1 / t219 ^ 2;
t297 = t271 * t257 - t274 * t262;
t213 = t291 * qJD(4) + t280 * t297;
t261 = t263 * qJD(2);
t289 = t275 * t262;
t239 = t274 * t261 + t271 * t289;
t201 = t219 * qJD(5) + t213 * t276 + t239 * t279;
t218 = t236 * t276 + t246 * t279;
t214 = t218 ^ 2;
t209 = t214 * t216 + 0.1e1;
t320 = t216 * t218;
t306 = qJD(5) * t218;
t202 = t213 * t279 - t239 * t276 - t306;
t324 = t202 * t215 * t216;
t325 = (t201 * t320 - t214 * t324) / t209 ^ 2;
t323 = t204 * t291;
t212 = t236 * qJD(4) + t277 * t297;
t322 = t212 * t204;
t321 = t215 * t276;
t319 = t218 * t279;
t318 = t222 * t291;
t317 = t223 * t291;
t314 = t246 * t277;
t313 = t246 * t280;
t231 = t291 ^ 2;
t200 = t231 * t204 + 0.1e1;
t305 = 0.2e1 * (-t231 * t329 - t291 * t322) / t200 ^ 2;
t304 = -0.2e1 * t325;
t302 = t218 * t324;
t300 = -0.2e1 * t232 * t316;
t299 = qJD(5) * t313 - t297;
t294 = t216 * t319 - t321;
t234 = t244 * t280 - t301;
t293 = -t234 * t248 + t252 * t315;
t288 = -qJD(4) * t314 + qJD(5) * t296 + t239 * t280;
t255 = t272 * t262;
t237 = t271 * t261 - t274 * t289;
t229 = -t251 * qJD(4) - t256 * t280;
t221 = t276 * t296 + t279 * t313;
t220 = t276 * t313 - t279 * t296;
t211 = -t232 * qJD(4) + t238 * t280;
t207 = 0.1e1 / t209;
t198 = 0.1e1 / t200;
t196 = t224 * t327;
t195 = t293 * t224;
t192 = (-t222 * t243 - t223 * t258) * t277 + t298 * t196;
t191 = t298 * t195 - t222 * t234 + t223 * t252;
t189 = t293 * t326 + (t252 * t300 - t211 * t248 + (t210 * t252 + t228 * t234 + t229 * t232) * t249) * t224;
t187 = t326 * t327 + (t292 * t307 + (-t258 * t300 - t237 * t248 + (-t210 * t258 + t228 * t243 - t232 * t255) * t249) * t277) * t224;
t1 = [0, t187, 0, t189, 0, 0; 0 (-t192 * t323 - t203 * t314) * t305 + ((t239 * t277 + t246 * t307) * t203 + (-t322 + t328) * t192 + (-t314 * t190 + (-t258 * t307 - t187 * t232 - t196 * t210 - t255 * t277 + (-t196 * t251 - t243 * t277) * t194) * t317 + (-t243 * t307 - t187 * t251 - t196 * t228 - t237 * t277 + (t196 * t232 + t258 * t277) * t194) * t318) * t204) * t198, 0 (-t191 * t323 - t203 * t236) * t305 + (t191 * t328 + t213 * t203 + (-t236 * t190 - t191 * t212 + (-t189 * t232 - t195 * t210 + t229 + (-t195 * t251 - t234) * t194) * t317 + (-t189 * t251 - t195 * t228 - t211 + (t195 * t232 - t252) * t194) * t318) * t204) * t198, 0, 0; 0, 0.2e1 * (-t215 * t220 + t221 * t320) * t325 + (0.2e1 * t221 * t302 + t299 * t215 * t279 + t288 * t321 + (t299 * t218 * t276 - t221 * t201 - t220 * t202 - t288 * t319) * t216) * t207, 0, -t294 * t291 * t304 + (t294 * t212 - ((-qJD(5) * t215 - 0.2e1 * t302) * t279 + (t201 * t279 + (t202 - t306) * t276) * t216) * t291) * t207, t304 + 0.2e1 * (t201 * t216 * t207 + (-t207 * t324 - t216 * t325) * t218) * t218, 0;];
JaD_rot  = t1;
