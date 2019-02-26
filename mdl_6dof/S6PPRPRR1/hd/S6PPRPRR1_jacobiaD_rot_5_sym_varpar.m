% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPRPRR1_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_jacobiaD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:39:48
% EndTime: 2019-02-26 19:39:49
% DurationCPUTime: 0.82s
% Computational Cost: add. (4637->79), mult. (13855->173), div. (281->12), fcn. (18065->17), ass. (0->88)
t255 = sin(pkin(13));
t260 = cos(pkin(13));
t266 = sin(qJ(3));
t268 = cos(qJ(3));
t281 = t255 * t268 + t266 * t260;
t298 = qJD(3) * t281;
t251 = t266 * t255 - t260 * t268;
t263 = cos(pkin(7));
t243 = t251 * t263;
t256 = sin(pkin(12));
t257 = sin(pkin(11));
t261 = cos(pkin(12));
t262 = cos(pkin(11));
t264 = cos(pkin(6));
t285 = t262 * t264;
t245 = -t256 * t257 + t261 * t285;
t246 = t256 * t285 + t257 * t261;
t259 = sin(pkin(6));
t258 = sin(pkin(7));
t278 = t251 * t258;
t276 = t259 * t278;
t221 = -t243 * t245 - t246 * t281 + t262 * t276;
t232 = (t243 * t261 + t256 * t281) * t259 + t264 * t278;
t207 = atan2(t221, t232);
t202 = sin(t207);
t203 = cos(t207);
t196 = t202 * t221 + t203 * t232;
t193 = 0.1e1 / t196;
t287 = t257 * t264;
t247 = -t256 * t262 - t261 * t287;
t288 = t257 * t259;
t234 = -t247 * t258 + t263 * t288;
t265 = sin(qJ(5));
t267 = cos(qJ(5));
t242 = t281 * t258;
t244 = t281 * t263;
t248 = -t256 * t287 + t261 * t262;
t275 = t242 * t288 + t244 * t247 - t248 * t251;
t213 = t234 * t265 + t267 * t275;
t209 = 0.1e1 / t213;
t228 = 0.1e1 / t232;
t194 = 0.1e1 / t196 ^ 2;
t210 = 0.1e1 / t213 ^ 2;
t229 = 0.1e1 / t232 ^ 2;
t218 = t221 ^ 2;
t206 = t218 * t229 + 0.1e1;
t204 = 0.1e1 / t206;
t239 = t258 * t298;
t241 = t263 * t298;
t249 = t251 * qJD(3);
t286 = t259 * t262;
t214 = t239 * t286 - t241 * t245 + t246 * t249;
t226 = t239 * t264 + (t241 * t261 - t249 * t256) * t259;
t292 = t221 * t229;
t187 = (t214 * t228 - t226 * t292) * t204;
t282 = -t202 * t232 + t203 * t221;
t184 = t282 * t187 + t202 * t214 + t203 * t226;
t297 = t184 * t193 * t194;
t212 = -t234 * t267 + t265 * t275;
t208 = t212 ^ 2;
t199 = t208 * t210 + 0.1e1;
t238 = t258 * t249;
t240 = t263 * t249;
t217 = -t238 * t288 - t240 * t247 - t248 * t298;
t200 = t213 * qJD(5) + t217 * t265;
t293 = t210 * t212;
t283 = qJD(5) * t212;
t201 = t217 * t267 - t283;
t294 = t201 * t209 * t210;
t296 = (t200 * t293 - t208 * t294) / t199 ^ 2;
t223 = -t247 * t243 - t248 * t281 - t257 * t276;
t295 = t194 * t223;
t231 = t242 * t264 + (t244 * t261 - t251 * t256) * t259;
t291 = t221 * t231;
t290 = t226 * t228 * t229;
t280 = -t209 * t265 + t267 * t293;
t220 = t242 * t286 - t244 * t245 + t246 * t251;
t279 = -t220 * t228 + t229 * t291;
t227 = -t238 * t264 + (-t240 * t261 - t256 * t298) * t259;
t219 = t223 ^ 2;
t216 = -t239 * t288 - t241 * t247 + t248 * t249;
t215 = -t238 * t286 + t240 * t245 + t246 * t298;
t197 = 0.1e1 / t199;
t191 = t194 * t219 + 0.1e1;
t188 = t279 * t204;
t185 = -t282 * t188 + t202 * t220 + t203 * t231;
t183 = 0.2e1 * t279 / t206 ^ 2 * (t214 * t292 - t218 * t290) + (0.2e1 * t290 * t291 + t215 * t228 + (-t214 * t231 - t220 * t226 - t221 * t227) * t229) * t204;
t1 = [0, 0, t183, 0, 0, 0; 0, 0, 0.2e1 * (-t185 * t295 - t193 * t275) / t191 ^ 2 * (t216 * t295 - t219 * t297) + (t217 * t193 + (-t184 * t275 + t185 * t216) * t194 + (-0.2e1 * t185 * t297 + ((t183 * t221 - t188 * t214 + t227 + (t188 * t232 + t220) * t187) * t203 + (-t183 * t232 + t188 * t226 + t215 + (t188 * t221 - t231) * t187) * t202) * t194) * t223) / t191, 0, 0, 0; 0, 0, 0.2e1 * t280 * t223 * t296 + (-t280 * t216 + ((qJD(5) * t209 + 0.2e1 * t212 * t294) * t267 + (-t200 * t267 + (-t201 + t283) * t265) * t210) * t223) * t197, 0, -0.2e1 * t296 + 0.2e1 * (t197 * t200 * t210 + (-t197 * t294 - t210 * t296) * t212) * t212, 0;];
JaD_rot  = t1;
