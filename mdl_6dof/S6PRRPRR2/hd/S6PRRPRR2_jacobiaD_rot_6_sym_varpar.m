% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRR2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:04:42
% EndTime: 2019-02-26 20:04:43
% DurationCPUTime: 0.99s
% Computational Cost: add. (6552->112), mult. (9671->227), div. (577->12), fcn. (12382->13), ass. (0->111)
t257 = sin(pkin(11));
t259 = cos(pkin(11));
t261 = sin(qJ(2));
t260 = cos(pkin(6));
t262 = cos(qJ(2));
t289 = t260 * t262;
t242 = -t257 * t261 + t259 * t289;
t238 = t242 * qJD(2);
t290 = t260 * t261;
t243 = t257 * t262 + t259 * t290;
t255 = qJ(3) + pkin(12);
t250 = sin(t255);
t258 = sin(pkin(6));
t293 = t258 * t259;
t280 = t250 * t293;
t251 = cos(t255);
t286 = qJD(3) * t251;
t209 = -qJD(3) * t280 + t238 * t250 + t243 * t286;
t227 = t243 * t250 + t251 * t293;
t225 = t227 ^ 2;
t292 = t258 * t261;
t236 = t250 * t292 - t260 * t251;
t234 = 0.1e1 / t236 ^ 2;
t219 = t225 * t234 + 0.1e1;
t217 = 0.1e1 / t219;
t237 = t260 * t250 + t251 * t292;
t287 = qJD(2) * t262;
t279 = t258 * t287;
t223 = t237 * qJD(3) + t250 * t279;
t233 = 0.1e1 / t236;
t299 = t227 * t234;
t189 = (-t209 * t233 + t223 * t299) * t217;
t220 = atan2(-t227, t236);
t215 = sin(t220);
t216 = cos(t220);
t274 = -t215 * t236 - t216 * t227;
t185 = t274 * t189 - t215 * t209 + t216 * t223;
t199 = -t215 * t227 + t216 * t236;
t196 = 0.1e1 / t199;
t197 = 0.1e1 / t199 ^ 2;
t313 = t185 * t196 * t197;
t281 = t257 * t290;
t245 = t259 * t262 - t281;
t294 = t257 * t258;
t271 = -t245 * t250 + t251 * t294;
t312 = -0.2e1 * t271 * t313;
t291 = t258 * t262;
t270 = -t233 * t242 + t291 * t299;
t311 = t250 * t270;
t300 = t223 * t233 * t234;
t310 = -0.2e1 * (t209 * t299 - t225 * t300) / t219 ^ 2;
t231 = t245 * t251 + t250 * t294;
t256 = qJ(5) + qJ(6);
t253 = cos(t256);
t244 = t257 * t289 + t259 * t261;
t252 = sin(t256);
t297 = t244 * t252;
t214 = t231 * t253 + t297;
t206 = 0.1e1 / t214;
t207 = 0.1e1 / t214 ^ 2;
t241 = -qJD(2) * t281 + t259 * t287;
t254 = qJD(5) + qJD(6);
t276 = t231 * t254 - t241;
t240 = t244 * qJD(2);
t212 = t271 * qJD(3) - t240 * t251;
t295 = t244 * t254;
t277 = t212 + t295;
t200 = t277 * t252 + t276 * t253;
t296 = t244 * t253;
t213 = t231 * t252 - t296;
t205 = t213 ^ 2;
t204 = t205 * t207 + 0.1e1;
t305 = t207 * t213;
t201 = -t276 * t252 + t277 * t253;
t307 = t201 * t206 * t207;
t309 = (t200 * t305 - t205 * t307) / t204 ^ 2;
t308 = t197 * t271;
t306 = t206 * t252;
t211 = t231 * qJD(3) - t240 * t250;
t304 = t211 * t197;
t303 = t213 * t253;
t302 = t215 * t271;
t301 = t216 * t271;
t298 = t244 * t250;
t288 = qJD(2) * t261;
t226 = t271 ^ 2;
t195 = t226 * t197 + 0.1e1;
t285 = 0.2e1 * (-t226 * t313 - t271 * t304) / t195 ^ 2;
t284 = -0.2e1 * t309;
t282 = t213 * t307;
t278 = -0.2e1 * t227 * t300;
t275 = t251 * t295 - t240;
t273 = t207 * t303 - t306;
t229 = t243 * t251 - t280;
t272 = -t229 * t233 + t237 * t299;
t269 = qJD(3) * t298 - t241 * t251 + t245 * t254;
t239 = t243 * qJD(2);
t224 = -t236 * qJD(3) + t251 * t279;
t222 = t245 * t252 - t251 * t296;
t221 = -t245 * t253 - t251 * t297;
t210 = -t227 * qJD(3) + t238 * t251;
t202 = 0.1e1 / t204;
t192 = 0.1e1 / t195;
t191 = t217 * t311;
t190 = t272 * t217;
t187 = (-t215 * t242 + t216 * t291) * t250 + t274 * t191;
t186 = t274 * t190 - t215 * t229 + t216 * t237;
t184 = t272 * t310 + (t237 * t278 - t210 * t233 + (t209 * t237 + t223 * t229 + t224 * t227) * t234) * t217;
t182 = t310 * t311 + (t270 * t286 + (t278 * t291 + t233 * t239 + (t223 * t242 + (t209 * t262 - t227 * t288) * t258) * t234) * t250) * t217;
t181 = t284 + 0.2e1 * (t200 * t207 * t202 + (-t202 * t307 - t207 * t309) * t213) * t213;
t1 = [0, t182, t184, 0, 0, 0; 0 (-t187 * t308 + t196 * t298) * t285 + ((-t241 * t250 - t244 * t286) * t196 + (-t304 + t312) * t187 + (t298 * t185 + (-t182 * t227 - t191 * t209 + (-t250 * t288 + t262 * t286) * t258 + (-t191 * t236 - t242 * t250) * t189) * t301 + (-t242 * t286 - t182 * t236 - t191 * t223 + t239 * t250 + (t191 * t227 - t250 * t291) * t189) * t302) * t197) * t192 (-t186 * t308 - t196 * t231) * t285 + (t186 * t312 + t212 * t196 + (-t231 * t185 - t186 * t211 + (-t184 * t227 - t190 * t209 + t224 + (-t190 * t236 - t229) * t189) * t301 + (-t184 * t236 - t190 * t223 - t210 + (t190 * t227 - t237) * t189) * t302) * t197) * t192, 0, 0, 0; 0, 0.2e1 * (-t206 * t221 + t222 * t305) * t309 + (0.2e1 * t222 * t282 - t275 * t206 * t253 + t269 * t306 + (-t275 * t213 * t252 - t222 * t200 - t221 * t201 - t269 * t303) * t207) * t202, -t273 * t271 * t284 + (t273 * t211 - ((-t206 * t254 - 0.2e1 * t282) * t253 + (t200 * t253 + (-t213 * t254 + t201) * t252) * t207) * t271) * t202, 0, t181, t181;];
JaD_rot  = t1;
