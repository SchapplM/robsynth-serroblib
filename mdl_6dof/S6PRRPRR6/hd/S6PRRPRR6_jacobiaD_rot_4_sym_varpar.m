% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRR6_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_jacobiaD_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:06
% EndTime: 2019-02-26 20:07:08
% DurationCPUTime: 1.11s
% Computational Cost: add. (4823->122), mult. (15462->242), div. (514->12), fcn. (19568->15), ass. (0->121)
t259 = sin(qJ(2));
t261 = cos(qJ(2));
t323 = cos(pkin(12));
t324 = cos(pkin(6));
t289 = t324 * t323;
t322 = sin(pkin(12));
t277 = -t322 * t259 + t261 * t289;
t242 = t277 * qJD(2);
t254 = sin(pkin(7));
t255 = sin(pkin(6));
t292 = t254 * t255 * t323;
t327 = -qJD(3) * t292 + t242;
t258 = sin(qJ(3));
t276 = -t259 * t289 - t322 * t261;
t271 = qJD(2) * t276;
t260 = cos(qJ(3));
t273 = t260 * t277;
t326 = qJD(3) * t273 + t258 * t271;
t257 = cos(pkin(7));
t272 = t277 * t258;
t268 = t257 * t272 - t260 * t276;
t306 = t257 * t260;
t197 = t268 * qJD(3) + t327 * t258 - t271 * t306;
t313 = t276 * t258;
t221 = -t257 * t273 + t260 * t292 - t313;
t219 = t221 ^ 2;
t302 = t260 * t261;
t305 = t258 * t259;
t282 = t257 * t302 - t305;
t297 = t254 * t324;
t233 = -t282 * t255 - t260 * t297;
t231 = 0.1e1 / t233 ^ 2;
t213 = t219 * t231 + 0.1e1;
t314 = t221 * t231;
t303 = t259 * t260;
t304 = t258 * t261;
t280 = t257 * t304 + t303;
t281 = t257 * t303 + t304;
t290 = qJD(3) * t297;
t217 = t258 * t290 + (t281 * qJD(2) + t280 * qJD(3)) * t255;
t230 = 0.1e1 / t233;
t315 = t217 * t230 * t231;
t325 = -0.2e1 * (t197 * t314 - t219 * t315) / t213 ^ 2;
t214 = atan2(-t221, t233);
t209 = sin(t214);
t210 = cos(t214);
t191 = -t209 * t221 + t210 * t233;
t188 = 0.1e1 / t191;
t288 = t324 * t322;
t275 = t259 * t288 - t323 * t261;
t274 = t323 * t259 + t261 * t288;
t296 = t255 * t322;
t291 = t254 * t296;
t278 = -t257 * t274 + t291;
t225 = t278 * t258 - t260 * t275;
t235 = t254 * t274 + t257 * t296;
t253 = sin(pkin(13));
t256 = cos(pkin(13));
t208 = t225 * t256 + t235 * t253;
t202 = 0.1e1 / t208;
t189 = 0.1e1 / t191 ^ 2;
t203 = 0.1e1 / t208 ^ 2;
t211 = 0.1e1 / t213;
t181 = (-t197 * t230 + t217 * t314) * t211;
t287 = -t209 * t233 - t210 * t221;
t177 = t287 * t181 - t197 * t209 + t210 * t217;
t321 = t177 * t188 * t189;
t312 = t275 * t258;
t224 = -t260 * t291 + t274 * t306 - t312;
t320 = t189 * t224;
t243 = t274 * qJD(2);
t244 = t275 * qJD(2);
t307 = t257 * t258;
t200 = t244 * t307 - t243 * t260 + (t278 * t260 + t312) * qJD(3);
t310 = t253 * t254;
t196 = t200 * t256 - t244 * t310;
t319 = t196 * t202 * t203;
t207 = t225 * t253 - t235 * t256;
t318 = t203 * t207;
t317 = t209 * t224;
t316 = t210 * t224;
t311 = t253 * t202;
t309 = t254 * t256;
t308 = t256 * t207;
t220 = t224 ^ 2;
t187 = t189 * t220 + 0.1e1;
t199 = t225 * qJD(3) - t243 * t258 - t244 * t306;
t301 = 0.2e1 * (t199 * t320 - t220 * t321) / t187 ^ 2;
t201 = t207 ^ 2;
t194 = t201 * t203 + 0.1e1;
t195 = t200 * t253 + t244 * t309;
t300 = 0.2e1 * (t195 * t318 - t201 * t319) / t194 ^ 2;
t299 = t207 * t319;
t298 = qJD(3) * t313;
t294 = -0.2e1 * t221 * t315;
t293 = 0.2e1 * t224 * t321;
t223 = -t258 * t292 + t268;
t234 = t280 * t255 + t258 * t297;
t284 = -t223 * t230 + t234 * t314;
t227 = -t276 * t306 + t272;
t241 = t281 * t255;
t283 = -t227 * t230 + t241 * t314;
t228 = -t258 * t274 - t275 * t306;
t229 = -t260 * t274 + t275 * t307;
t279 = -t257 * t305 + t302;
t226 = (t282 * qJD(2) + t279 * qJD(3)) * t255;
t218 = t260 * t290 + (t279 * qJD(2) + t282 * qJD(3)) * t255;
t216 = t229 * t256 - t275 * t310;
t215 = t229 * t253 + t275 * t309;
t206 = -t228 * qJD(3) + t243 * t307 + t244 * t260;
t205 = t242 * t306 + t257 * t298 + t326;
t198 = t326 * t257 + t327 * t260 + t298;
t192 = 0.1e1 / t194;
t185 = 0.1e1 / t187;
t183 = t283 * t211;
t182 = t284 * t211;
t179 = t287 * t183 - t209 * t227 + t210 * t241;
t178 = t287 * t182 - t209 * t223 + t210 * t234;
t176 = t283 * t325 + (t241 * t294 - t205 * t230 + (t197 * t241 + t217 * t227 + t221 * t226) * t231) * t211;
t175 = t284 * t325 + (t234 * t294 - t198 * t230 + (t197 * t234 + t217 * t223 + t218 * t221) * t231) * t211;
t1 = [0, t176, t175, 0, 0, 0; 0 (t179 * t320 - t188 * t228) * t301 + ((t229 * qJD(3) - t243 * t306 + t244 * t258) * t188 + t179 * t293 + (-t228 * t177 - t179 * t199 - (-t176 * t221 - t183 * t197 + t226 + (-t183 * t233 - t227) * t181) * t316 - (-t176 * t233 - t183 * t217 - t205 + (t183 * t221 - t241) * t181) * t317) * t189) * t185 (t178 * t320 - t188 * t225) * t301 + (t178 * t293 + t200 * t188 + (-t225 * t177 - t178 * t199 - (-t175 * t221 - t182 * t197 + t218 + (-t182 * t233 - t223) * t181) * t316 - (-t175 * t233 - t182 * t217 - t198 + (t182 * t221 - t234) * t181) * t317) * t189) * t185, 0, 0, 0; 0 (-t202 * t215 + t216 * t318) * t300 + ((t206 * t253 + t243 * t309) * t202 + 0.2e1 * t216 * t299 + (-t215 * t196 - (t206 * t256 - t243 * t310) * t207 - t216 * t195) * t203) * t192 (-t203 * t308 + t311) * t224 * t300 + (-0.2e1 * t224 * t256 * t299 - t199 * t311 + (t199 * t308 + (t195 * t256 + t196 * t253) * t224) * t203) * t192, 0, 0, 0;];
JaD_rot  = t1;
