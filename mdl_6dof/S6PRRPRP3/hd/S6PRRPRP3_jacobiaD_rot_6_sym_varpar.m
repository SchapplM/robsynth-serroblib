% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRP3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:02:23
% EndTime: 2019-02-26 20:02:25
% DurationCPUTime: 1.71s
% Computational Cost: add. (11295->157), mult. (21125->308), div. (784->12), fcn. (27225->13), ass. (0->123)
t260 = pkin(11) + qJ(5);
t258 = sin(t260);
t259 = cos(t260);
t262 = sin(qJ(3));
t264 = cos(qJ(3));
t263 = sin(qJ(2));
t265 = cos(qJ(2));
t321 = cos(pkin(10));
t322 = cos(pkin(6));
t286 = t322 * t321;
t320 = sin(pkin(10));
t274 = -t263 * t286 - t320 * t265;
t261 = sin(pkin(6));
t291 = t261 * t321;
t276 = t262 * t291 + t264 * t274;
t284 = t320 * t263 - t265 * t286;
t221 = t284 * t258 - t259 * t276;
t237 = t262 * t274 - t264 * t291;
t248 = t284 * qJD(2);
t225 = t237 * qJD(3) - t248 * t264;
t272 = qJD(2) * t274;
t197 = t221 * qJD(5) + t225 * t258 + t259 * t272;
t277 = t284 * t259;
t219 = -t258 * t276 - t277;
t214 = t219 ^ 2;
t306 = t261 * t263;
t255 = t322 * t262 + t264 * t306;
t305 = t261 * t265;
t234 = t255 * t258 + t259 * t305;
t232 = 0.1e1 / t234 ^ 2;
t206 = t214 * t232 + 0.1e1;
t204 = 0.1e1 / t206;
t235 = t255 * t259 - t258 * t305;
t254 = -t262 * t306 + t322 * t264;
t304 = qJD(2) * t261;
t292 = t265 * t304;
t242 = t254 * qJD(3) + t264 * t292;
t293 = t263 * t304;
t211 = t235 * qJD(5) + t242 * t258 - t259 * t293;
t231 = 0.1e1 / t234;
t311 = t219 * t232;
t184 = (-t197 * t231 + t211 * t311) * t204;
t207 = atan2(-t219, t234);
t202 = sin(t207);
t203 = cos(t207);
t283 = -t202 * t234 - t203 * t219;
t180 = t283 * t184 - t202 * t197 + t203 * t211;
t196 = -t202 * t219 + t203 * t234;
t193 = 0.1e1 / t196;
t194 = 0.1e1 / t196 ^ 2;
t326 = t180 * t193 * t194;
t285 = t322 * t320;
t273 = -t321 * t263 - t265 * t285;
t253 = -t263 * t285 + t321 * t265;
t290 = t261 * t320;
t275 = -t253 * t264 - t262 * t290;
t222 = -t258 * t275 + t259 * t273;
t325 = -0.2e1 * t222;
t288 = 0.2e1 * t222 * t326;
t279 = -t231 * t237 + t254 * t311;
t324 = t258 * t279;
t314 = t211 * t231 * t232;
t323 = -0.2e1 * (t197 * t311 - t214 * t314) / t206 ^ 2;
t223 = -t258 * t273 - t259 * t275;
t216 = 0.1e1 / t223;
t217 = 0.1e1 / t223 ^ 2;
t239 = -t253 * t262 + t264 * t290;
t236 = t239 ^ 2;
t313 = t217 * t236;
t210 = 0.1e1 + t313;
t249 = t273 * qJD(2);
t227 = t239 * qJD(3) + t249 * t264;
t250 = t253 * qJD(2);
t200 = -t222 * qJD(5) + t227 * t259 + t250 * t258;
t317 = t200 * t216 * t217;
t295 = t236 * t317;
t226 = t275 * qJD(3) - t249 * t262;
t310 = t226 * t239;
t319 = (t217 * t310 - t295) / t210 ^ 2;
t318 = t194 * t222;
t316 = t202 * t222;
t315 = t203 * t222;
t312 = t217 * t239;
t309 = t239 * t258;
t308 = t258 * t264;
t307 = t259 * t264;
t303 = qJD(2) * t264;
t302 = qJD(3) * t262;
t301 = qJD(5) * t258;
t300 = qJD(5) * t259;
t299 = qJD(5) * t264;
t215 = t222 ^ 2;
t192 = t194 * t215 + 0.1e1;
t199 = t223 * qJD(5) + t227 * t258 - t250 * t259;
t298 = 0.2e1 * (t199 * t318 - t215 * t326) / t192 ^ 2;
t297 = 0.2e1 * t319;
t294 = t239 * t317;
t287 = -0.2e1 * t219 * t314;
t281 = -t221 * t231 + t235 * t311;
t278 = t264 * t284;
t228 = -t258 * t278 + t259 * t274;
t243 = (-t259 * t263 + t265 * t308) * t261;
t280 = -t228 * t231 + t243 * t311;
t241 = -t255 * qJD(3) - t262 * t292;
t230 = t253 * t258 + t273 * t307;
t229 = -t253 * t259 + t273 * t308;
t224 = t276 * qJD(3) + t248 * t262;
t213 = ((-qJD(2) + t299) * t265 * t259 + (-t265 * t302 + (qJD(5) - t303) * t263) * t258) * t261;
t212 = -t234 * qJD(5) + t242 * t259 + t258 * t293;
t208 = 0.1e1 / t210;
t201 = t248 * t259 - t274 * t301 - t278 * t300 + (t274 * t303 + t284 * t302) * t258;
t198 = qJD(5) * t277 + t225 * t259 - t258 * t272 + t276 * t301;
t190 = 0.1e1 / t192;
t189 = t204 * t324;
t187 = t280 * t204;
t186 = t281 * t204;
t183 = (-t202 * t237 + t203 * t254) * t258 + t283 * t189;
t182 = t283 * t187 - t202 * t228 + t203 * t243;
t181 = t283 * t186 - t202 * t221 + t203 * t235;
t179 = t280 * t323 + (t243 * t287 - t201 * t231 + (t197 * t243 + t211 * t228 + t213 * t219) * t232) * t204;
t177 = t281 * t323 + (t235 * t287 - t198 * t231 + (t197 * t235 + t211 * t221 + t212 * t219) * t232) * t204;
t176 = t323 * t324 + (t279 * t300 + (t254 * t287 - t224 * t231 + (t197 * t254 + t211 * t237 + t219 * t241) * t232) * t258) * t204;
t1 = [0, t179, t176, 0, t177, 0; 0 (t182 * t318 - t193 * t229) * t298 + (t182 * t288 + (-t229 * t180 - t182 * t199 - (-t179 * t219 - t187 * t197 + t213 + (-t187 * t234 - t228) * t184) * t315 - (-t179 * t234 - t187 * t211 - t201 + (t187 * t219 - t243) * t184) * t316) * t194 + ((t273 * t299 - t249) * t259 + (qJD(5) * t253 - t250 * t264 - t273 * t302) * t258) * t193) * t190 (t183 * t318 - t193 * t309) * t298 + ((t226 * t258 + t239 * t300) * t193 + t183 * t288 + (-t183 * t199 - t309 * t180 - (t254 * t300 - t176 * t219 - t189 * t197 + t241 * t258 + (-t189 * t234 - t237 * t258) * t184) * t315 - (-t237 * t300 - t176 * t234 - t189 * t211 - t224 * t258 + (t189 * t219 - t254 * t258) * t184) * t316) * t194) * t190, 0 (t181 * t318 - t193 * t223) * t298 + (t181 * t288 + t200 * t193 + (-t223 * t180 - t181 * t199 - (-t177 * t219 - t186 * t197 + t212 + (-t186 * t234 - t221) * t184) * t315 - (-t177 * t234 - t186 * t211 - t198 + (t186 * t219 - t235) * t184) * t316) * t194) * t190, 0; 0 (t216 * t262 * t273 + t230 * t312) * t297 + (0.2e1 * t230 * t294 + (-qJD(3) * t264 * t273 + t250 * t262) * t216 + (-(t249 * t258 - t250 * t307 + t253 * t300) * t239 - t230 * t226 - (-t262 * t200 - (t258 * t299 + t259 * t302) * t239) * t273) * t217) * t208 (-t216 * t275 + t259 * t313) * t297 + (0.2e1 * t259 * t295 - t216 * t227 + (-t200 * t275 + t236 * t301 - 0.2e1 * t259 * t310) * t217) * t208, 0, t312 * t319 * t325 + (t294 * t325 + (t199 * t239 + t222 * t226) * t217) * t208, 0;];
JaD_rot  = t1;
