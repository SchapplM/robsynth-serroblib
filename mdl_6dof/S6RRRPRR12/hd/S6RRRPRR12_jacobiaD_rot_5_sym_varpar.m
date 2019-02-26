% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR12_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR12_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:22:26
% EndTime: 2019-02-26 22:22:27
% DurationCPUTime: 1.42s
% Computational Cost: add. (4943->149), mult. (13478->297), div. (726->12), fcn. (17045->13), ass. (0->127)
t262 = cos(pkin(6));
t264 = sin(qJ(2));
t335 = sin(qJ(1));
t297 = t335 * t264;
t287 = t262 * t297;
t291 = t335 * qJD(2);
t266 = cos(qJ(2));
t267 = cos(qJ(1));
t313 = t267 * t266;
t261 = sin(pkin(6));
t315 = t261 * t267;
t340 = -qJD(1) * t287 - t264 * t291 + (qJD(2) * t262 + qJD(1)) * t313 - qJD(3) * t315;
t263 = sin(qJ(3));
t265 = cos(qJ(3));
t296 = t335 * t266;
t314 = t267 * t264;
t279 = -t262 * t314 - t296;
t231 = -t263 * t279 + t265 * t315;
t317 = t261 * t264;
t242 = -t262 * t265 + t263 * t317;
t220 = atan2(-t231, t242);
t215 = sin(t220);
t216 = cos(t220);
t198 = -t215 * t231 + t216 * t242;
t196 = 0.1e1 / t198 ^ 2;
t247 = -t287 + t313;
t298 = t261 * t335;
t278 = -t247 * t263 + t265 * t298;
t228 = t278 ^ 2;
t194 = t228 * t196 + 0.1e1;
t277 = -t262 * t296 - t314;
t224 = t279 * qJD(1) + t277 * qJD(2);
t237 = t247 * t265 + t263 * t298;
t295 = qJD(1) * t315;
t202 = t237 * qJD(3) + t224 * t263 - t265 * t295;
t328 = t202 * t196;
t227 = t231 ^ 2;
t240 = 0.1e1 / t242 ^ 2;
t219 = t227 * t240 + 0.1e1;
t217 = 0.1e1 / t219;
t292 = t335 * qJD(1);
t286 = t261 * t292;
t310 = qJD(3) * t265;
t204 = t263 * t340 - t265 * t286 - t279 * t310;
t243 = t262 * t263 + t265 * t317;
t311 = qJD(2) * t266;
t294 = t261 * t311;
t229 = t243 * qJD(3) + t263 * t294;
t239 = 0.1e1 / t242;
t322 = t231 * t240;
t283 = -t204 * t239 + t229 * t322;
t186 = t283 * t217;
t284 = -t215 * t242 - t216 * t231;
t181 = t284 * t186 - t215 * t204 + t216 * t229;
t195 = 0.1e1 / t198;
t197 = t195 * t196;
t333 = t181 * t197;
t308 = 0.2e1 * (-t228 * t333 - t278 * t328) / t194 ^ 2;
t339 = t229 * t240;
t299 = t262 * t313;
t244 = -t297 + t299;
t316 = t261 * t266;
t280 = -t239 * t244 + t316 * t322;
t338 = t263 * t280;
t205 = (qJD(3) * t279 + t286) * t263 + t340 * t265;
t260 = pkin(12) + qJ(5);
t258 = sin(t260);
t259 = cos(t260);
t214 = t237 * t259 - t258 * t277;
t208 = 0.1e1 / t214;
t209 = 0.1e1 / t214 ^ 2;
t337 = -0.2e1 * t231;
t336 = -0.2e1 * t278;
t203 = t278 * qJD(3) + t224 * t265 + t263 * t295;
t223 = -qJD(1) * t299 - t267 * t311 + (t262 * t291 + t292) * t264;
t189 = t214 * qJD(5) + t203 * t258 + t223 * t259;
t213 = t237 * t258 + t259 * t277;
t207 = t213 ^ 2;
t201 = t207 * t209 + 0.1e1;
t327 = t209 * t213;
t309 = qJD(5) * t213;
t190 = t203 * t259 - t223 * t258 - t309;
t330 = t190 * t208 * t209;
t332 = (t189 * t327 - t207 * t330) / t201 ^ 2;
t324 = t239 * t339;
t331 = (t204 * t322 - t227 * t324) / t219 ^ 2;
t329 = t196 * t278;
t326 = t215 * t278;
t325 = t216 * t278;
t323 = t231 * t239;
t321 = t277 * t263;
t320 = t277 * t265;
t319 = t258 * t208;
t318 = t259 * t213;
t312 = qJD(2) * t264;
t307 = -0.2e1 * t332;
t306 = 0.2e1 * t332;
t305 = -0.2e1 * t331;
t304 = t197 * t336;
t303 = t239 * t331;
t302 = t196 * t326;
t301 = t196 * t325;
t300 = t213 * t330;
t290 = 0.2e1 * t300;
t289 = t324 * t337;
t233 = -t263 * t315 - t265 * t279;
t285 = -qJD(5) * t320 + t224;
t212 = -t233 * t259 + t244 * t258;
t211 = -t233 * t258 - t244 * t259;
t282 = t209 * t318 - t319;
t281 = -t233 * t239 + t243 * t322;
t275 = -t215 + (t216 * t323 + t215) * t217;
t274 = -qJD(3) * t321 + qJD(5) * t247 + t223 * t265;
t230 = -t242 * qJD(3) + t265 * t294;
t225 = t277 * qJD(1) + t279 * qJD(2);
t222 = t247 * t258 + t259 * t320;
t221 = -t247 * t259 + t258 * t320;
t199 = 0.1e1 / t201;
t192 = 0.1e1 / t194;
t191 = t217 * t338;
t188 = t281 * t217;
t185 = t275 * t278;
t183 = (-t215 * t244 + t216 * t316) * t263 + t284 * t191;
t182 = t284 * t188 - t215 * t233 + t216 * t243;
t180 = t281 * t305 + (t243 * t289 - t205 * t239 + (t204 * t243 + t229 * t233 + t230 * t231) * t240) * t217;
t178 = t305 * t338 + (t280 * t310 + (t289 * t316 - t225 * t239 + (t229 * t244 + (t204 * t266 - t231 * t312) * t261) * t240) * t263) * t217;
t1 = [t303 * t336 + (-t202 * t239 - t278 * t339) * t217, t178, t180, 0, 0, 0; t231 * t195 * t308 + (-t204 * t195 + (t181 * t231 + t185 * t202) * t196) * t192 - (-t185 * t196 * t308 + (-0.2e1 * t185 * t333 + (-t186 * t217 * t323 + t305) * t302 + (t303 * t337 - t186 + (t186 - t283) * t217) * t301 - t275 * t328) * t192) * t278 (-t183 * t329 - t195 * t321) * t308 + (-t183 * t328 + (t223 * t263 + t277 * t310) * t195 + (t183 * t304 - t196 * t321) * t181 + (-t178 * t231 - t191 * t204 + (-t263 * t312 + t266 * t310) * t261 + (-t191 * t242 - t244 * t263) * t186) * t301 + (-t244 * t310 - t178 * t242 - t191 * t229 - t225 * t263 + (t191 * t231 - t263 * t316) * t186) * t302) * t192 (-t182 * t329 - t195 * t237) * t308 + (t182 * t181 * t304 + t203 * t195 + (-t237 * t181 - t182 * t202 + (-t180 * t231 - t188 * t204 + t230 + (-t188 * t242 - t233) * t186) * t325 + (-t180 * t242 - t188 * t229 - t205 + (t188 * t231 - t243) * t186) * t326) * t196) * t192, 0, 0, 0; (-t208 * t211 + t212 * t327) * t306 + ((t212 * qJD(5) - t205 * t258 - t225 * t259) * t208 + t212 * t290 + (-t211 * t190 - (-t211 * qJD(5) - t205 * t259 + t225 * t258) * t213 - t212 * t189) * t209) * t199 (-t208 * t221 + t222 * t327) * t306 + (t222 * t290 - t285 * t208 * t259 + t274 * t319 + (-t285 * t213 * t258 - t222 * t189 - t221 * t190 - t274 * t318) * t209) * t199, -t282 * t278 * t307 + (t282 * t202 - ((-qJD(5) * t208 - 0.2e1 * t300) * t259 + (t189 * t259 + (t190 - t309) * t258) * t209) * t278) * t199, 0, t307 + 0.2e1 * (t189 * t209 * t199 + (-t199 * t330 - t209 * t332) * t213) * t213, 0;];
JaD_rot  = t1;
