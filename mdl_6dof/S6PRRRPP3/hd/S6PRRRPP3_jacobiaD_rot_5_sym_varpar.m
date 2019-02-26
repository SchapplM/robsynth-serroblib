% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRPP3_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:10:02
% EndTime: 2019-02-26 20:10:04
% DurationCPUTime: 1.64s
% Computational Cost: add. (7193->155), mult. (20987->306), div. (784->12), fcn. (27058->13), ass. (0->122)
t249 = sin(qJ(3));
t252 = cos(qJ(3));
t250 = sin(qJ(2));
t253 = cos(qJ(2));
t309 = cos(pkin(10));
t310 = cos(pkin(6));
t275 = t310 * t309;
t308 = sin(pkin(10));
t262 = -t250 * t275 - t308 * t253;
t247 = sin(pkin(6));
t280 = t247 * t309;
t220 = t249 * t262 - t252 * t280;
t273 = t308 * t250 - t253 * t275;
t236 = t273 * qJD(2);
t203 = qJD(3) * t220 - t236 * t252;
t248 = sin(qJ(4));
t251 = cos(qJ(4));
t264 = t249 * t280 + t252 * t262;
t211 = t248 * t273 - t251 * t264;
t260 = t262 * qJD(2);
t185 = qJD(4) * t211 + t203 * t248 + t251 * t260;
t265 = t273 * t251;
t209 = -t248 * t264 - t265;
t206 = t209 ^ 2;
t294 = t247 * t250;
t243 = t310 * t249 + t252 * t294;
t292 = t251 * t253;
t229 = t243 * t248 + t247 * t292;
t225 = 0.1e1 / t229 ^ 2;
t197 = t206 * t225 + 0.1e1;
t195 = 0.1e1 / t197;
t242 = -t249 * t294 + t310 * t252;
t291 = qJD(2) * t247;
t281 = t253 * t291;
t228 = qJD(3) * t242 + t252 * t281;
t293 = t248 * t253;
t230 = t243 * t251 - t247 * t293;
t282 = t250 * t291;
t199 = qJD(4) * t230 + t228 * t248 - t251 * t282;
t224 = 0.1e1 / t229;
t301 = t209 * t225;
t172 = (-t185 * t224 + t199 * t301) * t195;
t198 = atan2(-t209, t229);
t192 = sin(t198);
t193 = cos(t198);
t272 = -t192 * t229 - t193 * t209;
t168 = t172 * t272 - t192 * t185 + t193 * t199;
t184 = -t192 * t209 + t193 * t229;
t181 = 0.1e1 / t184;
t182 = 0.1e1 / t184 ^ 2;
t313 = t168 * t181 * t182;
t274 = t310 * t308;
t241 = -t250 * t274 + t309 * t253;
t279 = t247 * t308;
t223 = t241 * t252 + t249 * t279;
t261 = -t309 * t250 - t253 * t274;
t212 = t223 * t248 + t251 * t261;
t276 = 0.2e1 * t212 * t313;
t267 = -t220 * t224 + t242 * t301;
t312 = t248 * t267;
t303 = t199 * t224 * t225;
t311 = -0.2e1 * (t185 * t301 - t206 * t303) / t197 ^ 2;
t263 = -t241 * t249 + t252 * t279;
t217 = 0.1e1 / t263;
t218 = 0.1e1 / t263 ^ 2;
t307 = t182 * t212;
t237 = t261 * qJD(2);
t205 = qJD(3) * t263 + t237 * t252;
t213 = t223 * t251 - t248 * t261;
t238 = t241 * qJD(2);
t187 = qJD(4) * t213 + t205 * t248 - t238 * t251;
t306 = t187 * t182;
t305 = t192 * t212;
t304 = t193 * t212;
t204 = qJD(3) * t223 + t237 * t249;
t219 = t217 * t218;
t302 = t204 * t219;
t300 = t213 * t218;
t299 = t213 * t223;
t298 = t217 * t263;
t297 = t263 * t248;
t296 = t261 * t249;
t295 = t261 * t252;
t290 = qJD(2) * t252;
t289 = qJD(3) * t249;
t288 = qJD(4) * t248;
t287 = qJD(4) * t251;
t286 = t252 * qJD(4);
t207 = t212 ^ 2;
t180 = t207 * t182 + 0.1e1;
t285 = 0.2e1 * (-t207 * t313 + t212 * t306) / t180 ^ 2;
t188 = -qJD(4) * t212 + t205 * t251 + t238 * t248;
t208 = t213 ^ 2;
t194 = t208 * t218 + 0.1e1;
t284 = 0.2e1 * (t188 * t300 + t208 * t302) / t194 ^ 2;
t277 = -0.2e1 * t209 * t303;
t271 = qJD(4) * t241 - t238 * t252;
t269 = -t211 * t224 + t230 * t301;
t266 = t252 * t273;
t214 = -t248 * t266 + t251 * t262;
t232 = (-t250 * t251 + t252 * t293) * t247;
t268 = -t214 * t224 + t232 * t301;
t227 = -qJD(3) * t243 - t249 * t281;
t216 = t241 * t248 + t251 * t295;
t215 = -t241 * t251 + t248 * t295;
t202 = qJD(3) * t264 + t236 * t249;
t201 = ((-qJD(2) + t286) * t292 + (-t253 * t289 + (qJD(4) - t290) * t250) * t248) * t247;
t200 = -qJD(4) * t229 + t228 * t251 + t248 * t282;
t190 = 0.1e1 / t194;
t189 = t236 * t251 - t262 * t288 - t266 * t287 + (t262 * t290 + t273 * t289) * t248;
t186 = qJD(4) * t265 + t203 * t251 - t248 * t260 + t264 * t288;
t178 = 0.1e1 / t180;
t177 = t195 * t312;
t176 = t268 * t195;
t174 = t269 * t195;
t171 = (-t192 * t220 + t193 * t242) * t248 + t272 * t177;
t170 = t176 * t272 - t192 * t214 + t193 * t232;
t169 = t174 * t272 - t192 * t211 + t193 * t230;
t167 = t268 * t311 + (t232 * t277 - t189 * t224 + (t185 * t232 + t199 * t214 + t201 * t209) * t225) * t195;
t165 = t269 * t311 + (t230 * t277 - t186 * t224 + (t185 * t230 + t199 * t211 + t200 * t209) * t225) * t195;
t164 = t311 * t312 + (t267 * t287 + (t242 * t277 - t202 * t224 + (t185 * t242 + t199 * t220 + t209 * t227) * t225) * t248) * t195;
t1 = [0, t167, t164, t165, 0, 0; 0 (t170 * t307 - t181 * t215) * t285 + (t170 * t276 + (-t215 * t168 - t170 * t187 - (-t167 * t209 - t176 * t185 + t201 + (-t176 * t229 - t214) * t172) * t304 - (-t167 * t229 - t176 * t199 - t189 + (t176 * t209 - t232) * t172) * t305) * t182 + ((t261 * t286 - t237) * t251 + (-t261 * t289 + t271) * t248) * t181) * t178 (t171 * t307 - t181 * t297) * t285 + ((-t204 * t248 + t263 * t287) * t181 + (-t306 + t276) * t171 + (-t297 * t168 - (t242 * t287 - t164 * t209 - t177 * t185 + t227 * t248 + (-t177 * t229 - t220 * t248) * t172) * t304 - (-t220 * t287 - t164 * t229 - t177 * t199 - t202 * t248 + (t177 * t209 - t242 * t248) * t172) * t305) * t182) * t178 (t169 * t307 - t181 * t213) * t285 + (t169 * t276 + t188 * t181 + (-t213 * t168 - t169 * t187 - (-t165 * t209 - t174 * t185 + t200 + (-t174 * t229 - t211) * t172) * t304 - (-t165 * t229 - t174 * t199 - t186 + (t174 * t209 - t230) * t172) * t305) * t182) * t178, 0, 0; 0 (t216 * t217 + t296 * t300) * t284 + (-(t237 * t248 + t251 * t271) * t217 - (-(t248 * t286 + t251 * t289) * t217 + 0.2e1 * t249 * t213 * t302) * t261 + (-t188 * t296 - t216 * t204 + (-qJD(3) * t295 + t238 * t249) * t213) * t218) * t190 (t218 * t299 + t251 * t298) * t284 + (t288 * t298 + (-t188 * t223 - t205 * t213) * t218 + (-0.2e1 * t219 * t299 + (-t218 * t263 + t217) * t251) * t204) * t190, -t212 * t217 * t284 + (t204 * t212 * t218 + t187 * t217) * t190, 0, 0;];
JaD_rot  = t1;
