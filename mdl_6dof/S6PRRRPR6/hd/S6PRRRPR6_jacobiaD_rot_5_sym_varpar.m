% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRPR6_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:13:25
% EndTime: 2019-02-26 20:13:27
% DurationCPUTime: 1.64s
% Computational Cost: add. (7242->156), mult. (21125->307), div. (784->12), fcn. (27225->13), ass. (0->123)
t251 = sin(qJ(3));
t254 = cos(qJ(3));
t252 = sin(qJ(2));
t255 = cos(qJ(2));
t312 = cos(pkin(11));
t313 = cos(pkin(6));
t276 = t313 * t312;
t311 = sin(pkin(11));
t264 = -t252 * t276 - t311 * t255;
t249 = sin(pkin(6));
t281 = t249 * t312;
t223 = t251 * t264 - t254 * t281;
t274 = t311 * t252 - t255 * t276;
t239 = t274 * qJD(2);
t206 = qJD(3) * t223 - t239 * t254;
t250 = sin(qJ(4));
t253 = cos(qJ(4));
t266 = t251 * t281 + t254 * t264;
t216 = t274 * t250 - t253 * t266;
t262 = t264 * qJD(2);
t188 = qJD(4) * t216 + t206 * t250 + t253 * t262;
t267 = t274 * t253;
t214 = -t250 * t266 - t267;
t209 = t214 ^ 2;
t298 = t249 * t252;
t246 = t313 * t251 + t254 * t298;
t295 = t253 * t255;
t232 = t246 * t250 + t249 * t295;
t228 = 0.1e1 / t232 ^ 2;
t200 = t209 * t228 + 0.1e1;
t198 = 0.1e1 / t200;
t245 = -t251 * t298 + t313 * t254;
t294 = qJD(2) * t249;
t282 = t255 * t294;
t231 = qJD(3) * t245 + t254 * t282;
t297 = t250 * t255;
t233 = t246 * t253 - t249 * t297;
t283 = t252 * t294;
t202 = qJD(4) * t233 + t231 * t250 - t253 * t283;
t227 = 0.1e1 / t232;
t302 = t214 * t228;
t175 = (-t188 * t227 + t202 * t302) * t198;
t201 = atan2(-t214, t232);
t195 = sin(t201);
t196 = cos(t201);
t273 = -t195 * t232 - t196 * t214;
t171 = t273 * t175 - t195 * t188 + t196 * t202;
t187 = -t195 * t214 + t196 * t232;
t184 = 0.1e1 / t187;
t185 = 0.1e1 / t187 ^ 2;
t317 = t171 * t184 * t185;
t275 = t313 * t311;
t263 = -t312 * t252 - t255 * t275;
t244 = -t252 * t275 + t312 * t255;
t280 = t249 * t311;
t265 = -t244 * t254 - t251 * t280;
t217 = -t250 * t265 + t253 * t263;
t316 = -0.2e1 * t217;
t277 = 0.2e1 * t217 * t317;
t269 = -t223 * t227 + t245 * t302;
t315 = t250 * t269;
t304 = t202 * t227 * t228;
t314 = -0.2e1 * (t188 * t302 - t209 * t304) / t200 ^ 2;
t218 = -t250 * t263 - t253 * t265;
t211 = 0.1e1 / t218;
t212 = 0.1e1 / t218 ^ 2;
t225 = -t244 * t251 + t254 * t280;
t222 = t225 ^ 2;
t301 = t222 * t212;
t197 = 0.1e1 + t301;
t240 = t263 * qJD(2);
t207 = qJD(3) * t265 - t240 * t251;
t208 = qJD(3) * t225 + t240 * t254;
t241 = t244 * qJD(2);
t191 = -t217 * qJD(4) + t208 * t253 + t241 * t250;
t307 = t191 * t211 * t212;
t285 = t222 * t307;
t303 = t212 * t225;
t310 = (t207 * t303 - t285) / t197 ^ 2;
t309 = t185 * t217;
t190 = t218 * qJD(4) + t208 * t250 - t241 * t253;
t308 = t190 * t185;
t306 = t195 * t217;
t305 = t196 * t217;
t300 = t225 * t250;
t299 = t263 * t254;
t296 = t253 * t254;
t293 = qJD(2) * t254;
t292 = qJD(3) * t251;
t291 = qJD(4) * t250;
t290 = qJD(4) * t253;
t289 = t254 * qJD(4);
t210 = t217 ^ 2;
t183 = t210 * t185 + 0.1e1;
t288 = 0.2e1 * (-t210 * t317 + t217 * t308) / t183 ^ 2;
t287 = 0.2e1 * t310;
t284 = t225 * t307;
t278 = -0.2e1 * t214 * t304;
t271 = -t216 * t227 + t233 * t302;
t268 = t254 * t274;
t219 = -t250 * t268 + t253 * t264;
t235 = (-t252 * t253 + t254 * t297) * t249;
t270 = -t219 * t227 + t235 * t302;
t230 = -qJD(3) * t246 - t251 * t282;
t221 = t244 * t250 + t263 * t296;
t220 = -t244 * t253 + t250 * t299;
t205 = qJD(3) * t266 + t239 * t251;
t204 = ((-qJD(2) + t289) * t295 + (-t255 * t292 + (qJD(4) - t293) * t252) * t250) * t249;
t203 = -qJD(4) * t232 + t231 * t253 + t250 * t283;
t193 = 0.1e1 / t197;
t192 = t239 * t253 - t264 * t291 - t268 * t290 + (t264 * t293 + t274 * t292) * t250;
t189 = qJD(4) * t267 + t206 * t253 - t250 * t262 + t266 * t291;
t181 = 0.1e1 / t183;
t180 = t198 * t315;
t179 = t270 * t198;
t177 = t271 * t198;
t174 = (-t195 * t223 + t196 * t245) * t250 + t273 * t180;
t173 = t273 * t179 - t195 * t219 + t196 * t235;
t172 = t273 * t177 - t195 * t216 + t196 * t233;
t170 = t270 * t314 + (t235 * t278 - t192 * t227 + (t188 * t235 + t202 * t219 + t204 * t214) * t228) * t198;
t168 = t271 * t314 + (t233 * t278 - t189 * t227 + (t188 * t233 + t202 * t216 + t203 * t214) * t228) * t198;
t167 = t314 * t315 + (t269 * t290 + (t245 * t278 - t205 * t227 + (t188 * t245 + t202 * t223 + t214 * t230) * t228) * t250) * t198;
t1 = [0, t170, t167, t168, 0, 0; 0 (t173 * t309 - t184 * t220) * t288 + (t173 * t277 + (-t220 * t171 - t173 * t190 - (-t170 * t214 - t179 * t188 + t204 + (-t179 * t232 - t219) * t175) * t305 - (-t170 * t232 - t179 * t202 - t192 + (t179 * t214 - t235) * t175) * t306) * t185 + ((t263 * t289 - t240) * t253 + (qJD(4) * t244 - t241 * t254 - t263 * t292) * t250) * t184) * t181 (t174 * t309 - t184 * t300) * t288 + ((t207 * t250 + t225 * t290) * t184 + (-t308 + t277) * t174 + (-t300 * t171 - (t245 * t290 - t167 * t214 - t180 * t188 + t230 * t250 + (-t180 * t232 - t223 * t250) * t175) * t305 - (-t223 * t290 - t167 * t232 - t180 * t202 - t205 * t250 + (t180 * t214 - t245 * t250) * t175) * t306) * t185) * t181 (t172 * t309 - t184 * t218) * t288 + (t172 * t277 + t191 * t184 + (-t218 * t171 - t172 * t190 - (-t168 * t214 - t177 * t188 + t203 + (-t177 * t232 - t216) * t175) * t305 - (-t168 * t232 - t177 * t202 - t189 + (t177 * t214 - t233) * t175) * t306) * t185) * t181, 0, 0; 0 (t211 * t251 * t263 + t221 * t303) * t287 + (0.2e1 * t221 * t284 + (-qJD(3) * t299 + t241 * t251) * t211 + (-(t240 * t250 - t241 * t296 + t244 * t290) * t225 - t221 * t207 - (-t251 * t191 - (t250 * t289 + t253 * t292) * t225) * t263) * t212) * t193 (-t211 * t265 + t253 * t301) * t287 + (0.2e1 * t253 * t285 - t208 * t211 + (-0.2e1 * t207 * t225 * t253 - t191 * t265 + t222 * t291) * t212) * t193, t303 * t310 * t316 + (t284 * t316 + (t190 * t225 + t207 * t217) * t212) * t193, 0, 0;];
JaD_rot  = t1;
