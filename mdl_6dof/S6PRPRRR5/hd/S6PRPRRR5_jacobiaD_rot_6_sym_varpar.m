% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR5
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
% Datum: 2019-02-26 19:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRR5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:56:16
% EndTime: 2019-02-26 19:56:17
% DurationCPUTime: 1.10s
% Computational Cost: add. (9025->112), mult. (13312->233), div. (822->12), fcn. (17117->13), ass. (0->113)
t252 = qJ(4) + qJ(5);
t250 = cos(t252);
t251 = qJD(4) + qJD(5);
t249 = sin(t252);
t253 = sin(pkin(11));
t255 = cos(pkin(11));
t258 = sin(qJ(2));
t256 = cos(pkin(6));
t260 = cos(qJ(2));
t288 = t256 * t260;
t272 = -t253 * t258 + t255 * t288;
t269 = t272 * t249;
t289 = t256 * t258;
t245 = t253 * t260 + t255 * t289;
t254 = sin(pkin(6));
t292 = t254 * t255;
t277 = t245 * qJD(2) + t251 * t292;
t208 = t277 * t250 + t251 * t269;
t268 = t272 * t250;
t229 = t249 * t292 - t268;
t226 = t229 ^ 2;
t290 = t254 * t260;
t278 = t250 * t290;
t238 = t256 * t249 + t278;
t236 = 0.1e1 / t238 ^ 2;
t219 = t226 * t236 + 0.1e1;
t217 = 0.1e1 / t219;
t291 = t254 * t258;
t271 = qJD(2) * t291 - t251 * t256;
t279 = t249 * t290;
t224 = -t271 * t250 - t251 * t279;
t235 = 0.1e1 / t238;
t300 = t229 * t236;
t189 = (t208 * t235 - t224 * t300) * t217;
t220 = atan2(t229, t238);
t215 = sin(t220);
t216 = cos(t220);
t275 = -t215 * t238 + t216 * t229;
t185 = t275 * t189 + t215 * t208 + t216 * t224;
t199 = t215 * t229 + t216 * t238;
t196 = 0.1e1 / t199;
t197 = 0.1e1 / t199 ^ 2;
t314 = t185 * t196 * t197;
t246 = t253 * t288 + t255 * t258;
t293 = t253 * t254;
t227 = -t246 * t250 + t249 * t293;
t313 = 0.2e1 * t227 * t314;
t301 = t224 * t235 * t236;
t312 = (t208 * t300 - t226 * t301) / t219 ^ 2;
t281 = t229 * t291;
t270 = t235 * t245 + t236 * t281;
t311 = t250 * t270;
t228 = t246 * t249 + t250 * t293;
t259 = cos(qJ(6));
t280 = t253 * t289;
t247 = t255 * t260 - t280;
t257 = sin(qJ(6));
t297 = t247 * t257;
t214 = t228 * t259 + t297;
t210 = 0.1e1 / t214;
t211 = 0.1e1 / t214 ^ 2;
t287 = qJD(2) * t260;
t244 = -qJD(2) * t280 + t255 * t287;
t294 = t250 * t251;
t205 = t246 * t294 + (-t251 * t293 + t244) * t249;
t243 = t246 * qJD(2);
t200 = t214 * qJD(6) + t205 * t257 + t243 * t259;
t296 = t247 * t259;
t213 = t228 * t257 - t296;
t209 = t213 ^ 2;
t204 = t209 * t211 + 0.1e1;
t305 = t211 * t213;
t286 = qJD(6) * t213;
t201 = t205 * t259 - t243 * t257 - t286;
t308 = t201 * t210 * t211;
t310 = (t200 * t305 - t209 * t308) / t204 ^ 2;
t309 = t197 * t227;
t206 = t228 * t251 - t244 * t250;
t307 = t206 * t197;
t306 = t210 * t257;
t304 = t213 * t259;
t303 = t215 * t227;
t302 = t216 * t227;
t239 = t256 * t250 - t279;
t299 = t229 * t239;
t298 = t247 * t250;
t295 = t249 * t251;
t225 = t227 ^ 2;
t195 = t225 * t197 + 0.1e1;
t285 = 0.2e1 * (-t225 * t314 + t227 * t307) / t195 ^ 2;
t284 = -0.2e1 * t310;
t282 = t213 * t308;
t276 = qJD(6) * t247 * t249 + t244;
t274 = t211 * t304 - t306;
t230 = t250 * t292 + t269;
t273 = -t230 * t235 + t236 * t299;
t267 = -qJD(6) * t246 - t243 * t249 + t247 * t294;
t241 = t272 * qJD(2);
t223 = t271 * t249 - t251 * t278;
t222 = -t246 * t257 + t249 * t296;
t221 = t246 * t259 + t249 * t297;
t207 = t277 * t249 - t251 * t268;
t202 = 0.1e1 / t204;
t193 = 0.1e1 / t195;
t191 = t217 * t311;
t190 = t273 * t217;
t187 = (t215 * t245 - t216 * t291) * t250 + t275 * t191;
t186 = -t275 * t190 + t215 * t230 + t216 * t239;
t184 = 0.2e1 * t273 * t312 + (0.2e1 * t299 * t301 - t207 * t235 + (-t208 * t239 - t223 * t229 - t224 * t230) * t236) * t217;
t182 = -0.2e1 * t311 * t312 + (-t270 * t295 + (-0.2e1 * t281 * t301 + t235 * t241 + (-t224 * t245 + (t208 * t258 + t229 * t287) * t254) * t236) * t250) * t217;
t181 = t274 * t227 * t284 + (t274 * t206 + ((-qJD(6) * t210 - 0.2e1 * t282) * t259 + (t200 * t259 + (t201 - t286) * t257) * t211) * t227) * t202;
t180 = (t186 * t309 - t196 * t228) * t285 + (t186 * t313 + t205 * t196 + (-t228 * t185 - t186 * t206 - (t184 * t229 - t190 * t208 + t223 + (t190 * t238 + t230) * t189) * t302 - (-t184 * t238 + t190 * t224 - t207 + (t190 * t229 - t239) * t189) * t303) * t197) * t193;
t1 = [0, t182, 0, t184, t184, 0; 0 (t187 * t309 + t196 * t298) * t285 + ((t243 * t250 + t247 * t295) * t196 + (-t307 + t313) * t187 + (t298 * t185 - (t182 * t229 + t191 * t208 + (-t250 * t287 + t258 * t295) * t254 + (-t191 * t238 + t245 * t250) * t189) * t302 - (-t245 * t295 - t182 * t238 - t191 * t224 + t241 * t250 + (-t191 * t229 + t250 * t291) * t189) * t303) * t197) * t193, 0, t180, t180, 0; 0, 0.2e1 * (-t210 * t221 + t222 * t305) * t310 + (0.2e1 * t222 * t282 + t276 * t210 * t259 + t267 * t306 + (t276 * t213 * t257 - t222 * t200 - t221 * t201 - t267 * t304) * t211) * t202, 0, t181, t181, t284 + 0.2e1 * (t200 * t211 * t202 + (-t202 * t308 - t211 * t310) * t213) * t213;];
JaD_rot  = t1;
