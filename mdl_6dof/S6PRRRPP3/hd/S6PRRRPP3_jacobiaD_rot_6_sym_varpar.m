% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function JaD_rot = S6PRRRPP3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:10:02
% EndTime: 2019-02-26 20:10:03
% DurationCPUTime: 1.66s
% Computational Cost: add. (7193->157), mult. (20987->307), div. (784->12), fcn. (27058->13), ass. (0->120)
t247 = sin(pkin(10));
t249 = cos(pkin(10));
t256 = cos(qJ(2));
t250 = cos(pkin(6));
t253 = sin(qJ(2));
t289 = t250 * t253;
t239 = t247 * t256 + t249 * t289;
t252 = sin(qJ(3));
t248 = sin(pkin(6));
t255 = cos(qJ(3));
t290 = t248 * t255;
t222 = -t239 * t252 - t249 * t290;
t288 = t250 * t256;
t264 = -t247 * t253 + t249 * t288;
t234 = t264 * qJD(2);
t204 = t222 * qJD(3) + t234 * t255;
t235 = t239 * qJD(2);
t251 = sin(qJ(4));
t254 = cos(qJ(4));
t292 = t248 * t252;
t267 = -t239 * t255 + t249 * t292;
t263 = t267 * t251;
t187 = t235 * t251 + qJD(4) * t263 + (-qJD(4) * t264 + t204) * t254;
t211 = -t251 * t264 - t267 * t254;
t207 = t211 ^ 2;
t243 = t250 * t252 + t253 * t290;
t287 = t251 * t256;
t265 = -t243 * t254 + t248 * t287;
t226 = 0.1e1 / t265 ^ 2;
t198 = t207 * t226 + 0.1e1;
t196 = 0.1e1 / t198;
t291 = t248 * t253;
t242 = t250 * t255 - t252 * t291;
t284 = qJD(2) * t256;
t275 = t248 * t284;
t229 = t242 * qJD(3) + t255 * t275;
t285 = t254 * t256;
t230 = -t243 * t251 - t248 * t285;
t276 = qJD(2) * t291;
t201 = t230 * qJD(4) + t229 * t254 + t251 * t276;
t225 = 0.1e1 / t265;
t298 = t211 * t226;
t173 = (t187 * t225 + t201 * t298) * t196;
t199 = atan2(-t211, -t265);
t193 = sin(t199);
t194 = cos(t199);
t272 = t193 * t265 - t194 * t211;
t169 = t272 * t173 - t187 * t193 + t194 * t201;
t185 = -t193 * t211 - t194 * t265;
t182 = 0.1e1 / t185;
t183 = 0.1e1 / t185 ^ 2;
t307 = t169 * t182 * t183;
t277 = t247 * t289;
t241 = t249 * t256 - t277;
t224 = t241 * t255 + t247 * t292;
t240 = t247 * t288 + t249 * t253;
t214 = t224 * t254 + t240 * t251;
t274 = 0.2e1 * t214 * t307;
t268 = t222 * t225 + t242 * t298;
t306 = t254 * t268;
t300 = t201 * t225 * t226;
t305 = -0.2e1 * (t187 * t298 + t207 * t300) / t198 ^ 2;
t266 = -t241 * t252 + t247 * t290;
t219 = 0.1e1 / t266;
t220 = 0.1e1 / t266 ^ 2;
t304 = t183 * t214;
t236 = t240 * qJD(2);
t206 = t266 * qJD(3) - t236 * t255;
t213 = -t224 * t251 + t240 * t254;
t237 = -qJD(2) * t277 + t249 * t284;
t189 = t213 * qJD(4) + t206 * t254 + t237 * t251;
t303 = t189 * t183;
t302 = t193 * t214;
t301 = t194 * t214;
t205 = t224 * qJD(3) - t236 * t252;
t221 = t219 * t220;
t299 = t205 * t221;
t297 = t213 * t220;
t296 = t213 * t224;
t295 = t266 * t254;
t294 = t240 * t252;
t293 = t240 * t255;
t286 = t254 * t255;
t283 = qJD(3) * t252;
t282 = qJD(4) * t251;
t281 = qJD(4) * t255;
t209 = t214 ^ 2;
t181 = t183 * t209 + 0.1e1;
t280 = 0.2e1 * (-t209 * t307 + t214 * t303) / t181 ^ 2;
t188 = -t214 * qJD(4) - t206 * t251 + t237 * t254;
t208 = t213 ^ 2;
t195 = t208 * t220 + 0.1e1;
t279 = 0.2e1 * (t188 * t297 + t208 * t299) / t195 ^ 2;
t273 = 0.2e1 * t211 * t300;
t271 = -qJD(4) * t241 + t237 * t255;
t210 = t254 * t264 - t263;
t270 = -t210 * t225 + t230 * t298;
t215 = t239 * t251 + t264 * t286;
t232 = (t251 * t253 + t255 * t285) * t248;
t269 = t215 * t225 + t232 * t298;
t228 = -t243 * qJD(3) - t252 * t275;
t217 = -t240 * t286 + t241 * t251;
t216 = t241 * t254 + t251 * t293;
t203 = t267 * qJD(3) - t234 * t252;
t202 = ((qJD(2) - t281) * t287 + (-t256 * t283 + (-qJD(2) * t255 + qJD(4)) * t253) * t254) * t248;
t200 = t265 * qJD(4) - t229 * t251 + t254 * t276;
t191 = 0.1e1 / t195;
t190 = (-t264 * t281 + t234) * t251 + (qJD(4) * t239 - t235 * t255 - t264 * t283) * t254;
t186 = t211 * qJD(4) + t204 * t251 - t235 * t254;
t179 = 0.1e1 / t181;
t178 = t196 * t306;
t177 = t269 * t196;
t175 = t270 * t196;
t172 = (-t193 * t222 + t194 * t242) * t254 + t272 * t178;
t171 = t272 * t177 - t193 * t215 + t194 * t232;
t170 = t272 * t175 + t193 * t210 + t194 * t230;
t168 = t269 * t305 + (t232 * t273 + t190 * t225 + (t187 * t232 + t201 * t215 + t202 * t211) * t226) * t196;
t166 = t270 * t305 + (t230 * t273 - t186 * t225 + (t187 * t230 + t200 * t211 - t201 * t210) * t226) * t196;
t165 = t305 * t306 + (-t268 * t282 + (t242 * t273 + t203 * t225 + (t187 * t242 + t201 * t222 + t211 * t228) * t226) * t254) * t196;
t1 = [0, t168, t165, t166, 0, 0; 0 (t171 * t304 - t182 * t217) * t280 + (t171 * t274 + (-t217 * t169 - t171 * t189 - (-t168 * t211 - t177 * t187 + t202 + (t177 * t265 - t215) * t173) * t301 - (t168 * t265 - t177 * t201 - t190 + (t177 * t211 - t232) * t173) * t302) * t183 + ((t240 * t281 - t236) * t251 + (t240 * t283 - t271) * t254) * t182) * t179 (t172 * t304 - t182 * t295) * t280 + ((-t205 * t254 - t266 * t282) * t182 + (-t303 + t274) * t172 + (-t295 * t169 - (-t242 * t282 - t165 * t211 - t178 * t187 + t228 * t254 + (t178 * t265 - t222 * t254) * t173) * t301 - (t222 * t282 + t165 * t265 - t178 * t201 - t203 * t254 + (t178 * t211 - t242 * t254) * t173) * t302) * t183) * t179 (t170 * t304 - t182 * t213) * t280 + (t170 * t274 + t188 * t182 + (-t213 * t169 - t170 * t189 - (-t166 * t211 - t175 * t187 + t200 + (t175 * t265 + t210) * t173) * t301 - (t166 * t265 - t175 * t201 + t186 + (t175 * t211 - t230) * t173) * t302) * t183) * t179, 0, 0; 0 (t216 * t219 - t294 * t297) * t279 + (-(-t236 * t254 + t271 * t251) * t219 + (-(-t251 * t283 + t254 * t281) * t219 + 0.2e1 * t252 * t213 * t299) * t240 + (t188 * t294 - t216 * t205 + (qJD(3) * t293 + t237 * t252) * t213) * t220) * t191 (-t219 * t251 * t266 + t220 * t296) * t279 + (qJD(4) * t219 * t295 + (-t188 * t224 - t206 * t213) * t220 + (-0.2e1 * t221 * t296 + (t220 * t266 - t219) * t251) * t205) * t191, -t214 * t219 * t279 + (t205 * t214 * t220 + t189 * t219) * t191, 0, 0;];
JaD_rot  = t1;
