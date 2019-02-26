% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRP4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:41:10
% EndTime: 2019-02-26 22:41:11
% DurationCPUTime: 1.26s
% Computational Cost: add. (12642->127), mult. (10592->272), div. (1959->15), fcn. (13293->9), ass. (0->120)
t230 = qJ(2) + qJ(3);
t225 = cos(t230);
t229 = qJ(4) + qJ(5);
t222 = sin(t229);
t300 = sin(qJ(1));
t258 = t300 * t222;
t224 = cos(t229);
t231 = cos(qJ(1));
t278 = t231 * t224;
t204 = t225 * t278 + t258;
t198 = 0.1e1 / t204 ^ 2;
t223 = sin(t230);
t218 = t223 ^ 2;
t228 = t231 ^ 2;
t286 = t218 * t228;
t268 = t198 * t286;
t194 = 0.1e1 + t268;
t226 = qJD(4) + qJD(5);
t254 = qJD(1) * t300;
t227 = qJD(2) + qJD(3);
t281 = t227 * t231;
t261 = t223 * t281;
t241 = t225 * t254 + t261;
t256 = t300 * t226;
t279 = t231 * t222;
t183 = (-t225 * t226 + qJD(1)) * t279 + (t256 - t241) * t224;
t197 = 0.1e1 / t204;
t295 = t183 * t197 * t198;
t249 = t286 * t295;
t262 = t223 * t227 * t228;
t304 = (-t249 + (-t218 * t231 * t254 + t225 * t262) * t198) / t194 ^ 2;
t284 = t223 * t231;
t200 = t225 * t258 + t278;
t247 = t222 * t256;
t260 = t226 * t278;
t182 = t200 * qJD(1) + t222 * t261 - t225 * t260 - t247;
t257 = t300 * t224;
t203 = t225 * t279 - t257;
t215 = 0.1e1 / t222;
t216 = 0.1e1 / t222 ^ 2;
t219 = 0.1e1 / t223;
t220 = 0.1e1 / t223 ^ 2;
t282 = t225 * t227;
t263 = t220 * t282;
t283 = t224 * t226;
t288 = t215 * t219;
t303 = (t216 * t219 * t283 + t215 * t263) * t203 + t182 * t288;
t285 = t223 * t222;
t190 = atan2(-t200, t285);
t187 = cos(t190);
t186 = sin(t190);
t294 = t186 * t200;
t181 = t187 * t285 - t294;
t178 = 0.1e1 / t181;
t179 = 0.1e1 / t181 ^ 2;
t302 = -0.2e1 * t200;
t301 = 0.2e1 * t203;
t195 = t200 ^ 2;
t287 = t216 * t220;
t191 = t195 * t287 + 0.1e1;
t188 = 0.1e1 / t191;
t243 = t222 * t282 + t223 * t283;
t266 = t200 * t287;
t246 = t224 * t254;
t259 = t223 * t300;
t248 = t227 * t259;
t277 = qJD(1) * t231;
t184 = -t222 * t248 - t226 * t279 - t246 + (t222 * t277 + t224 * t256) * t225;
t269 = t184 * t288;
t170 = (t243 * t266 - t269) * t188;
t239 = -t170 * t200 + t243;
t165 = (-t170 * t285 - t184) * t186 + t239 * t187;
t180 = t178 * t179;
t299 = t165 * t180;
t217 = t215 * t216;
t221 = t219 / t218;
t264 = t220 * t283;
t298 = (t184 * t266 + (-t216 * t221 * t282 - t217 * t264) * t195) / t191 ^ 2;
t297 = t179 * t203;
t296 = t182 * t179;
t293 = t186 * t203;
t292 = t186 * t223;
t291 = t187 * t200;
t290 = t187 * t203;
t289 = t187 * t225;
t280 = t231 * t178;
t196 = t203 ^ 2;
t176 = t179 * t196 + 0.1e1;
t276 = 0.2e1 * (-t196 * t299 - t203 * t296) / t176 ^ 2;
t275 = -0.2e1 * t298;
t274 = 0.2e1 * t304;
t273 = t180 * t301;
t272 = t219 * t298;
t271 = t179 * t293;
t267 = t200 * t288;
t265 = t215 * t220 * t225;
t244 = t200 * t265 + t300;
t177 = t244 * t188;
t255 = t300 - t177;
t253 = t178 * t276;
t252 = t179 * t276;
t251 = t284 * t301;
t250 = t215 * t272;
t202 = t225 * t257 - t279;
t245 = t200 * t216 * t224 - t202 * t215;
t242 = t198 * t202 * t231 - t300 * t197;
t192 = 0.1e1 / t194;
t185 = t204 * qJD(1) - t224 * t248 - t225 * t247 - t260;
t174 = 0.1e1 / t176;
t173 = t245 * t219 * t188;
t169 = (-t186 + (t187 * t267 + t186) * t188) * t203;
t168 = -t177 * t291 + (t255 * t292 + t289) * t222;
t167 = t187 * t223 * t224 - t186 * t202 + (-t186 * t285 - t291) * t173;
t166 = t198 * t251 * t304 + (t251 * t295 + (t182 * t284 + (t223 * t254 - t225 * t281) * t203) * t198) * t192;
t164 = t244 * t275 + (t184 * t265 + t277 + (-t216 * t225 * t264 + (-0.2e1 * t221 * t225 ^ 2 - t219) * t227 * t215) * t200) * t188;
t162 = (t197 * t225 * t231 + t224 * t268) * t274 + (0.2e1 * t224 * t249 + t241 * t197 + ((t183 * t231 - 0.2e1 * t224 * t262) * t225 + (t222 * t226 * t228 + 0.2e1 * t231 * t246) * t218) * t198) * t192;
t161 = -0.2e1 * t245 * t272 + (-t245 * t263 + ((-t200 * t226 - t185) * t215 + (t217 * t283 * t302 + (t202 * t226 + t184) * t216) * t224) * t219) * t188;
t160 = t168 * t203 * t252 + (-(-t164 * t291 + (t170 * t294 - t184 * t187) * t177) * t297 + (-t223 * t280 - (-t177 * t292 + t186 * t259 + t289) * t297) * t283 + (t165 * t273 + t296) * t168) * t174 + (t253 * t284 + ((-t227 * t280 - (t255 * t227 - t170) * t271) * t225 + (t178 * t254 + (t231 * t165 - (-t164 + t277) * t293 - (t255 * t170 - t227) * t290) * t179) * t223) * t174) * t222;
t159 = (t167 * t297 - t178 * t204) * t276 + (t167 * t296 + t183 * t178 + (t167 * t273 - t179 * t204) * t165 - (-t226 * t285 + t224 * t282 - t161 * t200 - t173 * t184 + (-t173 * t285 - t202) * t170) * t179 * t290 - (-t185 + (-t161 * t222 - t170 * t224) * t223 - t239 * t173) * t271) * t174;
t1 = [t303 * t188 + t250 * t301, t164, t164, t161, t161, 0; t200 * t253 + (-t184 * t178 + (t165 * t200 + t169 * t182) * t179) * t174 + (t169 * t252 + (0.2e1 * t169 * t299 + (t182 * t188 - t182 - (-t170 * t188 * t267 + t275) * t203) * t179 * t186 + (-(t250 * t302 - t170) * t297 + (-(t170 + t269) * t203 + t303 * t200) * t179 * t188) * t187) * t174) * t203, t160, t160, t159, t159, 0; t242 * t223 * t274 + (-t242 * t282 + ((qJD(1) * t197 + 0.2e1 * t202 * t295) * t231 + (-t300 * t183 - t185 * t231 + t202 * t254) * t198) * t223) * t192, t162, t162, t166, t166, 0;];
JaD_rot  = t1;
