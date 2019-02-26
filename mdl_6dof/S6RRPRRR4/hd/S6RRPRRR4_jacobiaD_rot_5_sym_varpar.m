% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR4_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR4_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:55:46
% EndTime: 2019-02-26 21:55:47
% DurationCPUTime: 0.93s
% Computational Cost: add. (3975->97), mult. (10490->203), div. (466->12), fcn. (13487->13), ass. (0->100)
t245 = sin(pkin(6));
t299 = sin(pkin(12));
t300 = cos(pkin(12));
t302 = sin(qJ(2));
t303 = cos(qJ(2));
t258 = t303 * t299 + t302 * t300;
t231 = t258 * t245;
t222 = qJD(2) * t231;
t234 = t302 * t299 - t303 * t300;
t230 = t234 * t245;
t227 = 0.1e1 / t230 ^ 2;
t286 = t222 * t227;
t301 = cos(pkin(6));
t269 = t301 * t299;
t270 = t301 * t300;
t305 = t302 * t269 - t303 * t270;
t257 = t234 * qJD(2);
t246 = sin(qJ(1));
t247 = cos(qJ(1));
t214 = -t246 * t258 - t247 * t305;
t202 = atan2(t214, t230);
t197 = sin(t202);
t198 = cos(t202);
t211 = t214 ^ 2;
t201 = t211 * t227 + 0.1e1;
t199 = 0.1e1 / t201;
t226 = 0.1e1 / t230;
t288 = t214 * t226;
t304 = (t198 * t288 - t197) * t199 + t197;
t186 = t197 * t214 + t198 * t230;
t183 = 0.1e1 / t186;
t244 = qJ(4) + qJ(5);
t241 = sin(t244);
t242 = cos(t244);
t232 = t303 * t269 + t302 * t270;
t266 = t246 * t232 + t247 * t234;
t284 = t245 * t246;
t210 = t241 * t284 - t242 * t266;
t204 = 0.1e1 / t210;
t184 = 0.1e1 / t186 ^ 2;
t205 = 0.1e1 / t210 ^ 2;
t254 = t246 * t305;
t217 = -t247 * t258 + t254;
t212 = t217 ^ 2;
t182 = t212 * t184 + 0.1e1;
t225 = t232 * qJD(2);
t192 = t214 * qJD(1) - t246 * t225 - t247 * t257;
t292 = t192 * t184;
t282 = qJD(1) * t247;
t195 = qJD(1) * t254 - t247 * t225 + t246 * t257 - t258 * t282;
t265 = t195 * t226 - t214 * t286;
t177 = t265 * t199;
t268 = -t197 * t230 + t198 * t214;
t173 = t268 * t177 + t197 * t195 + t198 * t222;
t297 = t173 * t183 * t184;
t298 = (-t212 * t297 - t217 * t292) / t182 ^ 2;
t267 = -t247 * t232 + t246 * t234;
t287 = t214 * t231;
t263 = -t226 * t267 + t227 * t287;
t178 = t263 * t199;
t174 = -t268 * t178 + t197 * t267 + t198 * t231;
t296 = t174 * t217;
t243 = qJD(4) + qJD(5);
t261 = t243 * t266 + t245 * t282;
t224 = t305 * qJD(2);
t233 = t258 * qJD(2);
t193 = t267 * qJD(1) + t246 * t224 - t247 * t233;
t273 = t243 * t284 + t193;
t187 = t273 * t241 - t261 * t242;
t209 = -t241 * t266 - t242 * t284;
t203 = t209 ^ 2;
t191 = t203 * t205 + 0.1e1;
t289 = t205 * t209;
t188 = t261 * t241 + t273 * t242;
t293 = t188 * t204 * t205;
t295 = (t187 * t289 - t203 * t293) / t191 ^ 2;
t285 = t226 * t286;
t294 = (t214 * t227 * t195 - t211 * t285) / t201 ^ 2;
t291 = t197 * t217;
t290 = t198 * t217;
t283 = t245 * t247;
t281 = -0.2e1 * t298;
t280 = -0.2e1 * t297;
t279 = 0.2e1 * t295;
t278 = 0.2e1 * t294;
t277 = -0.2e1 * t226 * t294;
t276 = 0.2e1 * t209 * t293;
t256 = t266 * qJD(1) + t247 * t224 + t246 * t233;
t272 = t243 * t283 + t256;
t264 = -t241 * t204 + t242 * t289;
t262 = qJD(1) * t284 + t243 * t267;
t223 = t245 * t257;
t208 = t241 * t283 + t242 * t267;
t207 = t241 * t267 - t242 * t283;
t189 = 0.1e1 / t191;
t180 = 0.1e1 / t182;
t176 = t304 * t217;
t172 = t263 * t278 + (0.2e1 * t285 * t287 + t256 * t226 + (-t195 * t231 + t214 * t223 - t222 * t267) * t227) * t199;
t170 = -0.2e1 * t295 + 0.2e1 * (t187 * t205 * t189 + (-t189 * t293 - t205 * t295) * t209) * t209;
t1 = [t217 * t277 + (-t192 * t226 - t217 * t286) * t199, t172, 0, 0, 0, 0; t214 * t183 * t281 + (t195 * t183 + (-t173 * t214 - t176 * t192) * t184) * t180 + ((t176 * t280 - t304 * t292) * t180 + (t176 * t281 + ((-t177 * t199 * t288 + t278) * t291 + (t214 * t277 + t177 + (-t177 + t265) * t199) * t290) * t180) * t184) * t217, 0.2e1 * (t183 * t266 - t184 * t296) * t298 + (t193 * t183 + t280 * t296 + (t266 * t173 - t174 * t192 + (t172 * t214 - t178 * t195 - t223 + (t178 * t230 + t267) * t177) * t290 + (-t172 * t230 + t178 * t222 + t256 + (t178 * t214 - t231) * t177) * t291) * t184) * t180, 0, 0, 0, 0; (-t204 * t207 + t208 * t289) * t279 + ((t272 * t241 + t262 * t242) * t204 + t208 * t276 + (-t207 * t188 - (-t262 * t241 + t272 * t242) * t209 - t208 * t187) * t205) * t189, t264 * t217 * t279 + (t264 * t192 + ((t204 * t243 + t276) * t242 + (-t187 * t242 + (t209 * t243 - t188) * t241) * t205) * t217) * t189, 0, t170, t170, 0;];
JaD_rot  = t1;
