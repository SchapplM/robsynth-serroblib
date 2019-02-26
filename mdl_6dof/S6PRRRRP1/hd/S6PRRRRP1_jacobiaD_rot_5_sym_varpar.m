% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRRP1_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP1_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:15:13
% EndTime: 2019-02-26 20:15:14
% DurationCPUTime: 1.14s
% Computational Cost: add. (9025->114), mult. (13312->231), div. (822->12), fcn. (17117->13), ass. (0->113)
t259 = sin(pkin(11));
t261 = cos(pkin(11));
t264 = sin(qJ(2));
t262 = cos(pkin(6));
t266 = cos(qJ(2));
t293 = t262 * t266;
t247 = -t259 * t264 + t261 * t293;
t243 = t247 * qJD(2);
t294 = t262 * t264;
t248 = t259 * t266 + t261 * t294;
t258 = qJ(3) + qJ(4);
t255 = sin(t258);
t257 = qJD(3) + qJD(4);
t260 = sin(pkin(6));
t297 = t260 * t261;
t282 = t255 * t297;
t256 = cos(t258);
t299 = t256 * t257;
t210 = t243 * t255 + t248 * t299 - t257 * t282;
t232 = t248 * t255 + t256 * t297;
t230 = t232 ^ 2;
t296 = t260 * t264;
t284 = t255 * t296;
t240 = -t262 * t256 + t284;
t238 = 0.1e1 / t240 ^ 2;
t224 = t230 * t238 + 0.1e1;
t222 = 0.1e1 / t224;
t291 = qJD(2) * t266;
t275 = t257 * t262 + t260 * t291;
t283 = t256 * t296;
t228 = t255 * t275 + t257 * t283;
t237 = 0.1e1 / t240;
t304 = t232 * t238;
t194 = (-t210 * t237 + t228 * t304) * t222;
t225 = atan2(-t232, t240);
t220 = sin(t225);
t221 = cos(t225);
t278 = -t220 * t240 - t221 * t232;
t190 = t194 * t278 - t220 * t210 + t221 * t228;
t204 = -t220 * t232 + t221 * t240;
t201 = 0.1e1 / t204;
t202 = 0.1e1 / t204 ^ 2;
t318 = t190 * t201 * t202;
t285 = t259 * t294;
t250 = t261 * t266 - t285;
t298 = t259 * t260;
t235 = t250 * t255 - t256 * t298;
t317 = 0.2e1 * t235 * t318;
t295 = t260 * t266;
t274 = -t237 * t247 + t295 * t304;
t316 = t255 * t274;
t305 = t228 * t237 * t238;
t315 = -0.2e1 * (t210 * t304 - t230 * t305) / t224 ^ 2;
t236 = t250 * t256 + t255 * t298;
t265 = cos(qJ(5));
t249 = t259 * t293 + t261 * t264;
t263 = sin(qJ(5));
t302 = t249 * t263;
t219 = t236 * t265 + t302;
t215 = 0.1e1 / t219;
t216 = 0.1e1 / t219 ^ 2;
t245 = t249 * qJD(2);
t280 = t257 * t298 - t245;
t300 = t255 * t257;
t213 = -t250 * t300 + t256 * t280;
t246 = -qJD(2) * t285 + t261 * t291;
t205 = qJD(5) * t219 + t213 * t263 - t246 * t265;
t301 = t249 * t265;
t218 = t236 * t263 - t301;
t214 = t218 ^ 2;
t209 = t214 * t216 + 0.1e1;
t309 = t216 * t218;
t290 = qJD(5) * t218;
t206 = t213 * t265 + t246 * t263 - t290;
t312 = t206 * t215 * t216;
t314 = (t205 * t309 - t214 * t312) / t209 ^ 2;
t313 = t202 * t235;
t212 = t250 * t299 + t255 * t280;
t311 = t212 * t202;
t310 = t215 * t263;
t308 = t218 * t265;
t307 = t220 * t235;
t306 = t221 * t235;
t303 = t249 * t255;
t292 = qJD(2) * t264;
t231 = t235 ^ 2;
t200 = t231 * t202 + 0.1e1;
t289 = 0.2e1 * (-t231 * t318 + t235 * t311) / t200 ^ 2;
t288 = -0.2e1 * t314;
t286 = t218 * t312;
t281 = -0.2e1 * t232 * t305;
t279 = qJD(5) * t249 * t256 - t245;
t277 = t216 * t308 - t310;
t234 = t248 * t256 - t282;
t241 = t262 * t255 + t283;
t276 = -t234 * t237 + t241 * t304;
t273 = qJD(5) * t250 - t246 * t256 + t249 * t300;
t244 = t248 * qJD(2);
t229 = t256 * t275 - t257 * t284;
t227 = t250 * t263 - t256 * t301;
t226 = -t250 * t265 - t256 * t302;
t211 = -t248 * t300 + (-t257 * t297 + t243) * t256;
t207 = 0.1e1 / t209;
t198 = 0.1e1 / t200;
t196 = t222 * t316;
t195 = t276 * t222;
t192 = (-t220 * t247 + t221 * t295) * t255 + t278 * t196;
t191 = t195 * t278 - t220 * t234 + t221 * t241;
t189 = t276 * t315 + (t241 * t281 - t211 * t237 + (t210 * t241 + t228 * t234 + t229 * t232) * t238) * t222;
t187 = t315 * t316 + (t274 * t299 + (t281 * t295 + t237 * t244 + (t228 * t247 + (t210 * t266 - t232 * t292) * t260) * t238) * t255) * t222;
t186 = t277 * t235 * t288 + (t277 * t212 + ((-qJD(5) * t215 - 0.2e1 * t286) * t265 + (t205 * t265 + (t206 - t290) * t263) * t216) * t235) * t207;
t185 = (t191 * t313 - t201 * t236) * t289 + (t191 * t317 + t213 * t201 + (-t236 * t190 - t191 * t212 - (-t189 * t232 - t195 * t210 + t229 + (-t195 * t240 - t234) * t194) * t306 - (-t189 * t240 - t195 * t228 - t211 + (t195 * t232 - t241) * t194) * t307) * t202) * t198;
t1 = [0, t187, t189, t189, 0, 0; 0 (t192 * t313 + t201 * t303) * t289 + ((-t246 * t255 - t249 * t299) * t201 + (-t311 + t317) * t192 + (t303 * t190 - (-t187 * t232 - t196 * t210 + (-t255 * t292 + t266 * t299) * t260 + (-t196 * t240 - t247 * t255) * t194) * t306 - (-t247 * t299 - t187 * t240 - t196 * t228 + t244 * t255 + (t196 * t232 - t255 * t295) * t194) * t307) * t202) * t198, t185, t185, 0, 0; 0, 0.2e1 * (-t215 * t226 + t227 * t309) * t314 + (0.2e1 * t227 * t286 - t279 * t215 * t265 + t273 * t310 + (-t218 * t263 * t279 - t227 * t205 - t226 * t206 - t273 * t308) * t216) * t207, t186, t186, t288 + 0.2e1 * (t205 * t216 * t207 + (-t207 * t312 - t216 * t314) * t218) * t218, 0;];
JaD_rot  = t1;
