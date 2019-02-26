% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PPRRRR3
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPRRRR3_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobiaD_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:44:05
% EndTime: 2019-02-26 19:44:06
% DurationCPUTime: 0.88s
% Computational Cost: add. (4180->76), mult. (12611->185), div. (275->12), fcn. (16541->17), ass. (0->91)
t235 = sin(pkin(14));
t240 = cos(pkin(14));
t241 = cos(pkin(13));
t236 = sin(pkin(13));
t244 = cos(pkin(6));
t272 = t236 * t244;
t232 = -t235 * t272 + t240 * t241;
t246 = sin(qJ(3));
t248 = cos(qJ(3));
t231 = -t235 * t241 - t240 * t272;
t243 = cos(pkin(7));
t238 = sin(pkin(7));
t239 = sin(pkin(6));
t271 = t238 * t239;
t258 = t231 * t243 + t236 * t271;
t292 = t232 * t246 - t258 * t248;
t268 = t240 * t243;
t270 = t238 * t244;
t290 = (-t235 * t246 + t248 * t268) * t239 + t248 * t270;
t221 = t232 * t248 + t258 * t246;
t267 = t241 * t244;
t230 = t235 * t267 + t236 * t240;
t257 = -t235 * t236 + t240 * t267;
t255 = -t241 * t271 + t257 * t243;
t289 = -t230 * t246 + t255 * t248;
t237 = sin(pkin(8));
t242 = cos(pkin(8));
t269 = t239 * t243;
t209 = t289 * t237 - (-t257 * t238 - t241 * t269) * t242;
t222 = -t290 * t237 + (-t240 * t271 + t243 * t244) * t242;
t204 = atan2(t209, t222);
t199 = sin(t204);
t200 = cos(t204);
t186 = t199 * t209 + t200 * t222;
t183 = 0.1e1 / t186;
t245 = sin(qJ(4));
t247 = cos(qJ(4));
t228 = -t231 * t238 + t236 * t269;
t277 = t228 * t237;
t260 = -t242 * t292 + t277;
t198 = t221 * t247 + t260 * t245;
t194 = 0.1e1 / t198;
t216 = 0.1e1 / t222;
t234 = t237 ^ 2;
t184 = 0.1e1 / t186 ^ 2;
t195 = 0.1e1 / t198 ^ 2;
t217 = 0.1e1 / t222 ^ 2;
t210 = t228 * t242 + t237 * t292;
t208 = t210 ^ 2;
t182 = t184 * t208 + 0.1e1;
t215 = t221 * qJD(3);
t219 = -t230 * t248 - t255 * t246;
t213 = t219 * qJD(3);
t227 = -t246 * t270 + (-t235 * t248 - t246 * t268) * t239;
t225 = t227 * qJD(3);
t281 = t209 * t217;
t207 = t209 ^ 2;
t203 = t207 * t217 + 0.1e1;
t201 = 0.1e1 / t203;
t282 = t201 * t237;
t178 = (t213 * t216 + t225 * t281) * t282;
t261 = -t199 * t222 + t200 * t209;
t175 = (t199 * t213 - t200 * t225) * t237 + t261 * t178;
t287 = t175 * t183 * t184;
t288 = (t184 * t210 * t215 * t237 - t208 * t287) / t182 ^ 2;
t265 = t242 * t247;
t278 = t221 * t245;
t197 = -t247 * t277 + t265 * t292 + t278;
t193 = t197 ^ 2;
t190 = t193 * t195 + 0.1e1;
t214 = t292 * qJD(3);
t266 = t242 * t245;
t192 = -t215 * t266 - t214 * t247 + (t260 * t247 - t278) * qJD(4);
t284 = t192 * t194 * t195;
t191 = t198 * qJD(4) - t214 * t245 + t215 * t265;
t285 = t191 * t195;
t286 = (-t193 * t284 + t197 * t285) / t190 ^ 2;
t206 = -t221 * t266 - t247 * t292;
t283 = t197 * t206;
t280 = t209 * t227;
t279 = t216 * t217 * t225;
t259 = t216 * t219 + t217 * t280;
t205 = t221 * t265 - t245 * t292;
t224 = t290 * qJD(3);
t212 = t289 * qJD(3);
t188 = 0.1e1 / t190;
t180 = 0.1e1 / t182;
t179 = t259 * t282;
t176 = (t199 * t219 - t200 * t227) * t237 + t261 * t179;
t174 = -0.2e1 * t259 * t234 / t203 ^ 2 * (t207 * t279 + t213 * t281) + (0.2e1 * t234 * t279 * t280 - t212 * t216 * t237 + (-t209 * t224 * t237 + (t213 * t227 + t219 * t225) * t234) * t217) * t201;
t1 = [0, 0, t174, 0, 0, 0; 0, 0 (-(-t186 * t179 * t178 + t261 * t174) * t184 * t180 + 0.2e1 * (t180 * t287 + t184 * t288) * t176) * t210 + (-0.2e1 * t221 * t183 * t288 + (-t214 * t183 + (-t221 * t175 - t176 * t215 + (-(t178 * t219 + t179 * t213 + t224) * t200 - (t178 * t227 + t179 * t225 - t212) * t199) * t210) * t184) * t180) * t237, 0, 0, 0; 0, 0, 0.2e1 * (-t194 * t205 + t195 * t283) * t286 + ((t206 * qJD(4) - t214 * t265 - t215 * t245) * t194 + 0.2e1 * t283 * t284 + (-t205 * t192 - (-t205 * qJD(4) + t214 * t266 - t215 * t247) * t197 - t206 * t191) * t195) * t188, -0.2e1 * t286 + 0.2e1 * (t188 * t285 + (-t188 * t284 - t195 * t286) * t197) * t197, 0, 0;];
JaD_rot  = t1;
