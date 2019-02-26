% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:51:16
% EndTime: 2019-02-26 20:51:17
% DurationCPUTime: 1.20s
% Computational Cost: add. (7491->94), mult. (9702->187), div. (711->12), fcn. (12089->11), ass. (0->95)
t271 = pkin(10) + qJ(3);
t260 = sin(t271);
t261 = cos(t271);
t290 = sin(qJ(5));
t292 = cos(qJ(5));
t208 = t260 * t292 - t261 * t290;
t291 = sin(qJ(1));
t199 = t208 * t291;
t207 = t260 * t290 + t261 * t292;
t186 = atan2(t199, t207);
t181 = sin(t186);
t182 = cos(t186);
t201 = t207 * t291;
t255 = -t181 * t207 + t182 * t199;
t197 = t199 ^ 2;
t205 = 0.1e1 / t207 ^ 2;
t185 = t197 * t205 + 0.1e1;
t183 = 0.1e1 / t185;
t204 = 0.1e1 / t207;
t273 = t199 * t208;
t247 = -t201 * t204 - t205 * t273;
t300 = t247 * t183;
t158 = -t181 * t201 + t182 * t208 + t255 * t300;
t293 = cos(qJ(1));
t203 = t208 * t293;
t319 = t158 * t203;
t171 = t181 * t199 + t182 * t207;
t168 = 0.1e1 / t171;
t169 = 0.1e1 / t171 ^ 2;
t202 = t207 * t293;
t198 = t203 ^ 2;
t167 = t198 * t169 + 0.1e1;
t308 = qJD(3) - qJD(5);
t175 = -t199 * qJD(1) + t308 * t202;
t281 = t175 * t169;
t176 = -t203 * qJD(1) - t308 * t201;
t188 = t308 * t208;
t274 = t199 * t205;
t249 = -t176 * t204 + t188 * t274;
t161 = t249 * t183;
t156 = t161 * t255 - t181 * t176 - t182 * t188;
t289 = t156 * t168 * t169;
t270 = 0.2e1 * (-t198 * t289 + t203 * t281) / t167 ^ 2;
t318 = (t168 * t202 + t169 * t319) * t270;
t174 = qJD(1) * t201 + t308 * t203;
t226 = sin(qJ(6));
t227 = cos(qJ(6));
t196 = t202 * t227 - t226 * t291;
t263 = qJD(1) * t293;
t172 = qJD(6) * t196 - t174 * t226 + t227 * t263;
t245 = -t202 * t226 - t227 * t291;
t309 = qJD(6) * t245;
t173 = -t174 * t227 - t226 * t263 + t309;
t189 = t245 ^ 2;
t191 = 0.1e1 / t196 ^ 2;
t180 = t189 * t191 + 0.1e1;
t178 = 0.1e1 / t180;
t190 = 0.1e1 / t196;
t276 = t191 * t245;
t248 = -t226 * t190 - t227 * t276;
t283 = t173 * t190 * t191;
t295 = -0.2e1 * t245;
t259 = t283 * t295;
t317 = (t248 * t175 - ((t226 * (-t173 - t309) - t172 * t227) * t191 + (qJD(6) * t190 + t259) * t227) * t203) * t178;
t316 = -t202 * t156 + t158 * t175;
t269 = 0.2e1 * t289;
t315 = -t174 * t168 - t269 * t319;
t187 = t308 * t207;
t313 = (t207 * t300 + t201) * t161 + t300 * t176 - t187;
t177 = qJD(1) * t202 - t308 * t199;
t312 = -(-t199 * t300 - t208) * t161 - t300 * t188 + t177;
t302 = t188 * t205;
t277 = t204 * t302;
t310 = ((t176 * t208 - t187 * t199 - t188 * t201) * t205 - t177 * t204 - 0.2e1 * t273 * t277) * t183;
t275 = t199 * t204;
t298 = t183 * (-t182 * t275 + t181) - t181;
t294 = -0.2e1 * t203;
t288 = (-t172 * t276 - t189 * t283) / t180 ^ 2;
t287 = (-t176 * t274 + t197 * t277) / t185 ^ 2;
t280 = t178 * t191;
t279 = t181 * t203;
t278 = t182 * t203;
t268 = -0.2e1 * t288;
t267 = -0.2e1 * t287;
t266 = t191 * t288;
t265 = t204 * t287;
t264 = t172 * t280;
t262 = qJD(1) * t291;
t194 = -t201 * t227 - t226 * t293;
t246 = t201 * t226 - t227 * t293;
t165 = 0.1e1 / t167;
t160 = t298 * t203;
t154 = t247 * t267 + t310;
t153 = 0.2e1 * t247 * t287 - t310;
t1 = [t265 * t294 + (t175 * t204 + t203 * t302) * t183, 0, t153, 0, t154, 0; -t199 * t168 * t270 + (-t176 * t168 + (-t156 * t199 - t160 * t175) * t169) * t165 - ((-t160 * t269 + t298 * t281) * t165 + (-t160 * t270 + ((t161 * t183 * t275 + t267) * t279 + (0.2e1 * t199 * t265 - t161 + (t161 - t249) * t183) * t278) * t165) * t169) * t203, 0, t318 + (((t153 * t199 + t313) * t278 + (-t153 * t207 + t312) * t279 - t316) * t169 - t315) * t165, 0, -t318 + (((t154 * t199 - t313) * t278 + (-t154 * t207 - t312) * t279 + t316) * t169 + t315) * t165, 0; (t266 * t295 - t264) * t194 - (-t173 * t280 + t190 * t268) * t246 + ((qJD(6) * t194 - t177 * t226 - t227 * t262) * t190 + (qJD(6) * t246 - t177 * t227 + t226 * t262) * t276 + t194 * t259) * t178, 0, t248 * t288 * t294 + t317, 0, -t248 * t203 * t268 - t317, t268 + (t264 - (-t178 * t283 - t266) * t245) * t295;];
JaD_rot  = t1;
