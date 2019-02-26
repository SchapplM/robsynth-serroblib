% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRP4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:03:04
% EndTime: 2019-02-26 20:03:05
% DurationCPUTime: 0.91s
% Computational Cost: add. (3002->107), mult. (9085->220), div. (559->12), fcn. (11668->13), ass. (0->104)
t221 = cos(pkin(6));
t227 = cos(qJ(2));
t277 = cos(pkin(10));
t248 = t277 * t227;
t219 = sin(pkin(10));
t224 = sin(qJ(2));
t262 = t219 * t224;
t210 = t221 * t248 - t262;
t203 = t210 * qJD(2);
t226 = cos(qJ(3));
t223 = sin(qJ(3));
t249 = t277 * t224;
t261 = t219 * t227;
t236 = -t221 * t249 - t261;
t220 = sin(pkin(6));
t250 = t220 * t277;
t279 = t236 * t223 - t226 * t250;
t175 = t279 * qJD(3) + t203 * t226;
t196 = -t223 * t250 - t236 * t226;
t193 = t196 ^ 2;
t259 = t220 * t226;
t214 = t221 * t223 + t224 * t259;
t208 = 0.1e1 / t214 ^ 2;
t188 = t193 * t208 + 0.1e1;
t186 = 0.1e1 / t188;
t260 = t220 * t223;
t213 = t221 * t226 - t224 * t260;
t258 = t220 * t227;
t251 = qJD(2) * t258;
t201 = t213 * qJD(3) + t226 * t251;
t207 = 0.1e1 / t214;
t267 = t196 * t208;
t158 = (-t175 * t207 + t201 * t267) * t186;
t189 = atan2(-t196, t214);
t184 = sin(t189);
t185 = cos(t189);
t243 = -t184 * t214 - t185 * t196;
t154 = t243 * t158 - t184 * t175 + t185 * t201;
t168 = -t184 * t196 + t185 * t214;
t165 = 0.1e1 / t168;
t166 = 0.1e1 / t168 ^ 2;
t283 = t154 * t165 * t166;
t222 = sin(qJ(5));
t225 = cos(qJ(5));
t237 = -t221 * t261 - t249;
t212 = -t221 * t262 + t248;
t239 = -t212 * t223 + t219 * t259;
t242 = t222 * t237 - t225 * t239;
t282 = t242 * qJD(5);
t199 = t212 * t226 + t219 * t260;
t281 = 0.2e1 * t199 * t283;
t238 = -t207 * t210 + t258 * t267;
t280 = t226 * t238;
t266 = t201 * t207 * t208;
t278 = -0.2e1 * (t175 * t267 - t193 * t266) / t188 ^ 2;
t264 = t237 * t225;
t183 = -t222 * t239 - t264;
t179 = 0.1e1 / t183;
t180 = 0.1e1 / t183 ^ 2;
t205 = t237 * qJD(2);
t176 = t199 * qJD(3) + t205 * t223;
t206 = t212 * qJD(2);
t169 = t183 * qJD(5) - t176 * t225 + t206 * t222;
t178 = t242 ^ 2;
t173 = t178 * t180 + 0.1e1;
t271 = t180 * t242;
t170 = t176 * t222 + t206 * t225 + t282;
t274 = t170 * t179 * t180;
t276 = (-t169 * t271 - t178 * t274) / t173 ^ 2;
t275 = t166 * t199;
t177 = t239 * qJD(3) + t205 * t226;
t273 = t177 * t166;
t272 = t179 * t225;
t270 = t242 * t222;
t269 = t184 * t199;
t268 = t185 * t199;
t265 = t237 * t223;
t263 = t237 * t226;
t257 = qJD(2) * t224;
t256 = qJD(3) * t223;
t194 = t199 ^ 2;
t164 = t166 * t194 + 0.1e1;
t255 = 0.2e1 * (-t194 * t283 + t199 * t273) / t164 ^ 2;
t254 = 0.2e1 * t276;
t247 = -0.2e1 * t242 * t274;
t246 = -0.2e1 * t196 * t266;
t244 = qJD(5) * t265 + t205;
t241 = -t180 * t270 + t272;
t240 = -t207 * t279 + t213 * t267;
t235 = -qJD(3) * t263 + qJD(5) * t212 + t206 * t223;
t204 = t236 * qJD(2);
t200 = -t214 * qJD(3) - t223 * t251;
t191 = t212 * t225 + t222 * t265;
t190 = t212 * t222 - t223 * t264;
t174 = t196 * qJD(3) + t203 * t223;
t171 = 0.1e1 / t173;
t161 = 0.1e1 / t164;
t160 = t186 * t280;
t159 = t240 * t186;
t156 = (-t184 * t210 + t185 * t258) * t226 + t243 * t160;
t155 = t243 * t159 - t184 * t279 + t185 * t213;
t153 = t240 * t278 + (t213 * t246 + t174 * t207 + (t175 * t213 + t196 * t200 + t201 * t279) * t208) * t186;
t151 = t278 * t280 + (-t238 * t256 + (t246 * t258 - t204 * t207 + (t201 * t210 + (t175 * t227 - t196 * t257) * t220) * t208) * t226) * t186;
t1 = [0, t151, t153, 0, 0, 0; 0 (t156 * t275 - t165 * t263) * t255 + ((-t206 * t226 - t237 * t256) * t165 + (-t273 + t281) * t156 + (-t263 * t154 - (-t151 * t196 - t160 * t175 + (-t226 * t257 - t227 * t256) * t220 + (-t160 * t214 - t210 * t226) * t158) * t268 - (t210 * t256 - t151 * t214 - t160 * t201 - t204 * t226 + (t160 * t196 - t226 * t258) * t158) * t269) * t166) * t161 (t155 * t275 - t165 * t239) * t255 + (t155 * t281 - t176 * t165 + (-t239 * t154 - t155 * t177 - (-t153 * t196 - t159 * t175 + t200 + (-t159 * t214 - t279) * t158) * t268 - (-t153 * t214 - t159 * t201 + t174 + (t159 * t196 - t213) * t158) * t269) * t166) * t161, 0, 0, 0; 0 (-t179 * t190 - t191 * t271) * t254 + (t191 * t247 + t244 * t179 * t222 + t235 * t272 + (t225 * t242 * t244 - t191 * t169 - t190 * t170 - t235 * t270) * t180) * t171, t241 * t199 * t254 + (-t241 * t177 + ((qJD(5) * t179 + t247) * t222 + (-t169 * t222 + (t170 + t282) * t225) * t180) * t199) * t171, 0, -0.2e1 * t276 - 0.2e1 * (t169 * t171 * t180 - (-t171 * t274 - t180 * t276) * t242) * t242, 0;];
JaD_rot  = t1;
