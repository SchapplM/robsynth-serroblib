% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR4_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:39:33
% EndTime: 2019-02-26 21:39:34
% DurationCPUTime: 0.91s
% Computational Cost: add. (3615->96), mult. (9869->202), div. (448->12), fcn. (12731->13), ass. (0->97)
t227 = qJ(4) + pkin(12);
t225 = sin(t227);
t226 = cos(t227);
t282 = sin(pkin(11));
t284 = cos(pkin(6));
t252 = t284 * t282;
t283 = cos(pkin(11));
t253 = t284 * t283;
t285 = sin(qJ(2));
t286 = cos(qJ(2));
t216 = t286 * t252 + t285 * t253;
t218 = t285 * t282 - t286 * t283;
t229 = sin(qJ(1));
t230 = cos(qJ(1));
t249 = t229 * t216 + t218 * t230;
t228 = sin(pkin(6));
t267 = t228 * t229;
t244 = t225 * t249 + t226 * t267;
t289 = t244 * qJD(4);
t241 = t286 * t282 + t285 * t283;
t215 = t241 * t228;
t206 = qJD(2) * t215;
t214 = t218 * t228;
t211 = 0.1e1 / t214 ^ 2;
t269 = t206 * t211;
t288 = t285 * t252 - t286 * t253;
t240 = t218 * qJD(2);
t198 = -t229 * t241 - t230 * t288;
t186 = atan2(t198, t214);
t181 = sin(t186);
t182 = cos(t186);
t195 = t198 ^ 2;
t185 = t195 * t211 + 0.1e1;
t183 = 0.1e1 / t185;
t210 = 0.1e1 / t214;
t271 = t198 * t210;
t287 = (t182 * t271 - t181) * t183 + t181;
t170 = t181 * t198 + t182 * t214;
t167 = 0.1e1 / t170;
t194 = t225 * t267 - t226 * t249;
t188 = 0.1e1 / t194;
t168 = 0.1e1 / t170 ^ 2;
t189 = 0.1e1 / t194 ^ 2;
t237 = t229 * t288;
t201 = -t230 * t241 + t237;
t196 = t201 ^ 2;
t166 = t168 * t196 + 0.1e1;
t209 = t216 * qJD(2);
t176 = t198 * qJD(1) - t229 * t209 - t230 * t240;
t275 = t176 * t168;
t265 = qJD(1) * t230;
t179 = qJD(1) * t237 - t230 * t209 + t229 * t240 - t241 * t265;
t248 = t179 * t210 - t198 * t269;
t161 = t248 * t183;
t251 = -t181 * t214 + t182 * t198;
t157 = t251 * t161 + t181 * t179 + t182 * t206;
t280 = t157 * t167 * t168;
t281 = (-t196 * t280 - t201 * t275) / t166 ^ 2;
t250 = -t216 * t230 + t229 * t218;
t270 = t198 * t215;
t246 = -t210 * t250 + t211 * t270;
t162 = t246 * t183;
t158 = -t251 * t162 + t181 * t250 + t182 * t215;
t279 = t158 * t201;
t208 = t288 * qJD(2);
t217 = t241 * qJD(2);
t177 = t250 * qJD(1) + t229 * t208 - t217 * t230;
t259 = t228 * t265;
t171 = t194 * qJD(4) + t177 * t225 - t226 * t259;
t187 = t244 ^ 2;
t175 = t187 * t189 + 0.1e1;
t272 = t189 * t244;
t172 = t177 * t226 + t225 * t259 + t289;
t276 = t172 * t188 * t189;
t278 = (-t171 * t272 - t187 * t276) / t175 ^ 2;
t268 = t210 * t269;
t277 = (t179 * t198 * t211 - t195 * t268) / t185 ^ 2;
t274 = t181 * t201;
t273 = t182 * t201;
t266 = t228 * t230;
t264 = -0.2e1 * t281;
t263 = -0.2e1 * t280;
t262 = 0.2e1 * t278;
t261 = 0.2e1 * t277;
t260 = qJD(1) * t267;
t258 = -0.2e1 * t210 * t277;
t257 = -0.2e1 * t244 * t276;
t247 = -t225 * t188 - t226 * t272;
t245 = -t225 * t250 + t226 * t266;
t192 = t225 * t266 + t226 * t250;
t239 = t249 * qJD(1) + t208 * t230 + t229 * t217;
t207 = t228 * t240;
t173 = 0.1e1 / t175;
t164 = 0.1e1 / t166;
t160 = t287 * t201;
t156 = t246 * t261 + (0.2e1 * t268 * t270 + t239 * t210 + (-t179 * t215 + t198 * t207 - t206 * t250) * t211) * t183;
t1 = [t201 * t258 + (-t176 * t210 - t201 * t269) * t183, t156, 0, 0, 0, 0; t198 * t167 * t264 + (t179 * t167 + (-t157 * t198 - t160 * t176) * t168) * t164 + ((t160 * t263 - t287 * t275) * t164 + (t160 * t264 + ((-t161 * t183 * t271 + t261) * t274 + (t198 * t258 + t161 + (-t161 + t248) * t183) * t273) * t164) * t168) * t201, 0.2e1 * (t167 * t249 - t168 * t279) * t281 + (t177 * t167 + t263 * t279 + (t249 * t157 - t158 * t176 + (t156 * t198 - t162 * t179 - t207 + (t162 * t214 + t250) * t161) * t273 + (-t156 * t214 + t162 * t206 + t239 + (t162 * t198 - t215) * t161) * t274) * t168) * t164, 0, 0, 0, 0; (t188 * t245 - t192 * t272) * t262 + ((t192 * qJD(4) + t225 * t239 + t226 * t260) * t188 + t192 * t257 + (t245 * t172 + (t245 * qJD(4) - t225 * t260 + t226 * t239) * t244 - t192 * t171) * t189) * t173, t247 * t201 * t262 + (t247 * t176 + ((qJD(4) * t188 + t257) * t226 + (-t171 * t226 + (-t172 - t289) * t225) * t189) * t201) * t173, 0, -0.2e1 * t278 - 0.2e1 * (t171 * t189 * t173 - (-t173 * t276 - t189 * t278) * t244) * t244, 0, 0;];
JaD_rot  = t1;
