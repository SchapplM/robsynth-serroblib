% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRPR5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:48:46
% EndTime: 2019-02-26 19:48:47
% DurationCPUTime: 0.99s
% Computational Cost: add. (5804->108), mult. (9085->220), div. (559->12), fcn. (11668->13), ass. (0->105)
t224 = cos(pkin(6));
t228 = cos(qJ(2));
t278 = cos(pkin(10));
t249 = t278 * t228;
t222 = sin(pkin(10));
t226 = sin(qJ(2));
t262 = t222 * t226;
t212 = t224 * t249 - t262;
t208 = t212 * qJD(2);
t221 = pkin(11) + qJ(4);
t220 = cos(t221);
t219 = sin(t221);
t250 = t278 * t226;
t261 = t222 * t228;
t237 = -t224 * t250 - t261;
t223 = sin(pkin(6));
t251 = t223 * t278;
t280 = t237 * t219 - t220 * t251;
t175 = t280 * qJD(4) + t208 * t220;
t198 = -t219 * t251 - t220 * t237;
t195 = t198 ^ 2;
t260 = t223 * t226;
t207 = t219 * t224 + t220 * t260;
t204 = 0.1e1 / t207 ^ 2;
t188 = t195 * t204 + 0.1e1;
t186 = 0.1e1 / t188;
t206 = -t219 * t260 + t220 * t224;
t259 = t223 * t228;
t252 = qJD(2) * t259;
t194 = qJD(4) * t206 + t220 * t252;
t203 = 0.1e1 / t207;
t267 = t198 * t204;
t158 = (-t175 * t203 + t194 * t267) * t186;
t189 = atan2(-t198, t207);
t182 = sin(t189);
t183 = cos(t189);
t244 = -t182 * t207 - t183 * t198;
t154 = t158 * t244 - t182 * t175 + t183 * t194;
t168 = -t182 * t198 + t183 * t207;
t165 = 0.1e1 / t168;
t166 = 0.1e1 / t168 ^ 2;
t284 = t154 * t165 * t166;
t227 = cos(qJ(6));
t214 = -t224 * t262 + t249;
t263 = t222 * t223;
t240 = -t214 * t219 + t220 * t263;
t225 = sin(qJ(6));
t238 = -t224 * t261 - t250;
t265 = t238 * t225;
t243 = -t227 * t240 + t265;
t283 = qJD(6) * t243;
t201 = t214 * t220 + t219 * t263;
t282 = 0.2e1 * t201 * t284;
t239 = -t203 * t212 + t259 * t267;
t281 = t220 * t239;
t268 = t194 * t203 * t204;
t279 = -0.2e1 * (t175 * t267 - t195 * t268) / t188 ^ 2;
t264 = t238 * t227;
t185 = -t225 * t240 - t264;
t179 = 0.1e1 / t185;
t180 = 0.1e1 / t185 ^ 2;
t210 = t238 * qJD(2);
t176 = qJD(4) * t201 + t210 * t219;
t211 = t214 * qJD(2);
t169 = qJD(6) * t185 - t176 * t227 + t211 * t225;
t178 = t243 ^ 2;
t173 = t178 * t180 + 0.1e1;
t272 = t180 * t243;
t170 = t176 * t225 + t211 * t227 + t283;
t275 = t170 * t179 * t180;
t277 = (-t169 * t272 - t178 * t275) / t173 ^ 2;
t276 = t166 * t201;
t177 = qJD(4) * t240 + t210 * t220;
t274 = t177 * t166;
t273 = t179 * t227;
t271 = t182 * t201;
t270 = t183 * t201;
t269 = t243 * t225;
t266 = t238 * t220;
t258 = qJD(2) * t226;
t257 = qJD(4) * t219;
t196 = t201 ^ 2;
t164 = t166 * t196 + 0.1e1;
t256 = 0.2e1 * (-t196 * t284 + t201 * t274) / t164 ^ 2;
t255 = 0.2e1 * t277;
t248 = -0.2e1 * t243 * t275;
t247 = -0.2e1 * t198 * t268;
t245 = qJD(6) * t219 * t238 + t210;
t242 = -t269 * t180 + t273;
t241 = -t203 * t280 + t206 * t267;
t236 = -qJD(4) * t266 + qJD(6) * t214 + t211 * t219;
t209 = t237 * qJD(2);
t193 = -qJD(4) * t207 - t219 * t252;
t191 = t214 * t227 + t219 * t265;
t190 = t214 * t225 - t219 * t264;
t174 = qJD(4) * t198 + t208 * t219;
t171 = 0.1e1 / t173;
t161 = 0.1e1 / t164;
t160 = t186 * t281;
t159 = t241 * t186;
t156 = (-t182 * t212 + t183 * t259) * t220 + t244 * t160;
t155 = t159 * t244 - t182 * t280 + t183 * t206;
t153 = t241 * t279 + (t206 * t247 + t174 * t203 + (t175 * t206 + t193 * t198 + t194 * t280) * t204) * t186;
t151 = t279 * t281 + (-t239 * t257 + (t247 * t259 - t203 * t209 + (t194 * t212 + (t175 * t228 - t198 * t258) * t223) * t204) * t220) * t186;
t1 = [0, t151, 0, t153, 0, 0; 0 (t156 * t276 - t165 * t266) * t256 + ((-t211 * t220 - t238 * t257) * t165 + (-t274 + t282) * t156 + (-t266 * t154 - (-t151 * t198 - t160 * t175 + (-t220 * t258 - t228 * t257) * t223 + (-t160 * t207 - t212 * t220) * t158) * t270 - (t212 * t257 - t151 * t207 - t160 * t194 - t209 * t220 + (t160 * t198 - t220 * t259) * t158) * t271) * t166) * t161, 0 (t155 * t276 - t165 * t240) * t256 + (t155 * t282 - t176 * t165 + (-t240 * t154 - t155 * t177 - (-t153 * t198 - t159 * t175 + t193 + (-t159 * t207 - t280) * t158) * t270 - (-t153 * t207 - t159 * t194 + t174 + (t159 * t198 - t206) * t158) * t271) * t166) * t161, 0, 0; 0 (-t179 * t190 - t191 * t272) * t255 + (t191 * t248 + t245 * t179 * t225 + t236 * t273 + (t227 * t243 * t245 - t191 * t169 - t190 * t170 - t236 * t269) * t180) * t171, 0, t242 * t201 * t255 + (-t242 * t177 + ((qJD(6) * t179 + t248) * t225 + (-t169 * t225 + (t170 + t283) * t227) * t180) * t201) * t171, 0, -0.2e1 * t277 - 0.2e1 * (t169 * t171 * t180 - (-t171 * t275 - t180 * t277) * t243) * t243;];
JaD_rot  = t1;
