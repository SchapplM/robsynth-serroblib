% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRP3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:51:34
% EndTime: 2019-02-26 19:51:35
% DurationCPUTime: 0.94s
% Computational Cost: add. (5804->110), mult. (9085->227), div. (559->12), fcn. (11668->13), ass. (0->106)
t228 = sin(pkin(10));
t230 = cos(pkin(10));
t233 = sin(qJ(2));
t231 = cos(pkin(6));
t235 = cos(qJ(2));
t261 = t231 * t235;
t217 = -t228 * t233 + t230 * t261;
t213 = t217 * qJD(2);
t262 = t231 * t233;
t218 = t228 * t235 + t230 * t262;
t227 = pkin(11) + qJ(4);
t225 = sin(t227);
t229 = sin(pkin(6));
t265 = t229 * t230;
t251 = t225 * t265;
t226 = cos(t227);
t258 = qJD(4) * t226;
t180 = -qJD(4) * t251 + t213 * t225 + t218 * t258;
t202 = t218 * t225 + t226 * t265;
t200 = t202 ^ 2;
t264 = t229 * t233;
t210 = t225 * t264 - t231 * t226;
t208 = 0.1e1 / t210 ^ 2;
t194 = t200 * t208 + 0.1e1;
t192 = 0.1e1 / t194;
t211 = t231 * t225 + t226 * t264;
t259 = qJD(2) * t235;
t250 = t229 * t259;
t198 = t211 * qJD(4) + t225 * t250;
t207 = 0.1e1 / t210;
t270 = t202 * t208;
t164 = (-t180 * t207 + t198 * t270) * t192;
t195 = atan2(-t202, t210);
t188 = sin(t195);
t189 = cos(t195);
t247 = -t188 * t210 - t189 * t202;
t160 = t247 * t164 - t188 * t180 + t189 * t198;
t174 = -t188 * t202 + t189 * t210;
t171 = 0.1e1 / t174;
t172 = 0.1e1 / t174 ^ 2;
t284 = t160 * t171 * t172;
t252 = t228 * t262;
t220 = t230 * t235 - t252;
t266 = t228 * t229;
t244 = -t220 * t225 + t226 * t266;
t283 = -0.2e1 * t244 * t284;
t263 = t229 * t235;
t243 = -t207 * t217 + t263 * t270;
t282 = t225 * t243;
t271 = t198 * t207 * t208;
t281 = -0.2e1 * (t180 * t270 - t200 * t271) / t194 ^ 2;
t206 = t220 * t226 + t225 * t266;
t234 = cos(qJ(5));
t219 = t228 * t261 + t230 * t233;
t232 = sin(qJ(5));
t268 = t219 * t232;
t191 = t206 * t234 + t268;
t185 = 0.1e1 / t191;
t186 = 0.1e1 / t191 ^ 2;
t215 = t219 * qJD(2);
t183 = t244 * qJD(4) - t215 * t226;
t216 = -qJD(2) * t252 + t230 * t259;
t175 = t191 * qJD(5) + t183 * t232 - t216 * t234;
t267 = t219 * t234;
t190 = t206 * t232 - t267;
t184 = t190 ^ 2;
t179 = t184 * t186 + 0.1e1;
t275 = t186 * t190;
t257 = qJD(5) * t190;
t176 = t183 * t234 + t216 * t232 - t257;
t278 = t176 * t185 * t186;
t280 = (t175 * t275 - t184 * t278) / t179 ^ 2;
t279 = t172 * t244;
t182 = t206 * qJD(4) - t215 * t225;
t277 = t182 * t172;
t276 = t185 * t232;
t274 = t188 * t244;
t273 = t189 * t244;
t272 = t190 * t234;
t269 = t219 * t225;
t260 = qJD(2) * t233;
t201 = t244 ^ 2;
t170 = t201 * t172 + 0.1e1;
t256 = 0.2e1 * (-t201 * t284 - t244 * t277) / t170 ^ 2;
t255 = -0.2e1 * t280;
t253 = t190 * t278;
t249 = -0.2e1 * t202 * t271;
t248 = qJD(5) * t219 * t226 - t215;
t246 = t186 * t272 - t276;
t204 = t218 * t226 - t251;
t245 = -t204 * t207 + t211 * t270;
t242 = qJD(4) * t269 + qJD(5) * t220 - t216 * t226;
t214 = t218 * qJD(2);
t199 = -t210 * qJD(4) + t226 * t250;
t197 = t220 * t232 - t226 * t267;
t196 = -t220 * t234 - t226 * t268;
t181 = -t202 * qJD(4) + t213 * t226;
t177 = 0.1e1 / t179;
t167 = 0.1e1 / t170;
t166 = t192 * t282;
t165 = t245 * t192;
t162 = (-t188 * t217 + t189 * t263) * t225 + t247 * t166;
t161 = t247 * t165 - t188 * t204 + t189 * t211;
t159 = t245 * t281 + (t211 * t249 - t181 * t207 + (t180 * t211 + t198 * t204 + t199 * t202) * t208) * t192;
t157 = t281 * t282 + (t243 * t258 + (t249 * t263 + t207 * t214 + (t198 * t217 + (t180 * t235 - t202 * t260) * t229) * t208) * t225) * t192;
t1 = [0, t157, 0, t159, 0, 0; 0 (-t162 * t279 + t171 * t269) * t256 + ((-t216 * t225 - t219 * t258) * t171 + (-t277 + t283) * t162 + (t269 * t160 + (-t157 * t202 - t166 * t180 + (-t225 * t260 + t235 * t258) * t229 + (-t166 * t210 - t217 * t225) * t164) * t273 + (-t217 * t258 - t157 * t210 - t166 * t198 + t214 * t225 + (t166 * t202 - t225 * t263) * t164) * t274) * t172) * t167, 0 (-t161 * t279 - t171 * t206) * t256 + (t161 * t283 + t183 * t171 + (-t206 * t160 - t161 * t182 + (-t159 * t202 - t165 * t180 + t199 + (-t165 * t210 - t204) * t164) * t273 + (-t159 * t210 - t165 * t198 - t181 + (t165 * t202 - t211) * t164) * t274) * t172) * t167, 0, 0; 0, 0.2e1 * (-t185 * t196 + t197 * t275) * t280 + (0.2e1 * t197 * t253 - t248 * t185 * t234 + t242 * t276 + (-t248 * t190 * t232 - t197 * t175 - t196 * t176 - t242 * t272) * t186) * t177, 0, -t246 * t244 * t255 + (t246 * t182 - ((-qJD(5) * t185 - 0.2e1 * t253) * t234 + (t175 * t234 + (t176 - t257) * t232) * t186) * t244) * t177, t255 + 0.2e1 * (t175 * t186 * t177 + (-t177 * t278 - t186 * t280) * t190) * t190, 0;];
JaD_rot  = t1;
