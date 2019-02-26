% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRPR4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:48:15
% EndTime: 2019-02-26 19:48:16
% DurationCPUTime: 0.96s
% Computational Cost: add. (6115->111), mult. (9085->227), div. (559->12), fcn. (11668->13), ass. (0->106)
t230 = sin(pkin(10));
t232 = cos(pkin(10));
t234 = sin(qJ(2));
t233 = cos(pkin(6));
t235 = cos(qJ(2));
t261 = t233 * t235;
t216 = -t230 * t234 + t232 * t261;
t212 = t216 * qJD(2);
t262 = t233 * t234;
t217 = t230 * t235 + t232 * t262;
t229 = pkin(11) + qJ(4);
t225 = sin(t229);
t231 = sin(pkin(6));
t265 = t231 * t232;
t251 = t225 * t265;
t227 = cos(t229);
t258 = qJD(4) * t227;
t185 = -qJD(4) * t251 + t212 * t225 + t217 * t258;
t201 = t217 * t225 + t227 * t265;
t199 = t201 ^ 2;
t264 = t231 * t234;
t210 = t225 * t264 - t233 * t227;
t208 = 0.1e1 / t210 ^ 2;
t193 = t199 * t208 + 0.1e1;
t191 = 0.1e1 / t193;
t211 = t233 * t225 + t227 * t264;
t259 = qJD(2) * t235;
t250 = t231 * t259;
t197 = t211 * qJD(4) + t225 * t250;
t207 = 0.1e1 / t210;
t269 = t201 * t208;
t163 = (-t185 * t207 + t197 * t269) * t191;
t194 = atan2(-t201, t210);
t189 = sin(t194);
t190 = cos(t194);
t247 = -t189 * t210 - t190 * t201;
t159 = t247 * t163 - t189 * t185 + t190 * t197;
t173 = -t189 * t201 + t190 * t210;
t170 = 0.1e1 / t173;
t171 = 0.1e1 / t173 ^ 2;
t283 = t159 * t170 * t171;
t252 = t230 * t262;
t219 = t232 * t235 - t252;
t266 = t230 * t231;
t244 = -t219 * t225 + t227 * t266;
t282 = -0.2e1 * t244 * t283;
t263 = t231 * t235;
t243 = -t207 * t216 + t263 * t269;
t281 = t225 * t243;
t270 = t197 * t207 * t208;
t280 = -0.2e1 * (t185 * t269 - t199 * t270) / t193 ^ 2;
t205 = t219 * t227 + t225 * t266;
t218 = t230 * t261 + t232 * t234;
t228 = pkin(12) + qJ(6);
t224 = sin(t228);
t226 = cos(t228);
t184 = t205 * t226 + t218 * t224;
t180 = 0.1e1 / t184;
t181 = 0.1e1 / t184 ^ 2;
t214 = t218 * qJD(2);
t188 = t244 * qJD(4) - t214 * t227;
t215 = -qJD(2) * t252 + t232 * t259;
t174 = t184 * qJD(6) + t188 * t224 - t215 * t226;
t183 = t205 * t224 - t218 * t226;
t179 = t183 ^ 2;
t178 = t179 * t181 + 0.1e1;
t275 = t181 * t183;
t257 = qJD(6) * t183;
t175 = t188 * t226 + t215 * t224 - t257;
t277 = t175 * t180 * t181;
t279 = (t174 * t275 - t179 * t277) / t178 ^ 2;
t278 = t171 * t244;
t276 = t180 * t224;
t274 = t183 * t226;
t187 = t205 * qJD(4) - t214 * t225;
t273 = t187 * t171;
t272 = t189 * t244;
t271 = t190 * t244;
t268 = t218 * t225;
t267 = t218 * t227;
t260 = qJD(2) * t234;
t200 = t244 ^ 2;
t169 = t200 * t171 + 0.1e1;
t256 = 0.2e1 * (-t200 * t283 - t244 * t273) / t169 ^ 2;
t255 = -0.2e1 * t279;
t253 = t183 * t277;
t249 = -0.2e1 * t201 * t270;
t248 = qJD(6) * t267 - t214;
t246 = t181 * t274 - t276;
t203 = t217 * t227 - t251;
t245 = -t203 * t207 + t211 * t269;
t242 = qJD(4) * t268 + qJD(6) * t219 - t215 * t227;
t213 = t217 * qJD(2);
t198 = -t210 * qJD(4) + t227 * t250;
t196 = t219 * t224 - t226 * t267;
t195 = -t219 * t226 - t224 * t267;
t186 = -t201 * qJD(4) + t212 * t227;
t176 = 0.1e1 / t178;
t166 = 0.1e1 / t169;
t165 = t191 * t281;
t164 = t245 * t191;
t161 = (-t189 * t216 + t190 * t263) * t225 + t247 * t165;
t160 = t247 * t164 - t189 * t203 + t190 * t211;
t158 = t245 * t280 + (t211 * t249 - t186 * t207 + (t185 * t211 + t197 * t203 + t198 * t201) * t208) * t191;
t156 = t280 * t281 + (t243 * t258 + (t249 * t263 + t207 * t213 + (t197 * t216 + (t185 * t235 - t201 * t260) * t231) * t208) * t225) * t191;
t1 = [0, t156, 0, t158, 0, 0; 0 (-t161 * t278 + t170 * t268) * t256 + ((-t215 * t225 - t218 * t258) * t170 + (-t273 + t282) * t161 + (t268 * t159 + (-t156 * t201 - t165 * t185 + (-t225 * t260 + t235 * t258) * t231 + (-t165 * t210 - t216 * t225) * t163) * t271 + (-t216 * t258 - t156 * t210 - t165 * t197 + t213 * t225 + (t165 * t201 - t225 * t263) * t163) * t272) * t171) * t166, 0 (-t160 * t278 - t170 * t205) * t256 + (t160 * t282 + t188 * t170 + (-t205 * t159 - t160 * t187 + (-t158 * t201 - t164 * t185 + t198 + (-t164 * t210 - t203) * t163) * t271 + (-t158 * t210 - t164 * t197 - t186 + (t164 * t201 - t211) * t163) * t272) * t171) * t166, 0, 0; 0, 0.2e1 * (-t180 * t195 + t196 * t275) * t279 + (0.2e1 * t196 * t253 - t248 * t180 * t226 + t242 * t276 + (-t248 * t183 * t224 - t196 * t174 - t195 * t175 - t242 * t274) * t181) * t176, 0, -t246 * t244 * t255 + (t246 * t187 - ((-qJD(6) * t180 - 0.2e1 * t253) * t226 + (t174 * t226 + (t175 - t257) * t224) * t181) * t244) * t176, 0, t255 + 0.2e1 * (t174 * t181 * t176 + (-t176 * t277 - t181 * t279) * t183) * t183;];
JaD_rot  = t1;
