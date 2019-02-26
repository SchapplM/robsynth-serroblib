% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPPR3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:59:24
% EndTime: 2019-02-26 19:59:25
% DurationCPUTime: 0.94s
% Computational Cost: add. (3002->107), mult. (9085->220), div. (559->12), fcn. (11668->13), ass. (0->103)
t211 = cos(pkin(6));
t217 = cos(qJ(2));
t266 = cos(pkin(10));
t238 = t266 * t217;
t209 = sin(pkin(10));
t214 = sin(qJ(2));
t252 = t209 * t214;
t200 = t211 * t238 - t252;
t193 = t200 * qJD(2);
t216 = cos(qJ(3));
t213 = sin(qJ(3));
t239 = t266 * t214;
t251 = t209 * t217;
t226 = -t211 * t239 - t251;
t210 = sin(pkin(6));
t240 = t210 * t266;
t268 = t226 * t213 - t216 * t240;
t165 = t268 * qJD(3) + t193 * t216;
t185 = -t213 * t240 - t226 * t216;
t182 = t185 ^ 2;
t249 = t210 * t216;
t204 = t211 * t213 + t214 * t249;
t198 = 0.1e1 / t204 ^ 2;
t178 = t182 * t198 + 0.1e1;
t176 = 0.1e1 / t178;
t250 = t210 * t213;
t203 = t211 * t216 - t214 * t250;
t248 = t210 * t217;
t241 = qJD(2) * t248;
t190 = t203 * qJD(3) + t216 * t241;
t197 = 0.1e1 / t204;
t257 = t185 * t198;
t148 = (-t165 * t197 + t190 * t257) * t176;
t179 = atan2(-t185, t204);
t174 = sin(t179);
t175 = cos(t179);
t233 = -t174 * t204 - t175 * t185;
t144 = t233 * t148 - t174 * t165 + t175 * t190;
t158 = -t174 * t185 + t175 * t204;
t155 = 0.1e1 / t158;
t156 = 0.1e1 / t158 ^ 2;
t271 = t144 * t155 * t156;
t228 = t211 * t252 - t238;
t188 = t209 * t250 - t216 * t228;
t270 = 0.2e1 * t188 * t271;
t229 = -t197 * t200 + t248 * t257;
t269 = t216 * t229;
t256 = t190 * t197 * t198;
t267 = -0.2e1 * (t165 * t257 - t182 * t256) / t178 ^ 2;
t212 = sin(qJ(6));
t215 = cos(qJ(6));
t227 = t211 * t251 + t239;
t230 = t209 * t249 + t213 * t228;
t173 = -t212 * t227 - t215 * t230;
t169 = 0.1e1 / t173;
t170 = 0.1e1 / t173 ^ 2;
t265 = t156 * t188;
t195 = t227 * qJD(2);
t166 = t188 * qJD(3) - t195 * t213;
t196 = t228 * qJD(2);
t254 = t227 * t215;
t172 = -t212 * t230 + t254;
t245 = qJD(6) * t172;
t160 = t166 * t215 + t196 * t212 - t245;
t264 = t160 * t169 * t170;
t159 = t173 * qJD(6) + t166 * t212 - t196 * t215;
t168 = t172 ^ 2;
t163 = t168 * t170 + 0.1e1;
t261 = t170 * t172;
t263 = 0.1e1 / t163 ^ 2 * (t159 * t261 - t168 * t264);
t262 = t169 * t212;
t260 = t172 * t215;
t259 = t174 * t188;
t258 = t175 * t188;
t255 = t227 * t213;
t253 = t227 * t216;
t247 = qJD(2) * t214;
t246 = qJD(3) * t213;
t183 = t188 ^ 2;
t154 = t156 * t183 + 0.1e1;
t167 = t230 * qJD(3) - t195 * t216;
t244 = 0.2e1 * (t167 * t265 - t183 * t271) / t154 ^ 2;
t242 = 0.2e1 * t263;
t237 = 0.2e1 * t172 * t264;
t236 = -0.2e1 * t185 * t256;
t234 = -qJD(6) * t255 - t195;
t232 = t170 * t260 - t262;
t231 = -t197 * t268 + t203 * t257;
t225 = -qJD(3) * t253 + qJD(6) * t228 + t196 * t213;
t194 = t226 * qJD(2);
t189 = -t204 * qJD(3) - t213 * t241;
t181 = t212 * t228 - t213 * t254;
t180 = -t212 * t255 - t215 * t228;
t164 = t185 * qJD(3) + t193 * t213;
t161 = 0.1e1 / t163;
t151 = 0.1e1 / t154;
t150 = t176 * t269;
t149 = t231 * t176;
t146 = (-t174 * t200 + t175 * t248) * t216 + t233 * t150;
t145 = t233 * t149 - t174 * t268 + t175 * t203;
t143 = t231 * t267 + (t203 * t236 + t164 * t197 + (t165 * t203 + t185 * t189 + t190 * t268) * t198) * t176;
t141 = t267 * t269 + (-t229 * t246 + (t236 * t248 - t194 * t197 + (t190 * t200 + (t165 * t217 - t185 * t247) * t210) * t198) * t216) * t176;
t1 = [0, t141, t143, 0, 0, 0; 0 (t146 * t265 + t155 * t253) * t244 + ((t196 * t216 + t227 * t246) * t155 + t146 * t270 + (-t146 * t167 + t253 * t144 - (-t141 * t185 - t150 * t165 + (-t216 * t247 - t217 * t246) * t210 + (-t150 * t204 - t200 * t216) * t148) * t258 - (t200 * t246 - t141 * t204 - t150 * t190 - t194 * t216 + (t150 * t185 - t216 * t248) * t148) * t259) * t156) * t151 (t145 * t265 - t155 * t230) * t244 + (t145 * t270 - t166 * t155 + (-t230 * t144 - t145 * t167 - (-t143 * t185 - t149 * t165 + t189 + (-t149 * t204 - t268) * t148) * t258 - (-t143 * t204 - t149 * t190 + t164 + (t149 * t185 - t203) * t148) * t259) * t156) * t151, 0, 0, 0; 0 (-t169 * t180 + t181 * t261) * t242 + (t181 * t237 + t234 * t169 * t215 + t225 * t262 + (t234 * t172 * t212 - t181 * t159 - t180 * t160 - t225 * t260) * t170) * t161, t232 * t188 * t242 + (-t232 * t167 + ((qJD(6) * t169 + t237) * t215 + (-t159 * t215 + (-t160 + t245) * t212) * t170) * t188) * t161, 0, 0, -0.2e1 * t263 + 0.2e1 * (t159 * t170 * t161 + (-t161 * t264 - t170 * t263) * t172) * t172;];
JaD_rot  = t1;
