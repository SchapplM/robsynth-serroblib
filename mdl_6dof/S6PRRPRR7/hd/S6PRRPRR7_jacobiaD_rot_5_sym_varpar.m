% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRR7_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:37
% EndTime: 2019-02-26 20:07:38
% DurationCPUTime: 0.97s
% Computational Cost: add. (3002->107), mult. (9085->220), div. (559->12), fcn. (11668->13), ass. (0->103)
t212 = cos(pkin(6));
t218 = cos(qJ(2));
t267 = cos(pkin(11));
t239 = t267 * t218;
t210 = sin(pkin(11));
t215 = sin(qJ(2));
t253 = t210 * t215;
t201 = t212 * t239 - t253;
t194 = t201 * qJD(2);
t217 = cos(qJ(3));
t214 = sin(qJ(3));
t240 = t267 * t215;
t252 = t210 * t218;
t227 = -t212 * t240 - t252;
t211 = sin(pkin(6));
t241 = t211 * t267;
t269 = t227 * t214 - t217 * t241;
t166 = t269 * qJD(3) + t194 * t217;
t187 = -t214 * t241 - t227 * t217;
t184 = t187 ^ 2;
t250 = t211 * t217;
t205 = t212 * t214 + t215 * t250;
t199 = 0.1e1 / t205 ^ 2;
t179 = t184 * t199 + 0.1e1;
t177 = 0.1e1 / t179;
t251 = t211 * t214;
t204 = t212 * t217 - t215 * t251;
t249 = t211 * t218;
t242 = qJD(2) * t249;
t192 = t204 * qJD(3) + t217 * t242;
t198 = 0.1e1 / t205;
t258 = t187 * t199;
t149 = (-t166 * t198 + t192 * t258) * t177;
t180 = atan2(-t187, t205);
t175 = sin(t180);
t176 = cos(t180);
t234 = -t175 * t205 - t176 * t187;
t145 = t234 * t149 - t175 * t166 + t176 * t192;
t159 = -t175 * t187 + t176 * t205;
t156 = 0.1e1 / t159;
t157 = 0.1e1 / t159 ^ 2;
t273 = t145 * t156 * t157;
t213 = sin(qJ(5));
t216 = cos(qJ(5));
t228 = -t212 * t252 - t240;
t203 = -t212 * t253 + t239;
t230 = -t203 * t214 + t210 * t250;
t233 = t213 * t228 - t216 * t230;
t272 = t233 * qJD(5);
t190 = t203 * t217 + t210 * t251;
t271 = 0.2e1 * t190 * t273;
t229 = -t198 * t201 + t249 * t258;
t270 = t217 * t229;
t257 = t192 * t198 * t199;
t268 = -0.2e1 * (t166 * t258 - t184 * t257) / t179 ^ 2;
t255 = t228 * t216;
t174 = -t213 * t230 - t255;
t170 = 0.1e1 / t174;
t171 = 0.1e1 / t174 ^ 2;
t266 = t157 * t190;
t196 = t228 * qJD(2);
t167 = t190 * qJD(3) + t196 * t214;
t197 = t203 * qJD(2);
t161 = t167 * t213 + t197 * t216 + t272;
t265 = t161 * t170 * t171;
t160 = t174 * qJD(5) - t167 * t216 + t197 * t213;
t169 = t233 ^ 2;
t164 = t169 * t171 + 0.1e1;
t262 = t171 * t233;
t264 = 0.1e1 / t164 ^ 2 * (-t160 * t262 - t169 * t265);
t263 = t170 * t216;
t261 = t233 * t213;
t260 = t175 * t190;
t259 = t176 * t190;
t256 = t228 * t214;
t254 = t228 * t217;
t248 = qJD(2) * t215;
t247 = qJD(3) * t214;
t185 = t190 ^ 2;
t155 = t157 * t185 + 0.1e1;
t168 = t230 * qJD(3) + t196 * t217;
t246 = 0.2e1 * (t168 * t266 - t185 * t273) / t155 ^ 2;
t244 = 0.2e1 * t264;
t238 = -0.2e1 * t233 * t265;
t237 = -0.2e1 * t187 * t257;
t235 = qJD(5) * t256 + t196;
t232 = -t171 * t261 + t263;
t231 = -t198 * t269 + t204 * t258;
t226 = -qJD(3) * t254 + qJD(5) * t203 + t197 * t214;
t195 = t227 * qJD(2);
t191 = -t205 * qJD(3) - t214 * t242;
t182 = t203 * t216 + t213 * t256;
t181 = t203 * t213 - t214 * t255;
t165 = t187 * qJD(3) + t194 * t214;
t162 = 0.1e1 / t164;
t152 = 0.1e1 / t155;
t151 = t177 * t270;
t150 = t231 * t177;
t147 = (-t175 * t201 + t176 * t249) * t217 + t234 * t151;
t146 = t234 * t150 - t175 * t269 + t176 * t204;
t144 = t231 * t268 + (t204 * t237 + t165 * t198 + (t166 * t204 + t187 * t191 + t192 * t269) * t199) * t177;
t142 = t268 * t270 + (-t229 * t247 + (t237 * t249 - t195 * t198 + (t192 * t201 + (t166 * t218 - t187 * t248) * t211) * t199) * t217) * t177;
t1 = [0, t142, t144, 0, 0, 0; 0 (t147 * t266 - t156 * t254) * t246 + ((-t197 * t217 - t228 * t247) * t156 + t147 * t271 + (-t147 * t168 - t254 * t145 - (-t142 * t187 - t151 * t166 + (-t217 * t248 - t218 * t247) * t211 + (-t151 * t205 - t201 * t217) * t149) * t259 - (t201 * t247 - t142 * t205 - t151 * t192 - t195 * t217 + (t151 * t187 - t217 * t249) * t149) * t260) * t157) * t152 (t146 * t266 - t156 * t230) * t246 + (t146 * t271 - t167 * t156 + (-t230 * t145 - t146 * t168 - (-t144 * t187 - t150 * t166 + t191 + (-t150 * t205 - t269) * t149) * t259 - (-t144 * t205 - t150 * t192 + t165 + (t150 * t187 - t204) * t149) * t260) * t157) * t152, 0, 0, 0; 0 (-t170 * t181 - t182 * t262) * t244 + (t182 * t238 + t235 * t170 * t213 + t226 * t263 + (t216 * t233 * t235 - t182 * t160 - t181 * t161 - t226 * t261) * t171) * t162, t232 * t190 * t244 + (-t232 * t168 + ((qJD(5) * t170 + t238) * t213 + (-t160 * t213 + (t161 + t272) * t216) * t171) * t190) * t162, 0, -0.2e1 * t264 - 0.2e1 * (t160 * t171 * t162 - (-t162 * t265 - t171 * t264) * t233) * t233, 0;];
JaD_rot  = t1;
