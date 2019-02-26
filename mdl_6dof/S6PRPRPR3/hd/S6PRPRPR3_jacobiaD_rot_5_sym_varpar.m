% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR3
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
% Datum: 2019-02-26 19:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRPR3_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:47:44
% EndTime: 2019-02-26 19:47:45
% DurationCPUTime: 0.95s
% Computational Cost: add. (4524->94), mult. (13242->201), div. (524->12), fcn. (17302->13), ass. (0->93)
t213 = sin(pkin(11));
t216 = cos(pkin(11));
t220 = sin(qJ(2));
t222 = cos(qJ(2));
t206 = t220 * t213 - t222 * t216;
t218 = cos(pkin(6));
t230 = t206 * t218;
t200 = qJD(2) * t230;
t234 = t222 * t213 + t220 * t216;
t205 = t234 * qJD(2);
t214 = sin(pkin(10));
t217 = cos(pkin(10));
t179 = -t217 * t200 - t214 * t205;
t203 = t234 * t218;
t186 = t217 * t203 - t214 * t206;
t219 = sin(qJ(4));
t215 = sin(pkin(6));
t246 = t215 * t219;
t238 = t217 * t246;
t221 = cos(qJ(4));
t241 = qJD(4) * t221;
t155 = -qJD(4) * t238 + t179 * t219 + t186 * t241;
t245 = t215 * t221;
t173 = t186 * t219 + t217 * t245;
t170 = t173 ^ 2;
t202 = t234 * t215;
t194 = t202 * t219 - t218 * t221;
t192 = 0.1e1 / t194 ^ 2;
t166 = t170 * t192 + 0.1e1;
t164 = 0.1e1 / t166;
t195 = t202 * t221 + t218 * t219;
t201 = t206 * t215;
t199 = qJD(2) * t201;
t168 = t195 * qJD(4) - t199 * t219;
t191 = 0.1e1 / t194;
t249 = t173 * t192;
t143 = (-t155 * t191 + t168 * t249) * t164;
t167 = atan2(-t173, t194);
t162 = sin(t167);
t163 = cos(t167);
t236 = -t162 * t194 - t163 * t173;
t140 = t236 * t143 - t162 * t155 + t163 * t168;
t154 = -t162 * t173 + t163 * t194;
t151 = 0.1e1 / t154;
t152 = 0.1e1 / t154 ^ 2;
t260 = t140 * t151 * t152;
t181 = t214 * t200 - t217 * t205;
t235 = -t214 * t203 - t217 * t206;
t231 = t214 * t245 - t219 * t235;
t158 = t231 * qJD(4) + t181 * t221;
t177 = t214 * t246 + t221 * t235;
t172 = t177 ^ 2;
t188 = t214 * t230 - t217 * t234;
t183 = 0.1e1 / t188 ^ 2;
t161 = t172 * t183 + 0.1e1;
t204 = t206 * qJD(2);
t229 = t218 * t205;
t180 = t217 * t204 + t214 * t229;
t182 = 0.1e1 / t188;
t184 = t182 * t183;
t259 = 0.2e1 * (t177 * t183 * t158 - t172 * t184 * t180) / t161 ^ 2;
t258 = -0.2e1 * t231 * t260;
t185 = -t214 * t234 - t217 * t230;
t232 = -t185 * t191 - t201 * t249;
t257 = t219 * t232;
t250 = t168 * t191 * t192;
t256 = -0.2e1 * (t155 * t249 - t170 * t250) / t166 ^ 2;
t254 = t152 * t231;
t157 = t177 * qJD(4) + t181 * t219;
t253 = t157 * t152;
t252 = t162 * t231;
t251 = t163 * t231;
t248 = t177 * t235;
t247 = t188 * t219;
t171 = t231 ^ 2;
t150 = t171 * t152 + 0.1e1;
t240 = 0.2e1 * (-t171 * t260 - t231 * t253) / t150 ^ 2;
t237 = -0.2e1 * t173 * t250;
t175 = t186 * t221 - t238;
t233 = -t175 * t191 + t195 * t249;
t198 = t215 * t205;
t178 = t214 * t204 - t217 * t229;
t169 = -t194 * qJD(4) - t199 * t221;
t159 = 0.1e1 / t161;
t156 = -t173 * qJD(4) + t179 * t221;
t148 = 0.1e1 / t150;
t145 = t164 * t257;
t144 = t233 * t164;
t142 = (-t162 * t185 - t163 * t201) * t219 + t236 * t145;
t141 = t236 * t144 - t162 * t175 + t163 * t195;
t139 = t233 * t256 + (t195 * t237 - t156 * t191 + (t155 * t195 + t168 * t175 + t169 * t173) * t192) * t164;
t137 = t256 * t257 + (t232 * t241 + (-t201 * t237 - t178 * t191 + (-t155 * t201 + t168 * t185 - t173 * t198) * t192) * t219) * t164;
t1 = [0, t137, 0, t139, 0, 0; 0 (-t142 * t254 - t151 * t247) * t240 + ((t180 * t219 + t188 * t241) * t151 + (-t253 + t258) * t142 + (-t247 * t140 + (-t201 * t241 - t137 * t173 - t145 * t155 - t198 * t219 + (-t145 * t194 - t185 * t219) * t143) * t251 + (-t185 * t241 - t137 * t194 - t145 * t168 - t178 * t219 + (t145 * t173 + t201 * t219) * t143) * t252) * t152) * t148, 0 (-t141 * t254 - t151 * t177) * t240 + (t141 * t258 + t158 * t151 + (-t177 * t140 - t141 * t157 + (-t139 * t173 - t144 * t155 + t169 + (-t144 * t194 - t175) * t143) * t251 + (-t139 * t194 - t144 * t168 - t156 + (t144 * t173 - t195) * t143) * t252) * t152) * t148, 0, 0; 0 (t182 * t188 * t221 + t183 * t248) * t259 + (qJD(4) * t182 * t247 + (-t158 * t235 - t181 * t177) * t183 + (0.2e1 * t184 * t248 + (t183 * t188 - t182) * t221) * t180) * t159, 0, t231 * t182 * t259 + (t180 * t183 * t231 + t157 * t182) * t159, 0, 0;];
JaD_rot  = t1;
