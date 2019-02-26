% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPPR1_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:58:09
% EndTime: 2019-02-26 19:58:09
% DurationCPUTime: 0.90s
% Computational Cost: add. (5351->102), mult. (8196->213), div. (535->12), fcn. (10572->13), ass. (0->102)
t202 = sin(pkin(10));
t205 = cos(pkin(10));
t207 = sin(qJ(2));
t206 = cos(pkin(6));
t208 = cos(qJ(2));
t231 = t206 * t208;
t190 = -t202 * t207 + t205 * t231;
t186 = t190 * qJD(2);
t232 = t206 * t207;
t191 = t202 * t208 + t205 * t232;
t200 = qJ(3) + pkin(11);
t198 = sin(t200);
t203 = sin(pkin(6));
t235 = t203 * t205;
t222 = t198 * t235;
t199 = cos(t200);
t228 = qJD(3) * t199;
t154 = -qJD(3) * t222 + t186 * t198 + t191 * t228;
t176 = t191 * t198 + t199 * t235;
t174 = t176 ^ 2;
t234 = t203 * t207;
t184 = t198 * t234 - t206 * t199;
t182 = 0.1e1 / t184 ^ 2;
t168 = t174 * t182 + 0.1e1;
t166 = 0.1e1 / t168;
t185 = t206 * t198 + t199 * t234;
t229 = qJD(2) * t208;
t221 = t203 * t229;
t172 = t185 * qJD(3) + t198 * t221;
t181 = 0.1e1 / t184;
t240 = t176 * t182;
t138 = (-t154 * t181 + t172 * t240) * t166;
t169 = atan2(-t176, t184);
t164 = sin(t169);
t165 = cos(t169);
t219 = -t164 * t184 - t165 * t176;
t134 = t219 * t138 - t164 * t154 + t165 * t172;
t148 = -t164 * t176 + t165 * t184;
t145 = 0.1e1 / t148;
t146 = 0.1e1 / t148 ^ 2;
t254 = t134 * t145 * t146;
t223 = t202 * t232;
t193 = t205 * t208 - t223;
t236 = t202 * t203;
t180 = t193 * t199 + t198 * t236;
t201 = sin(pkin(12));
t192 = t202 * t231 + t205 * t207;
t204 = cos(pkin(12));
t237 = t192 * t204;
t162 = t180 * t201 - t237;
t188 = t192 * qJD(2);
t217 = -t193 * t198 + t199 * t236;
t157 = t217 * qJD(3) - t188 * t199;
t189 = -qJD(2) * t223 + t205 * t229;
t153 = t157 * t204 + t189 * t201;
t238 = t192 * t201;
t163 = t180 * t204 + t238;
t159 = 0.1e1 / t163;
t160 = 0.1e1 / t163 ^ 2;
t248 = t153 * t159 * t160;
t253 = 0.2e1 * t162 * t248;
t252 = -0.2e1 * t217 * t254;
t233 = t203 * t208;
t215 = -t181 * t190 + t233 * t240;
t251 = t198 * t215;
t241 = t172 * t181 * t182;
t250 = -0.2e1 * (t154 * t240 - t174 * t241) / t168 ^ 2;
t249 = t146 * t217;
t156 = t180 * qJD(3) - t188 * t198;
t247 = t156 * t146;
t246 = t159 * t201;
t245 = t160 * t162;
t244 = t162 * t204;
t243 = t164 * t217;
t242 = t165 * t217;
t239 = t192 * t198;
t230 = qJD(2) * t207;
t175 = t217 ^ 2;
t144 = t175 * t146 + 0.1e1;
t227 = 0.2e1 * (-t175 * t254 - t217 * t247) / t144 ^ 2;
t158 = t162 ^ 2;
t151 = t158 * t160 + 0.1e1;
t152 = t157 * t201 - t189 * t204;
t226 = 0.2e1 * (t152 * t245 - t158 * t248) / t151 ^ 2;
t220 = -0.2e1 * t176 * t241;
t178 = t191 * t199 - t222;
t218 = -t178 * t181 + t185 * t240;
t216 = qJD(3) * t239 - t189 * t199;
t187 = t191 * qJD(2);
t173 = -t184 * qJD(3) + t199 * t221;
t171 = t193 * t201 - t199 * t237;
t170 = -t193 * t204 - t199 * t238;
t155 = -t176 * qJD(3) + t186 * t199;
t149 = 0.1e1 / t151;
t141 = 0.1e1 / t144;
t140 = t166 * t251;
t139 = t218 * t166;
t136 = (-t164 * t190 + t165 * t233) * t198 + t219 * t140;
t135 = t219 * t139 - t164 * t178 + t165 * t185;
t133 = t218 * t250 + (t185 * t220 - t155 * t181 + (t154 * t185 + t172 * t178 + t173 * t176) * t182) * t166;
t131 = t250 * t251 + (t215 * t228 + (t220 * t233 + t181 * t187 + (t172 * t190 + (t154 * t208 - t176 * t230) * t203) * t182) * t198) * t166;
t1 = [0, t131, t133, 0, 0, 0; 0 (-t136 * t249 + t145 * t239) * t227 + ((-t189 * t198 - t192 * t228) * t145 + (-t247 + t252) * t136 + (t239 * t134 + (-t131 * t176 - t140 * t154 + (-t198 * t230 + t208 * t228) * t203 + (-t140 * t184 - t190 * t198) * t138) * t242 + (-t190 * t228 - t131 * t184 - t140 * t172 + t187 * t198 + (t140 * t176 - t198 * t233) * t138) * t243) * t146) * t141 (-t135 * t249 - t145 * t180) * t227 + (t135 * t252 + t157 * t145 + (-t180 * t134 - t135 * t156 + (-t133 * t176 - t139 * t154 + t173 + (-t139 * t184 - t178) * t138) * t242 + (-t133 * t184 - t139 * t172 - t155 + (t139 * t176 - t185) * t138) * t243) * t146) * t141, 0, 0, 0; 0 (-t159 * t170 + t171 * t245) * t226 + ((t188 * t204 + t216 * t201) * t159 + t171 * t253 + (-t170 * t153 - (-t188 * t201 + t216 * t204) * t162 - t171 * t152) * t160) * t149 -(-t160 * t244 + t246) * t217 * t226 + (t217 * t204 * t253 - t156 * t246 + (t156 * t244 - (t152 * t204 + t153 * t201) * t217) * t160) * t149, 0, 0, 0;];
JaD_rot  = t1;
