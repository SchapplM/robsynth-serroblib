% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function JaD_rot = S6PRPRPR4_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:48:16
% EndTime: 2019-02-26 19:48:17
% DurationCPUTime: 0.89s
% Computational Cost: add. (5351->102), mult. (8196->213), div. (535->12), fcn. (10572->13), ass. (0->102)
t199 = sin(pkin(10));
t202 = cos(pkin(10));
t204 = sin(qJ(2));
t203 = cos(pkin(6));
t205 = cos(qJ(2));
t228 = t203 * t205;
t187 = -t199 * t204 + t202 * t228;
t183 = t187 * qJD(2);
t229 = t203 * t204;
t188 = t199 * t205 + t202 * t229;
t197 = pkin(11) + qJ(4);
t195 = sin(t197);
t200 = sin(pkin(6));
t232 = t200 * t202;
t219 = t195 * t232;
t196 = cos(t197);
t225 = qJD(4) * t196;
t151 = -qJD(4) * t219 + t183 * t195 + t188 * t225;
t173 = t188 * t195 + t196 * t232;
t171 = t173 ^ 2;
t231 = t200 * t204;
t181 = t195 * t231 - t203 * t196;
t179 = 0.1e1 / t181 ^ 2;
t165 = t171 * t179 + 0.1e1;
t163 = 0.1e1 / t165;
t182 = t203 * t195 + t196 * t231;
t226 = qJD(2) * t205;
t218 = t200 * t226;
t169 = t182 * qJD(4) + t195 * t218;
t178 = 0.1e1 / t181;
t237 = t173 * t179;
t135 = (-t151 * t178 + t169 * t237) * t163;
t166 = atan2(-t173, t181);
t161 = sin(t166);
t162 = cos(t166);
t216 = -t161 * t181 - t162 * t173;
t131 = t216 * t135 - t161 * t151 + t162 * t169;
t145 = -t161 * t173 + t162 * t181;
t142 = 0.1e1 / t145;
t143 = 0.1e1 / t145 ^ 2;
t251 = t131 * t142 * t143;
t220 = t199 * t229;
t190 = t202 * t205 - t220;
t233 = t199 * t200;
t177 = t190 * t196 + t195 * t233;
t198 = sin(pkin(12));
t189 = t199 * t228 + t202 * t204;
t201 = cos(pkin(12));
t234 = t189 * t201;
t159 = t177 * t198 - t234;
t185 = t189 * qJD(2);
t214 = -t190 * t195 + t196 * t233;
t154 = t214 * qJD(4) - t185 * t196;
t186 = -qJD(2) * t220 + t202 * t226;
t150 = t154 * t201 + t186 * t198;
t235 = t189 * t198;
t160 = t177 * t201 + t235;
t156 = 0.1e1 / t160;
t157 = 0.1e1 / t160 ^ 2;
t245 = t150 * t156 * t157;
t250 = 0.2e1 * t159 * t245;
t249 = -0.2e1 * t214 * t251;
t230 = t200 * t205;
t212 = -t178 * t187 + t230 * t237;
t248 = t195 * t212;
t238 = t169 * t178 * t179;
t247 = -0.2e1 * (t151 * t237 - t171 * t238) / t165 ^ 2;
t246 = t143 * t214;
t153 = t177 * qJD(4) - t185 * t195;
t244 = t153 * t143;
t243 = t156 * t198;
t242 = t157 * t159;
t241 = t159 * t201;
t240 = t161 * t214;
t239 = t162 * t214;
t236 = t189 * t195;
t227 = qJD(2) * t204;
t172 = t214 ^ 2;
t141 = t172 * t143 + 0.1e1;
t224 = 0.2e1 * (-t172 * t251 - t214 * t244) / t141 ^ 2;
t155 = t159 ^ 2;
t148 = t155 * t157 + 0.1e1;
t149 = t154 * t198 - t186 * t201;
t223 = 0.2e1 * (t149 * t242 - t155 * t245) / t148 ^ 2;
t217 = -0.2e1 * t173 * t238;
t175 = t188 * t196 - t219;
t215 = -t175 * t178 + t182 * t237;
t213 = qJD(4) * t236 - t186 * t196;
t184 = t188 * qJD(2);
t170 = -t181 * qJD(4) + t196 * t218;
t168 = t190 * t198 - t196 * t234;
t167 = -t190 * t201 - t196 * t235;
t152 = -t173 * qJD(4) + t183 * t196;
t146 = 0.1e1 / t148;
t138 = 0.1e1 / t141;
t137 = t163 * t248;
t136 = t215 * t163;
t133 = (-t161 * t187 + t162 * t230) * t195 + t216 * t137;
t132 = t216 * t136 - t161 * t175 + t162 * t182;
t130 = t215 * t247 + (t182 * t217 - t152 * t178 + (t151 * t182 + t169 * t175 + t170 * t173) * t179) * t163;
t128 = t247 * t248 + (t212 * t225 + (t217 * t230 + t178 * t184 + (t169 * t187 + (t151 * t205 - t173 * t227) * t200) * t179) * t195) * t163;
t1 = [0, t128, 0, t130, 0, 0; 0 (-t133 * t246 + t142 * t236) * t224 + ((-t186 * t195 - t189 * t225) * t142 + (-t244 + t249) * t133 + (t236 * t131 + (-t128 * t173 - t137 * t151 + (-t195 * t227 + t205 * t225) * t200 + (-t137 * t181 - t187 * t195) * t135) * t239 + (-t187 * t225 - t128 * t181 - t137 * t169 + t184 * t195 + (t137 * t173 - t195 * t230) * t135) * t240) * t143) * t138, 0 (-t132 * t246 - t142 * t177) * t224 + (t132 * t249 + t154 * t142 + (-t177 * t131 - t132 * t153 + (-t130 * t173 - t136 * t151 + t170 + (-t136 * t181 - t175) * t135) * t239 + (-t130 * t181 - t136 * t169 - t152 + (t136 * t173 - t182) * t135) * t240) * t143) * t138, 0, 0; 0 (-t156 * t167 + t168 * t242) * t223 + ((t185 * t201 + t213 * t198) * t156 + t168 * t250 + (-t167 * t150 - (-t185 * t198 + t213 * t201) * t159 - t168 * t149) * t157) * t146, 0 -(-t157 * t241 + t243) * t214 * t223 + (t214 * t201 * t250 - t153 * t243 + (t153 * t241 - (t149 * t201 + t150 * t198) * t214) * t157) * t146, 0, 0;];
JaD_rot  = t1;
