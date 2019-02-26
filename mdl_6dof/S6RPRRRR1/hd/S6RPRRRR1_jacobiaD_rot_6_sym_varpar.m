% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRR1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:14:51
% EndTime: 2019-02-26 21:14:52
% DurationCPUTime: 0.84s
% Computational Cost: add. (9710->97), mult. (5101->206), div. (1026->12), fcn. (5942->9), ass. (0->96)
t185 = qJ(1) + pkin(11);
t181 = sin(t185);
t248 = 0.2e1 * t181;
t179 = t181 ^ 2;
t184 = qJ(3) + qJ(4) + qJ(5);
t177 = sin(t184);
t173 = t177 ^ 2;
t178 = cos(t184);
t175 = 0.1e1 / t178 ^ 2;
t234 = t173 * t175;
t168 = t179 * t234 + 0.1e1;
t166 = 0.1e1 / t168;
t174 = 0.1e1 / t178;
t182 = cos(t185);
t219 = qJD(1) * t182;
t208 = t177 * t219;
t183 = qJD(3) + qJD(4) + qJD(5);
t226 = t181 * t183;
t212 = t175 * t226;
t140 = (-(-t178 * t226 - t208) * t174 + t173 * t212) * t166;
t247 = t140 - t226;
t187 = cos(qJ(6));
t222 = t182 * t187;
t186 = sin(qJ(6));
t225 = t181 * t186;
t162 = t178 * t222 + t225;
t203 = qJD(6) * t178 - qJD(1);
t229 = t177 * t183;
t246 = t203 * t186 + t187 * t229;
t227 = t181 * t177;
t165 = atan2(-t227, -t178);
t164 = cos(t165);
t163 = sin(t165);
t213 = t163 * t227;
t150 = -t164 * t178 - t213;
t147 = 0.1e1 / t150;
t156 = 0.1e1 / t162;
t148 = 0.1e1 / t150 ^ 2;
t157 = 0.1e1 / t162 ^ 2;
t245 = t166 - 0.1e1;
t236 = t164 * t177;
t135 = (-t140 * t181 + t183) * t236 + (t247 * t178 - t208) * t163;
t244 = t135 * t147 * t148;
t202 = -qJD(1) * t178 + qJD(6);
t198 = t202 * t187;
t145 = t181 * t198 - t246 * t182;
t243 = t145 * t156 * t157;
t172 = t177 * t173;
t231 = t174 * t177;
t195 = (t172 * t174 * t175 + t231) * t183;
t232 = t173 * t181;
t200 = t219 * t232;
t242 = (t175 * t200 + t179 * t195) / t168 ^ 2;
t241 = t148 * t177;
t240 = t148 * t182;
t196 = t178 * t225 + t222;
t211 = t186 * t229;
t144 = t196 * qJD(1) - t162 * qJD(6) + t182 * t211;
t223 = t182 * t186;
t224 = t181 * t187;
t161 = t178 * t223 - t224;
t155 = t161 ^ 2;
t154 = t155 * t157 + 0.1e1;
t238 = t157 * t161;
t239 = 0.1e1 / t154 ^ 2 * (-t144 * t238 - t155 * t243);
t237 = t163 * t181;
t235 = t173 * t174;
t180 = t182 ^ 2;
t233 = t173 * t180;
t230 = t177 * t182;
t228 = t178 * t183;
t221 = t183 * t147;
t220 = qJD(1) * t181;
t143 = t148 * t233 + 0.1e1;
t218 = 0.2e1 * (-t233 * t244 + (t177 * t180 * t228 - t200) * t148) / t143 ^ 2;
t217 = 0.2e1 * t244;
t216 = -0.2e1 * t239;
t215 = t161 * t243;
t214 = t148 * t230;
t207 = 0.1e1 + t234;
t206 = t177 * t218;
t205 = -0.2e1 * t177 * t242;
t204 = t242 * t248;
t201 = t164 * t166 * t235;
t199 = t207 * t182;
t197 = -t156 * t186 + t187 * t238;
t160 = -t178 * t224 + t223;
t152 = 0.1e1 / t154;
t151 = t207 * t181 * t166;
t141 = 0.1e1 / t143;
t139 = (t245 * t177 * t163 - t181 * t201) * t182;
t137 = -t178 * t237 + t236 + (t163 * t178 - t164 * t227) * t151;
t136 = -t207 * t204 + (qJD(1) * t199 + t195 * t248) * t166;
t133 = t197 * t216 * t230 + (t197 * t182 * t228 + (-t197 * t220 + ((-qJD(6) * t156 - 0.2e1 * t215) * t187 + (-t144 * t187 + (-qJD(6) * t161 + t145) * t186) * t157) * t182) * t177) * t152;
t132 = (t137 * t241 - t147 * t178) * t182 * t218 + ((-t147 * t220 + (-t137 * t183 - t135) * t240) * t178 + (-t182 * t221 - (-t136 * t164 * t181 - t247 * t163 + (t140 * t237 - t163 * t183 - t164 * t219) * t151) * t214 + (t148 * t220 + t182 * t217) * t137 - ((t136 - t219) * t163 + ((-t151 * t181 + 0.1e1) * t183 + (t151 - t181) * t140) * t164) * t178 * t240) * t177) * t141;
t1 = [t182 * t174 * t205 + (t183 * t199 - t220 * t231) * t166, 0, t136, t136, t136, 0; (t147 * t206 + (-t178 * t221 + (qJD(1) * t139 + t135) * t241) * t141) * t181 + (t148 * t206 * t139 + (-((t205 - t228 + (t140 * t174 * t232 + t228) * t166) * t163 + (t204 * t235 - t140 * t177 + (-t172 * t212 + (t140 - 0.2e1 * t226) * t177) * t166) * t164) * t214 + (-t148 * t228 + t177 * t217) * t139 + (-t147 + ((-t179 + t180) * t201 + t245 * t213) * t148) * t177 * qJD(1)) * t141) * t182, 0, t132, t132, t132, 0; 0.2e1 * (t156 * t196 + t160 * t238) * t239 + (0.2e1 * t160 * t215 + (t160 * t144 + t196 * t145 + (-t246 * t181 - t182 * t198) * t161) * t157 + (t202 * t223 + (-t203 * t187 + t211) * t181) * t156) * t152, 0, t133, t133, t133, t216 + 0.2e1 * (-t144 * t157 * t152 + (-t152 * t243 - t157 * t239) * t161) * t161;];
JaD_rot  = t1;
