% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:30:46
% EndTime: 2019-02-26 22:30:47
% DurationCPUTime: 0.76s
% Computational Cost: add. (10386->97), mult. (5101->208), div. (1026->12), fcn. (5942->9), ass. (0->95)
t185 = sin(qJ(1));
t247 = 0.2e1 * t185;
t182 = t185 ^ 2;
t179 = qJ(2) + qJ(3) + qJ(4) + pkin(11);
t176 = sin(t179);
t172 = t176 ^ 2;
t177 = cos(t179);
t174 = 0.1e1 / t177 ^ 2;
t232 = t172 * t174;
t169 = t182 * t232 + 0.1e1;
t163 = 0.1e1 / t169;
t173 = 0.1e1 / t177;
t187 = cos(qJ(1));
t218 = qJD(1) * t187;
t208 = t176 * t218;
t181 = qJD(2) + qJD(3) + qJD(4);
t226 = t181 * t185;
t211 = t174 * t226;
t141 = (-(-t177 * t226 - t208) * t173 + t172 * t211) * t163;
t246 = t141 - t226;
t186 = cos(qJ(6));
t220 = t187 * t186;
t184 = sin(qJ(6));
t223 = t185 * t184;
t168 = t177 * t220 + t223;
t224 = t185 * t176;
t158 = atan2(-t224, -t177);
t157 = cos(t158);
t156 = sin(t158);
t213 = t156 * t224;
t149 = -t157 * t177 - t213;
t146 = 0.1e1 / t149;
t160 = 0.1e1 / t168;
t147 = 0.1e1 / t149 ^ 2;
t161 = 0.1e1 / t168 ^ 2;
t245 = t163 - 0.1e1;
t237 = t157 * t176;
t136 = (-t141 * t185 + t181) * t237 + (t246 * t177 - t208) * t156;
t244 = t136 * t146 * t147;
t196 = t177 * t223 + t220;
t225 = t181 * t187;
t209 = t176 * t225;
t150 = t196 * qJD(1) - t168 * qJD(6) + t184 * t209;
t221 = t187 * t184;
t222 = t185 * t186;
t167 = t177 * t221 - t222;
t159 = t167 ^ 2;
t155 = t159 * t161 + 0.1e1;
t235 = t161 * t167;
t202 = -qJD(1) * t177 + qJD(6);
t203 = qJD(6) * t177 - qJD(1);
t151 = -t203 * t221 + (t202 * t185 - t209) * t186;
t239 = t151 * t160 * t161;
t243 = (-t150 * t235 - t159 * t239) / t155 ^ 2;
t171 = t176 * t172;
t229 = t173 * t176;
t195 = t181 * (t171 * t173 * t174 + t229);
t230 = t172 * t185;
t200 = t218 * t230;
t242 = (t174 * t200 + t182 * t195) / t169 ^ 2;
t241 = t147 * t176;
t240 = t147 * t187;
t238 = t156 * t185;
t236 = t160 * t184;
t234 = t167 * t186;
t233 = t172 * t173;
t183 = t187 ^ 2;
t231 = t172 * t183;
t228 = t176 * t187;
t227 = t177 * t181;
t219 = qJD(1) * t185;
t144 = t147 * t231 + 0.1e1;
t217 = 0.2e1 * (-t231 * t244 + (t176 * t183 * t227 - t200) * t147) / t144 ^ 2;
t216 = 0.2e1 * t244;
t215 = -0.2e1 * t243;
t214 = t147 * t228;
t212 = t167 * t239;
t207 = 0.1e1 + t232;
t206 = t176 * t217;
t205 = -0.2e1 * t176 * t242;
t204 = t242 * t247;
t201 = t157 * t163 * t233;
t199 = t207 * t187;
t198 = t202 * t187;
t197 = t161 * t234 - t236;
t166 = -t177 * t222 + t221;
t153 = 0.1e1 / t155;
t152 = t207 * t185 * t163;
t142 = 0.1e1 / t144;
t140 = (t245 * t176 * t156 - t185 * t201) * t187;
t138 = -t177 * t238 + t237 + (t156 * t177 - t157 * t224) * t152;
t137 = -t207 * t204 + (qJD(1) * t199 + t195 * t247) * t163;
t134 = t197 * t215 * t228 + (t197 * t177 * t225 + (-t197 * t219 + ((-qJD(6) * t160 - 0.2e1 * t212) * t186 + (-t150 * t186 + (-qJD(6) * t167 + t151) * t184) * t161) * t187) * t176) * t153;
t133 = (t138 * t241 - t146 * t177) * t187 * t217 + ((-t146 * t219 + (-t138 * t181 - t136) * t240) * t177 + (-t146 * t225 - (-t137 * t157 * t185 - t246 * t156 + (t141 * t238 - t156 * t181 - t157 * t218) * t152) * t214 + (t147 * t219 + t187 * t216) * t138 - ((t137 - t218) * t156 + ((-t152 * t185 + 0.1e1) * t181 + (t152 - t185) * t141) * t157) * t177 * t240) * t176) * t142;
t1 = [t187 * t173 * t205 + (t181 * t199 - t219 * t229) * t163, t137, t137, t137, 0, 0; (t146 * t206 + (-t146 * t227 + (qJD(1) * t140 + t136) * t241) * t142) * t185 + (t147 * t206 * t140 + (-((t205 - t227 + (t141 * t173 * t230 + t227) * t163) * t156 + (t204 * t233 - t141 * t176 + (-t171 * t211 + (t141 - 0.2e1 * t226) * t176) * t163) * t157) * t214 + (-t147 * t227 + t176 * t216) * t140 + (-t146 + ((-t182 + t183) * t201 + t245 * t213) * t147) * t176 * qJD(1)) * t142) * t187, t133, t133, t133, 0, 0; 0.2e1 * (t160 * t196 + t166 * t235) * t243 + (0.2e1 * t166 * t212 - t203 * t160 * t222 + (t181 * t224 + t198) * t236 + (t166 * t150 + t196 * t151 - t198 * t234 - (t176 * t181 * t186 + t203 * t184) * t167 * t185) * t161) * t153, t134, t134, t134, 0, t215 + 0.2e1 * (-t150 * t161 * t153 + (-t153 * t239 - t161 * t243) * t167) * t167;];
JaD_rot  = t1;
