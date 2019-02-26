% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRR1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:47:03
% EndTime: 2019-02-26 22:47:04
% DurationCPUTime: 0.86s
% Computational Cost: add. (13377->97), mult. (6392->208), div. (1299->12), fcn. (7429->9), ass. (0->95)
t185 = sin(qJ(1));
t247 = 0.2e1 * t185;
t182 = t185 ^ 2;
t181 = qJ(2) + qJ(3) + qJ(4) + qJ(5);
t177 = sin(t181);
t172 = t177 ^ 2;
t178 = cos(t181);
t174 = 0.1e1 / t178 ^ 2;
t232 = t172 * t174;
t169 = t182 * t232 + 0.1e1;
t167 = 0.1e1 / t169;
t173 = 0.1e1 / t178;
t187 = cos(qJ(1));
t218 = qJD(1) * t187;
t208 = t177 * t218;
t179 = qJD(2) + qJD(3) + qJD(4) + qJD(5);
t226 = t179 * t185;
t211 = t174 * t226;
t141 = (-(-t178 * t226 - t208) * t173 + t172 * t211) * t167;
t246 = t141 - t226;
t186 = cos(qJ(6));
t220 = t186 * t187;
t184 = sin(qJ(6));
t222 = t185 * t184;
t166 = t178 * t220 + t222;
t223 = t185 * t177;
t158 = atan2(-t223, -t178);
t157 = cos(t158);
t156 = sin(t158);
t212 = t156 * t223;
t151 = -t157 * t178 - t212;
t148 = 0.1e1 / t151;
t160 = 0.1e1 / t166;
t149 = 0.1e1 / t151 ^ 2;
t161 = 0.1e1 / t166 ^ 2;
t245 = t167 - 0.1e1;
t237 = t157 * t177;
t136 = (-t141 * t185 + t179) * t237 + (t246 * t178 - t208) * t156;
t244 = t136 * t148 * t149;
t171 = t177 * t172;
t229 = t173 * t177;
t195 = t179 * (t171 * t173 * t174 + t229);
t230 = t172 * t185;
t200 = t218 * t230;
t243 = (t174 * t200 + t182 * t195) / t169 ^ 2;
t202 = -qJD(1) * t178 + qJD(6);
t203 = qJD(6) * t178 - qJD(1);
t225 = t179 * t187;
t210 = t177 * t225;
t224 = t184 * t187;
t147 = -t203 * t224 + (t202 * t185 - t210) * t186;
t242 = t147 * t160 * t161;
t241 = t149 * t177;
t240 = t149 * t187;
t196 = t178 * t222 + t220;
t146 = t196 * qJD(1) - t166 * qJD(6) + t184 * t210;
t221 = t185 * t186;
t165 = t178 * t224 - t221;
t159 = t165 ^ 2;
t155 = t159 * t161 + 0.1e1;
t235 = t161 * t165;
t239 = 0.1e1 / t155 ^ 2 * (-t146 * t235 - t159 * t242);
t238 = t156 * t185;
t236 = t160 * t184;
t234 = t165 * t186;
t233 = t172 * t173;
t183 = t187 ^ 2;
t231 = t172 * t183;
t228 = t177 * t187;
t227 = t178 * t179;
t219 = qJD(1) * t185;
t144 = t149 * t231 + 0.1e1;
t217 = 0.2e1 * (-t231 * t244 + (t177 * t183 * t227 - t200) * t149) / t144 ^ 2;
t216 = 0.2e1 * t244;
t215 = -0.2e1 * t239;
t214 = t165 * t242;
t213 = t149 * t228;
t207 = 0.1e1 + t232;
t206 = t177 * t217;
t205 = -0.2e1 * t177 * t243;
t204 = t243 * t247;
t201 = t157 * t167 * t233;
t199 = t207 * t187;
t198 = t202 * t187;
t197 = t161 * t234 - t236;
t164 = -t178 * t221 + t224;
t153 = 0.1e1 / t155;
t152 = t207 * t185 * t167;
t142 = 0.1e1 / t144;
t140 = (t245 * t177 * t156 - t185 * t201) * t187;
t138 = -t178 * t238 + t237 + (t156 * t178 - t157 * t223) * t152;
t137 = -t207 * t204 + (qJD(1) * t199 + t195 * t247) * t167;
t134 = t197 * t215 * t228 + (t197 * t178 * t225 + (-t197 * t219 + ((-qJD(6) * t160 - 0.2e1 * t214) * t186 + (-t146 * t186 + (-qJD(6) * t165 + t147) * t184) * t161) * t187) * t177) * t153;
t133 = (t138 * t241 - t148 * t178) * t187 * t217 + ((-t148 * t219 + (-t138 * t179 - t136) * t240) * t178 + (-t148 * t225 - (-t137 * t157 * t185 - t246 * t156 + (t141 * t238 - t156 * t179 - t157 * t218) * t152) * t213 + (t149 * t219 + t187 * t216) * t138 - ((t137 - t218) * t156 + ((-t152 * t185 + 0.1e1) * t179 + (t152 - t185) * t141) * t157) * t178 * t240) * t177) * t142;
t1 = [t187 * t173 * t205 + (t179 * t199 - t219 * t229) * t167, t137, t137, t137, t137, 0; (t148 * t206 + (-t148 * t227 + (qJD(1) * t140 + t136) * t241) * t142) * t185 + (t149 * t206 * t140 + (-((t205 - t227 + (t141 * t173 * t230 + t227) * t167) * t156 + (t204 * t233 - t141 * t177 + (-t171 * t211 + (t141 - 0.2e1 * t226) * t177) * t167) * t157) * t213 + (-t149 * t227 + t177 * t216) * t140 + (-t148 + ((-t182 + t183) * t201 + t245 * t212) * t149) * t177 * qJD(1)) * t142) * t187, t133, t133, t133, t133, 0; 0.2e1 * (t160 * t196 + t164 * t235) * t239 + (0.2e1 * t164 * t214 - t203 * t160 * t221 + (t179 * t223 + t198) * t236 + (t164 * t146 + t196 * t147 - t198 * t234 - (t177 * t179 * t186 + t203 * t184) * t165 * t185) * t161) * t153, t134, t134, t134, t134, t215 + 0.2e1 * (-t146 * t161 * t153 + (-t153 * t242 - t161 * t239) * t165) * t165;];
JaD_rot  = t1;
