% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRP1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:24:50
% EndTime: 2019-02-26 21:24:51
% DurationCPUTime: 0.99s
% Computational Cost: add. (6874->125), mult. (6168->272), div. (1114->15), fcn. (7752->9), ass. (0->115)
t176 = qJ(2) + pkin(9);
t174 = cos(t176);
t175 = pkin(10) + qJ(5);
t171 = sin(t175);
t249 = sin(qJ(1));
t210 = t249 * t171;
t173 = cos(t175);
t178 = cos(qJ(1));
t230 = t178 * t173;
t153 = t174 * t230 + t210;
t147 = 0.1e1 / t153 ^ 2;
t172 = sin(t176);
t167 = t172 ^ 2;
t177 = t178 ^ 2;
t234 = t167 * t177;
t215 = t147 * t234;
t143 = 0.1e1 + t215;
t202 = qJD(1) * t249;
t227 = qJD(2) * t178;
t206 = t172 * t227;
t188 = t174 * t202 + t206;
t201 = t249 * qJD(5);
t231 = t178 * t171;
t132 = (-qJD(5) * t174 + qJD(1)) * t231 + (t201 - t188) * t173;
t146 = 0.1e1 / t153;
t244 = t132 * t146 * t147;
t196 = t234 * t244;
t207 = qJD(2) * t172 * t177;
t252 = (-t196 + (-t167 * t178 * t202 + t174 * t207) * t147) / t143 ^ 2;
t232 = t172 * t178;
t149 = t174 * t210 + t230;
t193 = t171 * t201;
t224 = qJD(5) * t178;
t204 = t173 * t224;
t131 = t149 * qJD(1) + t171 * t206 - t174 * t204 - t193;
t209 = t249 * t173;
t152 = t174 * t231 - t209;
t164 = 0.1e1 / t171;
t165 = 0.1e1 / t171 ^ 2;
t168 = 0.1e1 / t172;
t169 = 0.1e1 / t172 ^ 2;
t228 = qJD(2) * t174;
t208 = t169 * t228;
t225 = qJD(5) * t173;
t237 = t164 * t168;
t251 = (t165 * t168 * t225 + t164 * t208) * t152 + t131 * t237;
t233 = t172 * t171;
t139 = atan2(-t149, t233);
t136 = cos(t139);
t135 = sin(t139);
t243 = t135 * t149;
t130 = t136 * t233 - t243;
t127 = 0.1e1 / t130;
t128 = 0.1e1 / t130 ^ 2;
t250 = 0.2e1 * t152;
t144 = t149 ^ 2;
t236 = t165 * t169;
t140 = t144 * t236 + 0.1e1;
t137 = 0.1e1 / t140;
t189 = t171 * t228 + t172 * t225;
t213 = t149 * t236;
t211 = t172 * t249;
t194 = qJD(2) * t211;
t195 = t173 * t202;
t229 = qJD(1) * t178;
t133 = t173 * t201 * t174 - t195 + (t229 * t174 - t194 - t224) * t171;
t216 = t133 * t237;
t119 = (t189 * t213 - t216) * t137;
t186 = -t119 * t149 + t189;
t115 = (-t119 * t233 - t133) * t135 + t186 * t136;
t129 = t127 * t128;
t248 = t115 * t129;
t166 = t164 * t165;
t170 = t168 / t167;
t205 = t169 * t225;
t247 = (t133 * t213 + (-t165 * t170 * t228 - t166 * t205) * t144) / t140 ^ 2;
t246 = t128 * t152;
t245 = t131 * t128;
t242 = t135 * t152;
t241 = t135 * t172;
t240 = t136 * t149;
t239 = t136 * t152;
t238 = t136 * t174;
t235 = t165 * t173;
t226 = qJD(5) * t171;
t145 = t152 ^ 2;
t125 = t128 * t145 + 0.1e1;
t223 = 0.2e1 * (-t145 * t248 - t152 * t245) / t125 ^ 2;
t222 = -0.2e1 * t247;
t221 = 0.2e1 * t252;
t220 = t129 * t250;
t219 = t168 * t247;
t218 = t128 * t242;
t214 = t149 * t237;
t212 = t164 * t169 * t174;
t191 = t149 * t212 + t249;
t126 = t191 * t137;
t203 = t249 - t126;
t200 = t127 * t223;
t199 = t128 * t223;
t198 = t232 * t250;
t197 = t164 * t219;
t151 = t174 * t209 - t231;
t192 = t149 * t235 - t151 * t164;
t190 = t147 * t151 * t178 - t249 * t146;
t141 = 0.1e1 / t143;
t134 = t153 * qJD(1) - t173 * t194 - t174 * t193 - t204;
t123 = 0.1e1 / t125;
t122 = t192 * t168 * t137;
t118 = (-t135 + (t136 * t214 + t135) * t137) * t152;
t117 = -t126 * t240 + (t203 * t241 + t238) * t171;
t116 = t136 * t172 * t173 - t135 * t151 + (-t135 * t233 - t240) * t122;
t114 = t191 * t222 + (t133 * t212 + t229 + (-t165 * t174 * t205 + (-0.2e1 * t170 * t174 ^ 2 - t168) * t164 * qJD(2)) * t149) * t137;
t112 = -0.2e1 * t192 * t219 + (-t192 * t208 + (t133 * t235 - t134 * t164 + (t151 * t235 + (-0.2e1 * t166 * t173 ^ 2 - t164) * t149) * qJD(5)) * t168) * t137;
t1 = [t251 * t137 + t197 * t250, t114, 0, 0, t112, 0; t149 * t200 + (-t133 * t127 + (t115 * t149 + t118 * t131) * t128) * t123 + (t118 * t199 + (0.2e1 * t118 * t248 + (t131 * t137 - t131 - (-t119 * t137 * t214 + t222) * t152) * t128 * t135 + (-(-0.2e1 * t149 * t197 - t119) * t246 + (-(t119 + t216) * t152 + t251 * t149) * t128 * t137) * t136) * t123) * t152, t117 * t152 * t199 + (-(-t114 * t240 + (t119 * t243 - t133 * t136) * t126) * t246 + (t115 * t220 + t245) * t117 + (-t127 * t232 - (-t126 * t241 + t135 * t211 + t238) * t246) * t225) * t123 + (t200 * t232 + ((-t127 * t227 - (t203 * qJD(2) - t119) * t218) * t174 + (t127 * t202 + (t178 * t115 - (-t114 + t229) * t242 - (t203 * t119 - qJD(2)) * t239) * t128) * t172) * t123) * t171, 0, 0 (t116 * t246 - t127 * t153) * t223 + (t116 * t245 + t132 * t127 + (t116 * t220 - t153 * t128) * t115 - (t173 * t228 - t172 * t226 - t112 * t149 - t122 * t133 + (-t122 * t233 - t151) * t119) * t128 * t239 - (-t134 + (-t112 * t171 - t119 * t173) * t172 - t186 * t122) * t218) * t123, 0; t190 * t172 * t221 + (-t190 * t228 + ((qJD(1) * t146 + 0.2e1 * t151 * t244) * t178 + (-t249 * t132 - t134 * t178 + t151 * t202) * t147) * t172) * t141 (t146 * t174 * t178 + t173 * t215) * t221 + (0.2e1 * t173 * t196 + t188 * t146 + ((t132 * t178 - 0.2e1 * t173 * t207) * t174 + (t177 * t226 + 0.2e1 * t178 * t195) * t167) * t147) * t141, 0, 0, t147 * t198 * t252 + (t198 * t244 + (t131 * t232 + (t172 * t202 - t174 * t227) * t152) * t147) * t141, 0;];
JaD_rot  = t1;
