% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRP4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRP4_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_jacobiaD_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:26:47
% EndTime: 2019-02-26 21:26:49
% DurationCPUTime: 0.87s
% Computational Cost: add. (926->105), mult. (3345->234), div. (468->12), fcn. (4023->11), ass. (0->103)
t155 = sin(qJ(2));
t146 = t155 ^ 2;
t158 = cos(qJ(2));
t149 = 0.1e1 / t158 ^ 2;
t204 = t146 * t149;
t156 = sin(qJ(1));
t147 = t156 ^ 2;
t144 = t147 * t204 + 0.1e1;
t148 = 0.1e1 / t158;
t201 = t148 * t155;
t222 = t155 * t204;
t167 = qJD(2) * (t148 * t222 + t201);
t159 = cos(qJ(1));
t195 = qJD(1) * t159;
t202 = t146 * t156;
t173 = t195 * t202;
t207 = (t147 * t167 + t149 * t173) / t144 ^ 2;
t223 = -0.2e1 * t207;
t177 = 0.1e1 + t204;
t221 = t156 * t177;
t154 = sin(qJ(5));
t157 = cos(qJ(5));
t191 = qJD(5) * t159;
t196 = qJD(1) * t156;
t220 = t154 * t196 - t157 * t191;
t219 = t154 * t191 + t157 * t196;
t152 = sin(pkin(9));
t153 = cos(pkin(9));
t197 = t158 * t159;
t137 = t152 * t197 - t156 * t153;
t138 = t156 * t152 + t153 * t197;
t121 = t137 * t154 + t138 * t157;
t115 = 0.1e1 / t121;
t199 = t156 * t155;
t143 = atan2(t199, t158);
t140 = cos(t143);
t139 = sin(t143);
t186 = t139 * t199;
t125 = t140 * t158 + t186;
t122 = 0.1e1 / t125;
t116 = 0.1e1 / t121 ^ 2;
t123 = 0.1e1 / t125 ^ 2;
t218 = 0.2e1 * t155;
t141 = 0.1e1 / t144;
t217 = t141 - 0.1e1;
t198 = t156 * t158;
t135 = -t152 * t198 - t153 * t159;
t192 = qJD(2) * t159;
t180 = t155 * t192;
t128 = t135 * qJD(1) - t152 * t180;
t136 = t152 * t159 - t153 * t198;
t129 = t136 * qJD(1) - t153 * t180;
t104 = t121 * qJD(5) - t128 * t157 + t129 * t154;
t170 = t137 * t157 - t138 * t154;
t114 = t170 ^ 2;
t108 = t114 * t116 + 0.1e1;
t210 = t116 * t170;
t105 = t170 * qJD(5) + t128 * t154 + t129 * t157;
t117 = t115 * t116;
t213 = t105 * t117;
t216 = 0.1e1 / t108 ^ 2 * (-t104 * t210 - t114 * t213);
t151 = t159 ^ 2;
t203 = t146 * t151;
t113 = t123 * t203 + 0.1e1;
t193 = qJD(2) * t158;
t182 = t155 * t195;
t194 = qJD(2) * t156;
t110 = ((t156 * t193 + t182) * t148 + t194 * t204) * t141;
t205 = t140 * t155;
t101 = (t110 * t156 - qJD(2)) * t205 + (t182 + (-t110 + t194) * t158) * t139;
t214 = t101 * t122 * t123;
t215 = (-t203 * t214 + (t151 * t155 * t193 - t173) * t123) / t113 ^ 2;
t212 = t110 * t139;
t211 = t110 * t155;
t168 = -t152 * t154 - t153 * t157;
t200 = t155 * t159;
t133 = t168 * t200;
t209 = t116 * t133;
t208 = t123 * t155;
t127 = t141 * t221;
t206 = t127 * t156;
t190 = 0.2e1 * t216;
t189 = -0.2e1 * t214;
t188 = -0.2e1 * t117 * t170;
t187 = t123 * t200;
t185 = t141 * t146 * t148;
t181 = t155 * t194;
t176 = -0.2e1 * t155 * t215;
t175 = t148 * t223;
t174 = t156 * t185;
t172 = t177 * t159;
t171 = t135 * t157 - t136 * t154;
t119 = t135 * t154 + t136 * t157;
t169 = t152 * t157 - t153 * t154;
t132 = t169 * t200;
t131 = -t138 * qJD(1) + t153 * t181;
t130 = -t137 * qJD(1) + t152 * t181;
t111 = 0.1e1 / t113;
t109 = (-t217 * t155 * t139 + t140 * t174) * t159;
t106 = 0.1e1 / t108;
t103 = t139 * t198 - t205 + (-t139 * t158 + t140 * t199) * t127;
t102 = t221 * t223 + (qJD(1) * t172 + 0.2e1 * t156 * t167) * t141;
t1 = [t175 * t200 + (qJD(2) * t172 - t196 * t201) * t141, t102, 0, 0, 0, 0; (t122 * t176 + (t122 * t193 + (-qJD(1) * t109 - t101) * t208) * t111) * t156 + (t123 * t176 * t109 + (((-t110 * t174 - t217 * t193 + t207 * t218) * t139 + (t175 * t202 + t211 + (-t211 + (t218 + t222) * t194) * t141) * t140) * t187 + (t123 * t193 + t155 * t189) * t109 + (t122 + ((-t147 + t151) * t140 * t185 + t217 * t186) * t123) * t155 * qJD(1)) * t111) * t159, 0.2e1 * (-t103 * t208 + t122 * t158) * t159 * t215 + ((t122 * t196 + (qJD(2) * t103 + t101) * t159 * t123) * t158 + (t122 * t192 + (t102 * t140 * t156 - t139 * t194 - t206 * t212 + t212 + (qJD(2) * t139 + t140 * t195) * t127) * t187 + (-t123 * t196 + t159 * t189) * t103 + ((-t102 + t195) * t139 + ((-0.1e1 + t206) * qJD(2) + (-t127 + t156) * t110) * t140) * t123 * t197) * t155) * t111, 0, 0, 0, 0; (t115 * t171 - t119 * t210) * t190 + ((t119 * qJD(5) - t130 * t157 + t131 * t154) * t115 + t119 * t105 * t188 + (t171 * t105 + (t171 * qJD(5) + t130 * t154 + t131 * t157) * t170 - t119 * t104) * t116) * t106 (-t115 * t132 - t170 * t209) * t190 + (-t104 * t209 + (-t116 * t132 + t133 * t188) * t105 + (t169 * t115 + t168 * t210) * t158 * t192 + ((t220 * t115 + t219 * t210) * t153 + (-t219 * t115 + t220 * t210) * t152) * t155) * t106, 0, 0, -0.2e1 * t216 - 0.2e1 * (t104 * t116 * t106 - (-t106 * t213 - t116 * t216) * t170) * t170, 0;];
JaD_rot  = t1;
