% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR11_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR11_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:43:36
% EndTime: 2019-02-26 21:43:37
% DurationCPUTime: 0.71s
% Computational Cost: add. (1820->94), mult. (2734->208), div. (498->12), fcn. (3199->9), ass. (0->96)
t154 = sin(qJ(2));
t148 = 0.1e1 / t154 ^ 2;
t156 = cos(qJ(2));
t152 = t156 ^ 2;
t202 = t148 * t152;
t220 = t156 * t202;
t155 = sin(qJ(1));
t177 = 0.1e1 + t202;
t219 = t155 * t177;
t146 = qJD(4) + qJD(6);
t173 = qJD(1) * t154 + t146;
t157 = cos(qJ(1));
t189 = qJD(2) * t157;
t218 = t173 * t155 - t156 * t189;
t190 = qJD(2) * t156;
t217 = t155 * t190 + t173 * t157;
t196 = t155 * t156;
t139 = atan2(-t196, t154);
t137 = cos(t139);
t136 = sin(t139);
t182 = t136 * t196;
t125 = t137 * t154 - t182;
t122 = 0.1e1 / t125;
t145 = qJ(4) + pkin(10) + qJ(6);
t143 = sin(t145);
t144 = cos(t145);
t197 = t155 * t144;
t198 = t154 * t157;
t133 = t143 * t198 + t197;
t129 = 0.1e1 / t133;
t147 = 0.1e1 / t154;
t123 = 0.1e1 / t125 ^ 2;
t130 = 0.1e1 / t133 ^ 2;
t150 = t155 ^ 2;
t142 = t150 * t202 + 0.1e1;
t140 = 0.1e1 / t142;
t216 = t140 - 0.1e1;
t174 = t146 * t154 + qJD(1);
t168 = t174 * t157;
t113 = t143 * t168 + t218 * t144;
t132 = t143 * t155 - t144 * t198;
t128 = t132 ^ 2;
t118 = t128 * t130 + 0.1e1;
t205 = t130 * t132;
t114 = -t218 * t143 + t144 * t168;
t213 = t114 * t129 * t130;
t215 = (t113 * t205 - t128 * t213) / t118 ^ 2;
t193 = qJD(1) * t157;
t180 = t156 * t193;
t191 = qJD(2) * t155;
t115 = ((t154 * t191 - t180) * t147 + t191 * t202) * t140;
t203 = t137 * t156;
t109 = (-t115 * t155 + qJD(2)) * t203 + (-t180 + (-t115 + t191) * t154) * t136;
t214 = t109 * t122 * t123;
t212 = t115 * t136;
t211 = t115 * t156;
t210 = t123 * t156;
t209 = t123 * t157;
t166 = qJD(2) * (-t156 - t220) * t147;
t200 = t152 * t155;
t171 = t193 * t200;
t208 = (t148 * t171 + t150 * t166) / t142 ^ 2;
t127 = t140 * t219;
t207 = t127 * t155;
t206 = t129 * t144;
t204 = t132 * t143;
t153 = t157 ^ 2;
t201 = t152 * t153;
t199 = t154 * t155;
t195 = qJD(1) * t155;
t194 = qJD(1) * t156;
t192 = qJD(2) * t154;
t121 = t123 * t201 + 0.1e1;
t188 = 0.2e1 * (-t201 * t214 + (-t153 * t154 * t190 - t171) * t123) / t121 ^ 2;
t187 = 0.2e1 * t215;
t186 = 0.2e1 * t214;
t185 = -0.2e1 * t208;
t184 = t156 * t209;
t183 = t156 * t208;
t181 = t147 * t200;
t176 = t156 * t188;
t175 = 0.2e1 * t132 * t213;
t172 = t140 * t181;
t170 = t177 * t157;
t169 = t174 * t155;
t167 = t130 * t204 + t206;
t165 = t167 * t157;
t135 = -t143 * t199 + t144 * t157;
t134 = t143 * t157 + t154 * t197;
t119 = 0.1e1 / t121;
t116 = 0.1e1 / t118;
t112 = (t216 * t156 * t136 + t137 * t172) * t157;
t111 = t136 * t199 + t203 + (-t136 * t154 - t137 * t196) * t127;
t110 = t185 * t219 + (qJD(1) * t170 + 0.2e1 * t155 * t166) * t140;
t106 = -0.2e1 * t215 + 0.2e1 * (t113 * t116 * t130 + (-t116 * t213 - t130 * t215) * t132) * t132;
t1 = [0.2e1 * t147 * t157 * t183 + (t147 * t155 * t194 + qJD(2) * t170) * t140, t110, 0, 0, 0, 0; (t122 * t176 + (t122 * t192 + (qJD(1) * t112 + t109) * t210) * t119) * t155 + (t123 * t176 * t112 + (-((-t115 * t172 - t216 * t192 - 0.2e1 * t183) * t136 + (t181 * t185 - t211 + (t211 + (-0.2e1 * t156 - t220) * t191) * t140) * t137) * t184 + (t123 * t192 + t156 * t186) * t112 + (-t122 + ((t150 - t153) * t152 * t147 * t140 * t137 + t216 * t182) * t123) * t194) * t119) * t157 (t111 * t210 + t122 * t154) * t157 * t188 + ((t122 * t195 + (qJD(2) * t111 + t109) * t209) * t154 + (-t122 * t189 - (-t110 * t137 * t155 + t136 * t191 + t207 * t212 - t212 + (-qJD(2) * t136 - t137 * t193) * t127) * t184 + (t123 * t195 + t157 * t186) * t111 - ((-t110 + t193) * t136 + ((-0.1e1 + t207) * qJD(2) + (-t127 + t155) * t115) * t137) * t123 * t198) * t156) * t119, 0, 0, 0, 0; (-t129 * t134 + t135 * t205) * t187 + (t135 * t175 - t129 * t143 * t169 + t217 * t206 + (t132 * t144 * t169 - t135 * t113 - t134 * t114 + t217 * t204) * t130) * t116, t156 * t165 * t187 + (t165 * t192 + (t167 * t195 + ((t129 * t146 + t175) * t143 + (-t113 * t143 + (-t132 * t146 + t114) * t144) * t130) * t157) * t156) * t116, 0, t106, 0, t106;];
JaD_rot  = t1;
