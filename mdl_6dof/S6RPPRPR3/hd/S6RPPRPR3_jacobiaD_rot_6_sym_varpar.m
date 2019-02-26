% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRPR3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:26:47
% EndTime: 2019-02-26 20:26:48
% DurationCPUTime: 0.71s
% Computational Cost: add. (3055->93), mult. (2519->204), div. (480->12), fcn. (2968->9), ass. (0->96)
t145 = qJ(4) + pkin(10);
t141 = sin(t145);
t135 = 0.1e1 / t141 ^ 2;
t143 = cos(t145);
t139 = t143 ^ 2;
t192 = t135 * t139;
t146 = qJ(1) + pkin(9);
t144 = cos(t146);
t210 = 0.2e1 * t144;
t209 = t143 * t192;
t140 = t144 ^ 2;
t132 = t140 * t192 + 0.1e1;
t130 = 0.1e1 / t132;
t134 = 0.1e1 / t141;
t142 = sin(t146);
t183 = qJD(1) * t143;
t172 = t142 * t183;
t179 = qJD(4) * t144;
t104 = ((t141 * t179 + t172) * t134 + t179 * t192) * t130;
t208 = -t104 + t179;
t187 = t144 * t143;
t129 = atan2(-t187, t141);
t127 = sin(t129);
t128 = cos(t129);
t114 = -t127 * t187 + t128 * t141;
t111 = 0.1e1 / t114;
t147 = sin(qJ(6));
t186 = t144 * t147;
t148 = cos(qJ(6));
t188 = t142 * t148;
t124 = t141 * t188 + t186;
t120 = 0.1e1 / t124;
t112 = 0.1e1 / t114 ^ 2;
t121 = 0.1e1 / t124 ^ 2;
t207 = t130 - 0.1e1;
t137 = t142 ^ 2;
t191 = t137 * t139;
t107 = t112 * t191 + 0.1e1;
t182 = qJD(1) * t144;
t163 = t139 * t142 * t182;
t180 = qJD(4) * t143;
t194 = t128 * t143;
t99 = (-t104 * t144 + qJD(4)) * t194 + (t208 * t141 + t172) * t127;
t205 = t111 * t112 * t99;
t206 = 0.1e1 / t107 ^ 2 * (-t191 * t205 + (-t137 * t141 * t180 + t163) * t112);
t167 = qJD(6) * t141 + qJD(1);
t157 = t147 * t180 + t167 * t148;
t166 = qJD(1) * t141 + qJD(6);
t160 = t144 * t166;
t108 = t157 * t142 + t147 * t160;
t185 = t144 * t148;
t189 = t142 * t147;
t123 = t141 * t189 - t185;
t119 = t123 ^ 2;
t118 = t119 * t121 + 0.1e1;
t197 = t121 * t123;
t156 = -t167 * t147 + t148 * t180;
t109 = t156 * t142 + t148 * t160;
t202 = t109 * t120 * t121;
t204 = (t108 * t197 - t119 * t202) / t118 ^ 2;
t203 = t104 * t143;
t158 = qJD(4) * (-t143 - t209) * t134;
t201 = (-t135 * t163 + t140 * t158) / t132 ^ 2;
t200 = t112 * t142;
t199 = t112 * t143;
t198 = t120 * t147;
t196 = t123 * t148;
t195 = t127 * t144;
t193 = t134 * t139;
t190 = t142 * t143;
t184 = qJD(1) * t142;
t181 = qJD(4) * t141;
t178 = -0.2e1 * t205;
t177 = 0.2e1 * t204;
t176 = t143 * t206;
t175 = t143 * t201;
t174 = t112 * t190;
t173 = t130 * t193;
t171 = t143 * t182;
t170 = 0.1e1 + t192;
t169 = 0.2e1 * t123 * t202;
t168 = t201 * t210;
t165 = t144 * t173;
t164 = t207 * t143 * t127;
t162 = t170 * t142;
t161 = t142 * t166;
t159 = t121 * t196 - t198;
t126 = t141 * t185 - t189;
t125 = t141 * t186 + t188;
t116 = 0.1e1 / t118;
t115 = t170 * t144 * t130;
t105 = 0.1e1 / t107;
t103 = (-t128 * t165 - t164) * t142;
t101 = t141 * t195 + t194 + (-t127 * t141 - t128 * t187) * t115;
t100 = -t170 * t168 + (-qJD(1) * t162 + t158 * t210) * t130;
t1 = [-0.2e1 * t134 * t142 * t175 + (-qJD(4) * t162 + t134 * t171) * t130, 0, 0, t100, 0, 0; (0.2e1 * t111 * t176 + (t111 * t181 + (qJD(1) * t103 + t99) * t199) * t105) * t144 + (-0.2e1 * t112 * t176 * t103 + (((t104 * t165 + t207 * t181 + 0.2e1 * t175) * t127 + (t168 * t193 + t203 + (-t203 + (0.2e1 * t143 + t209) * t179) * t130) * t128) * t174 + (-t112 * t181 + t143 * t178) * t103 + (t111 + ((t137 - t140) * t128 * t173 - t144 * t164) * t112) * t183) * t105) * t142, 0, 0, 0.2e1 * (-t101 * t199 - t111 * t141) * t142 * t206 + ((t111 * t182 + (-qJD(4) * t101 - t99) * t200) * t141 + (t142 * qJD(4) * t111 + (-t100 * t128 * t144 + t208 * t127 + (-qJD(4) * t127 + t104 * t195 + t128 * t184) * t115) * t174 + (t112 * t182 + t142 * t178) * t101 + ((-t100 - t184) * t127 + ((t115 * t144 - 0.1e1) * qJD(4) + (-t115 + t144) * t104) * t128) * t141 * t200) * t143) * t105, 0, 0; (-t120 * t125 + t126 * t197) * t177 + (t126 * t169 - t161 * t198 + t157 * t120 * t144 + (-t156 * t123 * t144 - t126 * t108 - t125 * t109 + t161 * t196) * t121) * t116, 0, 0, t159 * t177 * t190 + (-t159 * t171 + (t159 * t181 + ((qJD(6) * t120 + t169) * t148 + (-t108 * t148 + (qJD(6) * t123 - t109) * t147) * t121) * t143) * t142) * t116, 0, -0.2e1 * t204 + 0.2e1 * (t108 * t116 * t121 + (-t116 * t202 - t121 * t204) * t123) * t123;];
JaD_rot  = t1;
