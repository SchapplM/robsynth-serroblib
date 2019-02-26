% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPPR2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:39:46
% EndTime: 2019-02-26 20:39:47
% DurationCPUTime: 0.71s
% Computational Cost: add. (3055->93), mult. (2519->205), div. (480->12), fcn. (2968->9), ass. (0->93)
t145 = qJ(3) + pkin(10);
t141 = sin(t145);
t135 = 0.1e1 / t141 ^ 2;
t143 = cos(t145);
t139 = t143 ^ 2;
t193 = t135 * t139;
t209 = t143 * t193;
t146 = qJ(1) + pkin(9);
t142 = sin(t146);
t168 = 0.1e1 + t193;
t208 = t142 * t168;
t137 = t142 ^ 2;
t132 = t137 * t193 + 0.1e1;
t130 = 0.1e1 / t132;
t134 = 0.1e1 / t141;
t144 = cos(t146);
t183 = qJD(1) * t144;
t171 = t143 * t183;
t181 = qJD(3) * t142;
t104 = ((t141 * t181 - t171) * t134 + t181 * t193) * t130;
t207 = -t104 + t181;
t147 = sin(qJ(6));
t148 = cos(qJ(6));
t165 = qJD(6) * t141 + qJD(1);
t180 = qJD(3) * t143;
t206 = t165 * t147 - t148 * t180;
t205 = t147 * t180 + t165 * t148;
t190 = t142 * t143;
t129 = atan2(-t190, t141);
t128 = cos(t129);
t127 = sin(t129);
t173 = t127 * t190;
t115 = t128 * t141 - t173;
t111 = 0.1e1 / t115;
t187 = t144 * t147;
t188 = t142 * t148;
t124 = t141 * t187 + t188;
t120 = 0.1e1 / t124;
t112 = 0.1e1 / t115 ^ 2;
t121 = 0.1e1 / t124 ^ 2;
t204 = t130 - 0.1e1;
t194 = t128 * t143;
t99 = (-t104 * t142 + qJD(3)) * t194 + (t141 * t207 - t171) * t127;
t203 = t111 * t112 * t99;
t164 = qJD(1) * t141 + qJD(6);
t159 = t164 * t148;
t108 = t142 * t159 + t144 * t206;
t186 = t144 * t148;
t189 = t142 * t147;
t123 = -t141 * t186 + t189;
t119 = t123 ^ 2;
t118 = t119 * t121 + 0.1e1;
t196 = t121 * t123;
t160 = t164 * t147;
t109 = -t142 * t160 + t144 * t205;
t200 = t109 * t120 * t121;
t202 = (t108 * t196 - t119 * t200) / t118 ^ 2;
t201 = t104 * t143;
t157 = qJD(3) * (-t143 - t209) * t134;
t191 = t139 * t142;
t162 = t183 * t191;
t199 = (t135 * t162 + t137 * t157) / t132 ^ 2;
t198 = t112 * t143;
t197 = t112 * t144;
t195 = t127 * t142;
t140 = t144 ^ 2;
t192 = t139 * t140;
t185 = qJD(1) * t142;
t184 = qJD(1) * t143;
t182 = qJD(3) * t141;
t107 = t112 * t192 + 0.1e1;
t179 = 0.2e1 / t107 ^ 2 * (-t192 * t203 + (-t140 * t141 * t180 - t162) * t112);
t178 = 0.2e1 * t203;
t177 = 0.2e1 * t202;
t176 = -0.2e1 * t199;
t175 = t143 * t199;
t174 = t143 * t197;
t172 = t134 * t191;
t167 = t143 * t179;
t166 = 0.2e1 * t123 * t200;
t163 = t130 * t172;
t161 = t168 * t144;
t158 = t120 * t148 + t147 * t196;
t156 = t158 * t144;
t126 = -t141 * t189 + t186;
t125 = t141 * t188 + t187;
t116 = 0.1e1 / t118;
t114 = t130 * t208;
t105 = 0.1e1 / t107;
t103 = (t204 * t143 * t127 + t128 * t163) * t144;
t101 = t141 * t195 + t194 + (-t127 * t141 - t128 * t190) * t114;
t100 = t176 * t208 + (qJD(1) * t161 + 0.2e1 * t142 * t157) * t130;
t1 = [0.2e1 * t134 * t144 * t175 + (t134 * t142 * t184 + qJD(3) * t161) * t130, 0, t100, 0, 0, 0; (t111 * t167 + (t111 * t182 + (qJD(1) * t103 + t99) * t198) * t105) * t142 + (t112 * t167 * t103 + (-((-t104 * t163 - t204 * t182 - 0.2e1 * t175) * t127 + (t172 * t176 - t201 + (t201 + (-0.2e1 * t143 - t209) * t181) * t130) * t128) * t174 + (t112 * t182 + t143 * t178) * t103 + (-t111 + ((t137 - t140) * t139 * t134 * t130 * t128 + t204 * t173) * t112) * t184) * t105) * t144, 0 (t101 * t198 + t111 * t141) * t144 * t179 + ((t111 * t185 + (qJD(3) * t101 + t99) * t197) * t141 + (-t144 * qJD(3) * t111 - (-t100 * t128 * t142 + t207 * t127 + (-qJD(3) * t127 + t104 * t195 - t128 * t183) * t114) * t174 + (t112 * t185 + t144 * t178) * t101 - ((-t100 + t183) * t127 + ((t114 * t142 - 0.1e1) * qJD(3) + (-t114 + t142) * t104) * t128) * t141 * t197) * t143) * t105, 0, 0, 0; (-t120 * t125 + t126 * t196) * t177 + (t126 * t166 + (-t126 * t108 - t125 * t109 + (t142 * t205 + t144 * t160) * t123) * t121 + (-t206 * t142 + t144 * t159) * t120) * t116, 0, t143 * t156 * t177 + (t156 * t182 + (t158 * t185 + ((qJD(6) * t120 + t166) * t147 + (-t108 * t147 + (-qJD(6) * t123 + t109) * t148) * t121) * t144) * t143) * t116, 0, 0, -0.2e1 * t202 + 0.2e1 * (t108 * t116 * t121 + (-t116 * t200 - t121 * t202) * t123) * t123;];
JaD_rot  = t1;
