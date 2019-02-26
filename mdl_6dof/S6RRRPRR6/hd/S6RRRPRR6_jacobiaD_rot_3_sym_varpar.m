% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR6_jacobiaD_rot_3_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_jacobiaD_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR6_jacobiaD_rot_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_jacobiaD_rot_3_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:18:48
% EndTime: 2019-02-26 22:18:49
% DurationCPUTime: 0.68s
% Computational Cost: add. (1002->94), mult. (2519->211), div. (480->12), fcn. (2968->9), ass. (0->92)
t126 = sin(qJ(1));
t119 = t126 ^ 2;
t125 = sin(qJ(2));
t118 = t125 ^ 2;
t128 = cos(qJ(2));
t121 = 0.1e1 / t128 ^ 2;
t173 = t118 * t121;
t113 = t119 * t173 + 0.1e1;
t117 = t125 * t118;
t120 = 0.1e1 / t128;
t172 = t120 * t125;
t137 = qJD(2) * (t117 * t120 * t121 + t172);
t129 = cos(qJ(1));
t163 = qJD(1) * t129;
t151 = t126 * t163;
t181 = 0.1e1 / t113 ^ 2 * (t119 * t137 + t151 * t173);
t193 = -0.2e1 * t181;
t111 = 0.1e1 / t113;
t146 = 0.1e1 + t173;
t190 = t126 * t146;
t98 = t111 * t190;
t192 = t126 * t98 - 0.1e1;
t124 = sin(qJ(3));
t127 = cos(qJ(3));
t165 = t129 * t127;
t107 = t126 * t124 + t128 * t165;
t102 = 0.1e1 / t107 ^ 2;
t166 = t129 * t124;
t168 = t126 * t127;
t106 = t128 * t166 - t168;
t175 = t106 * t127;
t101 = 0.1e1 / t107;
t177 = t101 * t124;
t139 = t102 * t175 - t177;
t100 = t106 ^ 2;
t99 = t100 * t102 + 0.1e1;
t96 = 0.1e1 / t99;
t191 = t139 * t96;
t169 = t126 * t125;
t110 = atan2(-t169, -t128);
t109 = cos(t110);
t108 = sin(t110);
t154 = t108 * t169;
t94 = -t109 * t128 - t154;
t91 = 0.1e1 / t94;
t92 = 0.1e1 / t94 ^ 2;
t189 = t111 - 0.1e1;
t123 = t129 ^ 2;
t161 = qJD(2) * t128;
t155 = t92 * t161;
t150 = t125 * t163;
t162 = qJD(2) * t126;
t174 = t109 * t125;
t149 = t121 * t162;
t85 = (-(-t126 * t161 - t150) * t120 + t118 * t149) * t111;
t80 = (-t126 * t85 + qJD(2)) * t174 + (-t150 + (t85 - t162) * t128) * t108;
t187 = t80 * t91 * t92;
t90 = t123 * t118 * t92 + 0.1e1;
t188 = (t123 * t125 * t155 + (-t123 * t187 - t92 * t151) * t118) / t90 ^ 2;
t176 = t102 * t106;
t143 = -qJD(1) * t128 + qJD(3);
t144 = qJD(3) * t128 - qJD(1);
t160 = qJD(2) * t129;
t148 = t125 * t160;
t87 = -t144 * t166 + (t143 * t126 - t148) * t127;
t183 = t101 * t102 * t87;
t167 = t126 * t128;
t138 = t124 * t167 + t165;
t86 = t138 * qJD(1) - t107 * qJD(3) + t124 * t148;
t186 = (-t100 * t183 - t86 * t176) / t99 ^ 2;
t88 = 0.1e1 / t90;
t185 = t88 * t92;
t184 = t91 * t88;
t179 = t129 * t92;
t178 = qJD(2) * t98;
t171 = t125 * t129;
t164 = qJD(1) * t126;
t159 = 0.2e1 * t187;
t158 = -0.2e1 * t186;
t157 = t91 * t188;
t156 = t106 * t183;
t153 = t111 * t118 * t120;
t147 = 0.2e1 * t92 * t188;
t145 = t120 * t193;
t142 = t126 * t153;
t141 = t146 * t129;
t140 = t143 * t129;
t105 = -t127 * t167 + t166;
t84 = (t189 * t125 * t108 - t109 * t142) * t129;
t83 = -t192 * t174 + (-t126 + t98) * t128 * t108;
t81 = t190 * t193 + (qJD(1) * t141 + 0.2e1 * t126 * t137) * t111;
t1 = [t145 * t171 + (qJD(2) * t141 - t164 * t172) * t111, t81, 0, 0, 0, 0; (-t161 * t184 + (0.2e1 * t157 + (qJD(1) * t84 + t80) * t185) * t125) * t126 + (t84 * t147 * t125 + (-t84 * t155 + (t84 * t159 + ((0.2e1 * t125 * t181 - t85 * t142 - t189 * t161) * t108 + (t118 * t126 * t145 + t125 * t85 + (t117 * t149 - (t85 - 0.2e1 * t162) * t125) * t111) * t109) * t179) * t125 + (-t91 + (-(t119 - t123) * t109 * t153 + t189 * t154) * t92) * t125 * qJD(1)) * t88) * t129 (-t164 * t184 + (-0.2e1 * t157 + (-qJD(2) * t83 - t80) * t185) * t129) * t128 + (t83 * t129 * t147 + (-t91 * t160 - ((-t126 * t81 - t163 * t98) * t109 + (t192 * t85 + t162 - t178) * t108) * t92 * t171 + (t129 * t159 + t92 * t164) * t83 - ((t81 - t163) * t108 + (t85 * t98 + qJD(2) + (-t85 - t178) * t126) * t109) * t128 * t179) * t88) * t125, 0, 0, 0, 0; 0.2e1 * (t101 * t138 + t105 * t176) * t186 + (0.2e1 * t105 * t156 - t144 * t101 * t168 + (t125 * t162 + t140) * t177 + (t138 * t87 + t105 * t86 - t140 * t175 - (qJD(2) * t125 * t127 + t144 * t124) * t106 * t126) * t102) * t96, t128 * t160 * t191 + (-t164 * t191 + (t139 * t158 + ((-qJD(3) * t101 - 0.2e1 * t156) * t127 + (-t127 * t86 + (-qJD(3) * t106 + t87) * t124) * t102) * t96) * t129) * t125, t158 + 0.2e1 * (-t86 * t102 * t96 + (-t102 * t186 - t96 * t183) * t106) * t106, 0, 0, 0;];
JaD_rot  = t1;
