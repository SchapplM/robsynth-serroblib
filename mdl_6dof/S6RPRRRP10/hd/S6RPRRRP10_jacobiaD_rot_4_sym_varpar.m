% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRRP10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRP10_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP10_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_jacobiaD_rot_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:13:04
% EndTime: 2019-02-26 21:13:05
% DurationCPUTime: 0.65s
% Computational Cost: add. (813->90), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->96)
t132 = cos(qJ(1));
t195 = 0.2e1 * t132;
t167 = qJD(3) * t132;
t126 = t132 ^ 2;
t128 = sin(qJ(3));
t121 = 0.1e1 / t128 ^ 2;
t131 = cos(qJ(3));
t125 = t131 ^ 2;
t179 = t121 * t125;
t118 = t126 * t179 + 0.1e1;
t116 = 0.1e1 / t118;
t120 = 0.1e1 / t128;
t158 = t121 * t167;
t129 = sin(qJ(1));
t171 = qJD(1) * t131;
t160 = t129 * t171;
t90 = ((t128 * t167 + t160) * t120 + t125 * t158) * t116;
t152 = -t90 + t167;
t173 = t132 * t131;
t115 = atan2(-t173, t128);
t113 = sin(t115);
t114 = cos(t115);
t99 = -t113 * t173 + t114 * t128;
t96 = 0.1e1 / t99;
t127 = sin(qJ(4));
t175 = t132 * t127;
t130 = cos(qJ(4));
t177 = t129 * t130;
t110 = t128 * t177 + t175;
t106 = 0.1e1 / t110;
t107 = 0.1e1 / t110 ^ 2;
t97 = 0.1e1 / t99 ^ 2;
t194 = t116 - 0.1e1;
t123 = t129 ^ 2;
t170 = qJD(1) * t132;
t147 = t125 * t129 * t170;
t169 = qJD(3) * t128;
t163 = t97 * t169;
t153 = -t132 * t90 + qJD(3);
t181 = t114 * t131;
t85 = t153 * t181 + (t152 * t128 + t160) * t113;
t192 = t85 * t96 * t97;
t95 = t123 * t125 * t97 + 0.1e1;
t193 = (t97 * t147 + (-t125 * t192 - t131 * t163) * t123) / t95 ^ 2;
t93 = 0.1e1 / t95;
t191 = t93 * t97;
t190 = t96 * t93;
t174 = t132 * t130;
t178 = t129 * t127;
t109 = t128 * t178 - t174;
t105 = t109 ^ 2;
t104 = t105 * t107 + 0.1e1;
t184 = t107 * t109;
t150 = qJD(1) * t128 + qJD(4);
t143 = t150 * t132;
t151 = qJD(4) * t128 + qJD(1);
t145 = t151 * t127;
t168 = qJD(3) * t131;
t92 = t130 * t143 + (t130 * t168 - t145) * t129;
t188 = t106 * t107 * t92;
t144 = t151 * t130;
t91 = t129 * t144 + (t129 * t168 + t143) * t127;
t189 = 0.1e1 / t104 ^ 2 * (-t105 * t188 + t91 * t184);
t187 = t129 * t97;
t124 = t131 * t125;
t141 = qJD(3) * (-t121 * t124 - t131) * t120;
t186 = (-t121 * t147 + t126 * t141) / t118 ^ 2;
t185 = t106 * t127;
t183 = t109 * t130;
t182 = t113 * t132;
t180 = t120 * t125;
t176 = t129 * t131;
t172 = qJD(1) * t129;
t166 = -0.2e1 * t192;
t165 = 0.2e1 * t189;
t164 = t96 * t193;
t162 = t131 * t186;
t161 = t116 * t180;
t159 = t131 * t170;
t157 = -0.2e1 * t97 * t193;
t156 = 0.1e1 + t179;
t155 = 0.2e1 * t109 * t188;
t154 = t186 * t195;
t149 = t132 * t161;
t148 = t194 * t131 * t113;
t146 = t156 * t129;
t142 = t107 * t183 - t185;
t140 = -t150 * t129 + t131 * t167;
t112 = t128 * t174 - t178;
t111 = t128 * t175 + t177;
t103 = t156 * t132 * t116;
t101 = 0.1e1 / t104;
t89 = (-t114 * t149 - t148) * t129;
t88 = t128 * t182 + t181 + (-t113 * t128 - t114 * t173) * t103;
t86 = -t156 * t154 + (-qJD(1) * t146 + t141 * t195) * t116;
t1 = [-0.2e1 * t129 * t120 * t162 + (-qJD(3) * t146 + t120 * t159) * t116, 0, t86, 0, 0, 0; (t169 * t190 + (0.2e1 * t164 + (qJD(1) * t89 + t85) * t191) * t131) * t132 + (t89 * t157 * t131 + (-t89 * t163 + (t89 * t166 + ((t90 * t149 + t194 * t169 + 0.2e1 * t162) * t113 + (t154 * t180 + t131 * t90 + (t124 * t158 + (-t90 + 0.2e1 * t167) * t131) * t116) * t114) * t187) * t131 + (t96 + ((t123 - t126) * t114 * t161 - t132 * t148) * t97) * t171) * t93) * t129, 0 (t170 * t190 + (-0.2e1 * t164 + (-qJD(3) * t88 - t85) * t191) * t129) * t128 + (t88 * t129 * t157 + (t129 * qJD(3) * t96 + (-t114 * t132 * t86 + t152 * t113 + (-qJD(3) * t113 + t114 * t172 + t182 * t90) * t103) * t97 * t176 + (t129 * t166 + t97 * t170) * t88 + ((-t86 - t172) * t113 + (t152 * t103 - t153) * t114) * t128 * t187) * t93) * t131, 0, 0, 0; (-t106 * t111 + t112 * t184) * t165 + (t112 * t155 + t132 * t106 * t144 + t140 * t185 + (t132 * t109 * t145 - t111 * t92 - t112 * t91 - t140 * t183) * t107) * t101, 0, t142 * t165 * t176 + (-t142 * t159 + (t142 * t169 + ((qJD(4) * t106 + t155) * t130 + (-t130 * t91 + (qJD(4) * t109 - t92) * t127) * t107) * t131) * t129) * t101, -0.2e1 * t189 + 0.2e1 * (t91 * t107 * t101 + (-t101 * t188 - t107 * t189) * t109) * t109, 0, 0;];
JaD_rot  = t1;
