% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRP10
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRP10_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP10_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_jacobiaD_rot_5_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:48:37
% EndTime: 2019-02-26 20:48:38
% DurationCPUTime: 0.66s
% Computational Cost: add. (624->92), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->93)
t129 = cos(qJ(1));
t123 = t129 ^ 2;
t125 = sin(qJ(3));
t118 = t125 ^ 2;
t128 = cos(qJ(3));
t121 = 0.1e1 / t128 ^ 2;
t171 = t118 * t121;
t114 = t123 * t171 + 0.1e1;
t117 = t125 * t118;
t120 = 0.1e1 / t128;
t170 = t120 * t125;
t138 = qJD(3) * (t117 * t120 * t121 + t170);
t126 = sin(qJ(1));
t163 = qJD(1) * t129;
t153 = t126 * t163;
t180 = 0.1e1 / t114 ^ 2 * (t123 * t138 - t153 * t171);
t192 = -0.2e1 * t180;
t111 = 0.1e1 / t114;
t148 = 0.1e1 + t171;
t190 = t129 * t148;
t99 = t111 * t190;
t191 = -t129 * t99 + 0.1e1;
t144 = qJD(1) * t128 + qJD(5);
t160 = qJD(3) * t129;
t189 = t125 * t160 + t144 * t126;
t166 = t129 * t125;
t113 = atan2(t166, t128);
t109 = sin(t113);
t110 = cos(t113);
t95 = t109 * t166 + t110 * t128;
t92 = 0.1e1 / t95;
t124 = sin(qJ(5));
t127 = cos(qJ(5));
t165 = t129 * t127;
t168 = t126 * t128;
t108 = -t124 * t168 + t165;
t102 = 0.1e1 / t108;
t103 = 0.1e1 / t108 ^ 2;
t93 = 0.1e1 / t95 ^ 2;
t188 = t111 - 0.1e1;
t119 = t126 ^ 2;
t161 = qJD(3) * t128;
t157 = t93 * t161;
t164 = qJD(1) * t126;
t154 = t125 * t164;
t173 = t110 * t125;
t152 = t121 * t160;
t86 = ((t128 * t160 - t154) * t120 + t118 * t152) * t111;
t81 = (t129 * t86 - qJD(3)) * t173 + (-t154 + (-t86 + t160) * t128) * t109;
t186 = t81 * t92 * t93;
t91 = t119 * t118 * t93 + 0.1e1;
t187 = (t119 * t125 * t157 + (-t119 * t186 + t153 * t93) * t118) / t91 ^ 2;
t167 = t129 * t124;
t107 = t127 * t168 + t167;
t101 = t107 ^ 2;
t100 = t101 * t103 + 0.1e1;
t175 = t103 * t107;
t145 = qJD(5) * t128 + qJD(1);
t162 = qJD(3) * t126;
t169 = t126 * t127;
t88 = -t145 * t169 + (t125 * t162 - t129 * t144) * t124;
t182 = t102 * t103 * t88;
t155 = t128 * t165;
t87 = -qJD(1) * t155 - qJD(5) * t165 + (qJD(3) * t125 * t127 + t124 * t145) * t126;
t185 = (-t101 * t182 - t175 * t87) / t100 ^ 2;
t89 = 0.1e1 / t91;
t184 = t89 * t93;
t183 = t92 * t89;
t179 = t126 * t93;
t177 = qJD(3) * t99;
t176 = t102 * t127;
t174 = t107 * t124;
t172 = t118 * t120;
t159 = 0.2e1 * t186;
t158 = 0.2e1 * t185;
t156 = t129 * t172;
t150 = -0.2e1 * t92 * t187;
t149 = 0.2e1 * t93 * t187;
t147 = 0.2e1 * t107 * t182;
t146 = 0.2e1 * t125 * t180;
t143 = t111 * t156;
t142 = t188 * t125 * t109;
t141 = t148 * t126;
t140 = t145 * t129;
t139 = t103 * t174 + t176;
t97 = 0.1e1 / t100;
t137 = t139 * t97;
t106 = -t128 * t167 - t169;
t105 = -t126 * t124 + t155;
t85 = (-t110 * t143 + t142) * t126;
t84 = -t191 * t173 + (t129 - t99) * t128 * t109;
t82 = t190 * t192 + (-qJD(1) * t141 + 0.2e1 * t129 * t138) * t111;
t1 = [t126 * t120 * t146 + (-qJD(3) * t141 - t163 * t170) * t111, 0, t82, 0, 0, 0; (t161 * t183 + (t150 + (-qJD(1) * t85 - t81) * t184) * t125) * t129 + (t85 * t149 * t125 + (-t85 * t157 + (t85 * t159 + ((-t143 * t86 - t161 * t188 + t146) * t109 + (t156 * t192 + t125 * t86 + (t117 * t152 - (t86 - 0.2e1 * t160) * t125) * t111) * t110) * t179) * t125 + (-t92 + (-(t119 - t123) * t111 * t110 * t172 - t129 * t142) * t93) * t125 * qJD(1)) * t89) * t126, 0 (t163 * t183 + (t150 + (-qJD(3) * t84 - t81) * t184) * t126) * t128 + (t84 * t126 * t149 + (-t92 * t162 - ((t129 * t82 - t164 * t99) * t110 + (t191 * t86 - t160 + t177) * t109) * t125 * t179 + (t126 * t159 - t163 * t93) * t84) * t89 - ((-t82 - t164) * t109 + (-t86 * t99 - qJD(3) + (t86 + t177) * t129) * t110) * t168 * t184) * t125, 0, 0, 0; (-t102 * t105 + t106 * t175) * t158 + (t106 * t147 - t102 * t124 * t140 - t189 * t176 + (t107 * t127 * t140 - t105 * t88 + t106 * t87 - t189 * t174) * t103) * t97, 0, -t126 * t137 * t161 + (-t137 * t163 + (t139 * t158 + ((qJD(5) * t102 + t147) * t124 + (t124 * t87 + (-qJD(5) * t107 + t88) * t127) * t103) * t97) * t126) * t125, 0, -0.2e1 * t185 + 0.2e1 * (-t87 * t103 * t97 + (-t103 * t185 - t182 * t97) * t107) * t107, 0;];
JaD_rot  = t1;
