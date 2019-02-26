% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRRP3_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_jacobiaD_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:31:27
% EndTime: 2019-02-26 20:31:28
% DurationCPUTime: 0.67s
% Computational Cost: add. (1787->89), mult. (2519->202), div. (480->12), fcn. (2968->9), ass. (0->95)
t121 = qJ(1) + pkin(9);
t120 = cos(t121);
t193 = 0.2e1 * t120;
t119 = sin(t121);
t127 = sin(qJ(5));
t128 = sin(qJ(4));
t129 = cos(qJ(5));
t171 = t128 * t129;
t107 = t119 * t171 + t120 * t127;
t104 = 0.1e1 / t107 ^ 2;
t172 = t127 * t128;
t106 = t119 * t172 - t120 * t129;
t178 = t106 * t129;
t103 = 0.1e1 / t107;
t180 = t103 * t127;
t141 = t104 * t178 - t180;
t102 = t106 ^ 2;
t101 = t102 * t104 + 0.1e1;
t97 = 0.1e1 / t101;
t192 = t141 * t97;
t167 = qJD(4) * t120;
t118 = t120 ^ 2;
t123 = 0.1e1 / t128 ^ 2;
t130 = cos(qJ(4));
t126 = t130 ^ 2;
t173 = t123 * t126;
t115 = t118 * t173 + 0.1e1;
t113 = 0.1e1 / t115;
t122 = 0.1e1 / t128;
t156 = t123 * t167;
t168 = qJD(1) * t130;
t158 = t119 * t168;
t166 = qJD(4) * t128;
t87 = ((t120 * t166 + t158) * t122 + t126 * t156) * t113;
t150 = -t87 + t167;
t151 = -t120 * t87 + qJD(4);
t175 = t120 * t130;
t112 = atan2(-t175, t128);
t110 = sin(t112);
t111 = cos(t112);
t99 = -t110 * t175 + t111 * t128;
t94 = 0.1e1 / t99;
t95 = 0.1e1 / t99 ^ 2;
t191 = t113 - 0.1e1;
t117 = t119 ^ 2;
t169 = qJD(1) * t120;
t145 = t119 * t126 * t169;
t160 = t95 * t166;
t176 = t111 * t130;
t82 = t151 * t176 + (t128 * t150 + t158) * t110;
t189 = t82 * t94 * t95;
t92 = t117 * t126 * t95 + 0.1e1;
t190 = (t95 * t145 + (-t126 * t189 - t130 * t160) * t117) / t92 ^ 2;
t179 = t104 * t106;
t149 = qJD(5) * t128 + qJD(1);
t165 = qJD(4) * t130;
t138 = -t127 * t149 + t129 * t165;
t148 = qJD(1) * t128 + qJD(5);
t142 = t120 * t148;
t89 = t119 * t138 + t129 * t142;
t185 = t103 * t104 * t89;
t139 = t127 * t165 + t129 * t149;
t88 = t119 * t139 + t127 * t142;
t188 = (-t102 * t185 + t88 * t179) / t101 ^ 2;
t90 = 0.1e1 / t92;
t187 = t90 * t95;
t186 = t94 * t90;
t125 = t130 * t126;
t140 = qJD(4) * (-t123 * t125 - t130) * t122;
t184 = 0.1e1 / t115 ^ 2 * (t118 * t140 - t123 * t145);
t183 = t119 * t95;
t177 = t110 * t128;
t174 = t122 * t126;
t170 = qJD(1) * t119;
t164 = -0.2e1 * t189;
t163 = 0.2e1 * t188;
t162 = t94 * t190;
t161 = t130 * t184;
t159 = t113 * t174;
t157 = t120 * t168;
t155 = -0.2e1 * t95 * t190;
t154 = 0.1e1 + t173;
t153 = 0.2e1 * t106 * t185;
t152 = t184 * t193;
t147 = t120 * t159;
t146 = t191 * t130 * t110;
t144 = t154 * t119;
t143 = t119 * t148;
t109 = -t119 * t127 + t120 * t171;
t108 = t119 * t129 + t120 * t172;
t100 = t154 * t120 * t113;
t86 = (-t111 * t147 - t146) * t119;
t85 = t120 * t177 + t176 + (-t111 * t175 - t177) * t100;
t83 = -t154 * t152 + (-qJD(1) * t144 + t140 * t193) * t113;
t1 = [-0.2e1 * t119 * t122 * t161 + (-qJD(4) * t144 + t122 * t157) * t113, 0, 0, t83, 0, 0; (t166 * t186 + (0.2e1 * t162 + (qJD(1) * t86 + t82) * t187) * t130) * t120 + (t86 * t155 * t130 + (-t86 * t160 + (t86 * t164 + ((t147 * t87 + t166 * t191 + 0.2e1 * t161) * t110 + (t152 * t174 + t130 * t87 + (t125 * t156 + (-t87 + 0.2e1 * t167) * t130) * t113) * t111) * t183) * t130 + (t94 + ((t117 - t118) * t111 * t159 - t120 * t146) * t95) * t168) * t90) * t119, 0, 0 (t169 * t186 + (-0.2e1 * t162 + (-qJD(4) * t85 - t82) * t187) * t119) * t128 + (t85 * t119 * t155 + (t119 * qJD(4) * t94 + (t119 * t164 + t169 * t95) * t85 + (((t100 * t170 - t120 * t83) * t111 + (-t151 * t100 + t150) * t110) * t130 + ((-t83 - t170) * t110 + (t100 * t150 - t151) * t111) * t128) * t183) * t90) * t130, 0, 0; (-t103 * t108 + t109 * t179) * t163 + (t109 * t153 - t143 * t180 + t139 * t103 * t120 + (-t106 * t120 * t138 - t108 * t89 - t109 * t88 + t143 * t178) * t104) * t97, 0, 0, -t157 * t192 + (t166 * t192 + (t141 * t163 + ((qJD(5) * t103 + t153) * t129 + (-t129 * t88 + (qJD(5) * t106 - t89) * t127) * t104) * t97) * t130) * t119, -0.2e1 * t188 + 0.2e1 * (t104 * t88 * t97 + (-t104 * t188 - t185 * t97) * t106) * t106, 0;];
JaD_rot  = t1;
