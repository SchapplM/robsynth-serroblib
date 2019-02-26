% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRP4_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_jacobiaD_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:45:18
% EndTime: 2019-02-26 20:45:18
% DurationCPUTime: 0.72s
% Computational Cost: add. (1787->89), mult. (2519->200), div. (480->12), fcn. (2968->9), ass. (0->92)
t122 = qJ(1) + pkin(9);
t120 = sin(t122);
t171 = qJD(3) * t120;
t118 = t120 ^ 2;
t129 = sin(qJ(3));
t124 = 0.1e1 / t129 ^ 2;
t131 = cos(qJ(3));
t127 = t131 ^ 2;
t177 = t124 * t127;
t116 = t118 * t177 + 0.1e1;
t114 = 0.1e1 / t116;
t123 = 0.1e1 / t129;
t159 = t124 * t171;
t121 = cos(t122);
t172 = qJD(1) * t131;
t160 = t121 * t172;
t170 = qJD(3) * t129;
t88 = ((t120 * t170 - t160) * t123 + t127 * t159) * t114;
t149 = -t88 + t171;
t150 = -t120 * t88 + qJD(3);
t128 = sin(qJ(5));
t130 = cos(qJ(5));
t148 = qJD(5) * t129 + qJD(1);
t169 = qJD(3) * t131;
t193 = t148 * t128 - t130 * t169;
t192 = t128 * t169 + t148 * t130;
t179 = t120 * t131;
t113 = atan2(-t179, t129);
t112 = cos(t113);
t111 = sin(t113);
t163 = t111 * t179;
t101 = t112 * t129 - t163;
t95 = 0.1e1 / t101;
t176 = t128 * t129;
t108 = t120 * t130 + t121 * t176;
t104 = 0.1e1 / t108;
t96 = 0.1e1 / t101 ^ 2;
t105 = 0.1e1 / t108 ^ 2;
t191 = t114 - 0.1e1;
t180 = t112 * t131;
t83 = t150 * t180 + (t149 * t129 - t160) * t111;
t190 = t83 * t95 * t96;
t175 = t129 * t130;
t107 = t120 * t128 - t121 * t175;
t103 = t107 ^ 2;
t102 = t103 * t105 + 0.1e1;
t182 = t105 * t107;
t147 = qJD(1) * t129 + qJD(5);
t143 = t147 * t128;
t90 = -t120 * t143 + t192 * t121;
t186 = t104 * t105 * t90;
t142 = t147 * t130;
t89 = t120 * t142 + t193 * t121;
t189 = (-t103 * t186 + t89 * t182) / t102 ^ 2;
t119 = t121 ^ 2;
t93 = t119 * t127 * t96 + 0.1e1;
t91 = 0.1e1 / t93;
t188 = t91 * t96;
t187 = t95 * t91;
t184 = t121 * t96;
t181 = t111 * t129;
t178 = t123 * t127;
t174 = qJD(1) * t120;
t173 = qJD(1) * t121;
t145 = t120 * t127 * t173;
t164 = t96 * t170;
t168 = 0.2e1 * (-t96 * t145 + (-t127 * t190 - t131 * t164) * t119) / t93 ^ 2;
t167 = 0.2e1 * t190;
t166 = 0.2e1 * t189;
t126 = t131 * t127;
t140 = qJD(3) * (-t124 * t126 - t131) * t123;
t165 = 0.2e1 / t116 ^ 2 * (t118 * t140 + t124 * t145);
t162 = t114 * t178;
t161 = t120 * t172;
t156 = t95 * t168;
t155 = t96 * t168;
t154 = 0.1e1 + t177;
t153 = 0.2e1 * t107 * t186;
t152 = t120 * t165;
t151 = t131 * t165;
t146 = t120 * t162;
t144 = t154 * t121;
t141 = t104 * t130 + t128 * t182;
t98 = 0.1e1 / t102;
t139 = t141 * t98;
t110 = -t120 * t176 + t121 * t130;
t109 = t120 * t175 + t121 * t128;
t100 = t154 * t120 * t114;
t87 = (t191 * t131 * t111 + t112 * t146) * t121;
t86 = t120 * t181 + t180 + (-t112 * t179 - t181) * t100;
t84 = -t154 * t152 + (qJD(1) * t144 + 0.2e1 * t120 * t140) * t114;
t1 = [t121 * t123 * t151 + (qJD(3) * t144 + t123 * t161) * t114, 0, t84, 0, 0, 0; (t170 * t187 + (t156 + (qJD(1) * t87 + t83) * t188) * t131) * t120 + (t87 * t155 * t131 + (t87 * t164 + (t87 * t167 + ((t88 * t146 + t191 * t170 + t151) * t111 + (t152 * t178 + t131 * t88 + (t126 * t159 - (t88 - 0.2e1 * t171) * t131) * t114) * t112) * t184) * t131 + (-t95 + (-(-t118 + t119) * t112 * t162 + t191 * t163) * t96) * t172) * t91) * t121, 0 (t174 * t187 + (t156 + (qJD(3) * t86 + t83) * t188) * t121) * t129 + (t86 * t121 * t155 + (-t121 * qJD(3) * t95 + (t121 * t167 + t96 * t174) * t86 + (-((-t100 * t173 - t120 * t84) * t112 + (-t150 * t100 + t149) * t111) * t131 - ((-t84 + t173) * t111 + (t149 * t100 - t150) * t112) * t129) * t184) * t91) * t131, 0, 0, 0; (-t104 * t109 + t110 * t182) * t166 + (t110 * t153 + (-t109 * t90 - t110 * t89 + (t192 * t120 + t121 * t143) * t107) * t105 + (-t193 * t120 + t121 * t142) * t104) * t98, 0, t139 * t161 + (t139 * t170 + (t141 * t166 + ((qJD(5) * t104 + t153) * t128 + (-t128 * t89 + (-qJD(5) * t107 + t90) * t130) * t105) * t98) * t131) * t121, 0, -0.2e1 * t189 + 0.2e1 * (t105 * t89 * t98 + (-t105 * t189 - t98 * t186) * t107) * t107, 0;];
JaD_rot  = t1;
