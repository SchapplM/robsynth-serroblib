% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRR14V3
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR14V3_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14V3_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobiaD_rot_4_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:12:04
% EndTime: 2019-04-12 15:12:05
% DurationCPUTime: 0.75s
% Computational Cost: add. (1002->94), mult. (2519->212), div. (480->12), fcn. (2968->9), ass. (0->91)
t125 = sin(qJ(1));
t118 = t125 ^ 2;
t124 = sin(qJ(2));
t117 = t124 ^ 2;
t127 = cos(qJ(2));
t120 = 0.1e1 / t127 ^ 2;
t172 = t117 * t120;
t112 = t118 * t172 + 0.1e1;
t116 = t124 * t117;
t119 = 0.1e1 / t127;
t171 = t119 * t124;
t136 = qJD(2) * (t116 * t119 * t120 + t171);
t128 = cos(qJ(1));
t162 = qJD(1) * t128;
t149 = t125 * t162;
t179 = 0.1e1 / t112 ^ 2 * (t118 * t136 + t149 * t172);
t191 = -0.2e1 * t179;
t110 = 0.1e1 / t112;
t145 = 0.1e1 + t172;
t188 = t125 * t145;
t97 = t110 * t188;
t190 = t125 * t97 - 0.1e1;
t123 = sin(qJ(4));
t126 = cos(qJ(4));
t164 = t127 * t128;
t106 = t125 * t123 + t126 * t164;
t101 = 0.1e1 / t106 ^ 2;
t166 = t125 * t126;
t105 = t123 * t164 - t166;
t174 = t105 * t126;
t100 = 0.1e1 / t106;
t176 = t100 * t123;
t138 = t101 * t174 - t176;
t99 = t105 ^ 2;
t98 = t101 * t99 + 0.1e1;
t95 = 0.1e1 / t98;
t189 = t138 * t95;
t167 = t125 * t124;
t109 = atan2(-t167, -t127);
t108 = cos(t109);
t107 = sin(t109);
t153 = t107 * t167;
t93 = -t108 * t127 - t153;
t90 = 0.1e1 / t93;
t91 = 0.1e1 / t93 ^ 2;
t187 = t110 - 0.1e1;
t122 = t128 ^ 2;
t160 = qJD(2) * t127;
t154 = t91 * t160;
t150 = t124 * t162;
t161 = qJD(2) * t125;
t173 = t108 * t124;
t148 = t120 * t161;
t84 = (-(-t125 * t160 - t150) * t119 + t117 * t148) * t110;
t79 = (-t125 * t84 + qJD(2)) * t173 + (-t150 + (t84 - t161) * t127) * t107;
t185 = t79 * t90 * t91;
t89 = t117 * t122 * t91 + 0.1e1;
t186 = (t122 * t124 * t154 + (-t122 * t185 - t91 * t149) * t117) / t89 ^ 2;
t175 = t101 * t105;
t142 = -qJD(1) * t127 + qJD(4);
t143 = qJD(4) * t127 - qJD(1);
t159 = qJD(2) * t128;
t147 = t124 * t159;
t170 = t123 * t128;
t86 = -t143 * t170 + (t142 * t125 - t147) * t126;
t181 = t100 * t101 * t86;
t165 = t125 * t127;
t137 = t123 * t165 + t126 * t128;
t85 = t137 * qJD(1) - t106 * qJD(4) + t123 * t147;
t184 = (-t85 * t175 - t99 * t181) / t98 ^ 2;
t87 = 0.1e1 / t89;
t183 = t87 * t91;
t182 = t90 * t87;
t177 = qJD(2) * t97;
t169 = t124 * t128;
t163 = qJD(1) * t125;
t158 = 0.2e1 * t185;
t157 = -0.2e1 * t184;
t156 = t90 * t186;
t155 = t105 * t181;
t152 = t110 * t117 * t119;
t146 = 0.2e1 * t91 * t186;
t144 = t119 * t191;
t141 = t125 * t152;
t140 = t145 * t128;
t139 = t142 * t128;
t104 = -t126 * t165 + t170;
t83 = (t187 * t124 * t107 - t108 * t141) * t128;
t82 = -t190 * t173 + (-t125 + t97) * t127 * t107;
t80 = t188 * t191 + (qJD(1) * t140 + 0.2e1 * t125 * t136) * t110;
t1 = [t144 * t169 + (qJD(2) * t140 - t163 * t171) * t110, t80, 0, 0, 0, 0; (-t160 * t182 + (0.2e1 * t156 + (qJD(1) * t83 + t79) * t183) * t124) * t125 + (t83 * t146 * t124 + (-t83 * t154 + (t83 * t158 + ((0.2e1 * t124 * t179 - t84 * t141 - t187 * t160) * t107 + (t117 * t125 * t144 + t124 * t84 + (t116 * t148 - (t84 - 0.2e1 * t161) * t124) * t110) * t108) * t91 * t128) * t124 + (-t90 + (-(t118 - t122) * t108 * t152 + t187 * t153) * t91) * t124 * qJD(1)) * t87) * t128 (-t163 * t182 + (-0.2e1 * t156 + (-qJD(2) * t82 - t79) * t183) * t128) * t127 + (t82 * t128 * t146 + (-t90 * t159 - ((-t125 * t80 - t162 * t97) * t108 + (t190 * t84 + t161 - t177) * t107) * t91 * t169 + (t128 * t158 + t91 * t163) * t82) * t87 - ((t80 - t162) * t107 + (t84 * t97 + qJD(2) + (-t84 - t177) * t125) * t108) * t164 * t183) * t124, 0, 0, 0, 0; 0.2e1 * (t100 * t137 + t104 * t175) * t184 + (0.2e1 * t104 * t155 - t143 * t100 * t166 + (t124 * t161 + t139) * t176 + (t137 * t86 + t104 * t85 - t139 * t174 - (qJD(2) * t124 * t126 + t143 * t123) * t105 * t125) * t101) * t95, t127 * t159 * t189 + (-t163 * t189 + (t138 * t157 + ((-qJD(4) * t100 - 0.2e1 * t155) * t126 + (-t126 * t85 + (-qJD(4) * t105 + t86) * t123) * t101) * t95) * t128) * t124, 0, t157 + 0.2e1 * (-t101 * t85 * t95 + (-t101 * t184 - t95 * t181) * t105) * t105, 0, 0;];
JaD_rot  = t1;
