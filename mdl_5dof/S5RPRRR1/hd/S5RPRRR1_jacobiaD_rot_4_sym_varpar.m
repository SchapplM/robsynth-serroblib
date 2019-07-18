% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:26
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPRRR1_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_jacobiaD_rot_4_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_jacobiaD_rot_4_sym_varpar: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_jacobiaD_rot_4_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:26:25
% EndTime: 2019-07-18 13:26:26
% DurationCPUTime: 0.82s
% Computational Cost: add. (1002->94), mult. (2519->208), div. (480->12), fcn. (2968->9), ass. (0->94)
t121 = sin(qJ(4));
t123 = sin(qJ(1));
t124 = cos(qJ(4));
t125 = cos(qJ(3));
t126 = cos(qJ(1));
t163 = t125 * t126;
t103 = t121 * t163 - t123 * t124;
t192 = 0.2e1 * t103;
t116 = t123 ^ 2;
t122 = sin(qJ(3));
t115 = t122 ^ 2;
t118 = 0.1e1 / t125 ^ 2;
t170 = t115 * t118;
t110 = t116 * t170 + 0.1e1;
t114 = t122 * t115;
t117 = 0.1e1 / t125;
t169 = t117 * t122;
t134 = qJD(3) * (t114 * t117 * t118 + t169);
t161 = qJD(1) * t126;
t147 = t123 * t161;
t175 = 0.1e1 / t110 ^ 2 * (t116 * t134 + t147 * t170);
t191 = -0.2e1 * t175;
t108 = 0.1e1 / t110;
t142 = 0.1e1 + t170;
t187 = t123 * t142;
t95 = t108 * t187;
t190 = t123 * t95 - 0.1e1;
t140 = -qJD(1) * t125 + qJD(4);
t157 = qJD(3) * t126;
t145 = t122 * t157;
t168 = t121 * t126;
t156 = qJD(4) * t125;
t186 = qJD(1) - t156;
t84 = t186 * t168 + (t140 * t123 - t145) * t124;
t166 = t123 * t121;
t104 = t124 * t163 + t166;
t99 = 0.1e1 / t104 ^ 2;
t189 = t84 * t99;
t177 = t103 * t99;
t98 = 0.1e1 / t104;
t136 = -t121 * t98 + t124 * t177;
t97 = t103 ^ 2;
t96 = t97 * t99 + 0.1e1;
t93 = 0.1e1 / t96;
t188 = t136 * t93;
t165 = t123 * t122;
t107 = atan2(-t165, -t125);
t106 = cos(t107);
t105 = sin(t107);
t150 = t105 * t165;
t91 = -t106 * t125 - t150;
t88 = 0.1e1 / t91;
t89 = 0.1e1 / t91 ^ 2;
t185 = t108 - 0.1e1;
t120 = t126 ^ 2;
t158 = qJD(3) * t125;
t151 = t89 * t158;
t148 = t122 * t161;
t159 = qJD(3) * t123;
t171 = t106 * t122;
t146 = t118 * t159;
t82 = (-(-t123 * t158 - t148) * t117 + t115 * t146) * t108;
t77 = (-t123 * t82 + qJD(3)) * t171 + (-t148 + (t82 - t159) * t125) * t105;
t183 = t77 * t88 * t89;
t87 = t115 * t120 * t89 + 0.1e1;
t184 = (t120 * t122 * t151 + (-t120 * t183 - t89 * t147) * t115) / t87 ^ 2;
t178 = t98 * t189;
t164 = t123 * t125;
t135 = t121 * t164 + t124 * t126;
t144 = t124 * t156;
t83 = t135 * qJD(1) - qJD(4) * t166 + t121 * t145 - t126 * t144;
t182 = (-t83 * t177 - t97 * t178) / t96 ^ 2;
t181 = t83 * t99;
t85 = 0.1e1 / t87;
t180 = t85 * t89;
t179 = t88 * t85;
t172 = qJD(3) * t95;
t167 = t122 * t126;
t162 = qJD(1) * t123;
t160 = qJD(3) * t122;
t155 = 0.2e1 * t183;
t154 = -0.2e1 * t182;
t153 = t88 * t184;
t149 = t108 * t115 * t117;
t143 = 0.2e1 * t89 * t184;
t141 = t117 * t191;
t139 = t123 * t149;
t138 = t142 * t126;
t137 = t178 * t192 + t181;
t102 = -t124 * t164 + t168;
t81 = (t185 * t122 * t105 - t106 * t139) * t126;
t80 = -t190 * t171 + (-t123 + t95) * t125 * t105;
t78 = t187 * t191 + (qJD(1) * t138 + 0.2e1 * t123 * t134) * t108;
t1 = [t141 * t167 + (qJD(3) * t138 - t162 * t169) * t108, 0, t78, 0, 0; (-t158 * t179 + (0.2e1 * t153 + (qJD(1) * t81 + t77) * t180) * t122) * t123 + (t81 * t143 * t122 + (-t81 * t151 + (t81 * t155 + ((0.2e1 * t122 * t175 - t82 * t139 - t185 * t158) * t105 + (t115 * t123 * t141 + t122 * t82 + (t114 * t146 - (t82 - 0.2e1 * t159) * t122) * t108) * t106) * t89 * t126) * t122 + (-t88 + (-(t116 - t120) * t106 * t149 + t185 * t150) * t89) * t122 * qJD(1)) * t85) * t126, 0, (-t162 * t179 + (-0.2e1 * t153 + (-qJD(3) * t80 - t77) * t180) * t126) * t125 + (t80 * t126 * t143 + (-t88 * t157 - ((-t123 * t78 - t161 * t95) * t106 + (t190 * t82 + t159 - t172) * t105) * t89 * t167 + (t126 * t155 + t89 * t162) * t80) * t85 - ((t78 - t161) * t105 + (t82 * t95 + qJD(3) + (-t82 - t172) * t123) * t106) * t163 * t180) * t122, 0, 0; 0.2e1 * (t102 * t177 + t135 * t98) * t182 + (t135 * t189 + t137 * t102 + ((qJD(1) * t124 + t121 * t160 - t144) * t98 - (-t121 * t186 + t124 * t160) * t177) * t123 - t136 * t126 * t140) * t93, 0, t125 * t157 * t188 + (-t162 * t188 + (t136 * t154 + ((-qJD(4) * t103 + t84) * t99 * t121 + (-qJD(4) * t98 - t137) * t124) * t93) * t126) * t122, t154 + (-t93 * t181 + (-t93 * t178 - t99 * t182) * t103) * t192, 0;];
JaD_rot  = t1;
