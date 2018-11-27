% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% JaD_rot [3x7]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_rot = S7RRRRRRR1_jacobiaD_rot_3_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobiaD_rot_3_floatb_twist_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_jacobiaD_rot_3_floatb_twist_sym_varpar: qJD has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobiaD_rot_3_floatb_twist_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-26 21:20:53
% EndTime: 2018-11-26 21:20:54
% DurationCPUTime: 0.73s
% Computational Cost: add. (624->92), mult. (2519->212), div. (480->12), fcn. (2968->9), ass. (0->92)
t123 = sin(qJ(1));
t116 = t123 ^ 2;
t122 = sin(qJ(2));
t115 = t122 ^ 2;
t125 = cos(qJ(2));
t118 = 0.1e1 / t125 ^ 2;
t168 = t115 * t118;
t111 = t116 * t168 + 0.1e1;
t114 = t122 * t115;
t117 = 0.1e1 / t125;
t167 = t117 * t122;
t134 = qJD(2) * (t114 * t117 * t118 + t167);
t126 = cos(qJ(1));
t159 = qJD(1) * t126;
t146 = t123 * t159;
t175 = 0.1e1 / t111 ^ 2 * (t116 * t134 + t146 * t168);
t188 = -0.2e1 * t175;
t108 = 0.1e1 / t111;
t142 = 0.1e1 + t168;
t185 = t123 * t142;
t96 = t108 * t185;
t187 = -t123 * t96 + 0.1e1;
t121 = sin(qJ(3));
t124 = cos(qJ(3));
t161 = t125 * t126;
t148 = t124 * t161;
t105 = -t123 * t121 + t148;
t100 = 0.1e1 / t105 ^ 2;
t163 = t123 * t124;
t104 = t121 * t161 + t163;
t170 = t104 * t124;
t99 = 0.1e1 / t105;
t174 = t121 * t99;
t135 = t100 * t170 - t174;
t98 = t104 ^ 2;
t97 = t100 * t98 + 0.1e1;
t94 = 0.1e1 / t97;
t186 = t135 * t94;
t139 = qJD(1) * t125 + qJD(3);
t156 = qJD(2) * t126;
t184 = t122 * t156 + t139 * t123;
t164 = t123 * t122;
t110 = atan2(t164, t125);
t107 = cos(t110);
t106 = sin(t110);
t150 = t106 * t164;
t92 = t107 * t125 + t150;
t89 = 0.1e1 / t92;
t90 = 0.1e1 / t92 ^ 2;
t183 = t108 - 0.1e1;
t120 = t126 ^ 2;
t157 = qJD(2) * t125;
t151 = t90 * t157;
t147 = t122 * t159;
t158 = qJD(2) * t123;
t169 = t107 * t122;
t145 = t118 * t158;
t83 = ((t123 * t157 + t147) * t117 + t115 * t145) * t108;
t78 = (t123 * t83 - qJD(2)) * t169 + (t147 + (-t83 + t158) * t125) * t106;
t181 = t78 * t89 * t90;
t88 = t115 * t120 * t90 + 0.1e1;
t182 = (t120 * t122 * t151 + (-t120 * t181 - t90 * t146) * t115) / t88 ^ 2;
t171 = t100 * t104;
t140 = qJD(3) * t125 + qJD(1);
t166 = t121 * t126;
t85 = -t184 * t124 - t140 * t166;
t177 = t99 * t100 * t85;
t84 = -qJD(3) * t148 + t184 * t121 - t124 * t159;
t180 = (-t84 * t171 - t98 * t177) / t97 ^ 2;
t86 = 0.1e1 / t88;
t179 = t86 * t90;
t178 = t89 * t86;
t172 = qJD(2) * t96;
t165 = t122 * t126;
t162 = t123 * t125;
t160 = qJD(1) * t123;
t155 = -0.2e1 * t181;
t154 = -0.2e1 * t180;
t153 = t89 * t182;
t152 = t104 * t177;
t149 = t108 * t115 * t117;
t143 = -0.2e1 * t90 * t182;
t141 = t117 * t188;
t138 = t123 * t149;
t137 = t142 * t126;
t136 = t139 * t126;
t103 = -t124 * t162 - t166;
t102 = -t121 * t162 + t124 * t126;
t82 = (-t183 * t122 * t106 + t107 * t138) * t126;
t81 = -t187 * t169 + (t123 - t96) * t125 * t106;
t79 = t185 * t188 + (qJD(1) * t137 + 0.2e1 * t123 * t134) * t108;
t1 = [t141 * t165 + (qJD(2) * t137 - t160 * t167) * t108, t79, 0, 0, 0, 0, 0; (t157 * t178 + (-0.2e1 * t153 + (-qJD(1) * t82 - t78) * t179) * t122) * t123 + (t82 * t143 * t122 + (t82 * t151 + (t82 * t155 + ((0.2e1 * t122 * t175 - t83 * t138 - t183 * t157) * t106 + (t115 * t123 * t141 + t122 * t83 + (t114 * t145 + (-t83 + 0.2e1 * t158) * t122) * t108) * t107) * t90 * t126) * t122 + (t89 + ((-t116 + t120) * t107 * t149 + t183 * t150) * t90) * t122 * qJD(1)) * t86) * t126 (t160 * t178 + (0.2e1 * t153 + (qJD(2) * t81 + t78) * t179) * t126) * t125 + (t81 * t126 * t143 + (t89 * t156 + ((t123 * t79 + t159 * t96) * t107 + (t187 * t83 - t158 + t172) * t106) * t90 * t165 + (t126 * t155 - t90 * t160) * t81) * t86 + ((-t79 + t159) * t106 + (-t83 * t96 - qJD(2) + (t83 + t172) * t123) * t107) * t161 * t179) * t122, 0, 0, 0, 0, 0; 0.2e1 * (-t102 * t99 + t103 * t171) * t180 + (0.2e1 * t103 * t152 - t140 * t99 * t163 + (t122 * t158 - t136) * t174 + (-t102 * t85 + t103 * t84 + t136 * t170 - (qJD(2) * t122 * t124 + t140 * t121) * t104 * t123) * t100) * t94, t125 * t156 * t186 + (-t160 * t186 + (t135 * t154 + ((-qJD(3) * t99 - 0.2e1 * t152) * t124 + (-t124 * t84 + (-qJD(3) * t104 + t85) * t121) * t100) * t94) * t126) * t122, t154 + 0.2e1 * (-t100 * t84 * t94 + (-t100 * t180 - t94 * t177) * t104) * t104, 0, 0, 0, 0;];
JaD_rot  = t1;
