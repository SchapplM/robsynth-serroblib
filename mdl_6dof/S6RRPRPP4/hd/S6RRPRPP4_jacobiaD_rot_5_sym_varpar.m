% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPP4_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP4_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_jacobiaD_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:36:41
% EndTime: 2019-02-26 21:36:41
% DurationCPUTime: 0.61s
% Computational Cost: add. (1161->91), mult. (2519->207), div. (480->12), fcn. (2968->9), ass. (0->91)
t136 = sin(qJ(1));
t138 = cos(qJ(1));
t135 = sin(qJ(2));
t154 = qJD(1) * t135 + qJD(4);
t137 = cos(qJ(2));
t174 = qJD(2) * t137;
t198 = t154 * t136 - t138 * t174;
t197 = t136 * t174 + t154 * t138;
t181 = t136 * t137;
t120 = atan2(-t181, t135);
t119 = cos(t120);
t118 = sin(t120);
t168 = t118 * t181;
t107 = t119 * t135 - t168;
t104 = 0.1e1 / t107;
t127 = qJ(4) + pkin(9);
t125 = sin(t127);
t126 = cos(t127);
t182 = t136 * t126;
t183 = t135 * t138;
t115 = t125 * t183 + t182;
t111 = 0.1e1 / t115;
t128 = 0.1e1 / t135;
t105 = 0.1e1 / t107 ^ 2;
t112 = 0.1e1 / t115 ^ 2;
t129 = 0.1e1 / t135 ^ 2;
t131 = t136 ^ 2;
t133 = t137 ^ 2;
t186 = t129 * t133;
t124 = t131 * t186 + 0.1e1;
t122 = 0.1e1 / t124;
t196 = t122 - 0.1e1;
t114 = t125 * t136 - t126 * t183;
t110 = t114 ^ 2;
t103 = t110 * t112 + 0.1e1;
t190 = t112 * t114;
t155 = qJD(4) * t135 + qJD(1);
t149 = t155 * t138;
t96 = -t198 * t125 + t126 * t149;
t192 = t111 * t112 * t96;
t95 = t125 * t149 + t198 * t126;
t195 = 0.1e1 / t103 ^ 2 * (-t110 * t192 + t95 * t190);
t134 = t138 ^ 2;
t185 = t133 * t134;
t102 = t105 * t185 + 0.1e1;
t98 = 0.1e1 / t102;
t194 = t105 * t98;
t175 = qJD(2) * t136;
t165 = t129 * t175;
t177 = qJD(1) * t138;
t166 = t137 * t177;
t97 = ((t135 * t175 - t166) * t128 + t133 * t165) * t122;
t156 = -t97 + t175;
t157 = -t136 * t97 + qJD(2);
t188 = t119 * t137;
t91 = t157 * t188 + (t156 * t135 - t166) * t118;
t193 = t104 * t105 * t91;
t191 = t111 * t126;
t189 = t114 * t125;
t187 = t128 * t133;
t184 = t135 * t136;
t179 = qJD(1) * t136;
t178 = qJD(1) * t137;
t176 = qJD(2) * t135;
t152 = t133 * t136 * t177;
t173 = 0.2e1 * (-t185 * t193 + (-t134 * t135 * t174 - t152) * t105) / t102 ^ 2;
t172 = 0.2e1 * t195;
t171 = 0.2e1 * t193;
t132 = t137 * t133;
t147 = qJD(2) * (-t129 * t132 - t137) * t128;
t170 = 0.2e1 * (t129 * t152 + t131 * t147) / t124 ^ 2;
t169 = t98 * t176;
t167 = t122 * t187;
t162 = 0.1e1 + t186;
t161 = t104 * t173;
t160 = 0.2e1 * t114 * t192;
t159 = t136 * t170;
t158 = t137 * t170;
t153 = t136 * t167;
t151 = t162 * t138;
t150 = t155 * t136;
t148 = t112 * t189 + t191;
t146 = t148 * t138;
t117 = -t125 * t184 + t126 * t138;
t116 = t125 * t138 + t135 * t182;
t109 = t162 * t136 * t122;
t100 = 0.1e1 / t103;
t94 = (t196 * t137 * t118 + t119 * t153) * t138;
t93 = t118 * t184 + t188 + (-t118 * t135 - t119 * t181) * t109;
t92 = -t162 * t159 + (qJD(1) * t151 + 0.2e1 * t136 * t147) * t122;
t1 = [t128 * t138 * t158 + (t128 * t136 * t178 + qJD(2) * t151) * t122, t92, 0, 0, 0, 0; (t104 * t169 + (t161 + (qJD(1) * t94 + t91) * t194) * t137) * t136 + (t94 * t137 * t98 * t171 + (t94 * t169 + (t94 * t173 + ((t97 * t153 + t196 * t176 + t158) * t118 + (t159 * t187 + t137 * t97 + (t132 * t165 - (t97 - 0.2e1 * t175) * t137) * t122) * t119) * t98 * t138) * t137) * t105 + (-t104 + ((t131 - t134) * t119 * t167 + t196 * t168) * t105) * t98 * t178) * t138 (t104 * t98 * t179 + (t161 + (qJD(2) * t93 + t91) * t194) * t138) * t135 + (t93 * t138 * t105 * t173 + ((-qJD(2) * t104 + t93 * t171) * t138 + (t93 * t179 + (-(-t109 * t177 - t136 * t92) * t119 - ((t109 * t136 - 0.1e1) * t97 + (-t109 + t136) * qJD(2)) * t118) * t137 * t138) * t105) * t98 - ((-t92 + t177) * t118 + (t156 * t109 - t157) * t119) * t183 * t194) * t137, 0, 0, 0, 0; (-t111 * t116 + t117 * t190) * t172 + (t117 * t160 - t111 * t125 * t150 + t197 * t191 + (t114 * t126 * t150 - t116 * t96 - t117 * t95 + t197 * t189) * t112) * t100, t137 * t146 * t172 + (t146 * t176 + (t148 * t179 + ((qJD(4) * t111 + t160) * t125 + (-t125 * t95 + (-qJD(4) * t114 + t96) * t126) * t112) * t138) * t137) * t100, 0, -0.2e1 * t195 + 0.2e1 * (t100 * t112 * t95 + (-t100 * t192 - t112 * t195) * t114) * t114, 0, 0;];
JaD_rot  = t1;
