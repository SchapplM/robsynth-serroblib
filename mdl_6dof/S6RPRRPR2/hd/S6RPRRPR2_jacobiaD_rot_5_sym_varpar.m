% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPR2_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:01:26
% EndTime: 2019-02-26 21:01:26
% DurationCPUTime: 0.71s
% Computational Cost: add. (2324->95), mult. (2519->213), div. (480->12), fcn. (2968->9), ass. (0->94)
t135 = qJ(1) + pkin(10);
t131 = sin(t135);
t128 = t131 ^ 2;
t141 = sin(qJ(3));
t137 = t141 ^ 2;
t142 = cos(qJ(3));
t139 = 0.1e1 / t142 ^ 2;
t182 = t137 * t139;
t125 = t128 * t182 + 0.1e1;
t136 = t141 * t137;
t138 = 0.1e1 / t142;
t151 = qJD(3) * (t136 * t139 + t141) * t138;
t133 = cos(t135);
t180 = qJD(1) * t133;
t188 = t131 * t137;
t156 = t180 * t188;
t196 = (t128 * t151 + t139 * t156) / t125 ^ 2;
t204 = -0.2e1 * t196;
t163 = 0.1e1 + t182;
t203 = t131 * t163;
t134 = qJ(4) + pkin(11);
t130 = sin(t134);
t132 = cos(t134);
t183 = t133 * t142;
t118 = t131 * t130 + t132 * t183;
t187 = t131 * t141;
t122 = atan2(-t187, -t142);
t120 = cos(t122);
t119 = sin(t122);
t169 = t119 * t187;
t108 = -t120 * t142 - t169;
t105 = 0.1e1 / t108;
t112 = 0.1e1 / t118;
t106 = 0.1e1 / t108 ^ 2;
t113 = 0.1e1 / t118 ^ 2;
t123 = 0.1e1 / t125;
t202 = t123 - 0.1e1;
t129 = t133 ^ 2;
t191 = t129 * t137;
t101 = t106 * t191 + 0.1e1;
t176 = qJD(3) * t142;
t178 = qJD(3) * t131;
t165 = t139 * t178;
t179 = qJD(1) * t141;
t166 = t133 * t179;
t98 = (-(-t131 * t176 - t166) * t138 + t137 * t165) * t123;
t160 = t98 - t178;
t161 = -t131 * t98 + qJD(3);
t192 = t120 * t141;
t92 = t161 * t192 + (t160 * t142 - t166) * t119;
t198 = t105 * t106 * t92;
t201 = 0.1e1 / t101 ^ 2 * (-t191 * t198 + (t129 * t141 * t176 - t156) * t106);
t189 = t131 * t132;
t117 = t130 * t183 - t189;
t111 = t117 ^ 2;
t104 = t111 * t113 + 0.1e1;
t194 = t113 * t117;
t158 = -qJD(1) * t142 + qJD(4);
t159 = qJD(4) * t142 - qJD(1);
t177 = qJD(3) * t141;
t164 = t133 * t177;
t185 = t133 * t130;
t97 = -t159 * t185 + (t158 * t131 - t164) * t132;
t197 = t112 * t113 * t97;
t186 = t131 * t142;
t152 = t130 * t186 + t133 * t132;
t96 = t152 * qJD(1) - t118 * qJD(4) + t130 * t164;
t200 = 0.1e1 / t104 ^ 2 * (-t111 * t197 - t96 * t194);
t99 = 0.1e1 / t101;
t199 = t106 * t99;
t195 = t112 * t130;
t193 = t117 * t132;
t184 = t133 * t141;
t181 = qJD(1) * t131;
t175 = 0.2e1 * t201;
t174 = 0.2e1 * t200;
t173 = 0.2e1 * t198;
t172 = t105 * t201;
t171 = t117 * t197;
t170 = t99 * t176;
t168 = t123 * t137 * t138;
t162 = t138 * t204;
t157 = t131 * t168;
t155 = t163 * t133;
t154 = t158 * t133;
t153 = t113 * t193 - t195;
t150 = t153 * t141;
t116 = -t132 * t186 + t185;
t110 = t123 * t203;
t102 = 0.1e1 / t104;
t95 = (t202 * t141 * t119 - t120 * t157) * t133;
t94 = -t119 * t186 + t192 + (t119 * t142 - t120 * t187) * t110;
t93 = t203 * t204 + (qJD(1) * t155 + 0.2e1 * t131 * t151) * t123;
t1 = [t162 * t184 + (-t131 * t138 * t179 + qJD(3) * t155) * t123, 0, t93, 0, 0, 0; (-t105 * t170 + (0.2e1 * t172 + (qJD(1) * t95 + t92) * t199) * t141) * t131 + (t95 * t141 * t99 * t173 + (-t95 * t170 + (t95 * t175 + ((0.2e1 * t141 * t196 - t98 * t157 - t202 * t176) * t119 + (t162 * t188 + t141 * t98 + (t136 * t165 - (t98 - 0.2e1 * t178) * t141) * t123) * t120) * t99 * t133) * t141) * t106 + (-t105 + ((-t128 + t129) * t120 * t168 + t202 * t169) * t106) * t99 * t179) * t133, 0 (-t105 * t99 * t181 + (-0.2e1 * t172 + (-qJD(3) * t94 - t92) * t199) * t133) * t142 + (t94 * t133 * t106 * t175 + ((-qJD(3) * t105 + t94 * t173) * t133 + (t94 * t181 + (-(-t110 * t180 - t131 * t93) * t120 - ((t110 * t131 - 0.1e1) * t98 + (-t110 + t131) * qJD(3)) * t119) * t184) * t106) * t99 - ((t93 - t180) * t119 + (t160 * t110 + t161) * t120) * t183 * t199) * t141, 0, 0, 0; (t112 * t152 + t116 * t194) * t174 + (0.2e1 * t116 * t171 - t159 * t112 * t189 + (t131 * t177 + t154) * t195 + (t152 * t97 + t116 * t96 - t154 * t193 - (t159 * t130 + t132 * t177) * t117 * t131) * t113) * t102, 0, -t133 * t150 * t174 + (-t150 * t181 + (t153 * t176 + ((-qJD(4) * t112 - 0.2e1 * t171) * t132 + (-t132 * t96 + (-qJD(4) * t117 + t97) * t130) * t113) * t141) * t133) * t102, -0.2e1 * t200 + 0.2e1 * (-t102 * t113 * t96 + (-t102 * t197 - t113 * t200) * t117) * t117, 0, 0;];
JaD_rot  = t1;
