% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRRR10
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRR10_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobia_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:52:56
% EndTime: 2019-02-26 22:52:57
% DurationCPUTime: 0.39s
% Computational Cost: add. (1426->58), mult. (4126->139), div. (85->9), fcn. (5651->17), ass. (0->76)
t128 = cos(pkin(6));
t131 = sin(qJ(2));
t167 = sin(qJ(1));
t145 = t167 * t131;
t134 = cos(qJ(2));
t135 = cos(qJ(1));
t149 = t135 * t134;
t117 = -t128 * t149 + t145;
t144 = t167 * t134;
t150 = t135 * t131;
t118 = t128 * t150 + t144;
t127 = cos(pkin(7));
t130 = sin(qJ(3));
t133 = cos(qJ(3));
t124 = sin(pkin(7));
t125 = sin(pkin(6));
t155 = t125 * t135;
t147 = t124 * t155;
t106 = (t117 * t127 + t147) * t133 + t118 * t130;
t116 = -t117 * t124 + t127 * t155;
t123 = sin(pkin(8));
t126 = cos(pkin(8));
t98 = -t106 * t123 + t116 * t126;
t153 = t127 * t133;
t156 = t124 * t128;
t104 = -(t133 * t156 + (-t130 * t131 + t134 * t153) * t125) * t123 + (-t125 * t134 * t124 + t128 * t127) * t126;
t97 = atan2(t98, t104);
t94 = sin(t97);
t95 = cos(t97);
t88 = t95 * t104 + t94 * t98;
t87 = 0.1e1 / t88 ^ 2;
t120 = -t128 * t145 + t149;
t119 = -t128 * t144 - t150;
t146 = t125 * t167;
t137 = t119 * t127 + t124 * t146;
t108 = -t120 * t130 + t137 * t133;
t138 = t119 * t124 - t127 * t146;
t99 = t108 * t123 + t138 * t126;
t166 = t87 * t99;
t129 = sin(qJ(4));
t136 = t108 * t126 - t138 * t123;
t109 = t120 * t133 + t137 * t130;
t132 = cos(qJ(4));
t160 = t109 * t132;
t93 = t136 * t129 + t160;
t91 = 0.1e1 / t93 ^ 2;
t161 = t109 * t129;
t92 = -t136 * t132 + t161;
t165 = t91 * t92;
t164 = t95 * t98;
t163 = t99 ^ 2 * t87;
t103 = 0.1e1 / t104 ^ 2;
t162 = t103 * t98;
t157 = t124 * t126;
t154 = t127 * t130;
t152 = t130 * t134;
t151 = t131 * t133;
t148 = t92 ^ 2 * t91 + 0.1e1;
t143 = -t104 * t94 + t164;
t142 = t106 * t126 + t116 * t123;
t110 = -t119 * t130 - t120 * t153;
t141 = t120 * t123 * t124 + t110 * t126;
t139 = t117 * t154 - t118 * t133 + t130 * t147;
t115 = -t130 * t156 + (-t127 * t152 - t151) * t125;
t112 = (-(-t127 * t151 - t152) * t123 + t131 * t157) * t125;
t111 = t119 * t133 - t120 * t154;
t102 = 0.1e1 / t104;
t101 = (t117 * t130 - t118 * t153) * t123 - t118 * t157;
t96 = 0.1e1 / (t98 ^ 2 * t103 + 0.1e1);
t90 = 0.1e1 / t93;
t89 = 0.1e1 / t148;
t86 = 0.1e1 / t88;
t85 = 0.1e1 / (0.1e1 + t163);
t84 = (t102 * t139 + t115 * t162) * t96 * t123;
t83 = (t101 * t102 - t112 * t162) * t96;
t1 = [t99 * t102 * t96, t83, t84, 0, 0, 0; (t98 * t86 + (t94 + (t102 * t164 - t94) * t96) * t163) * t85 ((-t110 * t123 + t120 * t157) * t86 + (t94 * t101 + t95 * t112 + t143 * t83) * t166) * t85 (t109 * t123 * t86 + (t143 * t84 + (-t115 * t95 + t139 * t94) * t123) * t166) * t85, 0, 0, 0; ((t129 * t139 - t142 * t132) * t90 - (t142 * t129 + t132 * t139) * t165) * t89 ((t111 * t129 - t141 * t132) * t90 - (t111 * t132 + t141 * t129) * t165) * t89 ((t108 * t129 + t126 * t160) * t90 - (t108 * t132 - t126 * t161) * t165) * t89, t148 * t89, 0, 0;];
Ja_rot  = t1;
