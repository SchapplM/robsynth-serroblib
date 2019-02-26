% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR15
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR15_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:57
% EndTime: 2019-02-26 22:38:57
% DurationCPUTime: 0.44s
% Computational Cost: add. (1577->57), mult. (4547->136), div. (107->9), fcn. (6260->15), ass. (0->71)
t127 = cos(pkin(6));
t133 = cos(qJ(2));
t159 = sin(qJ(1));
t139 = t159 * t133;
t130 = sin(qJ(2));
t134 = cos(qJ(1));
t142 = t134 * t130;
t119 = t127 * t142 + t139;
t129 = sin(qJ(3));
t132 = cos(qJ(3));
t140 = t159 * t130;
t143 = t133 * t134;
t118 = -t127 * t143 + t140;
t124 = sin(pkin(7));
t126 = cos(pkin(7));
t125 = sin(pkin(6));
t148 = t125 * t134;
t136 = t118 * t126 + t124 * t148;
t106 = -t119 * t132 + t136 * t129;
t115 = -t118 * t124 + t126 * t148;
t128 = sin(qJ(4));
t131 = cos(qJ(4));
t161 = t106 * t128 - t115 * t131;
t160 = t106 * t131 + t115 * t128;
t147 = t126 * t129;
t150 = t124 * t127;
t114 = t129 * t150 + (t130 * t132 + t133 * t147) * t125;
t117 = -t124 * t125 * t133 + t126 * t127;
t100 = t114 * t128 - t117 * t131;
t90 = atan2(t161, t100);
t87 = sin(t90);
t88 = cos(t90);
t86 = t100 * t88 + t161 * t87;
t85 = 0.1e1 / t86 ^ 2;
t120 = -t127 * t139 - t142;
t121 = -t127 * t140 + t143;
t141 = t125 * t159;
t138 = t124 * t141;
t108 = t121 * t132 + (t120 * t126 + t138) * t129;
t135 = -t120 * t124 + t126 * t141;
t95 = t108 * t128 - t135 * t131;
t158 = t85 * t95;
t157 = t85 * t95 ^ 2;
t156 = t88 * t161;
t99 = 0.1e1 / t100 ^ 2;
t155 = t161 * t99;
t146 = t126 * t132;
t107 = -t120 * t146 + t121 * t129 - t132 * t138;
t103 = 0.1e1 / t107 ^ 2;
t96 = t108 * t131 + t135 * t128;
t154 = t103 * t96;
t149 = t124 * t131;
t145 = t129 * t130;
t144 = t132 * t133;
t137 = -t100 * t87 + t156;
t104 = -t119 * t129 - t136 * t132;
t113 = t132 * t150 + (t126 * t144 - t145) * t125;
t110 = ((-t126 * t145 + t144) * t128 - t130 * t149) * t125;
t109 = t120 * t132 - t121 * t147;
t102 = 0.1e1 / t107;
t101 = t114 * t131 + t117 * t128;
t98 = 0.1e1 / t100;
t97 = (-t118 * t132 - t119 * t147) * t128 - t119 * t149;
t91 = 0.1e1 / (t103 * t96 ^ 2 + 0.1e1);
t89 = 0.1e1 / (t161 ^ 2 * t99 + 0.1e1);
t84 = 0.1e1 / t86;
t83 = 0.1e1 / (0.1e1 + t157);
t82 = (-t104 * t98 - t113 * t155) * t89 * t128;
t81 = (-t110 * t155 - t97 * t98) * t89;
t80 = (-t101 * t155 + t160 * t98) * t89;
t1 = [-t95 * t98 * t89, t81, t82, t80, 0, 0; (t161 * t84 - (-t87 + (-t98 * t156 + t87) * t89) * t157) * t83 ((t109 * t128 - t121 * t149) * t84 - (t110 * t88 + t137 * t81 - t87 * t97) * t158) * t83 (-t107 * t128 * t84 - (t137 * t82 + (-t104 * t87 + t113 * t88) * t128) * t158) * t83 (t96 * t84 - (t101 * t88 + t137 * t80 + t160 * t87) * t158) * t83, 0, 0; (t160 * t102 - t104 * t154) * t91 ((t121 * t124 * t128 + t109 * t131) * t102 - (t120 * t129 + t121 * t146) * t154) * t91 (-t102 * t107 * t131 - t108 * t154) * t91, -t95 * t102 * t91, 0, 0;];
Ja_rot  = t1;
