% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR9_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_jacobia_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:05:29
% EndTime: 2019-02-26 21:05:29
% DurationCPUTime: 0.45s
% Computational Cost: add. (1859->50), mult. (4121->117), div. (85->9), fcn. (5660->17), ass. (0->71)
t127 = cos(pkin(6));
t125 = cos(pkin(12));
t159 = sin(qJ(1));
t140 = t159 * t125;
t122 = sin(pkin(12));
t132 = cos(qJ(1));
t145 = t132 * t122;
t116 = t127 * t145 + t140;
t129 = sin(qJ(3));
t131 = cos(qJ(3));
t141 = t159 * t122;
t144 = t132 * t125;
t115 = -t127 * t144 + t141;
t123 = sin(pkin(7));
t126 = cos(pkin(7));
t124 = sin(pkin(6));
t147 = t124 * t132;
t137 = t115 * t126 + t123 * t147;
t105 = -t116 * t131 + t137 * t129;
t112 = -t115 * t123 + t126 * t147;
t121 = qJ(4) + pkin(13);
t119 = sin(t121);
t120 = cos(t121);
t161 = t105 * t119 - t112 * t120;
t95 = t105 * t120 + t112 * t119;
t136 = t127 * t140 + t145;
t142 = t124 * t159;
t160 = -t123 * t142 + t136 * t126;
t146 = t125 * t126;
t148 = t123 * t127;
t111 = t129 * t148 + (t122 * t131 + t129 * t146) * t124;
t114 = -t124 * t125 * t123 + t127 * t126;
t100 = t111 * t119 - t114 * t120;
t89 = atan2(t161, t100);
t84 = sin(t89);
t85 = cos(t89);
t82 = t85 * t100 + t161 * t84;
t81 = 0.1e1 / t82 ^ 2;
t117 = -t127 * t141 + t144;
t107 = t117 * t131 - t160 * t129;
t133 = t136 * t123 + t126 * t142;
t96 = t107 * t119 - t133 * t120;
t158 = t81 * t96;
t157 = t85 * t161;
t130 = cos(qJ(6));
t106 = t117 * t129 + t160 * t131;
t128 = sin(qJ(6));
t153 = t106 * t128;
t97 = t107 * t120 + t133 * t119;
t91 = t97 * t130 + t153;
t88 = 0.1e1 / t91 ^ 2;
t152 = t106 * t130;
t90 = t97 * t128 - t152;
t156 = t88 * t90;
t99 = 0.1e1 / t100 ^ 2;
t155 = t161 * t99;
t154 = t96 ^ 2 * t81;
t143 = t90 ^ 2 * t88 + 0.1e1;
t138 = -t100 * t84 + t157;
t134 = -t116 * t129 - t137 * t131;
t110 = t131 * t148 + (-t122 * t129 + t131 * t146) * t124;
t101 = t111 * t120 + t114 * t119;
t98 = 0.1e1 / t100;
t87 = 0.1e1 / t91;
t86 = 0.1e1 / (t161 ^ 2 * t99 + 0.1e1);
t83 = 0.1e1 / t143;
t80 = 0.1e1 / t82;
t79 = 0.1e1 / (0.1e1 + t154);
t78 = (-t110 * t155 - t134 * t98) * t86 * t119;
t77 = (-t101 * t155 + t95 * t98) * t86;
t1 = [-t96 * t98 * t86, 0, t78, t77, 0, 0; (t161 * t80 - (-t84 + (-t98 * t157 + t84) * t86) * t154) * t79, 0 (-t106 * t119 * t80 - (t138 * t78 + (t110 * t85 - t134 * t84) * t119) * t158) * t79 (t97 * t80 - (t85 * t101 + t138 * t77 + t84 * t95) * t158) * t79, 0, 0; ((t95 * t128 - t130 * t134) * t87 - (t128 * t134 + t95 * t130) * t156) * t83, 0 ((-t107 * t130 - t120 * t153) * t87 - (t107 * t128 - t120 * t152) * t156) * t83 (-t128 * t87 + t130 * t156) * t96 * t83, 0, t143 * t83;];
Ja_rot  = t1;
