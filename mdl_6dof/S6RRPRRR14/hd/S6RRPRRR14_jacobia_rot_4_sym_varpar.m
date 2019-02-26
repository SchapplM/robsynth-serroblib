% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2019-01-03 10:25
% Revision: 5fdbc45bcf2cc60deefd7ac2d71d743ed41bf7e4 (2018-12-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR14_jacobia_rot_4_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobia_rot_4_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobia_rot_4_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-01-03 10:25:33
% EndTime: 2019-01-03 10:25:34
% DurationCPUTime: 0.30s
% Computational Cost: add. (955->48), mult. (2758->114), div. (55->9), fcn. (3778->17), ass. (0->66)
t117 = cos(pkin(6));
t122 = cos(qJ(2));
t123 = cos(qJ(1));
t131 = t123 * t122;
t119 = sin(qJ(2));
t120 = sin(qJ(1));
t134 = t120 * t119;
t105 = -t117 * t131 + t134;
t112 = sin(pkin(7));
t116 = cos(pkin(7));
t113 = sin(pkin(6));
t136 = t113 * t123;
t104 = -t105 * t112 + t116 * t136;
t111 = sin(pkin(8));
t115 = cos(pkin(8));
t132 = t123 * t119;
t133 = t120 * t122;
t106 = t117 * t132 + t133;
t110 = sin(pkin(14));
t114 = cos(pkin(14));
t127 = t105 * t116 + t112 * t136;
t95 = t106 * t110 + t127 * t114;
t88 = t104 * t115 - t95 * t111;
t118 = sin(qJ(4));
t121 = cos(qJ(4));
t107 = -t117 * t133 - t132;
t137 = t113 * t120;
t126 = -t107 * t112 + t116 * t137;
t125 = t107 * t116 + t112 * t137;
t108 = -t117 * t134 + t131;
t139 = t108 * t110;
t97 = t125 * t114 - t139;
t124 = t126 * t111 + t97 * t115;
t98 = t108 * t114 + t125 * t110;
t83 = t124 * t118 + t98 * t121;
t81 = 0.1e1 / t83 ^ 2;
t82 = t98 * t118 - t124 * t121;
t144 = t81 * t82;
t135 = t114 * t116;
t94 = -(t117 * t112 * t114 + (-t110 * t119 + t122 * t135) * t113) * t111 + (-t113 * t122 * t112 + t117 * t116) * t115;
t87 = atan2(t88, t94);
t85 = cos(t87);
t143 = t85 * t88;
t84 = sin(t87);
t78 = t84 * t88 + t85 * t94;
t77 = 0.1e1 / t78 ^ 2;
t89 = t97 * t111 - t126 * t115;
t142 = t89 ^ 2 * t77;
t138 = t112 * t115;
t130 = t82 ^ 2 * t81 + 0.1e1;
t129 = t104 * t111 + t115 * t95;
t99 = -t107 * t110 - t108 * t135;
t128 = t108 * t111 * t112 + t115 * t99;
t101 = (-(-t110 * t122 - t119 * t135) * t111 + t119 * t138) * t113;
t100 = t107 * t114 - t116 * t139;
t96 = -t106 * t114 + t127 * t110;
t93 = 0.1e1 / t94 ^ 2;
t92 = 0.1e1 / t94;
t91 = (t105 * t110 - t106 * t135) * t111 - t106 * t138;
t86 = 0.1e1 / (t88 ^ 2 * t93 + 0.1e1);
t80 = 0.1e1 / t83;
t79 = 0.1e1 / t130;
t76 = 0.1e1 / t78;
t75 = 0.1e1 / (0.1e1 + t142);
t74 = (-t101 * t88 * t93 + t91 * t92) * t86;
t1 = [t89 * t92 * t86, t74, 0, 0, 0, 0; (t88 * t76 + (t84 + (t92 * t143 - t84) * t86) * t142) * t75 ((t108 * t138 - t99 * t111) * t76 + (t85 * t101 + t84 * t91 + (-t84 * t94 + t143) * t74) * t89 * t77) * t75, 0, 0, 0, 0; ((t96 * t118 - t129 * t121) * t80 - (t129 * t118 + t96 * t121) * t144) * t79 ((t100 * t118 - t128 * t121) * t80 - (t100 * t121 + t128 * t118) * t144) * t79, 0, t130 * t79, 0, 0;];
Ja_rot  = t1;
