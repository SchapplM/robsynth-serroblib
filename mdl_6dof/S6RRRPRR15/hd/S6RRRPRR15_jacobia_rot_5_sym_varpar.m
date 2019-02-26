% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR15_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:24:25
% EndTime: 2019-02-26 22:24:25
% DurationCPUTime: 0.30s
% Computational Cost: add. (833->47), mult. (2499->112), div. (85->9), fcn. (3431->15), ass. (0->64)
t105 = cos(pkin(7));
t107 = sin(qJ(3));
t111 = cos(qJ(3));
t103 = sin(pkin(7));
t104 = sin(pkin(6));
t113 = cos(qJ(1));
t126 = t104 * t113;
t121 = t103 * t126;
t108 = sin(qJ(2));
t109 = sin(qJ(1));
t112 = cos(qJ(2));
t128 = cos(pkin(6));
t119 = t113 * t128;
t96 = t109 * t108 - t112 * t119;
t97 = t108 * t119 + t109 * t112;
t81 = (t105 * t96 + t121) * t111 + t97 * t107;
t125 = t105 * t107;
t114 = t107 * t121 - t97 * t111 + t96 * t125;
t118 = t128 * t103;
t92 = t107 * t118 + (t108 * t111 + t112 * t125) * t104;
t80 = atan2(t114, t92);
t77 = sin(t80);
t78 = cos(t80);
t71 = t114 * t77 + t78 * t92;
t70 = 0.1e1 / t71 ^ 2;
t127 = t104 * t109;
t120 = t109 * t128;
t98 = -t113 * t108 - t112 * t120;
t116 = t103 * t127 + t105 * t98;
t99 = -t108 * t120 + t113 * t112;
t129 = t99 * t111;
t86 = t116 * t107 + t129;
t136 = t70 * t86;
t106 = sin(qJ(5));
t110 = cos(qJ(5));
t85 = t99 * t107 - t116 * t111;
t94 = -t98 * t103 + t105 * t127;
t76 = t85 * t106 + t94 * t110;
t74 = 0.1e1 / t76 ^ 2;
t75 = t94 * t106 - t85 * t110;
t135 = t74 * t75;
t134 = t78 * t114;
t90 = 0.1e1 / t92 ^ 2;
t133 = t114 * t90;
t132 = t86 ^ 2 * t70;
t131 = t103 * t99;
t124 = t107 * t108;
t123 = t111 * t112;
t122 = t75 ^ 2 * t74 + 0.1e1;
t117 = -t77 * t92 + t134;
t95 = (-t105 * t124 + t123) * t104;
t93 = -t96 * t103 + t105 * t126;
t91 = t111 * t118 + (t105 * t123 - t124) * t104;
t89 = 0.1e1 / t92;
t88 = t105 * t129 + t98 * t107;
t87 = -t96 * t111 - t97 * t125;
t79 = 0.1e1 / (t114 ^ 2 * t90 + 0.1e1);
t73 = 0.1e1 / t76;
t72 = 0.1e1 / t122;
t69 = 0.1e1 / t71;
t68 = 0.1e1 / (0.1e1 + t132);
t67 = (-t95 * t133 - t87 * t89) * t79;
t66 = (-t91 * t133 + t81 * t89) * t79;
t1 = [-t86 * t89 * t79, t67, t66, 0, 0, 0; (t114 * t69 - (-t77 + (-t89 * t134 + t77) * t79) * t132) * t68 ((t98 * t111 - t99 * t125) * t69 - (t117 * t67 - t77 * t87 + t78 * t95) * t136) * t68 (-t85 * t69 - (t117 * t66 + t77 * t81 + t78 * t91) * t136) * t68, 0, 0, 0; ((t93 * t106 + t110 * t81) * t73 - (-t106 * t81 + t93 * t110) * t135) * t72 ((t106 * t131 - t88 * t110) * t73 - (t88 * t106 + t110 * t131) * t135) * t72 (-t106 * t135 - t110 * t73) * t86 * t72, 0, t122 * t72, 0;];
Ja_rot  = t1;
