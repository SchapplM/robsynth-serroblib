% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR13
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR13_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:37:35
% EndTime: 2019-02-26 22:37:35
% DurationCPUTime: 0.24s
% Computational Cost: add. (707->44), mult. (1952->108), div. (91->9), fcn. (2758->15), ass. (0->63)
t133 = sin(qJ(1));
t105 = cos(pkin(6));
t113 = cos(qJ(2));
t120 = t133 * t113;
t109 = sin(qJ(2));
t114 = cos(qJ(1));
t124 = t114 * t109;
t100 = t105 * t124 + t120;
t108 = sin(qJ(3));
t112 = cos(qJ(3));
t126 = sin(pkin(6));
t117 = t114 * t126;
t88 = t100 * t108 + t112 * t117;
t119 = t109 * t126;
t97 = t105 * t112 - t108 * t119;
t85 = atan2(t88, t97);
t82 = sin(t85);
t83 = cos(t85);
t73 = t82 * t88 + t83 * t97;
t72 = 0.1e1 / t73 ^ 2;
t121 = t133 * t109;
t123 = t114 * t113;
t102 = -t105 * t121 + t123;
t116 = t126 * t133;
t91 = t102 * t108 - t112 * t116;
t132 = t72 * t91;
t106 = sin(qJ(6));
t110 = cos(qJ(6));
t101 = t105 * t120 + t124;
t107 = sin(qJ(4));
t111 = cos(qJ(4));
t93 = t102 * t112 + t108 * t116;
t80 = -t101 * t111 + t93 * t107;
t81 = t101 * t107 + t93 * t111;
t77 = t80 * t106 + t81 * t110;
t75 = 0.1e1 / t77 ^ 2;
t76 = t81 * t106 - t80 * t110;
t131 = t75 * t76;
t130 = t76 ^ 2 * t75;
t129 = t83 * t88;
t96 = 0.1e1 / t97 ^ 2;
t128 = t88 * t96;
t127 = t91 ^ 2 * t72;
t125 = t101 * t112;
t122 = 0.1e1 + t130;
t118 = t113 * t126;
t89 = t100 * t112 - t108 * t117;
t115 = -t82 * t97 + t129;
t99 = t105 * t123 - t121;
t98 = -t105 * t108 - t112 * t119;
t95 = 0.1e1 / t97;
t87 = t102 * t107 - t111 * t125;
t86 = -t102 * t111 - t107 * t125;
t84 = 0.1e1 / (t88 ^ 2 * t96 + 0.1e1);
t79 = t99 * t107 - t111 * t89;
t78 = -t107 * t89 - t99 * t111;
t74 = 0.1e1 / t77;
t71 = 0.1e1 / t73;
t70 = 0.1e1 / (0.1e1 + t127);
t69 = (t118 * t128 + t95 * t99) * t84 * t108;
t68 = (-t98 * t128 + t89 * t95) * t84;
t67 = 0.1e1 / t122;
t1 = [t91 * t95 * t84, t69, t68, 0, 0, 0; (t88 * t71 + (t82 + (t95 * t129 - t82) * t84) * t127) * t70 (t101 * t108 * t71 + (t115 * t69 + (-t83 * t118 + t82 * t99) * t108) * t132) * t70 (-t93 * t71 + (t115 * t68 + t82 * t89 + t83 * t98) * t132) * t70, 0, 0, 0; ((t79 * t106 - t78 * t110) * t74 - (t78 * t106 + t79 * t110) * t131) * t67 ((t87 * t106 - t86 * t110) * t74 - (t86 * t106 + t87 * t110) * t131) * t67 ((-t106 * t111 + t107 * t110) * t74 - (-t106 * t107 - t110 * t111) * t131) * t67 * t91 (-t77 * t74 - t130) * t67, 0, t122 * t67;];
Ja_rot  = t1;
