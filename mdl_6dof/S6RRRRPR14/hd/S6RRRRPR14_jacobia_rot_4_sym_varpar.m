% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR14_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_jacobia_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:11
% EndTime: 2019-02-26 22:38:11
% DurationCPUTime: 0.29s
% Computational Cost: add. (833->51), mult. (2499->114), div. (85->9), fcn. (3431->15), ass. (0->68)
t136 = sin(qJ(1));
t109 = cos(qJ(3));
t101 = sin(pkin(7));
t102 = sin(pkin(6));
t111 = cos(qJ(1));
t127 = t102 * t111;
t119 = t101 * t127;
t103 = cos(pkin(7));
t125 = t103 * t109;
t106 = sin(qJ(3));
t104 = cos(pkin(6));
t110 = cos(qJ(2));
t116 = t136 * t110;
t107 = sin(qJ(2));
t122 = t111 * t107;
t96 = t104 * t122 + t116;
t129 = t96 * t106;
t117 = t136 * t107;
t121 = t111 * t110;
t95 = -t104 * t121 + t117;
t78 = t109 * t119 + t95 * t125 + t129;
t128 = t101 * t104;
t88 = -t109 * t128 + (t106 * t107 - t110 * t125) * t102;
t77 = atan2(-t78, t88);
t74 = sin(t77);
t75 = cos(t77);
t68 = -t74 * t78 + t75 * t88;
t67 = 0.1e1 / t68 ^ 2;
t113 = t104 * t116 + t122;
t112 = t113 * t109;
t118 = t102 * t136;
t115 = t101 * t118;
t97 = -t104 * t117 + t121;
t82 = t103 * t112 + t97 * t106 - t109 * t115;
t135 = t67 * t82;
t105 = sin(qJ(4));
t108 = cos(qJ(4));
t83 = t97 * t109 + (-t113 * t103 + t115) * t106;
t91 = t113 * t101 + t103 * t118;
t73 = t91 * t105 + t83 * t108;
t71 = 0.1e1 / t73 ^ 2;
t72 = t83 * t105 - t91 * t108;
t134 = t71 * t72;
t133 = t75 * t78;
t87 = 0.1e1 / t88 ^ 2;
t132 = t78 * t87;
t131 = t82 ^ 2 * t67;
t130 = t101 * t97;
t126 = t103 * t106;
t124 = t106 * t110;
t123 = t107 * t109;
t120 = t72 ^ 2 * t71 + 0.1e1;
t114 = -t74 * t88 - t133;
t81 = t106 * t119 - t96 * t109 + t95 * t126;
t94 = (t103 * t123 + t124) * t102;
t90 = -t95 * t101 + t103 * t127;
t89 = t106 * t128 + (t103 * t124 + t123) * t102;
t86 = 0.1e1 / t88;
t85 = -t97 * t126 - t112;
t84 = -t95 * t106 + t96 * t125;
t76 = 0.1e1 / (t78 ^ 2 * t87 + 0.1e1);
t70 = 0.1e1 / t73;
t69 = 0.1e1 / t120;
t66 = 0.1e1 / t68;
t65 = 0.1e1 / (0.1e1 + t131);
t64 = (t94 * t132 - t84 * t86) * t76;
t63 = (t89 * t132 + t81 * t86) * t76;
t1 = [-t82 * t86 * t76, t64, t63, 0, 0, 0; ((-t129 + (-t103 * t95 - t119) * t109) * t66 - (-t74 + (t86 * t133 + t74) * t76) * t131) * t65 ((-t113 * t106 + t97 * t125) * t66 - (t114 * t64 - t74 * t84 + t75 * t94) * t135) * t65 (t83 * t66 - (t114 * t63 + t74 * t81 + t75 * t89) * t135) * t65, 0, 0, 0; ((t81 * t105 - t90 * t108) * t70 - (t90 * t105 + t81 * t108) * t134) * t69 ((t85 * t105 - t108 * t130) * t70 - (t105 * t130 + t85 * t108) * t134) * t69 (-t105 * t70 + t108 * t134) * t82 * t69, t120 * t69, 0, 0;];
Ja_rot  = t1;
