% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR13_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_jacobia_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:04
% EndTime: 2019-02-26 22:23:04
% DurationCPUTime: 0.29s
% Computational Cost: add. (899->52), mult. (2499->114), div. (85->9), fcn. (3431->15), ass. (0->69)
t143 = sin(qJ(1));
t113 = cos(pkin(6));
t115 = sin(qJ(2));
t124 = t143 * t115;
t117 = cos(qJ(2));
t118 = cos(qJ(1));
t128 = t118 * t117;
t101 = -t113 * t128 + t124;
t116 = cos(qJ(3));
t110 = sin(pkin(7));
t111 = sin(pkin(6));
t134 = t111 * t118;
t126 = t110 * t134;
t112 = cos(pkin(7));
t132 = t112 * t116;
t123 = t143 * t117;
t129 = t118 * t115;
t102 = t113 * t129 + t123;
t114 = sin(qJ(3));
t137 = t102 * t114;
t84 = t101 * t132 + t116 * t126 + t137;
t135 = t110 * t113;
t94 = -t116 * t135 + (t114 * t115 - t117 * t132) * t111;
t83 = atan2(-t84, t94);
t80 = sin(t83);
t81 = cos(t83);
t74 = -t80 * t84 + t81 * t94;
t73 = 0.1e1 / t74 ^ 2;
t103 = -t113 * t124 + t128;
t120 = t113 * t123 + t129;
t119 = t120 * t116;
t125 = t111 * t143;
t122 = t110 * t125;
t88 = t103 * t114 + t112 * t119 - t116 * t122;
t142 = t73 * t88;
t109 = pkin(13) + qJ(5);
t107 = sin(t109);
t108 = cos(t109);
t89 = t103 * t116 + (-t112 * t120 + t122) * t114;
t97 = t110 * t120 + t112 * t125;
t79 = t107 * t97 + t108 * t89;
t77 = 0.1e1 / t79 ^ 2;
t78 = t107 * t89 - t108 * t97;
t141 = t77 * t78;
t140 = t81 * t84;
t93 = 0.1e1 / t94 ^ 2;
t139 = t84 * t93;
t138 = t88 ^ 2 * t73;
t136 = t103 * t110;
t133 = t112 * t114;
t131 = t114 * t117;
t130 = t115 * t116;
t127 = t77 * t78 ^ 2 + 0.1e1;
t121 = -t80 * t94 - t140;
t87 = t101 * t133 - t102 * t116 + t114 * t126;
t100 = (t112 * t130 + t131) * t111;
t96 = -t101 * t110 + t112 * t134;
t95 = t114 * t135 + (t112 * t131 + t130) * t111;
t92 = 0.1e1 / t94;
t91 = -t103 * t133 - t119;
t90 = -t101 * t114 + t102 * t132;
t82 = 0.1e1 / (t84 ^ 2 * t93 + 0.1e1);
t76 = 0.1e1 / t79;
t75 = 0.1e1 / t127;
t72 = 0.1e1 / t74;
t71 = 0.1e1 / (0.1e1 + t138);
t70 = (t100 * t139 - t90 * t92) * t82;
t69 = (t139 * t95 + t87 * t92) * t82;
t1 = [-t88 * t92 * t82, t70, t69, 0, 0, 0; ((-t137 + (-t101 * t112 - t126) * t116) * t72 - (-t80 + (t140 * t92 + t80) * t82) * t138) * t71 ((t103 * t132 - t114 * t120) * t72 - (t81 * t100 + t121 * t70 - t80 * t90) * t142) * t71 (t89 * t72 - (t121 * t69 + t80 * t87 + t81 * t95) * t142) * t71, 0, 0, 0; ((t107 * t87 - t108 * t96) * t76 - (t107 * t96 + t108 * t87) * t141) * t75 ((t107 * t91 - t108 * t136) * t76 - (t107 * t136 + t108 * t91) * t141) * t75 (-t107 * t76 + t108 * t141) * t88 * t75, 0, t127 * t75, 0;];
Ja_rot  = t1;
