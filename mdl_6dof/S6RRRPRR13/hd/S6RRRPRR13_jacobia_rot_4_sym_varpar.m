% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
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

function Ja_rot = S6RRRPRR13_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_jacobia_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:04
% EndTime: 2019-02-26 22:23:04
% DurationCPUTime: 0.25s
% Computational Cost: add. (781->51), mult. (2353->113), div. (80->9), fcn. (3235->15), ass. (0->67)
t132 = sin(qJ(1));
t106 = cos(qJ(3));
t100 = sin(pkin(6));
t108 = cos(qJ(1));
t123 = t100 * t108;
t99 = sin(pkin(7));
t116 = t99 * t123;
t102 = cos(pkin(7));
t121 = t102 * t106;
t104 = sin(qJ(3));
t103 = cos(pkin(6));
t107 = cos(qJ(2));
t113 = t132 * t107;
t105 = sin(qJ(2));
t118 = t108 * t105;
t93 = t103 * t118 + t113;
t124 = t93 * t104;
t114 = t132 * t105;
t117 = t108 * t107;
t92 = -t103 * t117 + t114;
t75 = t106 * t116 + t92 * t121 + t124;
t125 = t103 * t99;
t85 = -t106 * t125 + (t104 * t105 - t107 * t121) * t100;
t74 = atan2(-t75, t85);
t71 = sin(t74);
t72 = cos(t74);
t65 = -t71 * t75 + t72 * t85;
t64 = 0.1e1 / t65 ^ 2;
t110 = t103 * t113 + t118;
t109 = t110 * t106;
t115 = t100 * t132;
t112 = t99 * t115;
t94 = -t103 * t114 + t117;
t79 = t102 * t109 + t94 * t104 - t106 * t112;
t131 = t64 * t79;
t101 = cos(pkin(13));
t80 = t94 * t106 + (-t110 * t102 + t112) * t104;
t88 = t102 * t115 + t110 * t99;
t98 = sin(pkin(13));
t70 = t80 * t101 + t88 * t98;
t68 = 0.1e1 / t70 ^ 2;
t69 = -t88 * t101 + t80 * t98;
t130 = t68 * t69;
t129 = t72 * t75;
t84 = 0.1e1 / t85 ^ 2;
t128 = t75 * t84;
t127 = t79 ^ 2 * t64;
t126 = t94 * t99;
t122 = t102 * t104;
t120 = t104 * t107;
t119 = t105 * t106;
t111 = -t71 * t85 - t129;
t78 = t104 * t116 - t93 * t106 + t92 * t122;
t91 = (t102 * t119 + t120) * t100;
t87 = t102 * t123 - t92 * t99;
t86 = t104 * t125 + (t102 * t120 + t119) * t100;
t83 = 0.1e1 / t85;
t82 = -t94 * t122 - t109;
t81 = -t92 * t104 + t93 * t121;
t73 = 0.1e1 / (t75 ^ 2 * t84 + 0.1e1);
t67 = 0.1e1 / t70;
t66 = 0.1e1 / (t69 ^ 2 * t68 + 0.1e1);
t63 = 0.1e1 / t65;
t62 = 0.1e1 / (0.1e1 + t127);
t61 = (t91 * t128 - t81 * t83) * t73;
t60 = (t86 * t128 + t78 * t83) * t73;
t1 = [-t79 * t83 * t73, t61, t60, 0, 0, 0; ((-t124 + (-t102 * t92 - t116) * t106) * t63 - (-t71 + (t83 * t129 + t71) * t73) * t127) * t62 ((-t110 * t104 + t94 * t121) * t63 - (t111 * t61 - t71 * t81 + t72 * t91) * t131) * t62 (t80 * t63 - (t111 * t60 + t71 * t78 + t72 * t86) * t131) * t62, 0, 0, 0; ((-t87 * t101 + t78 * t98) * t67 - (t78 * t101 + t87 * t98) * t130) * t66 ((-t101 * t126 + t82 * t98) * t67 - (t82 * t101 + t98 * t126) * t130) * t66 (t101 * t130 - t98 * t67) * t79 * t66, 0, 0, 0;];
Ja_rot  = t1;
