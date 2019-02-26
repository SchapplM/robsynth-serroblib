% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PPRRPR2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:57
% EndTime: 2019-02-26 19:40:57
% DurationCPUTime: 0.25s
% Computational Cost: add. (1048->43), mult. (3005->105), div. (65->9), fcn. (4119->17), ass. (0->66)
t102 = sin(pkin(11));
t107 = cos(pkin(7));
t101 = sin(pkin(12));
t105 = cos(pkin(12));
t106 = cos(pkin(11));
t108 = cos(pkin(6));
t126 = t102 * t108;
t116 = t106 * t101 + t105 * t126;
t103 = sin(pkin(7));
t104 = sin(pkin(6));
t124 = t104 * t103;
t132 = -t102 * t124 + t116 * t107;
t110 = sin(qJ(4));
t113 = cos(qJ(4));
t111 = sin(qJ(3));
t114 = cos(qJ(3));
t121 = t106 * t108;
t97 = -t102 * t101 + t105 * t121;
t117 = -t106 * t124 + t107 * t97;
t98 = t101 * t121 + t102 * t105;
t84 = t117 * t111 + t98 * t114;
t123 = t104 * t107;
t93 = -t97 * t103 - t106 * t123;
t79 = t93 * t110 + t84 * t113;
t122 = t105 * t107;
t125 = t103 * t108;
t92 = t111 * t125 + (t101 * t114 + t111 * t122) * t104;
t96 = -t105 * t124 + t108 * t107;
t90 = t96 * t110 + t92 * t113;
t77 = atan2(-t79, t90);
t74 = sin(t77);
t75 = cos(t77);
t68 = -t74 * t79 + t75 * t90;
t67 = 0.1e1 / t68 ^ 2;
t99 = -t101 * t126 + t106 * t105;
t86 = -t132 * t111 + t99 * t114;
t94 = t102 * t123 + t116 * t103;
t82 = t94 * t110 + t86 * t113;
t131 = t67 * t82;
t109 = sin(qJ(6));
t112 = cos(qJ(6));
t85 = t99 * t111 + t132 * t114;
t127 = t85 * t112;
t81 = t86 * t110 - t94 * t113;
t73 = t81 * t109 + t127;
t71 = 0.1e1 / t73 ^ 2;
t128 = t85 * t109;
t72 = -t81 * t112 + t128;
t130 = t71 * t72;
t88 = 0.1e1 / t90 ^ 2;
t129 = t79 * t88;
t120 = t72 ^ 2 * t71 + 0.1e1;
t118 = -t74 * t90 - t75 * t79;
t91 = t114 * t125 + (-t101 * t111 + t114 * t122) * t104;
t89 = -t92 * t110 + t96 * t113;
t87 = 0.1e1 / t90;
t83 = -t98 * t111 + t117 * t114;
t78 = t84 * t110 - t93 * t113;
t76 = 0.1e1 / (t79 ^ 2 * t88 + 0.1e1);
t70 = 0.1e1 / t73;
t69 = 0.1e1 / t120;
t66 = 0.1e1 / t68;
t65 = 0.1e1 / (t82 ^ 2 * t67 + 0.1e1);
t64 = (t91 * t129 - t83 * t87) * t76 * t113;
t63 = (t89 * t129 + t78 * t87) * t76;
t1 = [0, 0, t64, t63, 0, 0; 0, 0 (-t85 * t113 * t66 - (t118 * t64 + (-t74 * t83 + t75 * t91) * t113) * t131) * t65 (-t81 * t66 - (t118 * t63 + t74 * t78 + t75 * t89) * t131) * t65, 0, 0; 0, 0 ((t86 * t109 + t110 * t127) * t70 - (-t110 * t128 + t86 * t112) * t130) * t69 (-t109 * t130 - t112 * t70) * t82 * t69, 0, t120 * t69;];
Ja_rot  = t1;
