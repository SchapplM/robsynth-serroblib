% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRP10_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:51:03
% EndTime: 2019-02-26 21:51:04
% DurationCPUTime: 0.34s
% Computational Cost: add. (1452->45), mult. (2519->109), div. (107->9), fcn. (3557->13), ass. (0->64)
t111 = sin(qJ(5));
t114 = cos(qJ(5));
t110 = cos(pkin(6));
t115 = cos(qJ(2));
t116 = cos(qJ(1));
t121 = t116 * t115;
t112 = sin(qJ(2));
t113 = sin(qJ(1));
t124 = t113 * t112;
t120 = -t110 * t121 + t124;
t122 = t116 * t112;
t123 = t113 * t115;
t102 = t110 * t122 + t123;
t108 = pkin(11) + qJ(4);
t106 = sin(t108);
t107 = cos(t108);
t109 = sin(pkin(6));
t126 = t109 * t116;
t94 = -t102 * t107 + t106 * t126;
t139 = t94 * t111 + t114 * t120;
t118 = t120 * t111;
t138 = t94 * t114 - t118;
t128 = t109 * t112;
t99 = t110 * t106 + t107 * t128;
t90 = t109 * t115 * t114 + t99 * t111;
t78 = atan2(t139, t90);
t74 = sin(t78);
t75 = cos(t78);
t73 = t139 * t74 + t75 * t90;
t72 = 0.1e1 / t73 ^ 2;
t103 = t110 * t123 + t122;
t129 = t103 * t114;
t104 = -t110 * t124 + t121;
t127 = t109 * t113;
t96 = t104 * t107 + t106 * t127;
t84 = t96 * t111 - t129;
t136 = t72 * t84;
t135 = t75 * t139;
t130 = t103 * t111;
t85 = t96 * t114 + t130;
t80 = 0.1e1 / t85 ^ 2;
t95 = -t104 * t106 + t107 * t127;
t134 = t80 * t95;
t88 = 0.1e1 / t90 ^ 2;
t133 = t139 * t88;
t132 = t84 ^ 2 * t72;
t131 = t95 ^ 2 * t80;
t125 = t111 * t115;
t119 = -t74 * t90 + t135;
t117 = t102 * t106 + t107 * t126;
t98 = -t106 * t128 + t110 * t107;
t97 = (t107 * t125 - t112 * t114) * t109;
t91 = -t109 * t125 + t99 * t114;
t87 = 0.1e1 / t90;
t86 = -t102 * t114 - t107 * t118;
t79 = 0.1e1 / t85;
t77 = 0.1e1 / (0.1e1 + t131);
t76 = 0.1e1 / (t139 ^ 2 * t88 + 0.1e1);
t71 = 0.1e1 / t73;
t70 = 0.1e1 / (0.1e1 + t132);
t69 = (t117 * t87 - t98 * t133) * t76 * t111;
t68 = (-t97 * t133 - t86 * t87) * t76;
t67 = (-t91 * t133 + t138 * t87) * t76;
t1 = [-t84 * t87 * t76, t68, 0, t69, t67, 0; (t139 * t71 - (-t74 + (-t87 * t135 + t74) * t76) * t132) * t70 ((-t104 * t114 - t107 * t130) * t71 - (t119 * t68 - t74 * t86 + t75 * t97) * t136) * t70, 0 (t95 * t111 * t71 - (t119 * t69 + (t117 * t74 + t75 * t98) * t111) * t136) * t70 (t85 * t71 - (t119 * t67 + t138 * t74 + t75 * t91) * t136) * t70, 0; (t117 * t79 - t138 * t134) * t77 (t103 * t106 * t79 - (t104 * t111 - t107 * t129) * t134) * t77, 0 (-t114 * t131 - t79 * t96) * t77, t84 * t77 * t134, 0;];
Ja_rot  = t1;
