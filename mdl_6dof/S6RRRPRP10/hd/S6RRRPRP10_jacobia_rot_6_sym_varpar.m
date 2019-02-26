% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRP10_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:14:24
% EndTime: 2019-02-26 22:14:24
% DurationCPUTime: 0.30s
% Computational Cost: add. (1427->45), mult. (2519->109), div. (107->9), fcn. (3557->13), ass. (0->64)
t108 = pkin(11) + qJ(5);
t106 = sin(t108);
t107 = cos(t108);
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
t111 = sin(qJ(3));
t114 = cos(qJ(3));
t109 = sin(pkin(6));
t125 = t109 * t116;
t94 = -t102 * t114 + t111 * t125;
t139 = t94 * t106 + t120 * t107;
t118 = t120 * t106;
t138 = t94 * t107 - t118;
t127 = t109 * t114;
t101 = t110 * t111 + t112 * t127;
t126 = t109 * t115;
t89 = t101 * t106 + t107 * t126;
t77 = atan2(t139, t89);
t74 = sin(t77);
t75 = cos(t77);
t73 = t139 * t74 + t75 * t89;
t72 = 0.1e1 / t73 ^ 2;
t103 = t110 * t123 + t122;
t130 = t103 * t107;
t104 = -t110 * t124 + t121;
t128 = t109 * t111;
t96 = t104 * t114 + t113 * t128;
t84 = t96 * t106 - t130;
t136 = t72 * t84;
t135 = t75 * t139;
t85 = t103 * t106 + t96 * t107;
t80 = 0.1e1 / t85 ^ 2;
t95 = -t104 * t111 + t113 * t127;
t134 = t80 * t95;
t88 = 0.1e1 / t89 ^ 2;
t133 = t139 * t88;
t132 = t84 ^ 2 * t72;
t131 = t95 ^ 2 * t80;
t129 = t106 * t114;
t119 = -t74 * t89 + t135;
t117 = t102 * t111 + t114 * t125;
t100 = t110 * t114 - t112 * t128;
t97 = (-t107 * t112 + t115 * t129) * t109;
t90 = t101 * t107 - t106 * t126;
t87 = 0.1e1 / t89;
t86 = -t102 * t107 - t114 * t118;
t79 = 0.1e1 / t85;
t78 = 0.1e1 / (0.1e1 + t131);
t76 = 0.1e1 / (t139 ^ 2 * t88 + 0.1e1);
t71 = 0.1e1 / t73;
t70 = 0.1e1 / (0.1e1 + t132);
t69 = (-t100 * t133 + t117 * t87) * t76 * t106;
t68 = (-t97 * t133 - t86 * t87) * t76;
t67 = (-t90 * t133 + t138 * t87) * t76;
t1 = [-t84 * t87 * t76, t68, t69, 0, t67, 0; (t139 * t71 - (-t74 + (-t87 * t135 + t74) * t76) * t132) * t70 ((-t103 * t129 - t104 * t107) * t71 - (t119 * t68 - t74 * t86 + t75 * t97) * t136) * t70 (t95 * t106 * t71 - (t119 * t69 + (t100 * t75 + t117 * t74) * t106) * t136) * t70, 0 (t85 * t71 - (t119 * t67 + t138 * t74 + t75 * t90) * t136) * t70, 0; (t117 * t79 - t138 * t134) * t78 (t103 * t111 * t79 - (t104 * t106 - t114 * t130) * t134) * t78 (-t107 * t131 - t79 * t96) * t78, 0, t84 * t78 * t134, 0;];
Ja_rot  = t1;
