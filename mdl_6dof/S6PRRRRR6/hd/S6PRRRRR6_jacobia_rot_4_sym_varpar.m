% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRR6_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_jacobia_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:21:51
% EndTime: 2019-02-26 20:21:51
% DurationCPUTime: 0.26s
% Computational Cost: add. (1034->50), mult. (3001->123), div. (65->9), fcn. (4101->17), ass. (0->69)
t103 = sin(pkin(8));
t104 = sin(pkin(7));
t106 = cos(pkin(14));
t107 = cos(pkin(8));
t111 = sin(qJ(3));
t114 = cos(qJ(3));
t108 = cos(pkin(7));
t105 = sin(pkin(6));
t130 = t104 * t105;
t102 = sin(pkin(14));
t112 = sin(qJ(2));
t109 = cos(pkin(6));
t115 = cos(qJ(2));
t124 = t109 * t115;
t97 = -t102 * t112 + t106 * t124;
t116 = t106 * t130 - t108 * t97;
t127 = t105 * t108;
t125 = t109 * t112;
t98 = t102 * t115 + t106 * t125;
t82 = (-t98 * t111 - t114 * t116) * t103 - (-t97 * t104 - t106 * t127) * t107;
t126 = t108 * t114;
t128 = t104 * t109;
t87 = -(t114 * t128 + (-t111 * t112 + t115 * t126) * t105) * t103 + (t109 * t108 - t115 * t130) * t107;
t81 = atan2(t82, t87);
t78 = sin(t81);
t79 = cos(t81);
t72 = t78 * t82 + t79 * t87;
t71 = 0.1e1 / t72 ^ 2;
t99 = -t102 * t124 - t106 * t112;
t117 = t102 * t130 + t108 * t99;
t100 = -t102 * t125 + t106 * t115;
t131 = t100 * t111;
t89 = t114 * t117 - t131;
t96 = t102 * t127 - t99 * t104;
t83 = -t89 * t103 + t96 * t107;
t136 = t71 * t83;
t110 = sin(qJ(4));
t119 = t103 * t96 + t107 * t89;
t113 = cos(qJ(4));
t90 = t100 * t114 + t111 * t117;
t132 = t90 * t113;
t77 = t110 * t119 + t132;
t75 = 0.1e1 / t77 ^ 2;
t133 = t90 * t110;
t76 = -t113 * t119 + t133;
t135 = t75 * t76;
t86 = 0.1e1 / t87 ^ 2;
t134 = t82 * t86;
t129 = t104 * t107;
t123 = t111 * t115;
t122 = t112 * t114;
t121 = t76 ^ 2 * t75 + 0.1e1;
t120 = -t78 * t87 + t79 * t82;
t91 = -t100 * t126 - t99 * t111;
t118 = t100 * t103 * t104 + t107 * t91;
t95 = -t111 * t128 + (-t108 * t123 - t122) * t105;
t93 = (-(-t108 * t122 - t123) * t103 + t112 * t129) * t105;
t92 = -t108 * t131 + t99 * t114;
t88 = t111 * t116 - t98 * t114;
t85 = 0.1e1 / t87;
t84 = (-t97 * t111 - t126 * t98) * t103 - t98 * t129;
t80 = 0.1e1 / (t82 ^ 2 * t86 + 0.1e1);
t74 = 0.1e1 / t77;
t73 = 0.1e1 / t121;
t70 = 0.1e1 / t72;
t69 = 0.1e1 / (t83 ^ 2 * t71 + 0.1e1);
t68 = (t134 * t95 + t85 * t88) * t80 * t103;
t67 = (-t134 * t93 + t84 * t85) * t80;
t1 = [0, t67, t68, 0, 0, 0; 0 ((t100 * t129 - t91 * t103) * t70 - (t120 * t67 + t78 * t84 + t79 * t93) * t136) * t69 (t90 * t103 * t70 - (t120 * t68 + (t78 * t88 - t79 * t95) * t103) * t136) * t69, 0, 0, 0; 0 ((t92 * t110 - t113 * t118) * t74 - (t110 * t118 + t92 * t113) * t135) * t73 ((t107 * t132 + t89 * t110) * t74 - (-t107 * t133 + t89 * t113) * t135) * t73, t121 * t73, 0, 0;];
Ja_rot  = t1;
