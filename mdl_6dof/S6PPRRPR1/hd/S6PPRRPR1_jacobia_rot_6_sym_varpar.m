% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PPRRPR1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_jacobia_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:19
% EndTime: 2019-02-26 19:40:19
% DurationCPUTime: 0.26s
% Computational Cost: add. (1096->44), mult. (3005->105), div. (65->9), fcn. (4119->17), ass. (0->67)
t108 = sin(pkin(11));
t113 = cos(pkin(7));
t107 = sin(pkin(12));
t111 = cos(pkin(12));
t112 = cos(pkin(11));
t114 = cos(pkin(6));
t133 = t108 * t114;
t123 = t112 * t107 + t111 * t133;
t109 = sin(pkin(7));
t110 = sin(pkin(6));
t131 = t110 * t109;
t139 = -t108 * t131 + t123 * t113;
t115 = sin(qJ(4));
t117 = cos(qJ(4));
t128 = t112 * t114;
t124 = -t108 * t107 + t111 * t128;
t130 = t110 * t113;
t121 = -t124 * t109 - t112 * t130;
t101 = t107 * t128 + t108 * t111;
t116 = sin(qJ(3));
t118 = cos(qJ(3));
t120 = -t112 * t131 + t124 * t113;
t88 = t101 * t118 + t120 * t116;
t82 = t88 * t115 - t121 * t117;
t100 = -t111 * t131 + t114 * t113;
t129 = t111 * t113;
t132 = t109 * t114;
t98 = t116 * t132 + (t107 * t118 + t116 * t129) * t110;
t93 = -t100 * t117 + t98 * t115;
t81 = atan2(-t82, t93);
t78 = sin(t81);
t79 = cos(t81);
t72 = -t78 * t82 + t79 * t93;
t71 = 0.1e1 / t72 ^ 2;
t119 = t108 * t130 + t123 * t109;
t102 = -t107 * t133 + t112 * t111;
t90 = t102 * t118 - t139 * t116;
t85 = t90 * t115 - t119 * t117;
t138 = t71 * t85;
t106 = pkin(13) + qJ(6);
t105 = cos(t106);
t104 = sin(t106);
t89 = t102 * t116 + t139 * t118;
t135 = t89 * t104;
t86 = t119 * t115 + t90 * t117;
t77 = t86 * t105 + t135;
t75 = 0.1e1 / t77 ^ 2;
t134 = t89 * t105;
t76 = t86 * t104 - t134;
t137 = t75 * t76;
t92 = 0.1e1 / t93 ^ 2;
t136 = t82 * t92;
t127 = t76 ^ 2 * t75 + 0.1e1;
t125 = -t78 * t93 - t79 * t82;
t97 = t118 * t132 + (-t107 * t116 + t118 * t129) * t110;
t94 = t100 * t115 + t98 * t117;
t91 = 0.1e1 / t93;
t87 = -t101 * t116 + t120 * t118;
t84 = t121 * t115 + t88 * t117;
t80 = 0.1e1 / (t82 ^ 2 * t92 + 0.1e1);
t74 = 0.1e1 / t77;
t73 = 0.1e1 / t127;
t70 = 0.1e1 / t72;
t69 = 0.1e1 / (t85 ^ 2 * t71 + 0.1e1);
t68 = (t97 * t136 - t87 * t91) * t80 * t115;
t67 = (t94 * t136 - t84 * t91) * t80;
t1 = [0, 0, t68, t67, 0, 0; 0, 0 (-t89 * t115 * t70 - (t125 * t68 + (-t78 * t87 + t79 * t97) * t115) * t138) * t69 (t86 * t70 - (t125 * t67 - t78 * t84 + t79 * t94) * t138) * t69, 0, 0; 0, 0 ((-t90 * t105 - t117 * t135) * t74 - (t90 * t104 - t117 * t134) * t137) * t73 (-t104 * t74 + t105 * t137) * t85 * t73, 0, t127 * t73;];
Ja_rot  = t1;
