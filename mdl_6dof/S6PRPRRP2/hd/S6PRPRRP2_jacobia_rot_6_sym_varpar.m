% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRP2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:51:04
% EndTime: 2019-02-26 19:51:04
% DurationCPUTime: 0.24s
% Computational Cost: add. (1374->42), mult. (3740->105), div. (87->9), fcn. (5214->15), ass. (0->62)
t112 = sin(qJ(5));
t115 = cos(qJ(5));
t107 = sin(pkin(10));
t110 = cos(pkin(10));
t106 = sin(pkin(11));
t109 = cos(pkin(11));
t114 = sin(qJ(2));
t117 = cos(qJ(2));
t104 = t114 * t106 - t117 * t109;
t111 = cos(pkin(6));
t119 = t104 * t111;
t120 = t117 * t106 + t114 * t109;
t118 = -t107 * t120 - t110 * t119;
t116 = cos(qJ(4));
t108 = sin(pkin(6));
t113 = sin(qJ(4));
t125 = t108 * t113;
t103 = t120 * t111;
t94 = t110 * t103 - t107 * t104;
t89 = -t110 * t125 + t94 * t116;
t77 = t89 * t112 + t118 * t115;
t101 = t104 * t108;
t102 = t120 * t108;
t99 = t102 * t116 + t111 * t113;
t85 = -t101 * t115 + t99 * t112;
t73 = atan2(-t77, t85);
t70 = sin(t73);
t71 = cos(t73);
t69 = -t70 * t77 + t71 * t85;
t68 = 0.1e1 / t69 ^ 2;
t96 = t107 * t119 - t110 * t120;
t126 = t96 * t115;
t121 = -t107 * t103 - t110 * t104;
t91 = t107 * t125 + t116 * t121;
t80 = t91 * t112 + t126;
t130 = t68 * t80;
t81 = -t96 * t112 + t91 * t115;
t76 = 0.1e1 / t81 ^ 2;
t124 = t108 * t116;
t90 = t107 * t124 - t113 * t121;
t129 = t76 * t90;
t84 = 0.1e1 / t85 ^ 2;
t128 = t77 * t84;
t127 = t90 ^ 2 * t76;
t123 = t112 * t116;
t122 = -t70 * t85 - t71 * t77;
t98 = -t102 * t113 + t111 * t116;
t92 = -t101 * t123 - t102 * t115;
t88 = -t110 * t124 - t94 * t113;
t86 = t101 * t112 + t99 * t115;
t83 = 0.1e1 / t85;
t82 = -t94 * t115 + t118 * t123;
t79 = -t118 * t112 + t89 * t115;
t75 = 0.1e1 / t81;
t74 = 0.1e1 / (0.1e1 + t127);
t72 = 0.1e1 / (t77 ^ 2 * t84 + 0.1e1);
t67 = 0.1e1 / t69;
t66 = 0.1e1 / (t80 ^ 2 * t68 + 0.1e1);
t65 = (t98 * t128 - t83 * t88) * t72 * t112;
t64 = (t92 * t128 - t82 * t83) * t72;
t63 = (t86 * t128 - t79 * t83) * t72;
t1 = [0, t64, 0, t65, t63, 0; 0 ((-t115 * t121 + t96 * t123) * t67 - (t122 * t64 - t70 * t82 + t71 * t92) * t130) * t66, 0 (t90 * t112 * t67 - (t122 * t65 + (-t70 * t88 + t71 * t98) * t112) * t130) * t66 (t81 * t67 - (t122 * t63 - t70 * t79 + t71 * t86) * t130) * t66, 0; 0 (-t96 * t113 * t75 - (t112 * t121 + t116 * t126) * t129) * t74, 0 (-t115 * t127 - t75 * t91) * t74, t80 * t74 * t129, 0;];
Ja_rot  = t1;
