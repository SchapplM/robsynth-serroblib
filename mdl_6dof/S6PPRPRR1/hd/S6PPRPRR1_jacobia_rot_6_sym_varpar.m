% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PPRPRR1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_jacobia_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:39:48
% EndTime: 2019-02-26 19:39:48
% DurationCPUTime: 0.33s
% Computational Cost: add. (1594->47), mult. (4461->113), div. (65->9), fcn. (6121->19), ass. (0->71)
t116 = sin(qJ(5));
t119 = cos(qJ(5));
t108 = sin(pkin(7));
t112 = cos(pkin(11));
t106 = sin(pkin(12));
t107 = sin(pkin(11));
t111 = cos(pkin(12));
t114 = cos(pkin(6));
t128 = t112 * t114;
t123 = -t107 * t106 + t111 * t128;
t109 = sin(pkin(6));
t113 = cos(pkin(7));
t129 = t109 * t113;
t121 = -t123 * t108 - t112 * t129;
t105 = sin(pkin(13));
t110 = cos(pkin(13));
t117 = sin(qJ(3));
t120 = cos(qJ(3));
t102 = t117 * t105 - t120 * t110;
t130 = t109 * t112;
t125 = t120 * t105 + t117 * t110;
t95 = t125 * t108;
t97 = t125 * t113;
t99 = t106 * t128 + t107 * t111;
t82 = -t99 * t102 + t123 * t97 - t95 * t130;
t76 = t82 * t116 - t121 * t119;
t91 = t114 * t95 + (-t102 * t106 + t111 * t97) * t109;
t98 = -t109 * t111 * t108 + t114 * t113;
t88 = t91 * t116 - t98 * t119;
t75 = atan2(-t76, t88);
t72 = sin(t75);
t73 = cos(t75);
t66 = -t72 * t76 + t73 * t88;
t65 = 0.1e1 / t66 ^ 2;
t131 = t107 * t114;
t100 = -t112 * t106 - t111 * t131;
t101 = -t106 * t131 + t112 * t111;
t132 = t107 * t109;
t122 = t100 * t97 - t101 * t102 + t95 * t132;
t124 = -t100 * t108 + t107 * t129;
t79 = t116 * t122 - t124 * t119;
t137 = t65 * t79;
t118 = cos(qJ(6));
t115 = sin(qJ(6));
t94 = t102 * t108;
t96 = t102 * t113;
t84 = -t100 * t96 - t101 * t125 - t94 * t132;
t134 = t84 * t115;
t80 = t124 * t116 + t119 * t122;
t71 = t80 * t118 - t134;
t69 = 0.1e1 / t71 ^ 2;
t133 = t84 * t118;
t70 = t80 * t115 + t133;
t136 = t69 * t70;
t87 = 0.1e1 / t88 ^ 2;
t135 = t76 * t87;
t127 = t70 ^ 2 * t69 + 0.1e1;
t126 = -t72 * t88 - t73 * t76;
t90 = -t114 * t94 + (-t106 * t125 - t111 * t96) * t109;
t89 = t98 * t116 + t91 * t119;
t86 = 0.1e1 / t88;
t81 = -t123 * t96 - t125 * t99 + t94 * t130;
t78 = t121 * t116 + t82 * t119;
t74 = 0.1e1 / (t76 ^ 2 * t87 + 0.1e1);
t68 = 0.1e1 / t71;
t67 = 0.1e1 / t127;
t64 = 0.1e1 / t66;
t63 = 0.1e1 / (t79 ^ 2 * t65 + 0.1e1);
t62 = (t90 * t135 - t81 * t86) * t74 * t116;
t61 = (t89 * t135 - t78 * t86) * t74;
t1 = [0, 0, t62, 0, t61, 0; 0, 0 (t84 * t116 * t64 - (t126 * t62 + (-t72 * t81 + t73 * t90) * t116) * t137) * t63, 0 (t80 * t64 - (t126 * t61 - t72 * t78 + t73 * t89) * t137) * t63, 0; 0, 0 ((-t118 * t122 + t119 * t134) * t68 - (t115 * t122 + t119 * t133) * t136) * t67, 0 (-t115 * t68 + t118 * t136) * t79 * t67, t127 * t67;];
Ja_rot  = t1;
