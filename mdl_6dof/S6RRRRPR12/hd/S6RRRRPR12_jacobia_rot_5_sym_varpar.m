% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR12
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
% Datum: 2019-02-26 22:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR12_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_jacobia_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:36:58
% EndTime: 2019-02-26 22:36:59
% DurationCPUTime: 0.29s
% Computational Cost: add. (899->52), mult. (2499->114), div. (85->9), fcn. (3431->15), ass. (0->69)
t144 = sin(qJ(1));
t114 = cos(pkin(6));
t116 = sin(qJ(2));
t125 = t144 * t116;
t118 = cos(qJ(2));
t119 = cos(qJ(1));
t129 = t119 * t118;
t102 = -t114 * t129 + t125;
t117 = cos(qJ(3));
t111 = sin(pkin(7));
t112 = sin(pkin(6));
t135 = t112 * t119;
t127 = t111 * t135;
t113 = cos(pkin(7));
t133 = t113 * t117;
t124 = t144 * t118;
t130 = t119 * t116;
t103 = t114 * t130 + t124;
t115 = sin(qJ(3));
t138 = t103 * t115;
t85 = t102 * t133 + t117 * t127 + t138;
t136 = t111 * t114;
t95 = -t117 * t136 + (t115 * t116 - t118 * t133) * t112;
t84 = atan2(-t85, t95);
t81 = sin(t84);
t82 = cos(t84);
t75 = -t81 * t85 + t82 * t95;
t74 = 0.1e1 / t75 ^ 2;
t104 = -t114 * t125 + t129;
t121 = t114 * t124 + t130;
t120 = t121 * t117;
t126 = t112 * t144;
t123 = t111 * t126;
t89 = t104 * t115 + t113 * t120 - t117 * t123;
t143 = t74 * t89;
t110 = qJ(4) + pkin(13);
t108 = sin(t110);
t109 = cos(t110);
t90 = t104 * t117 + (-t113 * t121 + t123) * t115;
t98 = t111 * t121 + t113 * t126;
t80 = t98 * t108 + t90 * t109;
t78 = 0.1e1 / t80 ^ 2;
t79 = t90 * t108 - t98 * t109;
t142 = t78 * t79;
t141 = t82 * t85;
t94 = 0.1e1 / t95 ^ 2;
t140 = t85 * t94;
t139 = t89 ^ 2 * t74;
t137 = t104 * t111;
t134 = t113 * t115;
t132 = t115 * t118;
t131 = t116 * t117;
t128 = t79 ^ 2 * t78 + 0.1e1;
t122 = -t81 * t95 - t141;
t88 = t102 * t134 - t103 * t117 + t115 * t127;
t101 = (t113 * t131 + t132) * t112;
t97 = -t102 * t111 + t113 * t135;
t96 = t115 * t136 + (t113 * t132 + t131) * t112;
t93 = 0.1e1 / t95;
t92 = -t104 * t134 - t120;
t91 = -t102 * t115 + t103 * t133;
t83 = 0.1e1 / (t85 ^ 2 * t94 + 0.1e1);
t77 = 0.1e1 / t80;
t76 = 0.1e1 / t128;
t73 = 0.1e1 / t75;
t72 = 0.1e1 / (0.1e1 + t139);
t71 = (t101 * t140 - t91 * t93) * t83;
t70 = (t140 * t96 + t88 * t93) * t83;
t1 = [-t89 * t93 * t83, t71, t70, 0, 0, 0; ((-t138 + (-t102 * t113 - t127) * t117) * t73 - (-t81 + (t141 * t93 + t81) * t83) * t139) * t72 ((t104 * t133 - t115 * t121) * t73 - (t82 * t101 + t122 * t71 - t81 * t91) * t143) * t72 (t90 * t73 - (t122 * t70 + t81 * t88 + t82 * t96) * t143) * t72, 0, 0, 0; ((t88 * t108 - t97 * t109) * t77 - (t97 * t108 + t88 * t109) * t142) * t76 ((t92 * t108 - t109 * t137) * t77 - (t108 * t137 + t92 * t109) * t142) * t76 (-t108 * t77 + t109 * t142) * t89 * t76, t128 * t76, 0, 0;];
Ja_rot  = t1;
