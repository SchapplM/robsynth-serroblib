% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PPRRRR3_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobia_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:43:55
% EndTime: 2019-02-26 19:43:56
% DurationCPUTime: 0.36s
% Computational Cost: add. (1774->52), mult. (5219->122), div. (65->9), fcn. (7059->19), ass. (0->74)
t108 = sin(pkin(8));
t109 = sin(pkin(7));
t113 = cos(pkin(8));
t114 = cos(pkin(7));
t118 = sin(qJ(3));
t121 = cos(qJ(3));
t106 = sin(pkin(14));
t107 = sin(pkin(13));
t111 = cos(pkin(14));
t112 = cos(pkin(13));
t115 = cos(pkin(6));
t133 = t112 * t115;
t127 = -t107 * t106 + t111 * t133;
t110 = sin(pkin(6));
t134 = t112 * t110;
t125 = -t109 * t134 + t127 * t114;
t126 = t106 * t133 + t107 * t111;
t123 = -t126 * t118 + t125 * t121;
t145 = t123 * t113 + (-t127 * t109 - t114 * t134) * t108;
t117 = sin(qJ(4));
t120 = cos(qJ(4));
t95 = t125 * t118 + t126 * t121;
t80 = t95 * t117 - t145 * t120;
t135 = t111 * t114;
t137 = t109 * t115;
t101 = t118 * t137 + (t106 * t121 + t118 * t135) * t110;
t100 = t121 * t137 + (-t106 * t118 + t121 * t135) * t110;
t136 = t110 * t109;
t129 = t100 * t113 + (-t111 * t136 + t115 * t114) * t108;
t89 = t101 * t117 - t129 * t120;
t79 = atan2(-t80, t89);
t76 = sin(t79);
t77 = cos(t79);
t70 = -t76 * t80 + t77 * t89;
t69 = 0.1e1 / t70 ^ 2;
t132 = t113 * t120;
t138 = t107 * t115;
t104 = -t112 * t106 - t111 * t138;
t102 = t107 * t110 * t114 - t104 * t109;
t139 = t102 * t108;
t105 = -t106 * t138 + t112 * t111;
t128 = t104 * t114 + t107 * t136;
t97 = t105 * t121 + t128 * t118;
t140 = t97 * t117;
t96 = -t105 * t118 + t128 * t121;
t83 = -t120 * t139 - t96 * t132 + t140;
t144 = t69 * t83;
t116 = sin(qJ(5));
t119 = cos(qJ(5));
t84 = t97 * t120 + (t113 * t96 + t139) * t117;
t91 = t102 * t113 - t96 * t108;
t75 = t91 * t116 + t84 * t119;
t73 = 0.1e1 / t75 ^ 2;
t74 = t84 * t116 - t91 * t119;
t143 = t73 * t74;
t88 = 0.1e1 / t89 ^ 2;
t142 = t80 * t88;
t141 = t108 * t97;
t131 = t74 ^ 2 * t73 + 0.1e1;
t130 = -t76 * t89 - t77 * t80;
t92 = t100 * t117 + t101 * t132;
t90 = t101 * t120 + t129 * t117;
t87 = 0.1e1 / t89;
t86 = -t113 * t140 + t96 * t120;
t85 = t123 * t117 + t95 * t132;
t82 = t145 * t117 + t95 * t120;
t78 = 0.1e1 / (t80 ^ 2 * t88 + 0.1e1);
t72 = 0.1e1 / t75;
t71 = 0.1e1 / t131;
t68 = 0.1e1 / t70;
t67 = 0.1e1 / (t83 ^ 2 * t69 + 0.1e1);
t66 = (t92 * t142 - t85 * t87) * t78;
t65 = (t90 * t142 - t82 * t87) * t78;
t1 = [0, 0, t66, t65, 0, 0; 0, 0 ((t96 * t117 + t97 * t132) * t68 - (t130 * t66 - t76 * t85 + t77 * t92) * t144) * t67 (t84 * t68 - (t130 * t65 - t76 * t82 + t77 * t90) * t144) * t67, 0, 0; 0, 0 ((t86 * t116 - t119 * t141) * t72 - (t116 * t141 + t86 * t119) * t143) * t71 (-t116 * t72 + t119 * t143) * t83 * t71, t131 * t71, 0;];
Ja_rot  = t1;
