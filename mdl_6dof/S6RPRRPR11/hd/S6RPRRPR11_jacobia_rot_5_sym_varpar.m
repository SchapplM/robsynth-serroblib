% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR11_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_jacobia_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:06:39
% EndTime: 2019-02-26 21:06:40
% DurationCPUTime: 0.46s
% Computational Cost: add. (1353->49), mult. (3877->116), div. (80->9), fcn. (5331->17), ass. (0->69)
t114 = sin(qJ(4));
t116 = cos(qJ(4));
t113 = cos(pkin(6));
t111 = cos(pkin(12));
t144 = sin(qJ(1));
t126 = t144 * t111;
t107 = sin(pkin(12));
t118 = cos(qJ(1));
t130 = t118 * t107;
t103 = t113 * t130 + t126;
t115 = sin(qJ(3));
t117 = cos(qJ(3));
t127 = t144 * t107;
t129 = t118 * t111;
t102 = -t113 * t129 + t127;
t108 = sin(pkin(7));
t112 = cos(pkin(7));
t109 = sin(pkin(6));
t132 = t109 * t118;
t123 = t102 * t112 + t108 * t132;
t92 = -t103 * t117 + t123 * t115;
t99 = -t102 * t108 + t112 * t132;
t146 = t92 * t114 - t99 * t116;
t82 = t99 * t114 + t92 * t116;
t122 = t113 * t126 + t130;
t128 = t109 * t144;
t145 = -t108 * t128 + t122 * t112;
t101 = -t109 * t111 * t108 + t113 * t112;
t131 = t111 * t112;
t133 = t108 * t113;
t98 = t115 * t133 + (t107 * t117 + t115 * t131) * t109;
t87 = -t101 * t116 + t98 * t114;
t78 = atan2(t146, t87);
t75 = sin(t78);
t76 = cos(t78);
t69 = t146 * t75 + t76 * t87;
t68 = 0.1e1 / t69 ^ 2;
t119 = t122 * t108 + t112 * t128;
t104 = -t113 * t127 + t129;
t94 = t104 * t117 - t145 * t115;
t83 = t94 * t114 - t119 * t116;
t143 = t68 * t83;
t110 = cos(pkin(13));
t106 = sin(pkin(13));
t93 = t104 * t115 + t145 * t117;
t138 = t93 * t106;
t84 = t119 * t114 + t94 * t116;
t74 = t84 * t110 + t138;
t72 = 0.1e1 / t74 ^ 2;
t137 = t93 * t110;
t73 = t84 * t106 - t137;
t142 = t72 * t73;
t141 = t76 * t146;
t86 = 0.1e1 / t87 ^ 2;
t140 = t146 * t86;
t139 = t83 ^ 2 * t68;
t124 = -t75 * t87 + t141;
t120 = -t103 * t115 - t123 * t117;
t97 = t117 * t133 + (-t107 * t115 + t117 * t131) * t109;
t88 = t101 * t114 + t98 * t116;
t85 = 0.1e1 / t87;
t77 = 0.1e1 / (t146 ^ 2 * t86 + 0.1e1);
t71 = 0.1e1 / t74;
t70 = 0.1e1 / (t73 ^ 2 * t72 + 0.1e1);
t67 = 0.1e1 / t69;
t66 = 0.1e1 / (0.1e1 + t139);
t65 = (-t120 * t85 - t97 * t140) * t77 * t114;
t64 = (-t88 * t140 + t82 * t85) * t77;
t1 = [-t83 * t85 * t77, 0, t65, t64, 0, 0; (t146 * t67 - (-t75 + (-t85 * t141 + t75) * t77) * t139) * t66, 0 (-t93 * t114 * t67 - (t124 * t65 + (-t120 * t75 + t76 * t97) * t114) * t143) * t66 (t84 * t67 - (t124 * t64 + t75 * t82 + t76 * t88) * t143) * t66, 0, 0; ((t82 * t106 - t110 * t120) * t71 - (t106 * t120 + t82 * t110) * t142) * t70, 0 ((-t94 * t110 - t116 * t138) * t71 - (t94 * t106 - t116 * t137) * t142) * t70 (-t106 * t71 + t110 * t142) * t83 * t70, 0, 0;];
Ja_rot  = t1;
