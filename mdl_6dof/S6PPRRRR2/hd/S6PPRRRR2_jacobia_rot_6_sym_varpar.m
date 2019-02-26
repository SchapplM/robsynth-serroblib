% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PPRRRR2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_jacobia_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:43:16
% EndTime: 2019-02-26 19:43:16
% DurationCPUTime: 0.29s
% Computational Cost: add. (1197->44), mult. (3249->105), div. (70->9), fcn. (4448->17), ass. (0->68)
t114 = sin(pkin(12));
t119 = cos(pkin(7));
t113 = sin(pkin(13));
t117 = cos(pkin(13));
t118 = cos(pkin(12));
t120 = cos(pkin(6));
t139 = t114 * t120;
t129 = t118 * t113 + t117 * t139;
t115 = sin(pkin(7));
t116 = sin(pkin(6));
t137 = t116 * t115;
t145 = -t114 * t137 + t129 * t119;
t121 = sin(qJ(4));
t123 = cos(qJ(4));
t134 = t118 * t120;
t130 = -t114 * t113 + t117 * t134;
t136 = t116 * t119;
t127 = -t130 * t115 - t118 * t136;
t107 = t113 * t134 + t114 * t117;
t122 = sin(qJ(3));
t124 = cos(qJ(3));
t126 = -t118 * t137 + t130 * t119;
t94 = t107 * t124 + t126 * t122;
t88 = t94 * t121 - t127 * t123;
t135 = t117 * t119;
t138 = t115 * t120;
t104 = t122 * t138 + (t113 * t124 + t122 * t135) * t116;
t106 = -t117 * t137 + t120 * t119;
t99 = t104 * t121 - t106 * t123;
t87 = atan2(-t88, t99);
t84 = sin(t87);
t85 = cos(t87);
t78 = -t84 * t88 + t85 * t99;
t77 = 0.1e1 / t78 ^ 2;
t125 = t114 * t136 + t129 * t115;
t108 = -t113 * t139 + t118 * t117;
t96 = t108 * t124 - t145 * t122;
t91 = t96 * t121 - t125 * t123;
t144 = t77 * t91;
t112 = qJ(5) + qJ(6);
t111 = cos(t112);
t110 = sin(t112);
t95 = t108 * t122 + t145 * t124;
t141 = t95 * t110;
t92 = t125 * t121 + t96 * t123;
t83 = t92 * t111 + t141;
t81 = 0.1e1 / t83 ^ 2;
t140 = t95 * t111;
t82 = t92 * t110 - t140;
t143 = t81 * t82;
t98 = 0.1e1 / t99 ^ 2;
t142 = t88 * t98;
t133 = t82 ^ 2 * t81 + 0.1e1;
t131 = -t84 * t99 - t85 * t88;
t103 = t124 * t138 + (-t113 * t122 + t124 * t135) * t116;
t100 = t104 * t123 + t106 * t121;
t97 = 0.1e1 / t99;
t93 = -t107 * t122 + t126 * t124;
t90 = t127 * t121 + t94 * t123;
t86 = 0.1e1 / (t88 ^ 2 * t98 + 0.1e1);
t80 = 0.1e1 / t83;
t79 = 0.1e1 / t133;
t76 = 0.1e1 / t78;
t75 = 0.1e1 / (t91 ^ 2 * t77 + 0.1e1);
t74 = (t103 * t142 - t93 * t97) * t86 * t121;
t73 = (t100 * t142 - t90 * t97) * t86;
t72 = t133 * t79;
t1 = [0, 0, t74, t73, 0, 0; 0, 0 (-t95 * t121 * t76 - (t131 * t74 + (t103 * t85 - t84 * t93) * t121) * t144) * t75 (t92 * t76 - (t85 * t100 + t131 * t73 - t84 * t90) * t144) * t75, 0, 0; 0, 0 ((-t96 * t111 - t123 * t141) * t80 - (t96 * t110 - t123 * t140) * t143) * t79 (-t110 * t80 + t111 * t143) * t91 * t79, t72, t72;];
Ja_rot  = t1;
