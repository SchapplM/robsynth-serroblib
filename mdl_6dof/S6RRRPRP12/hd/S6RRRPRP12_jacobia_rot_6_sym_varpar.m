% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP12
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRP12_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:15:31
% EndTime: 2019-02-26 22:15:32
% DurationCPUTime: 0.29s
% Computational Cost: add. (907->45), mult. (2525->110), div. (107->9), fcn. (3566->13), ass. (0->62)
t104 = sin(qJ(5));
t108 = cos(qJ(5));
t105 = sin(qJ(3));
t109 = cos(qJ(3));
t102 = sin(pkin(6));
t111 = cos(qJ(1));
t121 = t102 * t111;
t103 = cos(pkin(6));
t106 = sin(qJ(2));
t116 = t111 * t106;
t107 = sin(qJ(1));
t110 = cos(qJ(2));
t118 = t107 * t110;
t98 = t103 * t116 + t118;
t113 = t98 * t105 + t109 * t121;
t115 = t111 * t110;
t119 = t107 * t106;
t97 = -t103 * t115 + t119;
t80 = t113 * t104 + t97 * t108;
t133 = -t97 * t104 + t113 * t108;
t123 = t102 * t105;
t95 = -t103 * t109 + t106 * t123;
t88 = -t102 * t110 * t104 - t95 * t108;
t75 = atan2(t133, t88);
t71 = sin(t75);
t72 = cos(t75);
t70 = t133 * t71 + t72 * t88;
t69 = 0.1e1 / t70 ^ 2;
t100 = -t103 * t119 + t115;
t122 = t102 * t109;
t112 = t100 * t105 - t107 * t122;
t99 = t103 * t118 + t116;
t124 = t99 * t104;
t81 = -t112 * t108 + t124;
t132 = t69 * t81;
t131 = t72 * t133;
t82 = t112 * t104 + t99 * t108;
t77 = 0.1e1 / t82 ^ 2;
t92 = t100 * t109 + t107 * t123;
t130 = t77 * t92;
t87 = 0.1e1 / t88 ^ 2;
t129 = t133 * t87;
t128 = t81 ^ 2 * t69;
t127 = t92 ^ 2 * t77;
t120 = t105 * t108;
t117 = t108 * t110;
t114 = -t71 * t88 + t131;
t90 = -t105 * t121 + t98 * t109;
t96 = t103 * t105 + t106 * t122;
t94 = (t104 * t106 - t105 * t117) * t102;
t89 = -t102 * t117 + t95 * t104;
t86 = 0.1e1 / t88;
t83 = t98 * t104 + t97 * t120;
t76 = 0.1e1 / t82;
t74 = 0.1e1 / (0.1e1 + t127);
t73 = 0.1e1 / (t133 ^ 2 * t87 + 0.1e1);
t68 = 0.1e1 / t70;
t67 = 0.1e1 / (0.1e1 + t128);
t66 = (t96 * t129 + t86 * t90) * t73 * t108;
t65 = (-t94 * t129 - t83 * t86) * t73;
t64 = (-t89 * t129 - t80 * t86) * t73;
t1 = [-t81 * t86 * t73, t65, t66, 0, t64, 0; (t133 * t68 - (-t71 + (-t86 * t131 + t71) * t73) * t128) * t67 ((t100 * t104 + t99 * t120) * t68 - (t114 * t65 - t71 * t83 + t72 * t94) * t132) * t67 (-t92 * t108 * t68 - (t114 * t66 + (t71 * t90 - t72 * t96) * t108) * t132) * t67, 0 (t82 * t68 - (t114 * t64 - t71 * t80 + t72 * t89) * t132) * t67, 0; (-t80 * t130 + t90 * t76) * t74 (t99 * t109 * t76 + (t100 * t108 - t105 * t124) * t130) * t74 (t104 * t127 + t112 * t76) * t74, 0, -t81 * t74 * t130, 0;];
Ja_rot  = t1;
