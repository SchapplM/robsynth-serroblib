% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRP14_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:53:24
% EndTime: 2019-02-26 21:53:25
% DurationCPUTime: 0.32s
% Computational Cost: add. (905->44), mult. (2519->109), div. (107->9), fcn. (3557->13), ass. (0->63)
t107 = sin(qJ(5));
t111 = cos(qJ(5));
t106 = cos(pkin(6));
t109 = sin(qJ(2));
t114 = cos(qJ(1));
t119 = t114 * t109;
t110 = sin(qJ(1));
t113 = cos(qJ(2));
t120 = t110 * t113;
t116 = t106 * t119 + t120;
t118 = t114 * t113;
t121 = t110 * t109;
t101 = -t106 * t118 + t121;
t108 = sin(qJ(4));
t112 = cos(qJ(4));
t105 = sin(pkin(6));
t123 = t105 * t114;
t94 = -t101 * t108 + t112 * t123;
t136 = t94 * t107 + t116 * t111;
t115 = t116 * t107;
t135 = t94 * t111 - t115;
t124 = t105 * t113;
t100 = t106 * t112 - t108 * t124;
t89 = -t105 * t109 * t111 + t100 * t107;
t77 = atan2(t136, t89);
t73 = sin(t77);
t74 = cos(t77);
t72 = t136 * t73 + t74 * t89;
t71 = 0.1e1 / t72 ^ 2;
t103 = -t106 * t121 + t118;
t126 = t103 * t111;
t102 = t106 * t120 + t119;
t125 = t105 * t110;
t92 = t102 * t108 + t112 * t125;
t80 = t92 * t107 - t126;
t133 = t71 * t80;
t132 = t74 * t136;
t127 = t103 * t107;
t81 = t92 * t111 + t127;
t79 = 0.1e1 / t81 ^ 2;
t91 = t102 * t112 - t108 * t125;
t131 = t79 * t91;
t130 = t80 ^ 2 * t71;
t87 = 0.1e1 / t89 ^ 2;
t129 = t136 * t87;
t128 = t91 ^ 2 * t79;
t122 = t107 * t109;
t117 = -t73 * t89 + t132;
t93 = t101 * t112 + t108 * t123;
t99 = -t106 * t108 - t112 * t124;
t96 = (t108 * t122 - t111 * t113) * t105;
t90 = t100 * t111 + t105 * t122;
t86 = 0.1e1 / t89;
t85 = t101 * t111 + t108 * t115;
t78 = 0.1e1 / t81;
t76 = 0.1e1 / (0.1e1 + t128);
t75 = 0.1e1 / (t136 ^ 2 * t87 + 0.1e1);
t70 = 0.1e1 / t72;
t69 = 0.1e1 / (0.1e1 + t130);
t68 = (-t99 * t129 - t86 * t93) * t75 * t107;
t67 = (-t96 * t129 - t85 * t86) * t75;
t66 = (-t90 * t129 + t135 * t86) * t75;
t1 = [-t80 * t86 * t75, t67, 0, t68, t66, 0; (t136 * t70 - (-t73 + (-t86 * t132 + t73) * t75) * t130) * t69 ((t102 * t111 + t108 * t127) * t70 - (t117 * t67 - t73 * t85 + t74 * t96) * t133) * t69, 0 (t91 * t107 * t70 - (t117 * t68 + (-t73 * t93 + t74 * t99) * t107) * t133) * t69 (t81 * t70 - (t117 * t66 + t135 * t73 + t74 * t90) * t133) * t69, 0; (-t135 * t131 - t93 * t78) * t76 (t103 * t112 * t78 - (-t102 * t107 + t108 * t126) * t131) * t76, 0 (-t111 * t128 - t78 * t92) * t76, t80 * t76 * t131, 0;];
Ja_rot  = t1;
