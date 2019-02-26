% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPP1_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:02:47
% EndTime: 2019-02-26 22:02:48
% DurationCPUTime: 0.24s
% Computational Cost: add. (665->45), mult. (2053->108), div. (80->9), fcn. (2819->13), ass. (0->60)
t105 = sin(pkin(6));
t107 = cos(pkin(6));
t112 = cos(qJ(2));
t108 = sin(qJ(3));
t109 = sin(qJ(2));
t126 = t108 * t109;
t114 = t105 * t112 + t107 * t126;
t104 = sin(pkin(10));
t106 = cos(pkin(10));
t111 = cos(qJ(3));
t124 = t109 * t111;
t87 = t104 * t124 + t106 * t114;
t110 = sin(qJ(1));
t130 = t105 * t109;
t117 = t106 * t130;
t128 = t106 * t107;
t113 = cos(qJ(1));
t119 = t113 * t111;
t122 = t110 * t108;
t95 = t112 * t122 + t119;
t120 = t113 * t108;
t121 = t111 * t112;
t96 = t110 * t121 - t120;
t79 = t104 * t96 - t110 * t117 + t128 * t95;
t77 = atan2(-t79, t87);
t74 = sin(t77);
t75 = cos(t77);
t73 = -t74 * t79 + t75 * t87;
t72 = 0.1e1 / t73 ^ 2;
t98 = t112 * t119 + t122;
t131 = t98 * t104;
t97 = t110 * t111 - t112 * t120;
t81 = -t113 * t117 - t128 * t97 + t131;
t133 = t81 ^ 2 * t72;
t70 = 0.1e1 / (0.1e1 + t133);
t71 = 0.1e1 / t73;
t137 = t70 * t71;
t136 = t72 * t81;
t135 = t75 * t79;
t86 = 0.1e1 / t87 ^ 2;
t134 = t79 * t86;
t123 = t109 * t113;
t82 = t98 * t106 + (t105 * t123 + t107 * t97) * t104;
t91 = -t105 * t97 + t107 * t123;
t90 = 0.1e1 / t91 ^ 2;
t132 = t82 * t90;
t127 = t107 * t112;
t125 = t109 * t110;
t115 = -t74 * t87 - t135;
t94 = (-t104 * t108 + t111 * t128) * t109;
t89 = 0.1e1 / t91;
t88 = t104 * t121 + (t108 * t127 - t130) * t106;
t85 = 0.1e1 / t87;
t84 = t87 * t110;
t83 = -t104 * t95 + t128 * t96;
t78 = 0.1e1 / (t82 ^ 2 * t90 + 0.1e1);
t76 = 0.1e1 / (t79 ^ 2 * t86 + 0.1e1);
t69 = (t134 * t94 - t83 * t85) * t76;
t68 = (t134 * t88 + t84 * t85) * t76;
t1 = [-t81 * t85 * t76, t68, t69, 0, 0, 0; -t79 * t137 - (-t74 + (t135 * t85 + t74) * t76) * t70 * t133 -(t115 * t68 + t74 * t84 + t75 * t88) * t70 * t136 - t87 * t113 * t137 ((t104 * t97 + t128 * t98) * t71 - (t115 * t69 - t74 * t83 + t75 * t94) * t136) * t70, 0, 0, 0; ((-t96 * t106 + (-t105 * t125 + t107 * t95) * t104) * t89 - (-t105 * t95 - t107 * t125) * t132) * t78 ((t104 * t114 - t106 * t124) * t89 - (-t105 * t126 + t127) * t132) * t78 * t113 ((t106 * t97 - t107 * t131) * t89 - t98 * t105 * t132) * t78, 0, 0, 0;];
Ja_rot  = t1;
