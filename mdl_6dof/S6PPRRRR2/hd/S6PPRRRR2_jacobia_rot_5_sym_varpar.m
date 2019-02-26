% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Ja_rot = S6PPRRRR2_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_jacobia_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:43:16
% EndTime: 2019-02-26 19:43:16
% DurationCPUTime: 0.25s
% Computational Cost: add. (1048->42), mult. (3005->102), div. (65->9), fcn. (4119->17), ass. (0->66)
t101 = cos(pkin(7));
t100 = cos(pkin(12));
t102 = cos(pkin(6));
t96 = sin(pkin(12));
t123 = t102 * t96;
t95 = sin(pkin(13));
t99 = cos(pkin(13));
t115 = t100 * t95 + t99 * t123;
t97 = sin(pkin(7));
t98 = sin(pkin(6));
t124 = t98 * t97;
t129 = t115 * t101 - t96 * t124;
t104 = sin(qJ(4));
t107 = cos(qJ(4));
t119 = t100 * t102;
t114 = t99 * t119 - t96 * t95;
t120 = t98 * t101;
t109 = -t100 * t120 - t114 * t97;
t105 = sin(qJ(3));
t108 = cos(qJ(3));
t110 = -t100 * t124 + t114 * t101;
t92 = t95 * t119 + t96 * t99;
t79 = t110 * t105 + t92 * t108;
t73 = t79 * t104 - t109 * t107;
t113 = t102 * t97 + t99 * t120;
t125 = t95 * t98;
t89 = t113 * t105 + t108 * t125;
t91 = t102 * t101 - t99 * t124;
t84 = t89 * t104 - t91 * t107;
t72 = atan2(-t73, t84);
t69 = sin(t72);
t70 = cos(t72);
t63 = -t69 * t73 + t70 * t84;
t62 = 0.1e1 / t63 ^ 2;
t111 = t115 * t97 + t96 * t120;
t93 = t100 * t99 - t95 * t123;
t81 = -t129 * t105 + t93 * t108;
t76 = t81 * t104 - t111 * t107;
t128 = t62 * t76;
t106 = cos(qJ(5));
t103 = sin(qJ(5));
t80 = t93 * t105 + t129 * t108;
t122 = t80 * t103;
t77 = t111 * t104 + t81 * t107;
t68 = t77 * t106 + t122;
t66 = 0.1e1 / t68 ^ 2;
t121 = t80 * t106;
t67 = t77 * t103 - t121;
t127 = t66 * t67;
t83 = 0.1e1 / t84 ^ 2;
t126 = t73 * t83;
t117 = t67 ^ 2 * t66 + 0.1e1;
t116 = -t69 * t84 - t70 * t73;
t88 = -t105 * t125 + t113 * t108;
t85 = t91 * t104 + t89 * t107;
t82 = 0.1e1 / t84;
t78 = -t92 * t105 + t110 * t108;
t75 = t109 * t104 + t79 * t107;
t71 = 0.1e1 / (t73 ^ 2 * t83 + 0.1e1);
t65 = 0.1e1 / t68;
t64 = 0.1e1 / t117;
t61 = 0.1e1 / t63;
t60 = 0.1e1 / (t76 ^ 2 * t62 + 0.1e1);
t59 = (t88 * t126 - t78 * t82) * t71 * t104;
t58 = (t85 * t126 - t75 * t82) * t71;
t1 = [0, 0, t59, t58, 0, 0; 0, 0 (-t80 * t104 * t61 - (t116 * t59 + (-t69 * t78 + t70 * t88) * t104) * t128) * t60 (t77 * t61 - (t116 * t58 - t69 * t75 + t70 * t85) * t128) * t60, 0, 0; 0, 0 ((-t81 * t106 - t107 * t122) * t65 - (t81 * t103 - t107 * t121) * t127) * t64 (-t103 * t65 + t106 * t127) * t76 * t64, t117 * t64, 0;];
Ja_rot  = t1;
