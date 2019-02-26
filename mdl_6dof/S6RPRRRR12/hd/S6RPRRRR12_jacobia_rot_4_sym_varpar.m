% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRR12_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobia_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:08
% EndTime: 2019-02-26 21:21:08
% DurationCPUTime: 0.27s
% Computational Cost: add. (950->43), mult. (2744->104), div. (55->9), fcn. (3759->17), ass. (0->63)
t102 = cos(pkin(8));
t103 = cos(pkin(7));
t106 = sin(qJ(3));
t109 = cos(qJ(3));
t100 = sin(pkin(6));
t110 = cos(qJ(1));
t122 = t100 * t110;
t99 = sin(pkin(7));
t118 = t99 * t122;
t104 = cos(pkin(6));
t101 = cos(pkin(14));
t119 = t110 * t101;
t107 = sin(qJ(1));
t97 = sin(pkin(14));
t129 = t107 * t97;
t92 = -t104 * t119 + t129;
t120 = t107 * t101;
t128 = t110 * t97;
t93 = t104 * t128 + t120;
t84 = (t103 * t92 + t118) * t109 + t93 * t106;
t91 = t103 * t122 - t92 * t99;
t98 = sin(pkin(8));
t77 = t91 * t102 - t84 * t98;
t105 = sin(qJ(4));
t123 = t100 * t107;
t94 = -t104 * t120 - t128;
t114 = t103 * t123 - t94 * t99;
t115 = t103 * t94 + t99 * t123;
t95 = -t104 * t129 + t119;
t86 = -t95 * t106 + t115 * t109;
t111 = t86 * t102 + t114 * t98;
t108 = cos(qJ(4));
t87 = t115 * t106 + t95 * t109;
t126 = t87 * t108;
t72 = t111 * t105 + t126;
t70 = 0.1e1 / t72 ^ 2;
t127 = t87 * t105;
t71 = -t111 * t108 + t127;
t133 = t70 * t71;
t130 = t104 * t99;
t82 = -(t109 * t130 + (t101 * t103 * t109 - t106 * t97) * t100) * t98 + (-t100 * t101 * t99 + t104 * t103) * t102;
t76 = atan2(t77, t82);
t74 = cos(t76);
t132 = t74 * t77;
t73 = sin(t76);
t67 = t73 * t77 + t74 * t82;
t66 = 0.1e1 / t67 ^ 2;
t78 = -t114 * t102 + t86 * t98;
t131 = t78 ^ 2 * t66;
t121 = t103 * t106;
t117 = t71 ^ 2 * t70 + 0.1e1;
t116 = t102 * t84 + t91 * t98;
t112 = t106 * t118 - t93 * t109 + t92 * t121;
t90 = -t106 * t130 + (-t101 * t121 - t109 * t97) * t100;
t81 = 0.1e1 / t82 ^ 2;
t80 = 0.1e1 / t82;
t75 = 0.1e1 / (t77 ^ 2 * t81 + 0.1e1);
t69 = 0.1e1 / t72;
t68 = 0.1e1 / t117;
t65 = 0.1e1 / t67;
t64 = 0.1e1 / (0.1e1 + t131);
t63 = (t77 * t81 * t90 + t112 * t80) * t98 * t75;
t1 = [t78 * t80 * t75, 0, t63, 0, 0, 0; (t77 * t65 + (t73 + (t80 * t132 - t73) * t75) * t131) * t64, 0 (t87 * t98 * t65 + ((t112 * t73 - t74 * t90) * t98 + (-t73 * t82 + t132) * t63) * t78 * t66) * t64, 0, 0, 0; ((t105 * t112 - t116 * t108) * t69 - (t116 * t105 + t108 * t112) * t133) * t68, 0 ((t102 * t126 + t86 * t105) * t69 - (-t102 * t127 + t86 * t108) * t133) * t68, t117 * t68, 0, 0;];
Ja_rot  = t1;
