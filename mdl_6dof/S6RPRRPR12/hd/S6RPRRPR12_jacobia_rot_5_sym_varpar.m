% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR12
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR12_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:07:16
% EndTime: 2019-02-26 21:07:16
% DurationCPUTime: 0.36s
% Computational Cost: add. (1158->43), mult. (3332->103), div. (77->9), fcn. (4592->15), ass. (0->63)
t106 = sin(qJ(4));
t109 = cos(qJ(4));
t107 = sin(qJ(3));
t110 = cos(qJ(3));
t101 = sin(pkin(7));
t104 = cos(pkin(7));
t102 = sin(pkin(6));
t111 = cos(qJ(1));
t123 = t102 * t111;
t105 = cos(pkin(6));
t103 = cos(pkin(12));
t118 = t111 * t103;
t100 = sin(pkin(12));
t108 = sin(qJ(1));
t121 = t108 * t100;
t96 = -t105 * t118 + t121;
t115 = t101 * t123 + t104 * t96;
t119 = t111 * t100;
t120 = t108 * t103;
t97 = t105 * t119 + t120;
t86 = t115 * t107 - t97 * t110;
t93 = -t96 * t101 + t104 * t123;
t136 = t86 * t106 - t93 * t109;
t135 = t93 * t106 + t86 * t109;
t114 = t105 * t120 + t119;
t124 = t102 * t108;
t134 = -t101 * t124 + t114 * t104;
t122 = t103 * t104;
t125 = t101 * t105;
t92 = t107 * t125 + (t100 * t110 + t107 * t122) * t102;
t95 = -t102 * t103 * t101 + t105 * t104;
t80 = t92 * t106 - t95 * t109;
t71 = atan2(t136, t80);
t68 = sin(t71);
t69 = cos(t71);
t67 = t136 * t68 + t69 * t80;
t66 = 0.1e1 / t67 ^ 2;
t112 = t114 * t101 + t104 * t124;
t98 = -t105 * t121 + t118;
t88 = -t134 * t107 + t98 * t110;
t76 = t88 * t106 - t112 * t109;
t133 = t66 * t76;
t132 = t69 * t136;
t79 = 0.1e1 / t80 ^ 2;
t131 = t136 * t79;
t130 = t76 ^ 2 * t66;
t77 = t112 * t106 + t88 * t109;
t87 = t98 * t107 + t134 * t110;
t83 = 0.1e1 / t87 ^ 2;
t129 = t77 * t83;
t116 = -t68 * t80 + t132;
t84 = -t97 * t107 - t115 * t110;
t91 = t110 * t125 + (-t100 * t107 + t110 * t122) * t102;
t82 = 0.1e1 / t87;
t81 = t95 * t106 + t92 * t109;
t78 = 0.1e1 / t80;
t72 = 0.1e1 / (t77 ^ 2 * t83 + 0.1e1);
t70 = 0.1e1 / (t136 ^ 2 * t79 + 0.1e1);
t65 = 0.1e1 / t67;
t64 = 0.1e1 / (0.1e1 + t130);
t63 = (-t91 * t131 - t78 * t84) * t70 * t106;
t62 = (-t81 * t131 + t135 * t78) * t70;
t1 = [-t76 * t78 * t70, 0, t63, t62, 0, 0; (t136 * t65 - (-t68 + (-t78 * t132 + t68) * t70) * t130) * t64, 0 (-t87 * t106 * t65 - (t116 * t63 + (-t68 * t84 + t69 * t91) * t106) * t133) * t64 (t77 * t65 - (t116 * t62 + t135 * t68 + t69 * t81) * t133) * t64, 0, 0; (-t84 * t129 + t135 * t82) * t72, 0 (-t109 * t82 * t87 - t88 * t129) * t72, -t76 * t82 * t72, 0, 0;];
Ja_rot  = t1;
