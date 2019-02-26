% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRR3_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_jacobia_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:24
% EndTime: 2019-02-26 20:05:24
% DurationCPUTime: 0.24s
% Computational Cost: add. (1073->42), mult. (3064->103), div. (65->9), fcn. (4203->17), ass. (0->62)
t106 = cos(pkin(12));
t104 = sin(pkin(6));
t103 = sin(pkin(7));
t101 = sin(pkin(13));
t105 = cos(pkin(13));
t110 = sin(qJ(3));
t113 = cos(qJ(3));
t97 = t110 * t101 - t113 * t105;
t116 = t97 * t103;
t115 = t104 * t116;
t117 = t113 * t101 + t110 * t105;
t107 = cos(pkin(7));
t91 = t97 * t107;
t102 = sin(pkin(12));
t111 = sin(qJ(2));
t108 = cos(pkin(6));
t114 = cos(qJ(2));
t120 = t108 * t114;
t93 = -t102 * t111 + t106 * t120;
t121 = t108 * t111;
t94 = t102 * t114 + t106 * t121;
t77 = t106 * t115 - t117 * t94 - t93 * t91;
t84 = (t111 * t117 + t114 * t91) * t104 + t108 * t116;
t71 = atan2(t77, t84);
t68 = sin(t71);
t69 = cos(t71);
t66 = t68 * t77 + t69 * t84;
t65 = 0.1e1 / t66 ^ 2;
t95 = -t102 * t120 - t106 * t111;
t96 = -t102 * t121 + t106 * t114;
t78 = -t102 * t115 - t117 * t96 - t95 * t91;
t126 = t65 * t78;
t109 = sin(qJ(5));
t112 = cos(qJ(5));
t122 = t102 * t104;
t90 = t117 * t103;
t92 = t117 * t107;
t80 = t90 * t122 + t95 * t92 - t96 * t97;
t88 = -t95 * t103 + t107 * t122;
t75 = t88 * t109 + t80 * t112;
t73 = 0.1e1 / t75 ^ 2;
t74 = t80 * t109 - t88 * t112;
t125 = t73 * t74;
t82 = 0.1e1 / t84 ^ 2;
t124 = t77 * t82;
t123 = t103 * t96;
t119 = t74 ^ 2 * t73 + 0.1e1;
t118 = -t68 * t84 + t69 * t77;
t87 = (-t111 * t91 + t114 * t117) * t104;
t86 = -t96 * t92 - t95 * t97;
t85 = -t117 * t93 + t94 * t91;
t83 = t108 * t90 + (-t111 * t97 + t114 * t92) * t104;
t81 = 0.1e1 / t84;
t76 = t106 * t104 * t90 - t93 * t92 + t94 * t97;
t72 = 0.1e1 / t75;
t70 = 0.1e1 / (t77 ^ 2 * t82 + 0.1e1);
t67 = 0.1e1 / t119;
t64 = 0.1e1 / t66;
t63 = 0.1e1 / (t78 ^ 2 * t65 + 0.1e1);
t62 = (-t87 * t124 + t81 * t85) * t70;
t61 = (-t83 * t124 + t76 * t81) * t70;
t1 = [0, t62, t61, 0, 0, 0; 0 ((t117 * t95 - t96 * t91) * t64 + (t118 * t62 + t68 * t85 + t69 * t87) * t126) * t63 (t80 * t64 + (t118 * t61 + t68 * t76 + t69 * t83) * t126) * t63, 0, 0, 0; 0 ((t86 * t109 - t112 * t123) * t72 - (t109 * t123 + t86 * t112) * t125) * t67 (t109 * t72 - t112 * t125) * t78 * t67, 0, t119 * t67, 0;];
Ja_rot  = t1;
