% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRR4_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_jacobia_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:20:26
% EndTime: 2019-02-26 20:20:26
% DurationCPUTime: 0.19s
% Computational Cost: add. (721->41), mult. (1971->97), div. (70->9), fcn. (2694->15), ass. (0->62)
t105 = cos(pkin(13));
t106 = cos(pkin(7));
t102 = sin(pkin(13));
t109 = sin(qJ(2));
t107 = cos(pkin(6));
t111 = cos(qJ(2));
t120 = t107 * t111;
t113 = -t102 * t109 + t105 * t120;
t103 = sin(pkin(7));
t104 = sin(pkin(6));
t123 = t104 * t103;
t130 = -t105 * t123 + t113 * t106;
t108 = sin(qJ(3));
t110 = cos(qJ(3));
t121 = t107 * t109;
t94 = t102 * t111 + t105 * t121;
t79 = t94 * t108 - t130 * t110;
t122 = t106 * t110;
t124 = t103 * t107;
t88 = -t110 * t124 + (t108 * t109 - t111 * t122) * t104;
t78 = atan2(-t79, t88);
t75 = sin(t78);
t76 = cos(t78);
t69 = -t75 * t79 + t76 * t88;
t68 = 0.1e1 / t69 ^ 2;
t116 = t102 * t123;
t96 = -t102 * t121 + t105 * t111;
t125 = t96 * t108;
t95 = -t102 * t120 - t105 * t109;
t82 = -t110 * t116 - t95 * t122 + t125;
t129 = t68 * t82;
t101 = qJ(4) + qJ(5);
t100 = cos(t101);
t83 = t96 * t110 + (t106 * t95 + t116) * t108;
t90 = t102 * t104 * t106 - t95 * t103;
t99 = sin(t101);
t74 = t83 * t100 + t90 * t99;
t72 = 0.1e1 / t74 ^ 2;
t73 = -t90 * t100 + t83 * t99;
t128 = t72 * t73;
t87 = 0.1e1 / t88 ^ 2;
t127 = t79 * t87;
t126 = t103 * t96;
t119 = t108 * t111;
t118 = t109 * t110;
t117 = t73 ^ 2 * t72 + 0.1e1;
t114 = -t75 * t88 - t76 * t79;
t93 = (t106 * t118 + t119) * t104;
t89 = t108 * t124 + (t106 * t119 + t118) * t104;
t86 = 0.1e1 / t88;
t85 = -t106 * t125 + t95 * t110;
t84 = t113 * t108 + t94 * t122;
t81 = t130 * t108 + t94 * t110;
t77 = 0.1e1 / (t79 ^ 2 * t87 + 0.1e1);
t71 = 0.1e1 / t74;
t70 = 0.1e1 / t117;
t67 = 0.1e1 / t69;
t66 = 0.1e1 / (t82 ^ 2 * t68 + 0.1e1);
t65 = (t93 * t127 - t84 * t86) * t77;
t64 = (t89 * t127 - t81 * t86) * t77;
t63 = t117 * t70;
t1 = [0, t65, t64, 0, 0, 0; 0 ((t95 * t108 + t96 * t122) * t67 - (t114 * t65 - t75 * t84 + t76 * t93) * t129) * t66 (t83 * t67 - (t114 * t64 - t75 * t81 + t76 * t89) * t129) * t66, 0, 0, 0; 0 ((-t100 * t126 + t85 * t99) * t71 - (t85 * t100 + t99 * t126) * t128) * t70 (t100 * t128 - t99 * t71) * t82 * t70, t63, t63, 0;];
Ja_rot  = t1;
