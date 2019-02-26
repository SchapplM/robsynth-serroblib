% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPP9_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:30:16
% EndTime: 2019-02-26 22:30:16
% DurationCPUTime: 0.29s
% Computational Cost: add. (896->45), mult. (2497->112), div. (107->9), fcn. (3529->13), ass. (0->61)
t101 = sin(qJ(4));
t105 = cos(qJ(4));
t102 = sin(qJ(3));
t106 = cos(qJ(3));
t108 = cos(qJ(1));
t99 = sin(pkin(6));
t119 = t108 * t99;
t100 = cos(pkin(6));
t103 = sin(qJ(2));
t111 = t108 * t103;
t104 = sin(qJ(1));
t107 = cos(qJ(2));
t114 = t104 * t107;
t94 = t100 * t111 + t114;
t87 = t102 * t119 - t94 * t106;
t110 = t108 * t107;
t115 = t104 * t103;
t93 = -t100 * t110 + t115;
t73 = -t87 * t101 - t93 * t105;
t127 = -t93 * t101 + t87 * t105;
t120 = t106 * t99;
t92 = t100 * t102 + t103 * t120;
t84 = -t99 * t107 * t101 + t92 * t105;
t72 = atan2(t127, t84);
t68 = sin(t72);
t69 = cos(t72);
t67 = t127 * t68 + t69 * t84;
t66 = 0.1e1 / t67 ^ 2;
t95 = t100 * t114 + t111;
t116 = t95 * t101;
t121 = t102 * t99;
t96 = -t100 * t115 + t110;
t89 = t104 * t121 + t96 * t106;
t77 = t89 * t105 + t116;
t126 = t66 * t77;
t125 = t69 * t127;
t80 = 0.1e1 / t84 ^ 2;
t124 = t127 * t80;
t76 = -t89 * t101 + t95 * t105;
t88 = t96 * t102 - t104 * t120;
t82 = 0.1e1 / t88 ^ 2;
t123 = t76 * t82;
t122 = t77 ^ 2 * t66;
t113 = t105 * t106;
t112 = t105 * t107;
t109 = -t68 * t84 + t125;
t85 = -t94 * t102 - t106 * t119;
t91 = t100 * t106 - t103 * t121;
t90 = (t101 * t103 + t106 * t112) * t99;
t83 = -t92 * t101 - t99 * t112;
t81 = 0.1e1 / t88;
t79 = 0.1e1 / t84;
t78 = t94 * t101 - t93 * t113;
t71 = 0.1e1 / (t76 ^ 2 * t82 + 0.1e1);
t70 = 0.1e1 / (t127 ^ 2 * t80 + 0.1e1);
t65 = 0.1e1 / t67;
t64 = 0.1e1 / (0.1e1 + t122);
t63 = (-t91 * t124 - t79 * t85) * t70 * t105;
t62 = (-t90 * t124 - t78 * t79) * t70;
t61 = (-t83 * t124 + t73 * t79) * t70;
t1 = [-t77 * t79 * t70, t62, t63, t61, 0, 0; (t127 * t65 - (-t68 + (-t79 * t125 + t68) * t70) * t122) * t64 ((t96 * t101 - t95 * t113) * t65 - (t109 * t62 - t68 * t78 + t69 * t90) * t126) * t64 (-t88 * t105 * t65 - (t109 * t63 + (-t68 * t85 + t69 * t91) * t105) * t126) * t64 (t76 * t65 - (t109 * t61 + t68 * t73 + t69 * t83) * t126) * t64, 0, 0; (-t85 * t123 + t73 * t81) * t71 ((t96 * t105 + t106 * t116) * t81 + t95 * t102 * t123) * t71 (t101 * t81 * t88 - t89 * t123) * t71, -t77 * t81 * t71, 0, 0;];
Ja_rot  = t1;
