% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Ja_rot = S6RRRRPP9_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:30:06
% EndTime: 2019-02-26 22:30:06
% DurationCPUTime: 0.32s
% Computational Cost: add. (896->45), mult. (2497->111), div. (107->9), fcn. (3529->13), ass. (0->62)
t102 = sin(qJ(4));
t106 = cos(qJ(4));
t101 = cos(pkin(6));
t108 = cos(qJ(2));
t109 = cos(qJ(1));
t113 = t109 * t108;
t104 = sin(qJ(2));
t105 = sin(qJ(1));
t116 = t105 * t104;
t112 = -t101 * t113 + t116;
t103 = sin(qJ(3));
t107 = cos(qJ(3));
t100 = sin(pkin(6));
t118 = t100 * t109;
t114 = t109 * t104;
t115 = t105 * t108;
t95 = t101 * t114 + t115;
t87 = t103 * t118 - t95 * t107;
t130 = t87 * t102 + t112 * t106;
t110 = t112 * t102;
t129 = t87 * t106 - t110;
t119 = t100 * t107;
t94 = t101 * t103 + t104 * t119;
t83 = t100 * t108 * t106 + t94 * t102;
t72 = atan2(t130, t83);
t68 = sin(t72);
t69 = cos(t72);
t67 = t130 * t68 + t69 * t83;
t66 = 0.1e1 / t67 ^ 2;
t96 = t101 * t115 + t114;
t121 = t96 * t106;
t120 = t100 * t103;
t97 = -t101 * t116 + t113;
t89 = t105 * t120 + t97 * t107;
t76 = t89 * t102 - t121;
t127 = t66 * t76;
t126 = t69 * t130;
t80 = 0.1e1 / t83 ^ 2;
t125 = t130 * t80;
t124 = t76 ^ 2 * t66;
t122 = t96 * t102;
t77 = t89 * t106 + t122;
t88 = t97 * t103 - t105 * t119;
t82 = 0.1e1 / t88 ^ 2;
t123 = t77 * t82;
t117 = t102 * t108;
t111 = -t68 * t83 + t126;
t85 = -t95 * t103 - t107 * t118;
t93 = t101 * t107 - t104 * t120;
t90 = (-t104 * t106 + t107 * t117) * t100;
t84 = -t100 * t117 + t94 * t106;
t81 = 0.1e1 / t88;
t79 = 0.1e1 / t83;
t78 = -t95 * t106 - t107 * t110;
t71 = 0.1e1 / (t77 ^ 2 * t82 + 0.1e1);
t70 = 0.1e1 / (t130 ^ 2 * t80 + 0.1e1);
t65 = 0.1e1 / t67;
t64 = 0.1e1 / (0.1e1 + t124);
t63 = (-t93 * t125 - t79 * t85) * t70 * t102;
t62 = (-t90 * t125 - t78 * t79) * t70;
t61 = (-t84 * t125 + t129 * t79) * t70;
t1 = [-t76 * t79 * t70, t62, t63, t61, 0, 0; (t130 * t65 - (-t68 + (-t79 * t126 + t68) * t70) * t124) * t64 ((-t97 * t106 - t107 * t122) * t65 - (t111 * t62 - t68 * t78 + t69 * t90) * t127) * t64 (-t88 * t102 * t65 - (t111 * t63 + (-t68 * t85 + t69 * t93) * t102) * t127) * t64 (t77 * t65 - (t111 * t61 + t129 * t68 + t69 * t84) * t127) * t64, 0, 0; (-t85 * t123 + t129 * t81) * t71 ((t97 * t102 - t107 * t121) * t81 + t96 * t103 * t123) * t71 (-t106 * t81 * t88 - t89 * t123) * t71, -t76 * t81 * t71, 0, 0;];
Ja_rot  = t1;
