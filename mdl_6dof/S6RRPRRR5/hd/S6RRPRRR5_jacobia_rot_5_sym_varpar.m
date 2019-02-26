% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR5_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:56:19
% EndTime: 2019-02-26 21:56:19
% DurationCPUTime: 0.22s
% Computational Cost: add. (867->39), mult. (2363->96), div. (85->9), fcn. (3329->15), ass. (0->57)
t104 = sin(qJ(4));
t108 = cos(qJ(4));
t100 = sin(pkin(6));
t110 = cos(qJ(1));
t117 = t100 * t110;
t106 = sin(qJ(1));
t101 = cos(pkin(12));
t105 = sin(qJ(2));
t109 = cos(qJ(2));
t99 = sin(pkin(12));
t113 = t109 * t101 - t105 * t99;
t102 = cos(pkin(6));
t112 = t105 * t101 + t109 * t99;
t93 = t112 * t102;
t82 = t106 * t113 + t110 * t93;
t75 = t82 * t104 + t108 * t117;
t92 = t112 * t100;
t88 = -t102 * t108 + t92 * t104;
t74 = atan2(-t75, t88);
t71 = sin(t74);
t72 = cos(t74);
t65 = -t71 * t75 + t72 * t88;
t64 = 0.1e1 / t65 ^ 2;
t114 = -t106 * t93 + t110 * t113;
t118 = t100 * t106;
t79 = t104 * t114 - t108 * t118;
t125 = t64 * t79;
t107 = cos(qJ(5));
t103 = sin(qJ(5));
t111 = t113 * t102;
t84 = -t106 * t111 - t110 * t112;
t120 = t84 * t103;
t80 = t104 * t118 + t108 * t114;
t70 = t80 * t107 - t120;
t68 = 0.1e1 / t70 ^ 2;
t119 = t84 * t107;
t69 = t80 * t103 + t119;
t124 = t68 * t69;
t123 = t72 * t75;
t87 = 0.1e1 / t88 ^ 2;
t122 = t75 * t87;
t121 = t79 ^ 2 * t64;
t116 = t69 ^ 2 * t68 + 0.1e1;
t77 = -t104 * t117 + t82 * t108;
t115 = -t71 * t88 - t123;
t91 = t113 * t100;
t89 = t102 * t104 + t92 * t108;
t86 = 0.1e1 / t88;
t81 = -t106 * t112 + t110 * t111;
t73 = 0.1e1 / (t75 ^ 2 * t87 + 0.1e1);
t67 = 0.1e1 / t70;
t66 = 0.1e1 / t116;
t63 = 0.1e1 / t65;
t62 = 0.1e1 / (0.1e1 + t121);
t61 = (t91 * t122 - t81 * t86) * t73 * t104;
t60 = (t89 * t122 - t77 * t86) * t73;
t1 = [-t79 * t86 * t73, t61, 0, t60, 0, 0; (-t75 * t63 - (-t71 + (t86 * t123 + t71) * t73) * t121) * t62 (t84 * t104 * t63 - (t115 * t61 + (-t71 * t81 + t72 * t91) * t104) * t125) * t62, 0 (t80 * t63 - (t115 * t60 - t71 * t77 + t72 * t89) * t125) * t62, 0, 0; ((-t103 * t77 - t81 * t107) * t67 - (t81 * t103 - t107 * t77) * t124) * t66 ((-t107 * t114 + t108 * t120) * t67 - (t103 * t114 + t108 * t119) * t124) * t66, 0 (-t103 * t67 + t107 * t124) * t79 * t66, t116 * t66, 0;];
Ja_rot  = t1;
