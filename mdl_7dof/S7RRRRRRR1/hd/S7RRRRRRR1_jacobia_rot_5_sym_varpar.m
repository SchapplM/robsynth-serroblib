% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% Ja_rot [3x7]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S7RRRRRRR1_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobia_rot_5_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobia_rot_5_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:54:20
% EndTime: 2019-02-26 22:54:20
% DurationCPUTime: 0.24s
% Computational Cost: add. (609->45), mult. (1737->113), div. (115->9), fcn. (2479->13), ass. (0->60)
t102 = sin(qJ(1));
t101 = sin(qJ(2));
t104 = cos(qJ(4));
t117 = t101 * t104;
t100 = sin(qJ(3));
t107 = cos(qJ(1));
t112 = t107 * t100;
t105 = cos(qJ(3));
t106 = cos(qJ(2));
t113 = t105 * t106;
t90 = t102 * t113 + t112;
t99 = sin(qJ(4));
t77 = -t102 * t117 + t90 * t99;
t116 = t101 * t105;
t87 = t104 * t106 + t116 * t99;
t76 = atan2(-t77, t87);
t73 = sin(t76);
t74 = cos(t76);
t67 = -t73 * t77 + t74 * t87;
t66 = 0.1e1 / t67 ^ 2;
t115 = t101 * t107;
t111 = t107 * t105;
t114 = t102 * t100;
t93 = t106 * t111 - t114;
t81 = -t104 * t115 + t93 * t99;
t125 = t66 * t81;
t103 = cos(qJ(5));
t92 = -t102 * t105 - t106 * t112;
t98 = sin(qJ(5));
t120 = t92 * t98;
t82 = t104 * t93 + t115 * t99;
t72 = t103 * t82 + t120;
t70 = 0.1e1 / t72 ^ 2;
t119 = t92 * t103;
t71 = t82 * t98 - t119;
t124 = t70 * t71;
t123 = t74 * t77;
t86 = 0.1e1 / t87 ^ 2;
t122 = t77 * t86;
t121 = t81 ^ 2 * t66;
t118 = t100 * t101;
t110 = t70 * t71 ^ 2 + 0.1e1;
t109 = t101 * t112;
t108 = -t73 * t87 - t123;
t79 = t101 * t102 * t99 + t104 * t90;
t88 = t104 * t116 - t106 * t99;
t91 = t113 * t99 - t117;
t89 = t106 * t114 - t111;
t85 = 0.1e1 / t87;
t84 = t88 * t107;
t83 = t87 * t102;
t75 = 0.1e1 / (t77 ^ 2 * t86 + 0.1e1);
t69 = 0.1e1 / t72;
t68 = 0.1e1 / t110;
t65 = 0.1e1 / t67;
t64 = 0.1e1 / (0.1e1 + t121);
t63 = (-t118 * t122 + t85 * t89) * t99 * t75;
t62 = (t122 * t91 + t83 * t85) * t75;
t61 = (t122 * t88 - t79 * t85) * t75;
t1 = [-t81 * t85 * t75, t62, t63, t61, 0, 0, 0; (-t77 * t65 - (-t73 + (t123 * t85 + t73) * t75) * t121) * t64 (-(t108 * t62 + t73 * t83 + t74 * t91) * t125 - t87 * t65 * t107) * t64 (t92 * t99 * t65 - ((-t118 * t74 + t73 * t89) * t99 + t108 * t63) * t125) * t64 (t82 * t65 - (t108 * t61 - t73 * t79 + t74 * t88) * t125) * t64, 0, 0, 0; ((-t103 * t89 - t79 * t98) * t69 - (-t103 * t79 + t89 * t98) * t124) * t68 ((-t103 * t109 - t84 * t98) * t69 - (-t103 * t84 + t109 * t98) * t124) * t68 ((t103 * t93 + t104 * t120) * t69 - (t104 * t119 - t93 * t98) * t124) * t68 (t103 * t124 - t69 * t98) * t81 * t68, t110 * t68, 0, 0;];
Ja_rot  = t1;
