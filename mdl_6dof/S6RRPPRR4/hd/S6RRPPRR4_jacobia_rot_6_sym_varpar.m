% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:11
% EndTime: 2019-02-26 21:30:11
% DurationCPUTime: 0.27s
% Computational Cost: add. (867->38), mult. (2363->97), div. (85->9), fcn. (3329->15), ass. (0->56)
t105 = sin(qJ(5));
t109 = cos(qJ(5));
t101 = sin(pkin(6));
t111 = cos(qJ(1));
t119 = t101 * t111;
t107 = sin(qJ(1));
t103 = cos(pkin(6));
t100 = sin(pkin(11));
t102 = cos(pkin(11));
t106 = sin(qJ(2));
t110 = cos(qJ(2));
t98 = t106 * t100 - t110 * t102;
t113 = t98 * t103;
t115 = t110 * t100 + t106 * t102;
t86 = -t107 * t115 - t111 * t113;
t81 = t105 * t119 - t86 * t109;
t94 = t98 * t101;
t92 = t103 * t105 - t94 * t109;
t78 = atan2(t81, t92);
t75 = sin(t78);
t76 = cos(t78);
t69 = t75 * t81 + t76 * t92;
t68 = 0.1e1 / t69 ^ 2;
t112 = t107 * t113 - t111 * t115;
t120 = t101 * t107;
t79 = t105 * t120 + t112 * t109;
t126 = t68 * t79;
t104 = sin(qJ(6));
t108 = cos(qJ(6));
t96 = t115 * t103;
t116 = -t107 * t96 - t111 * t98;
t80 = -t112 * t105 + t109 * t120;
t74 = t104 * t116 + t80 * t108;
t72 = 0.1e1 / t74 ^ 2;
t73 = t80 * t104 - t108 * t116;
t125 = t72 * t73;
t124 = t76 * t81;
t123 = t79 ^ 2 * t68;
t91 = 0.1e1 / t92 ^ 2;
t122 = t81 * t91;
t121 = t116 * t105;
t118 = t73 ^ 2 * t72 + 0.1e1;
t117 = -t75 * t92 + t124;
t87 = t107 * t98 - t111 * t96;
t114 = t86 * t105 + t109 * t119;
t95 = t115 * t101;
t93 = t103 * t109 + t94 * t105;
t90 = 0.1e1 / t92;
t77 = 0.1e1 / (t81 ^ 2 * t91 + 0.1e1);
t71 = 0.1e1 / t74;
t70 = 0.1e1 / t118;
t67 = 0.1e1 / t69;
t66 = 0.1e1 / (0.1e1 + t123);
t65 = (t95 * t122 - t87 * t90) * t77 * t109;
t64 = (t114 * t90 - t93 * t122) * t77;
t1 = [-t79 * t90 * t77, t65, 0, 0, t64, 0; (t81 * t67 - (-t75 + (-t90 * t124 + t75) * t77) * t123) * t66 (-t116 * t109 * t67 - (t117 * t65 + (-t75 * t87 - t76 * t95) * t109) * t126) * t66, 0, 0 (t80 * t67 - (t114 * t75 + t117 * t64 + t76 * t93) * t126) * t66, 0; ((t104 * t114 - t87 * t108) * t71 - (t87 * t104 + t108 * t114) * t125) * t70 ((t104 * t121 - t112 * t108) * t71 - (t112 * t104 + t108 * t121) * t125) * t70, 0, 0 (-t104 * t71 + t108 * t125) * t79 * t70, t118 * t70;];
Ja_rot  = t1;
