% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRP5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:48:21
% EndTime: 2019-02-26 21:48:21
% DurationCPUTime: 0.25s
% Computational Cost: add. (867->39), mult. (2363->96), div. (85->9), fcn. (3329->15), ass. (0->57)
t107 = sin(qJ(4));
t111 = cos(qJ(4));
t103 = sin(pkin(6));
t113 = cos(qJ(1));
t119 = t103 * t113;
t109 = sin(qJ(1));
t105 = cos(pkin(6));
t102 = sin(pkin(11));
t104 = cos(pkin(11));
t108 = sin(qJ(2));
t112 = cos(qJ(2));
t115 = t112 * t102 + t108 * t104;
t96 = t115 * t105;
t97 = t108 * t102 - t112 * t104;
t85 = -t109 * t97 + t113 * t96;
t78 = t85 * t107 + t111 * t119;
t95 = t115 * t103;
t91 = -t105 * t111 + t95 * t107;
t77 = atan2(-t78, t91);
t74 = sin(t77);
t75 = cos(t77);
t68 = -t74 * t78 + t75 * t91;
t67 = 0.1e1 / t68 ^ 2;
t116 = -t109 * t96 - t113 * t97;
t120 = t103 * t109;
t82 = t107 * t116 - t111 * t120;
t127 = t67 * t82;
t110 = cos(qJ(5));
t106 = sin(qJ(5));
t114 = t97 * t105;
t87 = t109 * t114 - t113 * t115;
t122 = t87 * t106;
t83 = t107 * t120 + t111 * t116;
t73 = t83 * t110 - t122;
t71 = 0.1e1 / t73 ^ 2;
t121 = t87 * t110;
t72 = t83 * t106 + t121;
t126 = t71 * t72;
t125 = t75 * t78;
t90 = 0.1e1 / t91 ^ 2;
t124 = t78 * t90;
t123 = t82 ^ 2 * t67;
t118 = t72 ^ 2 * t71 + 0.1e1;
t80 = -t107 * t119 + t85 * t111;
t117 = -t74 * t91 - t125;
t94 = t97 * t103;
t92 = t105 * t107 + t95 * t111;
t89 = 0.1e1 / t91;
t84 = -t109 * t115 - t113 * t114;
t76 = 0.1e1 / (t78 ^ 2 * t90 + 0.1e1);
t70 = 0.1e1 / t73;
t69 = 0.1e1 / t118;
t66 = 0.1e1 / t68;
t65 = 0.1e1 / (0.1e1 + t123);
t64 = (-t94 * t124 - t84 * t89) * t76 * t107;
t63 = (t92 * t124 - t80 * t89) * t76;
t1 = [-t82 * t89 * t76, t64, 0, t63, 0, 0; (-t78 * t66 - (-t74 + (t89 * t125 + t74) * t76) * t123) * t65 (t87 * t107 * t66 - (t117 * t64 + (-t74 * t84 - t75 * t94) * t107) * t127) * t65, 0 (t83 * t66 - (t117 * t63 - t74 * t80 + t75 * t92) * t127) * t65, 0, 0; ((-t106 * t80 - t84 * t110) * t70 - (t84 * t106 - t110 * t80) * t126) * t69 ((-t110 * t116 + t111 * t122) * t70 - (t106 * t116 + t111 * t121) * t126) * t69, 0 (-t106 * t70 + t110 * t126) * t82 * t69, t118 * t69, 0;];
Ja_rot  = t1;
