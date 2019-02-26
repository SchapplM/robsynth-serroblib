% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR10_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:35:53
% EndTime: 2019-02-26 22:35:53
% DurationCPUTime: 0.22s
% Computational Cost: add. (1185->38), mult. (1738->90), div. (115->9), fcn. (2485->13), ass. (0->58)
t106 = cos(qJ(1));
t99 = sin(pkin(6));
t115 = t106 * t99;
t100 = cos(pkin(6));
t102 = sin(qJ(2));
t110 = t106 * t102;
t103 = sin(qJ(1));
t105 = cos(qJ(2));
t111 = t103 * t105;
t92 = t100 * t110 + t111;
t98 = qJ(3) + qJ(4);
t96 = sin(t98);
t97 = cos(t98);
t82 = -t96 * t115 + t92 * t97;
t118 = t102 * t99;
t90 = t100 * t96 + t97 * t118;
t80 = atan2(-t82, t90);
t75 = sin(t80);
t76 = cos(t80);
t71 = -t75 * t82 + t76 * t90;
t70 = 0.1e1 / t71 ^ 2;
t117 = t103 * t99;
t109 = t106 * t105;
t112 = t103 * t102;
t94 = -t100 * t112 + t109;
t86 = t96 * t117 + t94 * t97;
t123 = t70 * t86;
t101 = sin(qJ(6));
t104 = cos(qJ(6));
t93 = t100 * t111 + t110;
t113 = t93 * t104;
t85 = -t97 * t117 + t94 * t96;
t78 = t85 * t101 + t113;
t74 = 0.1e1 / t78 ^ 2;
t114 = t93 * t101;
t77 = -t85 * t104 + t114;
t122 = t74 * t77;
t121 = t76 * t82;
t88 = 0.1e1 / t90 ^ 2;
t120 = t82 * t88;
t119 = t86 ^ 2 * t70;
t116 = t105 * t99;
t108 = t77 ^ 2 * t74 + 0.1e1;
t107 = -t75 * t90 - t121;
t81 = t97 * t115 + t92 * t96;
t91 = t100 * t109 - t112;
t89 = t100 * t97 - t96 * t118;
t87 = 0.1e1 / t90;
t79 = 0.1e1 / (t82 ^ 2 * t88 + 0.1e1);
t73 = 0.1e1 / t78;
t72 = 0.1e1 / t108;
t69 = 0.1e1 / t71;
t68 = 0.1e1 / (0.1e1 + t119);
t67 = (t116 * t120 - t87 * t91) * t97 * t79;
t66 = (t89 * t120 + t81 * t87) * t79;
t65 = (-t101 * t122 - t104 * t73) * t86 * t72;
t64 = (-t85 * t69 - (t107 * t66 + t75 * t81 + t76 * t89) * t123) * t68;
t1 = [-t86 * t87 * t79, t67, t66, t66, 0, 0; (-t82 * t69 - (-t75 + (t87 * t121 + t75) * t79) * t119) * t68 (-t93 * t97 * t69 - ((t76 * t116 - t75 * t91) * t97 + t107 * t67) * t123) * t68, t64, t64, 0, 0; ((t91 * t101 + t104 * t81) * t73 - (-t101 * t81 + t91 * t104) * t122) * t72 ((t94 * t101 + t96 * t113) * t73 - (t94 * t104 - t96 * t114) * t122) * t72, t65, t65, 0, t108 * t72;];
Ja_rot  = t1;
