% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRP7_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:43:03
% EndTime: 2019-02-26 22:43:03
% DurationCPUTime: 0.24s
% Computational Cost: add. (1185->38), mult. (1738->90), div. (115->9), fcn. (2485->13), ass. (0->58)
t100 = sin(pkin(6));
t107 = cos(qJ(1));
t114 = t100 * t107;
t101 = cos(pkin(6));
t103 = sin(qJ(2));
t111 = t107 * t103;
t104 = sin(qJ(1));
t106 = cos(qJ(2));
t112 = t104 * t106;
t92 = t101 * t111 + t112;
t99 = qJ(3) + qJ(4);
t97 = sin(t99);
t98 = cos(t99);
t81 = t98 * t114 + t92 * t97;
t117 = t100 * t103;
t89 = -t101 * t98 + t97 * t117;
t80 = atan2(-t81, t89);
t75 = sin(t80);
t76 = cos(t80);
t71 = -t75 * t81 + t76 * t89;
t70 = 0.1e1 / t71 ^ 2;
t116 = t100 * t104;
t110 = t107 * t106;
t113 = t104 * t103;
t94 = -t101 * t113 + t110;
t85 = -t98 * t116 + t94 * t97;
t124 = t70 * t85;
t105 = cos(qJ(5));
t102 = sin(qJ(5));
t93 = t101 * t112 + t111;
t119 = t93 * t102;
t86 = t97 * t116 + t94 * t98;
t78 = t86 * t105 + t119;
t74 = 0.1e1 / t78 ^ 2;
t118 = t93 * t105;
t77 = t86 * t102 - t118;
t123 = t74 * t77;
t122 = t76 * t81;
t88 = 0.1e1 / t89 ^ 2;
t121 = t81 * t88;
t120 = t85 ^ 2 * t70;
t115 = t100 * t106;
t109 = t77 ^ 2 * t74 + 0.1e1;
t83 = -t97 * t114 + t92 * t98;
t108 = -t75 * t89 - t122;
t91 = t101 * t110 - t113;
t90 = t101 * t97 + t98 * t117;
t87 = 0.1e1 / t89;
t79 = 0.1e1 / (t81 ^ 2 * t88 + 0.1e1);
t73 = 0.1e1 / t78;
t72 = 0.1e1 / t109;
t69 = 0.1e1 / t71;
t68 = 0.1e1 / (0.1e1 + t120);
t67 = (t115 * t121 - t87 * t91) * t97 * t79;
t66 = (t90 * t121 - t83 * t87) * t79;
t65 = (-t102 * t73 + t105 * t123) * t85 * t72;
t64 = (t86 * t69 - (t108 * t66 - t75 * t83 + t76 * t90) * t124) * t68;
t1 = [-t85 * t87 * t79, t67, t66, t66, 0, 0; (-t81 * t69 - (-t75 + (t87 * t122 + t75) * t79) * t120) * t68 (-t93 * t97 * t69 - ((t76 * t115 - t75 * t91) * t97 + t108 * t67) * t124) * t68, t64, t64, 0, 0; ((-t102 * t83 - t91 * t105) * t73 - (t91 * t102 - t105 * t83) * t123) * t72 ((-t94 * t105 - t98 * t119) * t73 - (t94 * t102 - t98 * t118) * t123) * t72, t65, t65, t109 * t72, 0;];
Ja_rot  = t1;
