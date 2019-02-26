% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function Ja_rot = S6RRRRRP7_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:43:13
% EndTime: 2019-02-26 22:43:13
% DurationCPUTime: 0.20s
% Computational Cost: add. (1185->38), mult. (1738->90), div. (115->9), fcn. (2485->13), ass. (0->58)
t104 = qJ(3) + qJ(4);
t102 = sin(t104);
t103 = cos(t104);
t105 = sin(pkin(6));
t112 = cos(qJ(1));
t119 = t105 * t112;
t106 = cos(pkin(6));
t108 = sin(qJ(2));
t116 = t112 * t108;
t109 = sin(qJ(1));
t111 = cos(qJ(2));
t117 = t109 * t111;
t97 = t106 * t116 + t117;
t86 = t97 * t102 + t103 * t119;
t122 = t105 * t108;
t94 = t102 * t122 - t106 * t103;
t85 = atan2(-t86, t94);
t80 = sin(t85);
t81 = cos(t85);
t76 = -t80 * t86 + t81 * t94;
t75 = 0.1e1 / t76 ^ 2;
t121 = t105 * t109;
t115 = t112 * t111;
t118 = t109 * t108;
t99 = -t106 * t118 + t115;
t90 = t99 * t102 - t103 * t121;
t129 = t75 * t90;
t110 = cos(qJ(5));
t107 = sin(qJ(5));
t98 = t106 * t117 + t116;
t124 = t98 * t107;
t91 = t102 * t121 + t99 * t103;
t83 = t91 * t110 + t124;
t79 = 0.1e1 / t83 ^ 2;
t123 = t98 * t110;
t82 = t91 * t107 - t123;
t128 = t79 * t82;
t127 = t81 * t86;
t93 = 0.1e1 / t94 ^ 2;
t126 = t86 * t93;
t125 = t90 ^ 2 * t75;
t120 = t105 * t111;
t114 = t82 ^ 2 * t79 + 0.1e1;
t88 = -t102 * t119 + t97 * t103;
t113 = -t80 * t94 - t127;
t96 = t106 * t115 - t118;
t95 = t106 * t102 + t103 * t122;
t92 = 0.1e1 / t94;
t84 = 0.1e1 / (t86 ^ 2 * t93 + 0.1e1);
t78 = 0.1e1 / t83;
t77 = 0.1e1 / t114;
t74 = 0.1e1 / t76;
t73 = 0.1e1 / (0.1e1 + t125);
t72 = (t120 * t126 - t92 * t96) * t84 * t102;
t71 = (t95 * t126 - t88 * t92) * t84;
t70 = (-t107 * t78 + t110 * t128) * t90 * t77;
t69 = (t91 * t74 - (t113 * t71 - t80 * t88 + t81 * t95) * t129) * t73;
t1 = [-t90 * t92 * t84, t72, t71, t71, 0, 0; (-t86 * t74 - (-t80 + (t92 * t127 + t80) * t84) * t125) * t73 (-t98 * t102 * t74 - (t113 * t72 + (t81 * t120 - t80 * t96) * t102) * t129) * t73, t69, t69, 0, 0; ((-t107 * t88 - t96 * t110) * t78 - (t96 * t107 - t110 * t88) * t128) * t77 ((-t103 * t124 - t99 * t110) * t78 - (-t103 * t123 + t99 * t107) * t128) * t77, t70, t70, t114 * t77, 0;];
Ja_rot  = t1;
