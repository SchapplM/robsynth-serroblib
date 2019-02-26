% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR7_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:34:06
% EndTime: 2019-02-26 22:34:07
% DurationCPUTime: 0.24s
% Computational Cost: add. (1756->38), mult. (1738->90), div. (115->9), fcn. (2485->13), ass. (0->58)
t104 = qJ(3) + qJ(4) + pkin(12);
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
t81 = atan2(-t86, t94);
t78 = sin(t81);
t79 = cos(t81);
t76 = -t78 * t86 + t79 * t94;
t75 = 0.1e1 / t76 ^ 2;
t121 = t105 * t109;
t115 = t112 * t111;
t118 = t109 * t108;
t99 = -t106 * t118 + t115;
t90 = t99 * t102 - t103 * t121;
t129 = t75 * t90;
t128 = t79 * t86;
t110 = cos(qJ(6));
t107 = sin(qJ(6));
t98 = t106 * t117 + t116;
t124 = t98 * t107;
t91 = t102 * t121 + t99 * t103;
t85 = t91 * t110 + t124;
t83 = 0.1e1 / t85 ^ 2;
t123 = t98 * t110;
t84 = t91 * t107 - t123;
t127 = t83 * t84;
t93 = 0.1e1 / t94 ^ 2;
t126 = t86 * t93;
t125 = t90 ^ 2 * t75;
t120 = t105 * t111;
t114 = t84 ^ 2 * t83 + 0.1e1;
t88 = -t102 * t119 + t97 * t103;
t113 = -t78 * t94 - t128;
t96 = t106 * t115 - t118;
t95 = t106 * t102 + t103 * t122;
t92 = 0.1e1 / t94;
t82 = 0.1e1 / t85;
t80 = 0.1e1 / (t86 ^ 2 * t93 + 0.1e1);
t77 = 0.1e1 / t114;
t74 = 0.1e1 / t76;
t73 = 0.1e1 / (0.1e1 + t125);
t72 = (t120 * t126 - t92 * t96) * t80 * t102;
t71 = (t95 * t126 - t88 * t92) * t80;
t70 = (-t107 * t82 + t110 * t127) * t90 * t77;
t69 = (t91 * t74 - (t113 * t71 - t78 * t88 + t79 * t95) * t129) * t73;
t1 = [-t90 * t92 * t80, t72, t71, t71, 0, 0; (-t86 * t74 - (-t78 + (t92 * t128 + t78) * t80) * t125) * t73 (-t98 * t102 * t74 - (t113 * t72 + (t79 * t120 - t78 * t96) * t102) * t129) * t73, t69, t69, 0, 0; ((-t107 * t88 - t96 * t110) * t82 - (t96 * t107 - t110 * t88) * t127) * t77 ((-t103 * t124 - t99 * t110) * t82 - (-t103 * t123 + t99 * t107) * t127) * t77, t70, t70, 0, t114 * t77;];
Ja_rot  = t1;
