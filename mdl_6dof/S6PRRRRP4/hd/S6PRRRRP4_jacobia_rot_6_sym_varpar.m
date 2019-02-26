% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRP4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:16:54
% EndTime: 2019-02-26 20:16:54
% DurationCPUTime: 0.28s
% Computational Cost: add. (1511->41), mult. (2652->100), div. (114->9), fcn. (3738->13), ass. (0->60)
t111 = qJ(4) + qJ(5);
t109 = sin(t111);
t110 = cos(t111);
t112 = sin(pkin(11));
t114 = cos(pkin(11));
t117 = sin(qJ(2));
t115 = cos(pkin(6));
t119 = cos(qJ(2));
t123 = t115 * t119;
t120 = -t112 * t117 + t114 * t123;
t124 = t115 * t117;
t104 = t112 * t119 + t114 * t124;
t118 = cos(qJ(3));
t113 = sin(pkin(6));
t116 = sin(qJ(3));
t127 = t113 * t116;
t98 = t104 * t118 - t114 * t127;
t86 = t109 * t98 + t120 * t110;
t126 = t113 * t118;
t108 = t115 * t116 + t117 * t126;
t125 = t113 * t119;
t94 = t108 * t109 + t110 * t125;
t82 = atan2(-t86, t94);
t79 = sin(t82);
t80 = cos(t82);
t78 = -t79 * t86 + t80 * t94;
t77 = 0.1e1 / t78 ^ 2;
t106 = -t112 * t124 + t114 * t119;
t100 = t106 * t118 + t112 * t127;
t105 = t112 * t123 + t114 * t117;
t129 = t105 * t110;
t89 = t100 * t109 - t129;
t133 = t77 * t89;
t90 = t100 * t110 + t105 * t109;
t85 = 0.1e1 / t90 ^ 2;
t99 = -t106 * t116 + t112 * t126;
t132 = t85 * t99 ^ 2;
t131 = t85 * t99;
t93 = 0.1e1 / t94 ^ 2;
t130 = t86 * t93;
t128 = t109 * t118;
t83 = 0.1e1 / (0.1e1 + t132);
t122 = t89 * t83 * t131;
t121 = -t79 * t94 - t80 * t86;
t107 = t115 * t118 - t117 * t127;
t101 = (-t110 * t117 + t119 * t128) * t113;
t97 = -t104 * t116 - t114 * t126;
t95 = t108 * t110 - t109 * t125;
t92 = 0.1e1 / t94;
t91 = -t104 * t110 + t120 * t128;
t88 = -t120 * t109 + t98 * t110;
t84 = 0.1e1 / t90;
t81 = 0.1e1 / (t86 ^ 2 * t93 + 0.1e1);
t76 = 0.1e1 / t78;
t75 = 0.1e1 / (t77 * t89 ^ 2 + 0.1e1);
t74 = (t107 * t130 - t92 * t97) * t81 * t109;
t73 = (t101 * t130 - t91 * t92) * t81;
t72 = (t95 * t130 - t88 * t92) * t81;
t71 = (t90 * t76 - (t121 * t72 - t79 * t88 + t80 * t95) * t133) * t75;
t1 = [0, t73, t74, t72, t72, 0; 0 ((-t105 * t128 - t106 * t110) * t76 - (t80 * t101 + t121 * t73 - t79 * t91) * t133) * t75 (t99 * t109 * t76 - (t121 * t74 + (t107 * t80 - t79 * t97) * t109) * t133) * t75, t71, t71, 0; 0 (t105 * t116 * t84 - (t106 * t109 - t118 * t129) * t131) * t83 (-t100 * t84 - t110 * t132) * t83, t122, t122, 0;];
Ja_rot  = t1;
