% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRR6_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:50:04
% EndTime: 2019-02-26 22:50:04
% DurationCPUTime: 0.24s
% Computational Cost: add. (1326->39), mult. (1821->91), div. (120->9), fcn. (2597->13), ass. (0->59)
t116 = qJ(3) + qJ(4);
t112 = sin(t116);
t114 = cos(t116);
t118 = cos(pkin(6));
t117 = sin(pkin(6));
t119 = sin(qJ(2));
t132 = t117 * t119;
t103 = t112 * t132 - t118 * t114;
t122 = cos(qJ(1));
t126 = t122 * t119;
t120 = sin(qJ(1));
t121 = cos(qJ(2));
t127 = t120 * t121;
t106 = t118 * t126 + t127;
t129 = t117 * t122;
t95 = t106 * t112 + t114 * t129;
t94 = atan2(-t95, t103);
t91 = sin(t94);
t92 = cos(t94);
t85 = t92 * t103 - t91 * t95;
t84 = 0.1e1 / t85 ^ 2;
t125 = t122 * t121;
t128 = t120 * t119;
t108 = -t118 * t128 + t125;
t131 = t117 * t120;
t99 = t108 * t112 - t114 * t131;
t138 = t84 * t99;
t100 = t108 * t114 + t112 * t131;
t107 = t118 * t127 + t126;
t115 = qJ(5) + qJ(6);
t111 = sin(t115);
t113 = cos(t115);
t90 = t100 * t113 + t107 * t111;
t88 = 0.1e1 / t90 ^ 2;
t89 = t100 * t111 - t107 * t113;
t137 = t88 * t89;
t136 = t92 * t95;
t135 = t99 ^ 2 * t84;
t102 = 0.1e1 / t103 ^ 2;
t134 = t102 * t95;
t133 = t107 * t114;
t130 = t117 * t121;
t124 = t89 ^ 2 * t88 + 0.1e1;
t97 = t106 * t114 - t112 * t129;
t123 = -t103 * t91 - t136;
t105 = t118 * t125 - t128;
t104 = t118 * t112 + t114 * t132;
t101 = 0.1e1 / t103;
t93 = 0.1e1 / (t95 ^ 2 * t102 + 0.1e1);
t87 = 0.1e1 / t90;
t86 = 0.1e1 / t124;
t83 = 0.1e1 / t85;
t82 = 0.1e1 / (0.1e1 + t135);
t81 = (-t101 * t105 + t130 * t134) * t93 * t112;
t80 = (-t101 * t97 + t104 * t134) * t93;
t79 = t124 * t86;
t78 = (-t111 * t87 + t113 * t137) * t99 * t86;
t77 = (t100 * t83 - (t92 * t104 + t123 * t80 - t91 * t97) * t138) * t82;
t1 = [-t99 * t101 * t93, t81, t80, t80, 0, 0; (-t95 * t83 - (-t91 + (t101 * t136 + t91) * t93) * t135) * t82 (-t107 * t112 * t83 - (t123 * t81 + (-t105 * t91 + t92 * t130) * t112) * t138) * t82, t77, t77, 0, 0; ((-t105 * t113 - t111 * t97) * t87 - (t105 * t111 - t113 * t97) * t137) * t86 ((-t108 * t113 - t111 * t133) * t87 - (t108 * t111 - t113 * t133) * t137) * t86, t78, t78, t79, t79;];
Ja_rot  = t1;
