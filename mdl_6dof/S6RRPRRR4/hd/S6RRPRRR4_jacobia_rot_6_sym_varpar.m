% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR4
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
% Datum: 2019-02-26 21:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:55:46
% EndTime: 2019-02-26 21:55:46
% DurationCPUTime: 0.30s
% Computational Cost: add. (1726->40), mult. (3150->96), div. (115->9), fcn. (4438->15), ass. (0->60)
t123 = qJ(4) + qJ(5);
t121 = sin(t123);
t122 = cos(t123);
t127 = cos(pkin(6));
t124 = sin(pkin(12));
t126 = cos(pkin(12));
t129 = sin(qJ(2));
t132 = cos(qJ(2));
t135 = t132 * t124 + t129 * t126;
t115 = t135 * t127;
t116 = t129 * t124 - t132 * t126;
t130 = sin(qJ(1));
t133 = cos(qJ(1));
t136 = -t130 * t115 - t133 * t116;
t125 = sin(pkin(6));
t140 = t125 * t130;
t102 = t121 * t140 + t122 * t136;
t131 = cos(qJ(6));
t134 = t116 * t127;
t106 = t130 * t134 - t133 * t135;
t128 = sin(qJ(6));
t142 = t106 * t128;
t92 = t102 * t131 - t142;
t90 = 0.1e1 / t92 ^ 2;
t141 = t106 * t131;
t91 = t102 * t128 + t141;
t147 = t90 * t91;
t114 = t135 * t125;
t110 = t114 * t121 - t127 * t122;
t104 = t133 * t115 - t130 * t116;
t139 = t125 * t133;
t97 = t104 * t121 + t122 * t139;
t96 = atan2(-t97, t110);
t94 = cos(t96);
t146 = t94 * t97;
t101 = t121 * t136 - t122 * t140;
t93 = sin(t96);
t87 = t94 * t110 - t93 * t97;
t86 = 0.1e1 / t87 ^ 2;
t145 = t101 * t86;
t144 = t101 ^ 2 * t86;
t109 = 0.1e1 / t110 ^ 2;
t143 = t109 * t97;
t138 = t91 ^ 2 * t90 + 0.1e1;
t99 = t104 * t122 - t121 * t139;
t137 = -t110 * t93 - t146;
t113 = t116 * t125;
t111 = t114 * t122 + t127 * t121;
t108 = 0.1e1 / t110;
t103 = -t130 * t135 - t133 * t134;
t95 = 0.1e1 / (t97 ^ 2 * t109 + 0.1e1);
t89 = 0.1e1 / t92;
t88 = 0.1e1 / t138;
t85 = 0.1e1 / t87;
t84 = 0.1e1 / (0.1e1 + t144);
t83 = (-t103 * t108 - t113 * t143) * t95 * t121;
t82 = (-t108 * t99 + t111 * t143) * t95;
t81 = (-t128 * t89 + t131 * t147) * t88 * t101;
t80 = (t102 * t85 - (t94 * t111 + t137 * t82 - t93 * t99) * t145) * t84;
t1 = [-t101 * t108 * t95, t83, 0, t82, t82, 0; (-t97 * t85 - (-t93 + (t108 * t146 + t93) * t95) * t144) * t84 (t106 * t121 * t85 - (t137 * t83 + (-t103 * t93 - t113 * t94) * t121) * t145) * t84, 0, t80, t80, 0; ((-t103 * t131 - t128 * t99) * t89 - (t103 * t128 - t131 * t99) * t147) * t88 ((t122 * t142 - t131 * t136) * t89 - (t122 * t141 + t128 * t136) * t147) * t88, 0, t81, t81, t138 * t88;];
Ja_rot  = t1;
