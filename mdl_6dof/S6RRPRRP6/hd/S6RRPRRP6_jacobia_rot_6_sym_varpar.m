% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP6
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
% Datum: 2019-02-26 21:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRP6_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:48:50
% EndTime: 2019-02-26 21:48:50
% DurationCPUTime: 0.45s
% Computational Cost: add. (1741->46), mult. (4726->114), div. (107->9), fcn. (6600->15), ass. (0->65)
t126 = cos(pkin(6));
t123 = sin(pkin(11));
t125 = cos(pkin(11));
t129 = sin(qJ(2));
t133 = cos(qJ(2));
t137 = t133 * t123 + t129 * t125;
t118 = t137 * t126;
t121 = t129 * t123 - t133 * t125;
t130 = sin(qJ(1));
t134 = cos(qJ(1));
t109 = t134 * t118 - t130 * t121;
t128 = sin(qJ(4));
t132 = cos(qJ(4));
t124 = sin(pkin(6));
t143 = t124 * t134;
t103 = -t109 * t132 + t128 * t143;
t127 = sin(qJ(5));
t131 = cos(qJ(5));
t135 = t121 * t126;
t141 = -t130 * t137 - t134 * t135;
t154 = t103 * t127 - t141 * t131;
t140 = t141 * t127;
t153 = t103 * t131 + t140;
t117 = t137 * t124;
t114 = t117 * t132 + t126 * t128;
t116 = t121 * t124;
t98 = t114 * t127 - t116 * t131;
t86 = atan2(t154, t98);
t83 = sin(t86);
t84 = cos(t86);
t82 = t154 * t83 + t84 * t98;
t81 = 0.1e1 / t82 ^ 2;
t138 = -t130 * t118 - t134 * t121;
t144 = t124 * t130;
t105 = t128 * t144 + t132 * t138;
t111 = t130 * t135 - t134 * t137;
t145 = t111 * t131;
t93 = t105 * t127 + t145;
t151 = t81 * t93;
t150 = t84 * t154;
t97 = 0.1e1 / t98 ^ 2;
t149 = t154 * t97;
t148 = t93 ^ 2 * t81;
t104 = -t128 * t138 + t132 * t144;
t94 = t105 * t131 - t111 * t127;
t89 = 0.1e1 / t94 ^ 2;
t147 = t104 ^ 2 * t89;
t146 = t104 * t89;
t142 = t127 * t132;
t139 = -t83 * t98 + t150;
t136 = t109 * t128 + t132 * t143;
t113 = -t117 * t128 + t126 * t132;
t106 = -t116 * t142 - t117 * t131;
t99 = t114 * t131 + t116 * t127;
t96 = 0.1e1 / t98;
t95 = -t109 * t131 + t132 * t140;
t88 = 0.1e1 / t94;
t87 = 0.1e1 / (0.1e1 + t147);
t85 = 0.1e1 / (t154 ^ 2 * t97 + 0.1e1);
t80 = 0.1e1 / t82;
t79 = 0.1e1 / (0.1e1 + t148);
t78 = (-t113 * t149 + t136 * t96) * t85 * t127;
t77 = (-t106 * t149 - t95 * t96) * t85;
t76 = (-t99 * t149 + t153 * t96) * t85;
t1 = [-t93 * t96 * t85, t77, 0, t78, t76, 0; (t154 * t80 - (-t83 + (-t96 * t150 + t83) * t85) * t148) * t79 ((t111 * t142 - t131 * t138) * t80 - (t84 * t106 + t139 * t77 - t83 * t95) * t151) * t79, 0 (t104 * t127 * t80 - (t139 * t78 + (t113 * t84 + t136 * t83) * t127) * t151) * t79 (t94 * t80 - (t139 * t76 + t153 * t83 + t84 * t99) * t151) * t79, 0; (t136 * t88 - t153 * t146) * t87 (-t111 * t128 * t88 - (t127 * t138 + t132 * t145) * t146) * t87, 0 (-t105 * t88 - t131 * t147) * t87, t93 * t87 * t146, 0;];
Ja_rot  = t1;
