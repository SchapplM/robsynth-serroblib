% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPR8_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:37
% EndTime: 2019-02-26 20:14:37
% DurationCPUTime: 0.35s
% Computational Cost: add. (1524->59), mult. (4378->142), div. (95->9), fcn. (6002->17), ass. (0->76)
t126 = sin(pkin(6));
t132 = sin(qJ(3));
t133 = sin(qJ(2));
t136 = cos(qJ(3));
t137 = cos(qJ(2));
t128 = cos(pkin(7));
t147 = t128 * t132;
t125 = sin(pkin(7));
t129 = cos(pkin(6));
t150 = t125 * t129;
t113 = t132 * t150 + (t133 * t136 + t137 * t147) * t126;
t151 = t125 * t126;
t117 = t129 * t128 - t137 * t151;
t131 = sin(qJ(4));
t135 = cos(qJ(4));
t108 = t113 * t135 + t117 * t131;
t124 = sin(pkin(12));
t127 = cos(pkin(12));
t145 = t129 * t133;
t119 = t124 * t137 + t127 * t145;
t144 = t129 * t137;
t118 = -t124 * t133 + t127 * t144;
t138 = t118 * t128 - t127 * t151;
t104 = t119 * t136 + t138 * t132;
t148 = t126 * t128;
t114 = -t118 * t125 - t127 * t148;
t95 = t104 * t135 + t114 * t131;
t93 = atan2(-t95, t108);
t90 = sin(t93);
t91 = cos(t93);
t84 = t108 * t91 - t90 * t95;
t83 = 0.1e1 / t84 ^ 2;
t120 = -t124 * t144 - t127 * t133;
t121 = -t124 * t145 + t127 * t137;
t140 = t124 * t151;
t106 = t121 * t136 + (t120 * t128 + t140) * t132;
t115 = -t120 * t125 + t124 * t148;
t98 = t106 * t135 + t115 * t131;
t156 = t83 * t98;
t130 = sin(qJ(6));
t146 = t128 * t136;
t105 = -t120 * t146 + t121 * t132 - t136 * t140;
t134 = cos(qJ(6));
t152 = t105 * t134;
t97 = t106 * t131 - t115 * t135;
t89 = t130 * t97 + t152;
t87 = 0.1e1 / t89 ^ 2;
t153 = t105 * t130;
t88 = -t134 * t97 + t153;
t155 = t87 * t88;
t102 = 0.1e1 / t108 ^ 2;
t154 = t102 * t95;
t149 = t125 * t131;
t143 = t132 * t133;
t142 = t136 * t137;
t141 = t87 * t88 ^ 2 + 0.1e1;
t139 = -t108 * t90 - t91 * t95;
t112 = t136 * t150 + (t128 * t142 - t143) * t126;
t111 = ((-t128 * t143 + t142) * t135 + t133 * t149) * t126;
t110 = t120 * t136 - t121 * t147;
t109 = t120 * t132 + t121 * t146;
t107 = -t113 * t131 + t117 * t135;
t103 = -t119 * t132 + t138 * t136;
t101 = 0.1e1 / t108;
t100 = -t121 * t125 * t135 + t110 * t131;
t99 = (t118 * t136 - t119 * t147) * t135 + t119 * t149;
t94 = t104 * t131 - t114 * t135;
t92 = 0.1e1 / (t102 * t95 ^ 2 + 0.1e1);
t86 = 0.1e1 / t89;
t85 = 0.1e1 / t141;
t82 = 0.1e1 / t84;
t81 = 0.1e1 / (t83 * t98 ^ 2 + 0.1e1);
t80 = (-t101 * t103 + t112 * t154) * t92 * t135;
t79 = (-t101 * t99 + t111 * t154) * t92;
t78 = (t101 * t94 + t107 * t154) * t92;
t1 = [0, t79, t80, t78, 0, 0; 0 ((t110 * t135 + t121 * t149) * t82 - (t91 * t111 + t139 * t79 - t90 * t99) * t156) * t81 (-t105 * t135 * t82 - (t139 * t80 + (-t103 * t90 + t112 * t91) * t135) * t156) * t81 (-t97 * t82 - (t91 * t107 + t139 * t78 + t90 * t94) * t156) * t81, 0, 0; 0 ((-t100 * t134 + t109 * t130) * t86 - (t100 * t130 + t109 * t134) * t155) * t85 ((t106 * t130 + t131 * t152) * t86 - (t106 * t134 - t131 * t153) * t155) * t85 (-t130 * t155 - t134 * t86) * t98 * t85, 0, t141 * t85;];
Ja_rot  = t1;
