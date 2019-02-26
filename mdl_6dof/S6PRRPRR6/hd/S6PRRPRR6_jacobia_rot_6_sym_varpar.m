% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR6
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRR6_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_jacobia_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:05
% EndTime: 2019-02-26 20:07:06
% DurationCPUTime: 0.35s
% Computational Cost: add. (1983->60), mult. (4378->142), div. (95->9), fcn. (6002->17), ass. (0->77)
t129 = sin(pkin(6));
t134 = sin(qJ(3));
t135 = sin(qJ(2));
t137 = cos(qJ(3));
t138 = cos(qJ(2));
t131 = cos(pkin(7));
t151 = t131 * t134;
t128 = sin(pkin(7));
t132 = cos(pkin(6));
t153 = t128 * t132;
t116 = t134 * t153 + (t135 * t137 + t138 * t151) * t129;
t154 = t128 * t129;
t118 = t132 * t131 - t138 * t154;
t126 = pkin(13) + qJ(5);
t124 = sin(t126);
t125 = cos(t126);
t104 = t116 * t124 - t118 * t125;
t127 = sin(pkin(12));
t130 = cos(pkin(12));
t149 = t132 * t135;
t119 = t127 * t138 + t130 * t149;
t148 = t132 * t138;
t141 = -t127 * t135 + t130 * t148;
t139 = -t130 * t154 + t141 * t131;
t107 = t119 * t137 + t139 * t134;
t152 = t129 * t131;
t140 = -t141 * t128 - t130 * t152;
t95 = t107 * t124 - t140 * t125;
t94 = atan2(-t95, t104);
t89 = sin(t94);
t90 = cos(t94);
t85 = t104 * t90 - t89 * t95;
t84 = 0.1e1 / t85 ^ 2;
t120 = -t127 * t148 - t130 * t135;
t121 = -t127 * t149 + t130 * t138;
t144 = t127 * t154;
t109 = t121 * t137 + (t120 * t131 + t144) * t134;
t142 = -t120 * t128 + t127 * t152;
t98 = t109 * t124 - t142 * t125;
t160 = t84 * t98;
t136 = cos(qJ(6));
t150 = t131 * t137;
t108 = -t120 * t150 + t121 * t134 - t137 * t144;
t133 = sin(qJ(6));
t157 = t108 * t133;
t99 = t109 * t125 + t142 * t124;
t92 = t136 * t99 + t157;
t88 = 0.1e1 / t92 ^ 2;
t156 = t108 * t136;
t91 = t133 * t99 - t156;
t159 = t88 * t91;
t103 = 0.1e1 / t104 ^ 2;
t158 = t103 * t95;
t155 = t125 * t128;
t147 = t134 * t135;
t146 = t137 * t138;
t145 = t88 * t91 ^ 2 + 0.1e1;
t143 = -t104 * t89 - t90 * t95;
t115 = t137 * t153 + (t131 * t146 - t147) * t129;
t112 = ((-t131 * t147 + t146) * t124 - t135 * t155) * t129;
t111 = t120 * t137 - t121 * t151;
t110 = t120 * t134 + t121 * t150;
t106 = -t119 * t134 + t139 * t137;
t105 = t116 * t125 + t118 * t124;
t102 = 0.1e1 / t104;
t101 = t121 * t124 * t128 + t111 * t125;
t100 = (-t119 * t151 + t141 * t137) * t124 - t119 * t155;
t97 = t107 * t125 + t140 * t124;
t93 = 0.1e1 / (t103 * t95 ^ 2 + 0.1e1);
t87 = 0.1e1 / t92;
t86 = 0.1e1 / t145;
t83 = 0.1e1 / t85;
t82 = 0.1e1 / (t84 * t98 ^ 2 + 0.1e1);
t81 = (-t102 * t106 + t115 * t158) * t93 * t124;
t80 = (-t100 * t102 + t112 * t158) * t93;
t79 = (-t102 * t97 + t105 * t158) * t93;
t1 = [0, t80, t81, 0, t79, 0; 0 ((t111 * t124 - t121 * t155) * t83 - (-t89 * t100 + t90 * t112 + t143 * t80) * t160) * t82 (-t108 * t124 * t83 - (t143 * t81 + (-t106 * t89 + t115 * t90) * t124) * t160) * t82, 0 (t99 * t83 - (t90 * t105 + t143 * t79 - t89 * t97) * t160) * t82, 0; 0 ((t101 * t133 - t110 * t136) * t87 - (t101 * t136 + t110 * t133) * t159) * t86 ((-t109 * t136 - t125 * t157) * t87 - (t109 * t133 - t125 * t156) * t159) * t86, 0 (-t133 * t87 + t136 * t159) * t98 * t86, t145 * t86;];
Ja_rot  = t1;
