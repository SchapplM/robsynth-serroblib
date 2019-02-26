% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPR5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_jacobia_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:12:55
% EndTime: 2019-02-26 20:12:55
% DurationCPUTime: 0.34s
% Computational Cost: add. (1983->60), mult. (4378->142), div. (95->9), fcn. (6002->17), ass. (0->77)
t130 = sin(pkin(6));
t135 = sin(qJ(3));
t136 = sin(qJ(2));
t138 = cos(qJ(3));
t139 = cos(qJ(2));
t132 = cos(pkin(7));
t152 = t132 * t135;
t129 = sin(pkin(7));
t133 = cos(pkin(6));
t155 = t129 * t133;
t117 = t135 * t155 + (t136 * t138 + t139 * t152) * t130;
t154 = t130 * t129;
t119 = t133 * t132 - t139 * t154;
t127 = qJ(4) + pkin(13);
t125 = sin(t127);
t126 = cos(t127);
t105 = t117 * t125 - t119 * t126;
t128 = sin(pkin(12));
t131 = cos(pkin(12));
t150 = t133 * t136;
t120 = t128 * t139 + t131 * t150;
t149 = t133 * t139;
t142 = -t128 * t136 + t131 * t149;
t140 = -t131 * t154 + t142 * t132;
t108 = t120 * t138 + t140 * t135;
t153 = t130 * t132;
t141 = -t142 * t129 - t131 * t153;
t96 = t108 * t125 - t141 * t126;
t95 = atan2(-t96, t105);
t90 = sin(t95);
t91 = cos(t95);
t86 = t105 * t91 - t90 * t96;
t85 = 0.1e1 / t86 ^ 2;
t121 = -t128 * t149 - t131 * t136;
t122 = -t128 * t150 + t131 * t139;
t145 = t128 * t154;
t110 = t122 * t138 + (t121 * t132 + t145) * t135;
t143 = -t121 * t129 + t128 * t153;
t99 = t110 * t125 - t143 * t126;
t161 = t85 * t99;
t100 = t110 * t126 + t143 * t125;
t137 = cos(qJ(6));
t151 = t132 * t138;
t109 = -t121 * t151 + t122 * t135 - t138 * t145;
t134 = sin(qJ(6));
t158 = t109 * t134;
t93 = t100 * t137 + t158;
t89 = 0.1e1 / t93 ^ 2;
t157 = t109 * t137;
t92 = t100 * t134 - t157;
t160 = t89 * t92;
t104 = 0.1e1 / t105 ^ 2;
t159 = t104 * t96;
t156 = t126 * t129;
t148 = t135 * t136;
t147 = t138 * t139;
t146 = t89 * t92 ^ 2 + 0.1e1;
t144 = -t105 * t90 - t91 * t96;
t116 = t138 * t155 + (t132 * t147 - t148) * t130;
t113 = ((-t132 * t148 + t147) * t125 - t136 * t156) * t130;
t112 = t121 * t138 - t122 * t152;
t111 = t121 * t135 + t122 * t151;
t107 = -t120 * t135 + t140 * t138;
t106 = t117 * t126 + t119 * t125;
t103 = 0.1e1 / t105;
t102 = t122 * t125 * t129 + t112 * t126;
t101 = (-t120 * t152 + t142 * t138) * t125 - t120 * t156;
t98 = t108 * t126 + t141 * t125;
t94 = 0.1e1 / (t104 * t96 ^ 2 + 0.1e1);
t88 = 0.1e1 / t93;
t87 = 0.1e1 / t146;
t84 = 0.1e1 / t86;
t83 = 0.1e1 / (t85 * t99 ^ 2 + 0.1e1);
t82 = (-t103 * t107 + t116 * t159) * t94 * t125;
t81 = (-t101 * t103 + t113 * t159) * t94;
t80 = (-t103 * t98 + t106 * t159) * t94;
t1 = [0, t81, t82, t80, 0, 0; 0 ((t112 * t125 - t122 * t156) * t84 - (-t101 * t90 + t113 * t91 + t144 * t81) * t161) * t83 (-t109 * t125 * t84 - (t144 * t82 + (-t107 * t90 + t116 * t91) * t125) * t161) * t83 (t100 * t84 - (t106 * t91 + t144 * t80 - t90 * t98) * t161) * t83, 0, 0; 0 ((t102 * t134 - t111 * t137) * t88 - (t102 * t137 + t111 * t134) * t160) * t87 ((-t110 * t137 - t126 * t158) * t88 - (t110 * t134 - t126 * t157) * t160) * t87 (-t134 * t88 + t137 * t160) * t99 * t87, 0, t146 * t87;];
Ja_rot  = t1;
