% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRR8_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:07
% EndTime: 2019-02-26 20:08:07
% DurationCPUTime: 0.36s
% Computational Cost: add. (1524->59), mult. (4378->141), div. (95->9), fcn. (6002->17), ass. (0->76)
t121 = sin(pkin(6));
t127 = sin(qJ(3));
t128 = sin(qJ(2));
t131 = cos(qJ(3));
t132 = cos(qJ(2));
t123 = cos(pkin(7));
t144 = t123 * t131;
t120 = sin(pkin(7));
t124 = cos(pkin(6));
t147 = t120 * t124;
t109 = -t131 * t147 + (t128 * t127 - t132 * t144) * t121;
t148 = t120 * t121;
t114 = t124 * t123 - t132 * t148;
t126 = sin(qJ(5));
t130 = cos(qJ(5));
t104 = -t109 * t130 + t114 * t126;
t122 = cos(pkin(12));
t119 = sin(pkin(12));
t142 = t124 * t132;
t136 = -t119 * t128 + t122 * t142;
t145 = t121 * t123;
t111 = -t136 * t120 - t122 * t145;
t134 = -t122 * t148 + t136 * t123;
t143 = t124 * t128;
t135 = t119 * t132 + t122 * t143;
t133 = t135 * t127 - t134 * t131;
t90 = t111 * t126 - t133 * t130;
t89 = atan2(-t90, t104);
t86 = sin(t89);
t87 = cos(t89);
t80 = t87 * t104 - t86 * t90;
t79 = 0.1e1 / t80 ^ 2;
t115 = -t119 * t142 - t122 * t128;
t138 = t119 * t148;
t116 = -t119 * t143 + t122 * t132;
t149 = t116 * t127;
t102 = -t115 * t144 - t131 * t138 + t149;
t112 = -t115 * t120 + t119 * t145;
t93 = -t102 * t130 + t112 * t126;
t154 = t79 * t93;
t129 = cos(qJ(6));
t103 = t116 * t131 + (t115 * t123 + t138) * t127;
t125 = sin(qJ(6));
t151 = t103 * t125;
t94 = t102 * t126 + t112 * t130;
t85 = t94 * t129 + t151;
t83 = 0.1e1 / t85 ^ 2;
t150 = t103 * t129;
t84 = t94 * t125 - t150;
t153 = t83 * t84;
t100 = 0.1e1 / t104 ^ 2;
t152 = t100 * t90;
t146 = t120 * t126;
t141 = t128 * t131;
t140 = t132 * t127;
t139 = t84 ^ 2 * t83 + 0.1e1;
t137 = -t104 * t86 - t87 * t90;
t110 = t127 * t147 + (t123 * t140 + t141) * t121;
t108 = (t128 * t146 - (t123 * t141 + t140) * t130) * t121;
t107 = t115 * t131 - t123 * t149;
t106 = t115 * t127 + t116 * t144;
t105 = t109 * t126 + t114 * t130;
t101 = t134 * t127 + t135 * t131;
t99 = 0.1e1 / t104;
t96 = t116 * t120 * t130 + t106 * t126;
t95 = -t135 * t146 + (t136 * t127 + t135 * t144) * t130;
t92 = t111 * t130 + t133 * t126;
t88 = 0.1e1 / (t90 ^ 2 * t100 + 0.1e1);
t82 = 0.1e1 / t85;
t81 = 0.1e1 / t139;
t78 = 0.1e1 / t80;
t77 = 0.1e1 / (t93 ^ 2 * t79 + 0.1e1);
t76 = (t101 * t99 - t110 * t152) * t88 * t130;
t75 = (t108 * t152 + t95 * t99) * t88;
t74 = (t105 * t152 - t92 * t99) * t88;
t1 = [0, t75, t76, 0, t74, 0; 0 ((-t106 * t130 + t116 * t146) * t78 - (t87 * t108 + t137 * t75 + t86 * t95) * t154) * t77 (-t103 * t130 * t78 - (t137 * t76 + (t101 * t86 - t110 * t87) * t130) * t154) * t77, 0 (t94 * t78 - (t87 * t105 + t137 * t74 - t86 * t92) * t154) * t77, 0; 0 ((-t107 * t129 + t96 * t125) * t82 - (t107 * t125 + t96 * t129) * t153) * t81 ((t102 * t129 + t126 * t151) * t82 - (-t102 * t125 + t126 * t150) * t153) * t81, 0 (-t125 * t82 + t129 * t153) * t93 * t81, t139 * t81;];
Ja_rot  = t1;
