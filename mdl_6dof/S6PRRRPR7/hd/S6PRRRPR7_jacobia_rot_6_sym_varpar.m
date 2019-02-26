% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR7
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
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPR7_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_jacobia_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:04
% EndTime: 2019-02-26 20:14:04
% DurationCPUTime: 0.32s
% Computational Cost: add. (1590->60), mult. (4378->143), div. (95->9), fcn. (6002->17), ass. (0->76)
t130 = sin(pkin(12));
t133 = cos(pkin(12));
t138 = sin(qJ(2));
t135 = cos(pkin(6));
t141 = cos(qJ(2));
t151 = t135 * t141;
t123 = -t130 * t151 - t133 * t138;
t152 = t135 * t138;
t124 = -t130 * t152 + t133 * t141;
t134 = cos(pkin(7));
t137 = sin(qJ(3));
t140 = cos(qJ(3));
t131 = sin(pkin(7));
t132 = sin(pkin(6));
t156 = t132 * t131;
t147 = t130 * t156;
t110 = t124 * t140 + (t123 * t134 + t147) * t137;
t136 = sin(qJ(4));
t139 = cos(qJ(4));
t155 = t132 * t134;
t145 = -t123 * t131 + t130 * t155;
t102 = t110 * t139 + t145 * t136;
t153 = t134 * t140;
t109 = -t123 * t153 + t124 * t137 - t140 * t147;
t129 = pkin(13) + qJ(6);
t127 = sin(t129);
t128 = cos(t129);
t93 = t102 * t128 + t109 * t127;
t91 = 0.1e1 / t93 ^ 2;
t92 = t102 * t127 - t109 * t128;
t162 = t91 * t92;
t101 = t110 * t136 - t145 * t139;
t154 = t134 * t137;
t158 = t131 * t135;
t119 = t137 * t158 + (t138 * t140 + t141 * t154) * t132;
t121 = t135 * t134 - t141 * t156;
t111 = t119 * t136 - t121 * t139;
t122 = t130 * t141 + t133 * t152;
t144 = -t130 * t138 + t133 * t151;
t142 = -t133 * t156 + t144 * t134;
t108 = t122 * t140 + t142 * t137;
t143 = -t144 * t131 - t133 * t155;
t98 = t108 * t136 - t143 * t139;
t97 = atan2(-t98, t111);
t94 = sin(t97);
t95 = cos(t97);
t88 = t95 * t111 - t94 * t98;
t87 = 0.1e1 / t88 ^ 2;
t161 = t101 * t87;
t106 = 0.1e1 / t111 ^ 2;
t160 = t106 * t98;
t159 = t109 * t139;
t157 = t131 * t139;
t150 = t137 * t138;
t149 = t140 * t141;
t148 = t92 ^ 2 * t91 + 0.1e1;
t146 = -t111 * t94 - t95 * t98;
t118 = t140 * t158 + (t134 * t149 - t150) * t132;
t115 = ((-t134 * t150 + t149) * t136 - t138 * t157) * t132;
t114 = t123 * t140 - t124 * t154;
t113 = t123 * t137 + t124 * t153;
t112 = t119 * t139 + t121 * t136;
t107 = -t122 * t137 + t142 * t140;
t105 = 0.1e1 / t111;
t104 = t124 * t131 * t136 + t114 * t139;
t103 = (-t122 * t154 + t144 * t140) * t136 - t122 * t157;
t100 = t108 * t139 + t143 * t136;
t96 = 0.1e1 / (t98 ^ 2 * t106 + 0.1e1);
t90 = 0.1e1 / t93;
t89 = 0.1e1 / t148;
t86 = 0.1e1 / t88;
t85 = 0.1e1 / (t101 ^ 2 * t87 + 0.1e1);
t84 = (-t105 * t107 + t118 * t160) * t96 * t136;
t83 = (-t103 * t105 + t115 * t160) * t96;
t82 = (-t100 * t105 + t112 * t160) * t96;
t1 = [0, t83, t84, t82, 0, 0; 0 ((t114 * t136 - t124 * t157) * t86 - (-t94 * t103 + t95 * t115 + t146 * t83) * t161) * t85 (-t109 * t136 * t86 - (t146 * t84 + (-t107 * t94 + t118 * t95) * t136) * t161) * t85 (t102 * t86 - (-t94 * t100 + t95 * t112 + t146 * t82) * t161) * t85, 0, 0; 0 ((t104 * t127 - t113 * t128) * t90 - (t104 * t128 + t113 * t127) * t162) * t89 ((-t110 * t128 - t127 * t159) * t90 - (t110 * t127 - t128 * t159) * t162) * t89 (-t127 * t90 + t128 * t162) * t89 * t101, 0, t148 * t89;];
Ja_rot  = t1;
