% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPPRRR1
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PPPRRR1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_jacobia_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:38:44
% EndTime: 2019-02-26 19:38:44
% DurationCPUTime: 0.48s
% Computational Cost: add. (2922->55), mult. (8425->130), div. (65->9), fcn. (11413->21), ass. (0->81)
t128 = sin(pkin(8));
t134 = cos(pkin(8));
t126 = sin(pkin(13));
t132 = cos(pkin(13));
t133 = cos(pkin(12));
t127 = sin(pkin(12));
t136 = cos(pkin(6));
t162 = t127 * t136;
t124 = -t126 * t162 + t133 * t132;
t125 = sin(pkin(14));
t131 = cos(pkin(14));
t123 = -t133 * t126 - t132 * t162;
t135 = cos(pkin(7));
t129 = sin(pkin(7));
t130 = sin(pkin(6));
t161 = t129 * t130;
t150 = t123 * t135 + t127 * t161;
t147 = t124 * t125 - t150 * t131;
t159 = t130 * t135;
t151 = -t123 * t129 + t127 * t159;
t167 = -t151 * t128 + t147 * t134;
t158 = t132 * t135;
t160 = t129 * t136;
t119 = t130 * t126 * t131 + (t130 * t158 + t160) * t125;
t139 = sin(qJ(4));
t142 = cos(qJ(4));
t118 = t131 * t160 + (-t125 * t126 + t131 * t158) * t130;
t120 = -t132 * t161 + t136 * t135;
t154 = t118 * t134 + t120 * t128;
t112 = t119 * t142 + t154 * t139;
t116 = -t118 * t128 + t120 * t134;
t138 = sin(qJ(5));
t141 = cos(qJ(5));
t103 = t112 * t138 - t116 * t141;
t157 = t133 * t136;
t122 = t126 * t157 + t127 * t132;
t121 = -t127 * t126 + t132 * t157;
t152 = t121 * t135 - t133 * t161;
t114 = t122 * t131 + t152 * t125;
t148 = -t122 * t125 + t152 * t131;
t153 = -t121 * t129 - t133 * t159;
t144 = t153 * t128 + t148 * t134;
t106 = t114 * t142 + t144 * t139;
t145 = -t148 * t128 + t153 * t134;
t96 = t106 * t138 - t145 * t141;
t95 = atan2(-t96, t103);
t92 = sin(t95);
t93 = cos(t95);
t86 = t93 * t103 - t92 * t96;
t85 = 0.1e1 / t86 ^ 2;
t115 = t124 * t131 + t150 * t125;
t108 = t115 * t142 - t167 * t139;
t143 = t147 * t128 + t151 * t134;
t99 = t108 * t138 - t143 * t141;
t166 = t85 * t99;
t100 = t108 * t141 + t143 * t138;
t107 = t115 * t139 + t167 * t142;
t137 = sin(qJ(6));
t140 = cos(qJ(6));
t91 = t100 * t140 + t107 * t137;
t89 = 0.1e1 / t91 ^ 2;
t90 = t100 * t137 - t107 * t140;
t165 = t89 * t90;
t102 = 0.1e1 / t103 ^ 2;
t164 = t102 * t96;
t163 = t107 * t141;
t156 = t90 ^ 2 * t89 + 0.1e1;
t155 = -t103 * t92 - t93 * t96;
t111 = -t119 * t139 + t154 * t142;
t105 = -t114 * t139 + t144 * t142;
t104 = t112 * t141 + t116 * t138;
t101 = 0.1e1 / t103;
t98 = t106 * t141 + t145 * t138;
t94 = 0.1e1 / (t96 ^ 2 * t102 + 0.1e1);
t88 = 0.1e1 / t91;
t87 = 0.1e1 / t156;
t84 = 0.1e1 / t86;
t83 = 0.1e1 / (t99 ^ 2 * t85 + 0.1e1);
t82 = (-t101 * t105 + t111 * t164) * t94 * t138;
t81 = (-t101 * t98 + t104 * t164) * t94;
t1 = [0, 0, 0, t82, t81, 0; 0, 0, 0 (-t107 * t138 * t84 - (t155 * t82 + (-t105 * t92 + t111 * t93) * t138) * t166) * t83 (t100 * t84 - (t93 * t104 + t155 * t81 - t92 * t98) * t166) * t83, 0; 0, 0, 0 ((-t108 * t140 - t137 * t163) * t88 - (t108 * t137 - t140 * t163) * t165) * t87 (-t137 * t88 + t140 * t165) * t99 * t87, t156 * t87;];
Ja_rot  = t1;
