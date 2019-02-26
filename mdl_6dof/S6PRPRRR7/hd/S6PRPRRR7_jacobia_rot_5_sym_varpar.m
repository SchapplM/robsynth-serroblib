% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRR7_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobia_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:57:29
% EndTime: 2019-02-26 19:57:30
% DurationCPUTime: 0.40s
% Computational Cost: add. (1772->61), mult. (5223->145), div. (65->9), fcn. (7061->19), ass. (0->77)
t127 = sin(pkin(13));
t132 = cos(pkin(13));
t141 = cos(qJ(2));
t135 = cos(pkin(6));
t138 = sin(qJ(2));
t153 = t135 * t138;
t122 = t127 * t141 + t132 * t153;
t126 = sin(pkin(14));
t128 = sin(pkin(8));
t129 = sin(pkin(7));
t131 = cos(pkin(14));
t133 = cos(pkin(8));
t134 = cos(pkin(7));
t152 = t135 * t141;
t145 = -t127 * t138 + t132 * t152;
t130 = sin(pkin(6));
t157 = t132 * t130;
t144 = -t129 * t157 + t145 * t134;
t167 = (-t122 * t126 + t144 * t131) * t133 + (-t145 * t129 - t134 * t157) * t128;
t124 = -t127 * t153 + t132 * t141;
t123 = -t127 * t152 - t132 * t138;
t158 = t130 * t129;
t146 = t123 * t134 + t127 * t158;
t112 = -t124 * t126 + t146 * t131;
t120 = t127 * t130 * t134 - t123 * t129;
t166 = t112 * t133 + t120 * t128;
t154 = t134 * t141;
t159 = t129 * t135;
t119 = t130 * t138 * t131 + (t130 * t154 + t159) * t126;
t137 = sin(qJ(4));
t140 = cos(qJ(4));
t148 = (t131 * t159 + (-t138 * t126 + t131 * t154) * t130) * t133 + (t135 * t134 - t141 * t158) * t128;
t104 = t119 * t137 - t148 * t140;
t111 = t122 * t131 + t144 * t126;
t95 = t111 * t137 - t167 * t140;
t94 = atan2(-t95, t104);
t91 = sin(t94);
t92 = cos(t94);
t85 = t92 * t104 - t91 * t95;
t84 = 0.1e1 / t85 ^ 2;
t113 = t124 * t131 + t146 * t126;
t98 = t113 * t137 - t166 * t140;
t165 = t84 * t98;
t106 = -t112 * t128 + t120 * t133;
t136 = sin(qJ(5));
t139 = cos(qJ(5));
t99 = t113 * t140 + t166 * t137;
t90 = t106 * t136 + t99 * t139;
t88 = 0.1e1 / t90 ^ 2;
t89 = -t106 * t139 + t99 * t136;
t164 = t88 * t89;
t103 = 0.1e1 / t104 ^ 2;
t163 = t103 * t95;
t160 = t129 * t128;
t156 = t134 * t126;
t155 = t134 * t131;
t151 = t138 * t134;
t150 = t89 ^ 2 * t88 + 0.1e1;
t149 = -t104 * t91 - t92 * t95;
t114 = -t123 * t126 - t124 * t155;
t147 = t114 * t133 + t124 * t160;
t115 = t123 * t131 - t124 * t156;
t108 = ((-t126 * t151 + t131 * t141) * t137 + (-(-t126 * t141 - t131 * t151) * t133 - t138 * t160) * t140) * t130;
t107 = t124 * t129 * t133 - t114 * t128;
t105 = t119 * t140 + t148 * t137;
t102 = 0.1e1 / t104;
t101 = t115 * t140 + t147 * t137;
t100 = (-t122 * t156 + t145 * t131) * t137 + (-(-t122 * t155 - t145 * t126) * t133 - t122 * t160) * t140;
t97 = t111 * t140 + t167 * t137;
t93 = 0.1e1 / (t95 ^ 2 * t103 + 0.1e1);
t87 = 0.1e1 / t90;
t86 = 0.1e1 / t150;
t83 = 0.1e1 / t85;
t82 = 0.1e1 / (t98 ^ 2 * t84 + 0.1e1);
t81 = (-t100 * t102 + t108 * t163) * t93;
t80 = (-t102 * t97 + t105 * t163) * t93;
t1 = [0, t81, 0, t80, 0, 0; 0 ((t115 * t137 - t147 * t140) * t83 - (-t91 * t100 + t92 * t108 + t149 * t81) * t165) * t82, 0 (t99 * t83 - (t92 * t105 + t149 * t80 - t91 * t97) * t165) * t82, 0, 0; 0 ((t101 * t136 - t107 * t139) * t87 - (t101 * t139 + t107 * t136) * t164) * t86, 0 (-t136 * t87 + t139 * t164) * t98 * t86, t150 * t86, 0;];
Ja_rot  = t1;
