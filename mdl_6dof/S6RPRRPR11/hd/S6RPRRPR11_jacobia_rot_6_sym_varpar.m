% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR11_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_jacobia_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:06:44
% EndTime: 2019-02-26 21:06:44
% DurationCPUTime: 0.39s
% Computational Cost: add. (1506->50), mult. (4121->118), div. (85->9), fcn. (5660->17), ass. (0->70)
t133 = cos(pkin(6));
t131 = cos(pkin(12));
t164 = sin(qJ(1));
t146 = t164 * t131;
t128 = sin(pkin(12));
t138 = cos(qJ(1));
t151 = t138 * t128;
t122 = t133 * t151 + t146;
t135 = sin(qJ(3));
t137 = cos(qJ(3));
t147 = t164 * t128;
t150 = t138 * t131;
t121 = -t133 * t150 + t147;
t129 = sin(pkin(7));
t132 = cos(pkin(7));
t130 = sin(pkin(6));
t153 = t130 * t138;
t143 = t121 * t132 + t129 * t153;
t111 = -t122 * t137 + t143 * t135;
t118 = -t121 * t129 + t132 * t153;
t134 = sin(qJ(4));
t136 = cos(qJ(4));
t166 = t111 * t134 - t118 * t136;
t101 = t111 * t136 + t118 * t134;
t142 = t133 * t146 + t151;
t148 = t130 * t164;
t165 = -t129 * t148 + t142 * t132;
t123 = -t133 * t147 + t150;
t113 = t123 * t137 - t165 * t135;
t139 = t142 * t129 + t132 * t148;
t103 = t113 * t136 + t139 * t134;
t112 = t123 * t135 + t165 * t137;
t127 = pkin(13) + qJ(6);
t125 = sin(t127);
t126 = cos(t127);
t93 = t103 * t126 + t112 * t125;
t91 = 0.1e1 / t93 ^ 2;
t92 = t103 * t125 - t112 * t126;
t163 = t91 * t92;
t152 = t131 * t132;
t154 = t129 * t133;
t117 = t135 * t154 + (t128 * t137 + t135 * t152) * t130;
t120 = -t130 * t131 * t129 + t133 * t132;
t106 = t117 * t134 - t120 * t136;
t97 = atan2(t166, t106);
t95 = cos(t97);
t162 = t95 * t166;
t102 = t113 * t134 - t139 * t136;
t94 = sin(t97);
t88 = t95 * t106 + t166 * t94;
t87 = 0.1e1 / t88 ^ 2;
t161 = t102 * t87;
t160 = t102 ^ 2 * t87;
t105 = 0.1e1 / t106 ^ 2;
t159 = t105 * t166;
t158 = t112 * t136;
t149 = t92 ^ 2 * t91 + 0.1e1;
t144 = -t106 * t94 + t162;
t140 = -t122 * t135 - t143 * t137;
t116 = t137 * t154 + (-t128 * t135 + t137 * t152) * t130;
t107 = t117 * t136 + t120 * t134;
t104 = 0.1e1 / t106;
t96 = 0.1e1 / (t105 * t166 ^ 2 + 0.1e1);
t90 = 0.1e1 / t93;
t89 = 0.1e1 / t149;
t86 = 0.1e1 / t88;
t85 = 0.1e1 / (0.1e1 + t160);
t84 = (-t104 * t140 - t116 * t159) * t96 * t134;
t83 = (t101 * t104 - t107 * t159) * t96;
t1 = [-t102 * t104 * t96, 0, t84, t83, 0, 0; (t166 * t86 - (-t94 + (-t104 * t162 + t94) * t96) * t160) * t85, 0 (-t112 * t134 * t86 - (t144 * t84 + (t116 * t95 - t140 * t94) * t134) * t161) * t85 (t103 * t86 - (t101 * t94 + t95 * t107 + t144 * t83) * t161) * t85, 0, 0; ((t101 * t125 - t126 * t140) * t90 - (t101 * t126 + t125 * t140) * t163) * t89, 0 ((-t113 * t126 - t125 * t158) * t90 - (t113 * t125 - t126 * t158) * t163) * t89 (-t125 * t90 + t126 * t163) * t89 * t102, 0, t149 * t89;];
Ja_rot  = t1;
