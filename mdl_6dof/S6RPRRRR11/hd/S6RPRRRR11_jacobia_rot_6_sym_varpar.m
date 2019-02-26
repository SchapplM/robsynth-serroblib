% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRR11_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_jacobia_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:20:26
% EndTime: 2019-02-26 21:20:26
% DurationCPUTime: 0.42s
% Computational Cost: add. (1607->50), mult. (4365->118), div. (90->9), fcn. (5989->17), ass. (0->71)
t134 = cos(pkin(6));
t132 = cos(pkin(13));
t165 = sin(qJ(1));
t147 = t165 * t132;
t129 = sin(pkin(13));
t139 = cos(qJ(1));
t152 = t139 * t129;
t123 = t134 * t152 + t147;
t136 = sin(qJ(3));
t138 = cos(qJ(3));
t148 = t165 * t129;
t151 = t139 * t132;
t122 = -t134 * t151 + t148;
t130 = sin(pkin(7));
t133 = cos(pkin(7));
t131 = sin(pkin(6));
t154 = t131 * t139;
t144 = t122 * t133 + t130 * t154;
t112 = -t123 * t138 + t144 * t136;
t119 = -t122 * t130 + t133 * t154;
t135 = sin(qJ(4));
t137 = cos(qJ(4));
t167 = t112 * t135 - t119 * t137;
t102 = t112 * t137 + t119 * t135;
t143 = t134 * t147 + t152;
t149 = t131 * t165;
t166 = -t130 * t149 + t143 * t133;
t124 = -t134 * t148 + t151;
t114 = t124 * t138 - t166 * t136;
t140 = t143 * t130 + t133 * t149;
t104 = t114 * t137 + t140 * t135;
t113 = t124 * t136 + t166 * t138;
t128 = qJ(5) + qJ(6);
t126 = sin(t128);
t127 = cos(t128);
t94 = t104 * t127 + t113 * t126;
t92 = 0.1e1 / t94 ^ 2;
t93 = t104 * t126 - t113 * t127;
t164 = t92 * t93;
t153 = t132 * t133;
t155 = t130 * t134;
t118 = t136 * t155 + (t129 * t138 + t136 * t153) * t131;
t121 = -t130 * t131 * t132 + t133 * t134;
t107 = t118 * t135 - t121 * t137;
t98 = atan2(t167, t107);
t96 = cos(t98);
t163 = t96 * t167;
t103 = t114 * t135 - t140 * t137;
t95 = sin(t98);
t89 = t107 * t96 + t167 * t95;
t88 = 0.1e1 / t89 ^ 2;
t162 = t103 * t88;
t161 = t103 ^ 2 * t88;
t106 = 0.1e1 / t107 ^ 2;
t160 = t106 * t167;
t159 = t113 * t137;
t150 = t92 * t93 ^ 2 + 0.1e1;
t145 = -t107 * t95 + t163;
t141 = -t123 * t136 - t144 * t138;
t117 = t138 * t155 + (-t129 * t136 + t138 * t153) * t131;
t108 = t118 * t137 + t121 * t135;
t105 = 0.1e1 / t107;
t97 = 0.1e1 / (t106 * t167 ^ 2 + 0.1e1);
t91 = 0.1e1 / t94;
t90 = 0.1e1 / t150;
t87 = 0.1e1 / t89;
t86 = 0.1e1 / (0.1e1 + t161);
t85 = (-t105 * t141 - t117 * t160) * t97 * t135;
t84 = (t102 * t105 - t108 * t160) * t97;
t83 = t150 * t90;
t1 = [-t103 * t105 * t97, 0, t85, t84, 0, 0; (t167 * t87 - (-t95 + (-t105 * t163 + t95) * t97) * t161) * t86, 0 (-t113 * t135 * t87 - (t145 * t85 + (t117 * t96 - t141 * t95) * t135) * t162) * t86 (t104 * t87 - (t102 * t95 + t96 * t108 + t145 * t84) * t162) * t86, 0, 0; ((t102 * t126 - t127 * t141) * t91 - (t102 * t127 + t126 * t141) * t164) * t90, 0 ((-t114 * t127 - t126 * t159) * t91 - (t114 * t126 - t127 * t159) * t164) * t90 (-t126 * t91 + t127 * t164) * t90 * t103, t83, t83;];
Ja_rot  = t1;
