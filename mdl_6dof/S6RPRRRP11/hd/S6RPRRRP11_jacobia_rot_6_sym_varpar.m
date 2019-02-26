% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRP11_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:13:32
% EndTime: 2019-02-26 21:13:32
% DurationCPUTime: 0.41s
% Computational Cost: add. (1440->49), mult. (4121->118), div. (85->9), fcn. (5660->17), ass. (0->69)
t130 = cos(pkin(6));
t128 = cos(pkin(12));
t163 = sin(qJ(1));
t145 = t163 * t128;
t125 = sin(pkin(12));
t137 = cos(qJ(1));
t150 = t137 * t125;
t122 = t130 * t150 + t145;
t133 = sin(qJ(3));
t136 = cos(qJ(3));
t146 = t163 * t125;
t149 = t137 * t128;
t121 = -t130 * t149 + t146;
t126 = sin(pkin(7));
t129 = cos(pkin(7));
t127 = sin(pkin(6));
t152 = t127 * t137;
t142 = t121 * t129 + t126 * t152;
t111 = -t122 * t136 + t142 * t133;
t118 = -t121 * t126 + t129 * t152;
t132 = sin(qJ(4));
t135 = cos(qJ(4));
t165 = t111 * t132 - t118 * t135;
t101 = t111 * t135 + t118 * t132;
t141 = t130 * t145 + t150;
t147 = t127 * t163;
t164 = -t126 * t147 + t141 * t129;
t123 = -t130 * t146 + t149;
t113 = t123 * t136 - t164 * t133;
t138 = t141 * t126 + t129 * t147;
t103 = t113 * t135 + t138 * t132;
t112 = t123 * t133 + t164 * t136;
t131 = sin(qJ(5));
t134 = cos(qJ(5));
t93 = t103 * t134 + t112 * t131;
t91 = 0.1e1 / t93 ^ 2;
t92 = t103 * t131 - t112 * t134;
t162 = t91 * t92;
t151 = t128 * t129;
t153 = t126 * t130;
t117 = t133 * t153 + (t125 * t136 + t133 * t151) * t127;
t120 = -t127 * t128 * t126 + t130 * t129;
t106 = t117 * t132 - t120 * t135;
t97 = atan2(t165, t106);
t95 = cos(t97);
t161 = t95 * t165;
t102 = t113 * t132 - t138 * t135;
t94 = sin(t97);
t88 = t95 * t106 + t165 * t94;
t87 = 0.1e1 / t88 ^ 2;
t160 = t102 * t87;
t159 = t102 ^ 2 * t87;
t105 = 0.1e1 / t106 ^ 2;
t158 = t105 * t165;
t157 = t112 * t135;
t148 = t92 ^ 2 * t91 + 0.1e1;
t143 = -t106 * t94 + t161;
t139 = -t122 * t133 - t142 * t136;
t116 = t136 * t153 + (-t125 * t133 + t136 * t151) * t127;
t107 = t117 * t135 + t120 * t132;
t104 = 0.1e1 / t106;
t96 = 0.1e1 / (t105 * t165 ^ 2 + 0.1e1);
t90 = 0.1e1 / t93;
t89 = 0.1e1 / t148;
t86 = 0.1e1 / t88;
t85 = 0.1e1 / (0.1e1 + t159);
t84 = (-t104 * t139 - t116 * t158) * t96 * t132;
t83 = (t101 * t104 - t107 * t158) * t96;
t1 = [-t102 * t104 * t96, 0, t84, t83, 0, 0; (t165 * t86 - (-t94 + (-t104 * t161 + t94) * t96) * t159) * t85, 0 (-t112 * t132 * t86 - (t143 * t84 + (t116 * t95 - t139 * t94) * t132) * t160) * t85 (t103 * t86 - (t101 * t94 + t95 * t107 + t143 * t83) * t160) * t85, 0, 0; ((t101 * t131 - t134 * t139) * t90 - (t101 * t134 + t131 * t139) * t162) * t89, 0 ((-t113 * t134 - t131 * t157) * t90 - (t113 * t131 - t134 * t157) * t162) * t89 (-t131 * t90 + t134 * t162) * t89 * t102, t148 * t89, 0;];
Ja_rot  = t1;
