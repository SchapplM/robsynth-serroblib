% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR15_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:47
% EndTime: 2019-02-26 22:38:48
% DurationCPUTime: 0.60s
% Computational Cost: add. (1916->67), mult. (5494->159), div. (115->9), fcn. (7543->17), ass. (0->76)
t144 = sin(qJ(2));
t145 = sin(qJ(1));
t149 = cos(qJ(2));
t150 = cos(qJ(1));
t169 = cos(pkin(6));
t154 = t150 * t169;
t132 = t144 * t154 + t145 * t149;
t143 = sin(qJ(3));
t148 = cos(qJ(3));
t131 = t145 * t144 - t149 * t154;
t138 = sin(pkin(7));
t140 = cos(pkin(7));
t139 = sin(pkin(6));
t162 = t139 * t150;
t152 = t131 * t140 + t138 * t162;
t118 = -t132 * t148 + t152 * t143;
t126 = -t131 * t138 + t140 * t162;
t142 = sin(qJ(4));
t147 = cos(qJ(4));
t106 = t118 * t142 - t126 * t147;
t175 = t118 * t147 + t126 * t142;
t155 = t145 * t169;
t133 = -t150 * t144 - t149 * t155;
t134 = -t144 * t155 + t150 * t149;
t163 = t139 * t145;
t156 = t138 * t163;
t120 = t134 * t148 + (t133 * t140 + t156) * t143;
t128 = -t133 * t138 + t140 * t163;
t107 = t120 * t142 - t128 * t147;
t141 = sin(qJ(6));
t160 = t140 * t148;
t119 = -t133 * t160 + t134 * t143 - t148 * t156;
t146 = cos(qJ(6));
t166 = t119 * t146;
t98 = t107 * t141 + t166;
t96 = 0.1e1 / t98 ^ 2;
t167 = t119 * t141;
t97 = -t107 * t146 + t167;
t172 = t96 * t97;
t108 = t120 * t147 + t128 * t142;
t153 = t169 * t138;
t161 = t140 * t143;
t125 = t143 * t153 + (t144 * t148 + t149 * t161) * t139;
t130 = -t139 * t149 * t138 + t169 * t140;
t114 = t125 * t147 + t130 * t142;
t102 = atan2(t175, t114);
t100 = cos(t102);
t99 = sin(t102);
t93 = t100 * t114 + t175 * t99;
t92 = 0.1e1 / t93 ^ 2;
t171 = t108 * t92;
t170 = t108 ^ 2 * t92;
t112 = 0.1e1 / t114 ^ 2;
t168 = t175 * t112;
t164 = t138 * t142;
t159 = t143 * t144;
t158 = t148 * t149;
t157 = t97 ^ 2 * t96 + 0.1e1;
t151 = -t132 * t143 - t152 * t148;
t124 = t148 * t153 + (t140 * t158 - t159) * t139;
t123 = ((-t140 * t159 + t158) * t147 + t144 * t164) * t139;
t122 = t133 * t148 - t134 * t161;
t121 = t133 * t143 + t134 * t160;
t113 = -t125 * t142 + t130 * t147;
t111 = 0.1e1 / t114;
t110 = -t134 * t138 * t147 + t122 * t142;
t109 = (-t131 * t148 - t132 * t161) * t147 + t132 * t164;
t101 = 0.1e1 / (t112 * t175 ^ 2 + 0.1e1);
t95 = 0.1e1 / t98;
t94 = 0.1e1 / t157;
t91 = 0.1e1 / t93;
t90 = 0.1e1 / (0.1e1 + t170);
t89 = (-t111 * t151 - t124 * t168) * t147 * t101;
t88 = (-t109 * t111 - t123 * t168) * t101;
t87 = (-t106 * t111 - t113 * t168) * t101;
t1 = [-t108 * t111 * t101, t88, t89, t87, 0, 0; (t175 * t91 - (-t99 + (-t100 * t111 * t175 + t99) * t101) * t170) * t90 ((t122 * t147 + t134 * t164) * t91 - ((-t114 * t88 - t109) * t99 + (t175 * t88 + t123) * t100) * t171) * t90 (-t119 * t147 * t91 - ((-t114 * t89 - t147 * t151) * t99 + (t124 * t147 + t175 * t89) * t100) * t171) * t90 (-t107 * t91 - ((-t114 * t87 - t106) * t99 + (t175 * t87 + t113) * t100) * t171) * t90, 0, 0; ((-t106 * t146 + t141 * t151) * t95 - (t106 * t141 + t146 * t151) * t172) * t94 ((-t110 * t146 + t121 * t141) * t95 - (t110 * t141 + t121 * t146) * t172) * t94 ((t120 * t141 + t142 * t166) * t95 - (t120 * t146 - t142 * t167) * t172) * t94 (-t141 * t172 - t146 * t95) * t94 * t108, 0, t157 * t94;];
Ja_rot  = t1;
