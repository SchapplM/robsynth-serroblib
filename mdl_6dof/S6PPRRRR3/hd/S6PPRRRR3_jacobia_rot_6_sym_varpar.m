% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PPRRRR3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobia_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:43:56
% EndTime: 2019-02-26 19:43:56
% DurationCPUTime: 0.63s
% Computational Cost: add. (4253->74), mult. (12284->171), div. (95->9), fcn. (16637->21), ass. (0->89)
t147 = sin(pkin(14));
t148 = sin(pkin(13));
t152 = cos(pkin(14));
t153 = cos(pkin(13));
t156 = cos(pkin(6));
t179 = t153 * t156;
t142 = -t148 * t147 + t152 * t179;
t143 = t147 * t179 + t148 * t152;
t160 = sin(qJ(3));
t164 = cos(qJ(3));
t150 = sin(pkin(7));
t151 = sin(pkin(6));
t182 = t150 * t151;
t175 = t153 * t182;
t155 = cos(pkin(7));
t176 = t155 * t160;
t134 = t142 * t176 + t143 * t164 - t160 * t175;
t159 = sin(qJ(4));
t163 = cos(qJ(4));
t149 = sin(pkin(8));
t154 = cos(pkin(8));
t167 = -t143 * t160 + (t142 * t155 - t175) * t164;
t180 = t151 * t155;
t172 = -t142 * t150 - t153 * t180;
t165 = t149 * t172 + t154 * t167;
t121 = t134 * t163 + t159 * t165;
t158 = sin(qJ(5));
t162 = cos(qJ(5));
t166 = -t149 * t167 + t154 * t172;
t111 = t121 * t158 - t162 * t166;
t181 = t150 * t156;
t140 = t160 * t181 + (t147 * t164 + t152 * t176) * t151;
t139 = t164 * t181 + (t152 * t155 * t164 - t147 * t160) * t151;
t141 = -t152 * t182 + t156 * t155;
t173 = t139 * t154 + t141 * t149;
t132 = t140 * t163 + t159 * t173;
t137 = -t139 * t149 + t141 * t154;
t124 = t132 * t158 - t137 * t162;
t110 = atan2(-t111, t124);
t107 = sin(t110);
t108 = cos(t110);
t101 = -t107 * t111 + t108 * t124;
t100 = 0.1e1 / t101 ^ 2;
t184 = t148 * t156;
t145 = -t147 * t184 + t153 * t152;
t144 = -t153 * t147 - t152 * t184;
t170 = t144 * t155 + t148 * t182;
t135 = -t145 * t160 + t164 * t170;
t136 = t145 * t164 + t160 * t170;
t171 = -t144 * t150 + t148 * t180;
t169 = t171 * t149;
t123 = t136 * t163 + (t135 * t154 + t169) * t159;
t168 = -t135 * t149 + t154 * t171;
t114 = t123 * t158 - t162 * t168;
t188 = t100 * t114;
t115 = t123 * t162 + t158 * t168;
t177 = t154 * t163;
t122 = -t135 * t177 + t136 * t159 - t163 * t169;
t157 = sin(qJ(6));
t161 = cos(qJ(6));
t106 = t115 * t161 + t122 * t157;
t104 = 0.1e1 / t106 ^ 2;
t105 = t115 * t157 - t122 * t161;
t187 = t104 * t105;
t119 = 0.1e1 / t124 ^ 2;
t186 = t111 * t119;
t185 = t122 * t162;
t183 = t149 * t162;
t178 = t154 * t159;
t174 = t105 ^ 2 * t104 + 0.1e1;
t131 = -t140 * t159 + t163 * t173;
t128 = (t139 * t163 - t140 * t178) * t158 - t140 * t183;
t127 = t135 * t163 - t136 * t178;
t126 = t135 * t159 + t136 * t177;
t125 = t132 * t162 + t137 * t158;
t120 = -t134 * t159 + t163 * t165;
t118 = 0.1e1 / t124;
t117 = t136 * t149 * t158 + t127 * t162;
t116 = (-t134 * t178 + t163 * t167) * t158 - t134 * t183;
t113 = t121 * t162 + t158 * t166;
t109 = 0.1e1 / (t111 ^ 2 * t119 + 0.1e1);
t103 = 0.1e1 / t106;
t102 = 0.1e1 / t174;
t99 = 0.1e1 / t101;
t98 = 0.1e1 / (t114 ^ 2 * t100 + 0.1e1);
t97 = (-t118 * t120 + t131 * t186) * t158 * t109;
t96 = (-t116 * t118 + t128 * t186) * t109;
t95 = (-t113 * t118 + t125 * t186) * t109;
t1 = [0, 0, t96, t97, t95, 0; 0, 0 ((t127 * t158 - t136 * t183) * t99 - ((-t111 * t96 + t128) * t108 + (-t124 * t96 - t116) * t107) * t188) * t98 (-t122 * t158 * t99 - ((-t111 * t97 + t131 * t158) * t108 + (-t120 * t158 - t124 * t97) * t107) * t188) * t98 (t115 * t99 - ((-t111 * t95 + t125) * t108 + (-t124 * t95 - t113) * t107) * t188) * t98, 0; 0, 0 ((t117 * t157 - t126 * t161) * t103 - (t117 * t161 + t126 * t157) * t187) * t102 ((-t123 * t161 - t157 * t185) * t103 - (t123 * t157 - t161 * t185) * t187) * t102 (-t157 * t103 + t161 * t187) * t114 * t102, t174 * t102;];
Ja_rot  = t1;
