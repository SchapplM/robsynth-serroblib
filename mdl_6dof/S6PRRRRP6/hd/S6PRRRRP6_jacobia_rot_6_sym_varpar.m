% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRP6_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:18:02
% EndTime: 2019-02-26 20:18:03
% DurationCPUTime: 0.60s
% Computational Cost: add. (2847->74), mult. (8194->173), div. (117->9), fcn. (11175->17), ass. (0->85)
t153 = sin(pkin(12));
t156 = cos(pkin(12));
t166 = cos(qJ(2));
t158 = cos(pkin(6));
t162 = sin(qJ(2));
t177 = t158 * t162;
t148 = t153 * t166 + t156 * t177;
t161 = sin(qJ(3));
t165 = cos(qJ(3));
t176 = t158 * t166;
t147 = -t153 * t162 + t156 * t176;
t157 = cos(pkin(7));
t154 = sin(pkin(7));
t155 = sin(pkin(6));
t183 = t154 * t155;
t168 = -t147 * t157 + t156 * t183;
t167 = t148 * t161 + t168 * t165;
t134 = t148 * t165 - t168 * t161;
t180 = t155 * t157;
t143 = -t147 * t154 - t156 * t180;
t160 = sin(qJ(4));
t164 = cos(qJ(4));
t127 = t134 * t164 + t143 * t160;
t159 = sin(qJ(5));
t163 = cos(qJ(5));
t114 = t127 * t159 - t167 * t163;
t172 = t162 * t165;
t173 = t161 * t166;
t182 = t154 * t158;
t142 = t161 * t182 + (t157 * t173 + t172) * t155;
t146 = t158 * t157 - t166 * t183;
t138 = t142 * t164 + t146 * t160;
t171 = t165 * t166;
t174 = t161 * t162;
t141 = -t165 * t182 + (-t157 * t171 + t174) * t155;
t124 = t138 * t159 - t141 * t163;
t111 = atan2(-t114, t124);
t107 = sin(t111);
t108 = cos(t111);
t106 = -t107 * t114 + t108 * t124;
t105 = 0.1e1 / t106 ^ 2;
t149 = -t153 * t176 - t156 * t162;
t150 = -t153 * t177 + t156 * t166;
t170 = t153 * t183;
t136 = t150 * t165 + (t149 * t157 + t170) * t161;
t144 = -t149 * t154 + t153 * t180;
t129 = t136 * t164 + t144 * t160;
t178 = t157 * t165;
t135 = -t149 * t178 + t150 * t161 - t165 * t170;
t185 = t135 * t163;
t117 = t129 * t159 - t185;
t189 = t105 * t117;
t118 = t129 * t163 + t135 * t159;
t113 = 0.1e1 / t118 ^ 2;
t128 = -t136 * t160 + t144 * t164;
t188 = t113 * t128;
t122 = 0.1e1 / t124 ^ 2;
t187 = t114 * t122;
t186 = t128 ^ 2 * t113;
t181 = t154 * t160;
t179 = t157 * t161;
t175 = t159 * t164;
t169 = -t107 * t124 - t108 * t114;
t140 = t149 * t165 - t150 * t179;
t139 = t149 * t161 + t150 * t178;
t137 = -t142 * t160 + t146 * t164;
t132 = (((-t157 * t174 + t171) * t164 + t162 * t181) * t159 - (t157 * t172 + t173) * t163) * t155;
t131 = t140 * t164 + t150 * t181;
t130 = -t141 * t175 - t142 * t163;
t126 = -t134 * t160 + t143 * t164;
t125 = t138 * t163 + t141 * t159;
t121 = 0.1e1 / t124;
t120 = ((t147 * t165 - t148 * t179) * t164 + t148 * t181) * t159 - (t147 * t161 + t148 * t178) * t163;
t119 = -t134 * t163 - t167 * t175;
t116 = t127 * t163 + t167 * t159;
t112 = 0.1e1 / t118;
t110 = 0.1e1 / (0.1e1 + t186);
t109 = 0.1e1 / (t114 ^ 2 * t122 + 0.1e1);
t104 = 0.1e1 / t106;
t103 = 0.1e1 / (t117 ^ 2 * t105 + 0.1e1);
t102 = (-t121 * t126 + t137 * t187) * t159 * t109;
t101 = (-t120 * t121 + t132 * t187) * t109;
t100 = (-t119 * t121 + t130 * t187) * t109;
t99 = (-t116 * t121 + t125 * t187) * t109;
t1 = [0, t101, t100, t102, t99, 0; 0 ((t131 * t159 - t139 * t163) * t104 - (t169 * t101 - t107 * t120 + t108 * t132) * t189) * t103 ((-t135 * t175 - t136 * t163) * t104 - (t169 * t100 - t107 * t119 + t108 * t130) * t189) * t103 (t128 * t159 * t104 - ((-t107 * t126 + t108 * t137) * t159 + t169 * t102) * t189) * t103 (t118 * t104 - ((-t114 * t99 + t125) * t108 + (-t124 * t99 - t116) * t107) * t189) * t103, 0; 0 ((t150 * t154 * t164 - t140 * t160) * t112 - (t131 * t163 + t139 * t159) * t188) * t110 (t135 * t160 * t112 - (t136 * t159 - t164 * t185) * t188) * t110 (-t112 * t129 - t163 * t186) * t110, t117 * t110 * t188, 0;];
Ja_rot  = t1;
