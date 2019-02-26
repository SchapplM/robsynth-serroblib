% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRR6_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_jacobia_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:22:01
% EndTime: 2019-02-26 20:22:02
% DurationCPUTime: 0.53s
% Computational Cost: add. (2578->75), mult. (7606->169), div. (95->9), fcn. (10284->19), ass. (0->89)
t145 = sin(pkin(8));
t146 = sin(pkin(7));
t148 = cos(pkin(8));
t149 = cos(pkin(7));
t150 = cos(pkin(6));
t154 = sin(qJ(2));
t188 = cos(pkin(14));
t168 = t188 * t154;
t144 = sin(pkin(14));
t158 = cos(qJ(2));
t181 = t144 * t158;
t140 = t150 * t168 + t181;
t153 = sin(qJ(3));
t157 = cos(qJ(3));
t167 = t188 * t158;
t182 = t144 * t154;
t163 = -t150 * t167 + t182;
t147 = sin(pkin(6));
t169 = t147 * t188;
t162 = -t146 * t169 - t149 * t163;
t160 = -t140 * t153 + t157 * t162;
t190 = t160 * t148 + (t146 * t163 - t149 * t169) * t145;
t142 = -t150 * t182 + t167;
t141 = -t150 * t181 - t168;
t179 = t146 * t147;
t164 = t141 * t149 + t144 * t179;
t130 = -t142 * t153 + t157 * t164;
t156 = cos(qJ(4));
t177 = t148 * t156;
t138 = t144 * t147 * t149 - t141 * t146;
t183 = t138 * t145;
t131 = t142 * t157 + t153 * t164;
t152 = sin(qJ(4));
t184 = t131 * t152;
t113 = -t130 * t177 - t156 * t183 + t184;
t129 = t140 * t157 + t153 * t162;
t110 = t129 * t152 - t190 * t156;
t172 = t154 * t157;
t174 = t153 * t158;
t178 = t146 * t150;
t137 = t153 * t178 + (t149 * t174 + t172) * t147;
t171 = t157 * t158;
t173 = t154 * t153;
t136 = t157 * t178 + (t149 * t171 - t173) * t147;
t166 = t136 * t148 + (t150 * t149 - t158 * t179) * t145;
t121 = t137 * t152 - t156 * t166;
t109 = atan2(-t110, t121);
t106 = sin(t109);
t107 = cos(t109);
t100 = -t106 * t110 + t107 * t121;
t99 = 0.1e1 / t100 ^ 2;
t189 = t113 * t99;
t114 = t131 * t156 + (t130 * t148 + t183) * t152;
t123 = -t130 * t145 + t138 * t148;
t151 = sin(qJ(5));
t155 = cos(qJ(5));
t105 = t114 * t155 + t123 * t151;
t103 = 0.1e1 / t105 ^ 2;
t104 = t114 * t151 - t123 * t155;
t187 = t103 * t104;
t120 = 0.1e1 / t121 ^ 2;
t186 = t110 * t120;
t185 = t131 * t145;
t180 = t146 * t145;
t176 = t149 * t153;
t175 = t149 * t157;
t170 = t104 ^ 2 * t103 + 0.1e1;
t132 = -t141 * t153 - t142 * t175;
t165 = t132 * t148 + t142 * t180;
t133 = t141 * t157 - t142 * t176;
t126 = ((-t149 * t173 + t171) * t152 + (-(-t149 * t172 - t174) * t148 - t154 * t180) * t156) * t147;
t125 = t142 * t146 * t148 - t132 * t145;
t124 = t136 * t152 + t137 * t177;
t122 = t137 * t156 + t152 * t166;
t119 = 0.1e1 / t121;
t118 = t130 * t156 - t148 * t184;
t117 = t129 * t177 + t152 * t160;
t116 = t133 * t156 + t152 * t165;
t115 = (-t140 * t176 - t157 * t163) * t152 + (-(-t140 * t175 + t153 * t163) * t148 - t140 * t180) * t156;
t112 = t129 * t156 + t190 * t152;
t108 = 0.1e1 / (t110 ^ 2 * t120 + 0.1e1);
t102 = 0.1e1 / t105;
t101 = 0.1e1 / t170;
t98 = 0.1e1 / t100;
t97 = 0.1e1 / (t113 ^ 2 * t99 + 0.1e1);
t96 = (-t115 * t119 + t126 * t186) * t108;
t95 = (-t117 * t119 + t124 * t186) * t108;
t94 = (-t112 * t119 + t122 * t186) * t108;
t1 = [0, t96, t95, t94, 0, 0; 0 ((t133 * t152 - t156 * t165) * t98 - ((-t110 * t96 + t126) * t107 + (-t121 * t96 - t115) * t106) * t189) * t97 ((t130 * t152 + t131 * t177) * t98 - ((-t110 * t95 + t124) * t107 + (-t121 * t95 - t117) * t106) * t189) * t97 (t114 * t98 - ((-t110 * t94 + t122) * t107 + (-t121 * t94 - t112) * t106) * t189) * t97, 0, 0; 0 ((t116 * t151 - t125 * t155) * t102 - (t116 * t155 + t125 * t151) * t187) * t101 ((t118 * t151 - t155 * t185) * t102 - (t118 * t155 + t151 * t185) * t187) * t101 (-t151 * t102 + t155 * t187) * t113 * t101, t170 * t101, 0;];
Ja_rot  = t1;
