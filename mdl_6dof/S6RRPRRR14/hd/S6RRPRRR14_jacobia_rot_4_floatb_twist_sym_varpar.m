% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR14_jacobia_rot_4_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobia_rot_4_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobia_rot_4_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:21
% EndTime: 2018-12-10 18:38:21
% DurationCPUTime: 0.40s
% Computational Cost: add. (3541->72), mult. (3682->127), div. (56->9), fcn. (3874->29), ass. (0->85)
t180 = pkin(8) + qJ(4);
t191 = sin(t180) / 0.2e1;
t158 = pkin(6) + qJ(2);
t190 = cos(t158) / 0.2e1;
t182 = pkin(6) - qJ(2);
t178 = cos(t182);
t143 = t178 / 0.2e1 + t190;
t168 = sin(qJ(2));
t169 = sin(qJ(1));
t172 = cos(qJ(1));
t129 = -t169 * t143 - t172 * t168;
t176 = sin(t182);
t179 = sin(t158) / 0.2e1;
t139 = t179 - t176 / 0.2e1;
t171 = cos(qJ(2));
t131 = t169 * t139 - t172 * t171;
t156 = pkin(7) + pkin(14);
t144 = sin(t156) / 0.2e1;
t157 = pkin(7) - pkin(14);
t152 = sin(t157);
t132 = t144 + t152 / 0.2e1;
t145 = cos(t157) / 0.2e1;
t153 = cos(t156);
t134 = t145 + t153 / 0.2e1;
t159 = sin(pkin(14));
t162 = sin(pkin(6));
t184 = t162 * t169;
t118 = t129 * t134 + t131 * t159 + t132 * t184;
t160 = sin(pkin(8));
t164 = cos(pkin(8));
t161 = sin(pkin(7));
t165 = cos(pkin(7));
t173 = -t129 * t161 + t165 * t184;
t109 = t118 * t160 - t173 * t164;
t126 = -t172 * t143 + t169 * t168;
t127 = t172 * t139 + t169 * t171;
t183 = t162 * t172;
t116 = t126 * t134 + t127 * t159 + t132 * t183;
t125 = -t126 * t161 + t165 * t183;
t108 = -t116 * t160 + t125 * t164;
t138 = t179 + t176 / 0.2e1;
t142 = t190 - t178 / 0.2e1;
t166 = cos(pkin(6));
t114 = -(t166 * t132 + t138 * t134 + t142 * t159) * t160 + (-t138 * t161 + t166 * t165) * t164;
t107 = atan2(t108, t114);
t104 = sin(t107);
t105 = cos(t107);
t98 = t104 * t108 + t105 * t114;
t97 = 0.1e1 / t98 ^ 2;
t189 = t109 ^ 2 * t97;
t133 = t144 - t152 / 0.2e1;
t135 = t145 - t153 / 0.2e1;
t163 = cos(pkin(14));
t119 = t129 * t133 - t131 * t163 + t135 * t184;
t181 = pkin(8) - qJ(4);
t175 = sin(t181);
t137 = t191 - t175 / 0.2e1;
t174 = cos(t180) / 0.2e1;
t177 = cos(t181);
t140 = t174 - t177 / 0.2e1;
t170 = cos(qJ(4));
t103 = t118 * t137 + t119 * t170 - t173 * t140;
t101 = 0.1e1 / t103 ^ 2;
t136 = t191 + t175 / 0.2e1;
t141 = t177 / 0.2e1 + t174;
t167 = sin(qJ(4));
t102 = -t118 * t141 + t119 * t167 - t173 * t136;
t188 = t101 * t102;
t187 = t102 ^ 2 * t101;
t186 = t131 * t161;
t185 = t161 * t164;
t122 = t129 * t163 + t131 * t133;
t121 = -t129 * t159 + t131 * t134;
t120 = -(t142 * t134 - t138 * t159) * t160 - t142 * t185;
t117 = t126 * t133 - t127 * t163 + t135 * t183;
t113 = 0.1e1 / t114 ^ 2;
t112 = 0.1e1 / t114;
t111 = (t126 * t159 - t127 * t134) * t160 - t127 * t185;
t106 = 0.1e1 / (t108 ^ 2 * t113 + 0.1e1);
t100 = 0.1e1 / t103;
t99 = 0.1e1 / (0.1e1 + t187);
t96 = 0.1e1 / t98;
t95 = 0.1e1 / (0.1e1 + t189);
t94 = (-t108 * t113 * t120 + t111 * t112) * t106;
t1 = [t109 * t112 * t106, t94, 0, 0, 0, 0; (t108 * t96 + (t104 + (t105 * t108 * t112 - t104) * t106) * t189) * t95 ((-t121 * t160 - t131 * t185) * t96 + ((t108 * t94 + t120) * t105 + (-t114 * t94 + t111) * t104) * t109 * t97) * t95, 0, 0, 0, 0; ((-t116 * t141 + t117 * t167 - t125 * t136) * t100 - (t116 * t137 + t117 * t170 - t125 * t140) * t188) * t99 ((-t121 * t141 + t122 * t167 + t136 * t186) * t100 - (t121 * t137 + t122 * t170 + t140 * t186) * t188) * t99, 0 (t103 * t100 + t187) * t99, 0, 0;];
Ja_rot  = t1;
