% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR14_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobia_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:55:45
% EndTime: 2019-02-26 22:55:46
% DurationCPUTime: 0.70s
% Computational Cost: add. (2431->68), mult. (7151->158), div. (85->9), fcn. (9681->19), ass. (0->84)
t161 = cos(pkin(6));
t167 = cos(qJ(2));
t203 = sin(qJ(1));
t183 = t203 * t167;
t164 = sin(qJ(2));
t168 = cos(qJ(1));
t186 = t168 * t164;
t151 = t161 * t186 + t183;
t154 = sin(pkin(14));
t158 = cos(pkin(14));
t184 = t203 * t164;
t187 = t167 * t168;
t150 = -t161 * t187 + t184;
t156 = sin(pkin(7));
t160 = cos(pkin(7));
t157 = sin(pkin(6));
t191 = t157 * t168;
t177 = t150 * t160 + t156 * t191;
t140 = -t151 * t158 + t177 * t154;
t163 = sin(qJ(4));
t166 = cos(qJ(4));
t139 = t151 * t154 + t177 * t158;
t148 = -t150 * t156 + t160 * t191;
t155 = sin(pkin(8));
t159 = cos(pkin(8));
t204 = -t139 * t159 - t148 * t155;
t122 = -t140 * t166 + t204 * t163;
t209 = t140 * t163 + t204 * t166;
t152 = -t161 * t184 + t187;
t176 = t161 * t183 + t186;
t185 = t157 * t203;
t172 = t156 * t185 - t176 * t160;
t170 = t152 * t154 - t172 * t158;
t173 = t176 * t156 + t160 * t185;
t205 = -t173 * t155 + t170 * t159;
t189 = t160 * t167;
t192 = t156 * t161;
t147 = t157 * t158 * t164 + (t157 * t189 + t192) * t154;
t179 = (t158 * t192 + (-t154 * t164 + t158 * t189) * t157) * t159 + (-t156 * t157 * t167 + t160 * t161) * t155;
t130 = t147 * t163 - t179 * t166;
t119 = atan2(t209, t130);
t112 = sin(t119);
t113 = cos(t119);
t110 = t112 * t209 + t113 * t130;
t109 = 0.1e1 / t110 ^ 2;
t141 = t152 * t158 + t172 * t154;
t124 = t141 * t163 + t205 * t166;
t202 = t109 * t124;
t201 = t109 * t124 ^ 2;
t200 = t113 * t209;
t125 = t141 * t166 - t205 * t163;
t133 = t170 * t155 + t173 * t159;
t162 = sin(qJ(5));
t165 = cos(qJ(5));
t118 = t125 * t165 + t133 * t162;
t115 = 0.1e1 / t118 ^ 2;
t117 = t125 * t162 - t133 * t165;
t199 = t115 * t117;
t129 = 0.1e1 / t130 ^ 2;
t198 = t209 * t129;
t194 = t154 * t160;
t193 = t156 * t155;
t190 = t158 * t160;
t188 = t164 * t160;
t182 = t115 * t117 ^ 2 + 0.1e1;
t181 = -t112 * t130 + t200;
t142 = -t152 * t190 + t176 * t154;
t178 = t142 * t159 + t152 * t193;
t143 = -t152 * t194 - t176 * t158;
t135 = t152 * t156 * t159 - t142 * t155;
t134 = ((-t154 * t188 + t158 * t167) * t163 + (-(-t154 * t167 - t158 * t188) * t159 - t164 * t193) * t166) * t157;
t132 = -t139 * t155 + t148 * t159;
t131 = t147 * t166 + t179 * t163;
t128 = 0.1e1 / t130;
t127 = t143 * t166 + t178 * t163;
t126 = (-t150 * t158 - t151 * t194) * t163 + (-(t150 * t154 - t151 * t190) * t159 - t151 * t193) * t166;
t116 = 0.1e1 / (t129 * t209 ^ 2 + 0.1e1);
t114 = 0.1e1 / t118;
t111 = 0.1e1 / t182;
t108 = 0.1e1 / t110;
t107 = 0.1e1 / (0.1e1 + t201);
t106 = (-t126 * t128 - t134 * t198) * t116;
t105 = (-t122 * t128 - t131 * t198) * t116;
t1 = [-t124 * t128 * t116, t106, 0, t105, 0, 0; (t209 * t108 - (-t112 + (-t128 * t200 + t112) * t116) * t201) * t107 ((t143 * t163 - t178 * t166) * t108 - (t181 * t106 - t112 * t126 + t113 * t134) * t202) * t107, 0 (t125 * t108 - (t181 * t105 - t112 * t122 + t113 * t131) * t202) * t107, 0, 0; ((-t122 * t162 - t132 * t165) * t114 - (-t122 * t165 + t132 * t162) * t199) * t111 ((t127 * t162 - t135 * t165) * t114 - (t127 * t165 + t135 * t162) * t199) * t111, 0 (-t114 * t162 + t165 * t199) * t124 * t111, t182 * t111, 0;];
Ja_rot  = t1;
