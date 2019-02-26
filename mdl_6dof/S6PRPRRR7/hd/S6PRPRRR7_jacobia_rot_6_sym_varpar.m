% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function Ja_rot = S6PRPRRR7_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobia_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:57:30
% EndTime: 2019-02-26 19:57:31
% DurationCPUTime: 0.76s
% Computational Cost: add. (4256->82), mult. (12302->195), div. (95->9), fcn. (16658->21), ass. (0->98)
t166 = sin(pkin(8));
t171 = cos(pkin(8));
t165 = sin(pkin(13));
t170 = cos(pkin(13));
t181 = cos(qJ(2));
t173 = cos(pkin(6));
t177 = sin(qJ(2));
t198 = t173 * t177;
t162 = -t165 * t198 + t170 * t181;
t164 = sin(pkin(14));
t169 = cos(pkin(14));
t197 = t173 * t181;
t161 = -t165 * t197 - t170 * t177;
t172 = cos(pkin(7));
t167 = sin(pkin(7));
t168 = sin(pkin(6));
t206 = t167 * t168;
t189 = t161 * t172 + t165 * t206;
t186 = t162 * t164 - t189 * t169;
t202 = t168 * t172;
t190 = -t161 * t167 + t165 * t202;
t213 = -t190 * t166 + t186 * t171;
t160 = t165 * t181 + t170 * t198;
t159 = -t165 * t177 + t170 * t197;
t191 = t159 * t172 - t170 * t206;
t149 = t160 * t169 + t191 * t164;
t176 = sin(qJ(4));
t180 = cos(qJ(4));
t187 = -t160 * t164 + t191 * t169;
t192 = -t159 * t167 - t170 * t202;
t183 = t192 * t166 + t187 * t171;
t136 = t149 * t180 + t183 * t176;
t175 = sin(qJ(5));
t179 = cos(qJ(5));
t184 = -t187 * t166 + t192 * t171;
t124 = t136 * t175 - t184 * t179;
t199 = t172 * t181;
t204 = t167 * t173;
t156 = t168 * t177 * t169 + (t168 * t199 + t204) * t164;
t155 = t169 * t204 + (-t164 * t177 + t169 * t199) * t168;
t158 = t173 * t172 - t181 * t206;
t194 = t155 * t171 + t158 * t166;
t145 = t156 * t180 + t194 * t176;
t148 = -t155 * t166 + t158 * t171;
t133 = t145 * t175 - t148 * t179;
t123 = atan2(-t124, t133);
t120 = sin(t123);
t121 = cos(t123);
t114 = -t120 * t124 + t121 * t133;
t113 = 0.1e1 / t114 ^ 2;
t150 = t162 * t169 + t189 * t164;
t138 = t150 * t180 - t176 * t213;
t182 = t186 * t166 + t190 * t171;
t127 = t138 * t175 - t182 * t179;
t212 = t113 * t127;
t128 = t138 * t179 + t182 * t175;
t137 = t150 * t176 + t180 * t213;
t174 = sin(qJ(6));
t178 = cos(qJ(6));
t119 = t128 * t178 + t137 * t174;
t117 = 0.1e1 / t119 ^ 2;
t118 = t128 * t174 - t137 * t178;
t211 = t117 * t118;
t132 = 0.1e1 / t133 ^ 2;
t210 = t124 * t132;
t209 = t137 * t179;
t208 = t164 * t172;
t207 = t167 * t166;
t205 = t167 * t171;
t203 = t167 * t177;
t201 = t169 * t172;
t200 = t172 * t177;
t196 = t117 * t118 ^ 2 + 0.1e1;
t195 = -t120 * t133 - t121 * t124;
t152 = -t161 * t164 - t162 * t201;
t193 = t152 * t171 + t162 * t207;
t153 = t161 * t169 - t162 * t208;
t151 = -t159 * t164 - t160 * t201;
t146 = -t152 * t166 + t162 * t205;
t144 = -t156 * t176 + t194 * t180;
t141 = ((t171 * t175 * t176 + t166 * t179) * (-t164 * t181 - t169 * t200) + ((-t164 * t200 + t169 * t181) * t180 + t166 * t176 * t203) * t175 - t171 * t179 * t203) * t168;
t140 = t153 * t180 + t193 * t176;
t139 = t153 * t176 - t193 * t180;
t135 = -t149 * t176 + t183 * t180;
t134 = t145 * t179 + t148 * t175;
t131 = 0.1e1 / t133;
t130 = t140 * t179 + t146 * t175;
t129 = ((t159 * t169 - t160 * t208) * t180 + (t151 * t171 + t160 * t207) * t176) * t175 - (-t151 * t166 + t160 * t205) * t179;
t126 = t136 * t179 + t184 * t175;
t122 = 0.1e1 / (t124 ^ 2 * t132 + 0.1e1);
t116 = 0.1e1 / t119;
t115 = 0.1e1 / t196;
t112 = 0.1e1 / t114;
t111 = 0.1e1 / (t113 * t127 ^ 2 + 0.1e1);
t110 = (-t131 * t135 + t144 * t210) * t175 * t122;
t109 = (-t129 * t131 + t141 * t210) * t122;
t108 = (-t126 * t131 + t134 * t210) * t122;
t1 = [0, t109, 0, t110, t108, 0; 0 ((t140 * t175 - t146 * t179) * t112 - (t195 * t109 - t120 * t129 + t121 * t141) * t212) * t111, 0 (-t137 * t175 * t112 - ((-t120 * t135 + t121 * t144) * t175 + t195 * t110) * t212) * t111 (t128 * t112 - (t195 * t108 - t120 * t126 + t121 * t134) * t212) * t111, 0; 0 ((t130 * t174 - t139 * t178) * t116 - (t130 * t178 + t139 * t174) * t211) * t115, 0 ((-t138 * t178 - t174 * t209) * t116 - (t138 * t174 - t178 * t209) * t211) * t115 (-t174 * t116 + t178 * t211) * t127 * t115, t196 * t115;];
Ja_rot  = t1;
