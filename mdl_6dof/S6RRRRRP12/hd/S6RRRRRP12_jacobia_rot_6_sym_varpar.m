% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRP12_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:46:22
% EndTime: 2019-02-26 22:46:23
% DurationCPUTime: 0.92s
% Computational Cost: add. (3414->77), mult. (9812->180), div. (137->9), fcn. (13393->17), ass. (0->89)
t174 = cos(pkin(6));
t179 = sin(qJ(1));
t183 = cos(qJ(2));
t191 = t179 * t183;
t178 = sin(qJ(2));
t184 = cos(qJ(1));
t192 = t178 * t184;
t165 = t174 * t192 + t191;
t177 = sin(qJ(3));
t182 = cos(qJ(3));
t189 = t183 * t184;
t194 = t178 * t179;
t164 = -t174 * t189 + t194;
t171 = sin(pkin(7));
t173 = cos(pkin(7));
t172 = sin(pkin(6));
t200 = t172 * t184;
t186 = t164 * t173 + t171 * t200;
t152 = -t165 * t182 + t186 * t177;
t159 = -t164 * t171 + t173 * t200;
t176 = sin(qJ(4));
t181 = cos(qJ(4));
t140 = t152 * t181 + t159 * t176;
t175 = sin(qJ(5));
t180 = cos(qJ(5));
t185 = t165 * t177 + t186 * t182;
t218 = t140 * t175 + t185 * t180;
t217 = t140 * t180 - t185 * t175;
t138 = t152 * t176 - t159 * t181;
t193 = t178 * t182;
t195 = t177 * t183;
t203 = t171 * t174;
t158 = t177 * t203 + (t173 * t195 + t193) * t172;
t163 = -t171 * t172 * t183 + t173 * t174;
t149 = t158 * t181 + t163 * t176;
t190 = t182 * t183;
t196 = t177 * t178;
t157 = -t182 * t203 + (-t173 * t190 + t196) * t172;
t135 = t149 * t175 - t157 * t180;
t122 = atan2(t218, t135);
t119 = sin(t122);
t120 = cos(t122);
t118 = t119 * t218 + t120 * t135;
t117 = 0.1e1 / t118 ^ 2;
t166 = -t174 * t191 - t192;
t167 = -t174 * t194 + t189;
t201 = t172 * t179;
t188 = t171 * t201;
t154 = t167 * t182 + (t166 * t173 + t188) * t177;
t161 = -t166 * t171 + t173 * t201;
t142 = t154 * t181 + t161 * t176;
t198 = t173 * t182;
t153 = -t166 * t198 + t167 * t177 - t182 * t188;
t206 = t153 * t180;
t129 = t142 * t175 - t206;
t212 = t117 * t129;
t211 = t117 * t129 ^ 2;
t210 = t120 * t218;
t130 = t142 * t180 + t153 * t175;
t125 = 0.1e1 / t130 ^ 2;
t141 = -t154 * t176 + t161 * t181;
t209 = t125 * t141 ^ 2;
t208 = t125 * t141;
t134 = 0.1e1 / t135 ^ 2;
t207 = t218 * t134;
t202 = t171 * t176;
t199 = t173 * t177;
t197 = t175 * t181;
t187 = -t119 * t135 + t210;
t156 = t166 * t182 - t167 * t199;
t155 = t166 * t177 + t167 * t198;
t148 = -t158 * t176 + t163 * t181;
t145 = (((-t173 * t196 + t190) * t181 + t178 * t202) * t175 - (t173 * t193 + t195) * t180) * t172;
t144 = t156 * t181 + t167 * t202;
t143 = -t157 * t197 - t158 * t180;
t136 = t149 * t180 + t157 * t175;
t133 = 0.1e1 / t135;
t132 = ((-t164 * t182 - t165 * t199) * t181 + t165 * t202) * t175 - (-t164 * t177 + t165 * t198) * t180;
t131 = t152 * t180 - t185 * t197;
t124 = 0.1e1 / t130;
t123 = 0.1e1 / (0.1e1 + t209);
t121 = 0.1e1 / (t134 * t218 ^ 2 + 0.1e1);
t116 = 0.1e1 / t118;
t115 = 0.1e1 / (0.1e1 + t211);
t114 = (-t133 * t138 - t148 * t207) * t175 * t121;
t113 = (-t132 * t133 - t145 * t207) * t121;
t112 = (-t131 * t133 - t143 * t207) * t121;
t111 = (t133 * t217 - t136 * t207) * t121;
t1 = [-t129 * t133 * t121, t113, t112, t114, t111, 0; (t218 * t116 - (-t119 + (-t133 * t210 + t119) * t121) * t211) * t115 ((t144 * t175 - t155 * t180) * t116 - (t187 * t113 - t119 * t132 + t120 * t145) * t212) * t115 ((-t153 * t197 - t154 * t180) * t116 - (t187 * t112 - t119 * t131 + t120 * t143) * t212) * t115 (t141 * t175 * t116 - ((-t119 * t138 + t120 * t148) * t175 + t187 * t114) * t212) * t115 (t130 * t116 - (t187 * t111 + t119 * t217 + t120 * t136) * t212) * t115, 0; (-t138 * t124 - t217 * t208) * t123 ((t167 * t171 * t181 - t156 * t176) * t124 - (t144 * t180 + t155 * t175) * t208) * t123 (t153 * t176 * t124 - (t154 * t175 - t181 * t206) * t208) * t123 (-t124 * t142 - t180 * t209) * t123, t129 * t123 * t208, 0;];
Ja_rot  = t1;
