% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRR8_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_jacobia_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:51:23
% EndTime: 2019-02-26 22:51:24
% DurationCPUTime: 0.66s
% Computational Cost: add. (3129->66), mult. (6894->153), div. (145->9), fcn. (9467->17), ass. (0->84)
t167 = cos(pkin(6));
t173 = cos(qJ(2));
t205 = sin(qJ(1));
t183 = t205 * t173;
t170 = sin(qJ(2));
t174 = cos(qJ(1));
t187 = t174 * t170;
t157 = t167 * t187 + t183;
t169 = sin(qJ(3));
t172 = cos(qJ(3));
t184 = t205 * t170;
t186 = t174 * t173;
t156 = -t167 * t186 + t184;
t164 = sin(pkin(7));
t166 = cos(pkin(7));
t165 = sin(pkin(6));
t191 = t165 * t174;
t179 = t156 * t166 + t164 * t191;
t143 = -t157 * t172 + t179 * t169;
t153 = -t156 * t164 + t166 * t191;
t163 = qJ(4) + qJ(5);
t161 = sin(t163);
t162 = cos(t163);
t206 = t143 * t161 - t153 * t162;
t131 = t143 * t162 + t153 * t161;
t190 = t166 * t169;
t192 = t164 * t167;
t152 = t169 * t192 + (t170 * t172 + t173 * t190) * t165;
t155 = -t165 * t173 * t164 + t167 * t166;
t138 = t152 * t161 - t155 * t162;
t127 = atan2(t206, t138);
t120 = sin(t127);
t121 = cos(t127);
t118 = t120 * t206 + t121 * t138;
t117 = 0.1e1 / t118 ^ 2;
t178 = t167 * t183 + t187;
t185 = t165 * t205;
t181 = t164 * t185;
t158 = -t167 * t184 + t186;
t194 = t158 * t172;
t145 = t194 + (-t178 * t166 + t181) * t169;
t175 = t178 * t164 + t166 * t185;
t132 = t145 * t161 - t175 * t162;
t204 = t117 * t132;
t203 = t121 * t206;
t133 = t145 * t162 + t175 * t161;
t171 = cos(qJ(6));
t177 = t178 * t172;
t144 = t158 * t169 + t166 * t177 - t172 * t181;
t168 = sin(qJ(6));
t199 = t144 * t168;
t126 = t133 * t171 + t199;
t123 = 0.1e1 / t126 ^ 2;
t198 = t144 * t171;
t125 = t133 * t168 - t198;
t202 = t123 * t125;
t137 = 0.1e1 / t138 ^ 2;
t201 = t206 * t137;
t200 = t132 ^ 2 * t117;
t193 = t164 * t162;
t189 = t169 * t170;
t188 = t172 * t173;
t182 = t125 ^ 2 * t123 + 0.1e1;
t180 = -t120 * t138 + t203;
t176 = -t157 * t169 - t179 * t172;
t151 = t172 * t192 + (t166 * t188 - t189) * t165;
t148 = ((-t166 * t189 + t188) * t161 - t170 * t193) * t165;
t147 = -t158 * t190 - t177;
t146 = t166 * t194 - t178 * t169;
t139 = t152 * t162 + t155 * t161;
t136 = 0.1e1 / t138;
t135 = t158 * t164 * t161 + t147 * t162;
t134 = (-t156 * t172 - t157 * t190) * t161 - t157 * t193;
t124 = 0.1e1 / (t137 * t206 ^ 2 + 0.1e1);
t122 = 0.1e1 / t126;
t119 = 0.1e1 / t182;
t116 = 0.1e1 / t118;
t115 = 0.1e1 / (0.1e1 + t200);
t114 = (-t136 * t176 - t151 * t201) * t161 * t124;
t113 = (-t134 * t136 - t148 * t201) * t124;
t112 = (t131 * t136 - t139 * t201) * t124;
t111 = (-t168 * t122 + t171 * t202) * t132 * t119;
t110 = (t133 * t116 - (t180 * t112 + t120 * t131 + t121 * t139) * t204) * t115;
t1 = [-t132 * t136 * t124, t113, t114, t112, t112, 0; (t206 * t116 - (-t120 + (-t136 * t203 + t120) * t124) * t200) * t115 ((t147 * t161 - t158 * t193) * t116 - (t180 * t113 - t120 * t134 + t121 * t148) * t204) * t115 (-t144 * t161 * t116 - ((-t120 * t176 + t121 * t151) * t161 + t180 * t114) * t204) * t115, t110, t110, 0; ((t131 * t168 - t171 * t176) * t122 - (t131 * t171 + t168 * t176) * t202) * t119 ((t135 * t168 - t146 * t171) * t122 - (t135 * t171 + t146 * t168) * t202) * t119 ((-t145 * t171 - t162 * t199) * t122 - (t145 * t168 - t162 * t198) * t202) * t119, t111, t111, t182 * t119;];
Ja_rot  = t1;
