% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRR10_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobia_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:52:47
% EndTime: 2019-02-26 22:52:48
% DurationCPUTime: 0.91s
% Computational Cost: add. (3237->82), mult. (9534->184), div. (115->9), fcn. (12904->19), ass. (0->95)
t177 = cos(pkin(6));
t181 = sin(qJ(2));
t226 = sin(qJ(1));
t202 = t226 * t181;
t185 = cos(qJ(2));
t186 = cos(qJ(1));
t206 = t185 * t186;
t166 = -t177 * t206 + t202;
t201 = t226 * t185;
t205 = t186 * t181;
t167 = t177 * t205 + t201;
t180 = sin(qJ(3));
t184 = cos(qJ(3));
t173 = sin(pkin(7));
t174 = sin(pkin(6));
t214 = t174 * t186;
t204 = t173 * t214;
t176 = cos(pkin(7));
t212 = t176 * t180;
t156 = t166 * t212 - t167 * t184 + t180 * t204;
t179 = sin(qJ(4));
t183 = cos(qJ(4));
t155 = (t166 * t176 + t204) * t184 + t167 * t180;
t164 = -t166 * t173 + t176 * t214;
t172 = sin(pkin(8));
t175 = cos(pkin(8));
t198 = t155 * t175 + t164 * t172;
t136 = t156 * t183 + t179 * t198;
t133 = -t156 * t179 + t183 * t198;
t208 = t181 * t184;
t209 = t180 * t185;
t215 = t173 * t177;
t163 = t180 * t215 + (t176 * t209 + t208) * t174;
t207 = t184 * t185;
t210 = t180 * t181;
t162 = t184 * t215 + (t176 * t207 - t210) * t174;
t197 = t162 * t175 + (-t173 * t174 * t185 + t176 * t177) * t172;
t145 = t163 * t179 - t183 * t197;
t132 = atan2(-t133, t145);
t127 = sin(t132);
t128 = cos(t132);
t123 = -t127 * t133 + t128 * t145;
t122 = 0.1e1 / t123 ^ 2;
t169 = -t177 * t202 + t206;
t168 = -t177 * t201 - t205;
t203 = t174 * t226;
t193 = t168 * t176 + t173 * t203;
t189 = t169 * t180 - t184 * t193;
t187 = t189 * t183;
t194 = -t168 * t173 + t176 * t203;
t191 = t194 * t172;
t157 = t169 * t184 + t180 * t193;
t219 = t157 * t179;
t137 = t175 * t187 - t183 * t191 + t219;
t225 = t122 * t137;
t224 = t122 * t137 ^ 2;
t138 = t157 * t183 + (-t175 * t189 + t191) * t179;
t148 = t172 * t189 + t175 * t194;
t178 = sin(qJ(5));
t182 = cos(qJ(5));
t130 = t138 * t182 + t148 * t178;
t126 = 0.1e1 / t130 ^ 2;
t129 = t138 * t178 - t148 * t182;
t223 = t126 * t129;
t222 = t128 * t133;
t144 = 0.1e1 / t145 ^ 2;
t221 = t133 * t144;
t220 = t157 * t172;
t216 = t173 * t172;
t213 = t175 * t183;
t211 = t176 * t184;
t200 = t126 * t129 ^ 2 + 0.1e1;
t199 = -t127 * t145 - t222;
t158 = -t168 * t180 - t169 * t211;
t196 = t158 * t175 + t169 * t216;
t159 = t168 * t184 - t169 * t212;
t151 = t169 * t173 * t175 - t158 * t172;
t150 = ((-t176 * t210 + t207) * t179 + (-(-t176 * t208 - t209) * t175 - t181 * t216) * t183) * t174;
t149 = t162 * t179 + t163 * t213;
t147 = -t155 * t172 + t164 * t175;
t146 = t163 * t183 + t179 * t197;
t143 = 0.1e1 / t145;
t142 = -t175 * t219 - t187;
t141 = -t155 * t179 - t156 * t213;
t140 = t159 * t183 + t179 * t196;
t139 = (-t166 * t184 - t167 * t212) * t179 + (-(t166 * t180 - t167 * t211) * t175 - t167 * t216) * t183;
t131 = 0.1e1 / (t133 ^ 2 * t144 + 0.1e1);
t125 = 0.1e1 / t130;
t124 = 0.1e1 / t200;
t121 = 0.1e1 / t123;
t120 = 0.1e1 / (0.1e1 + t224);
t119 = (-t139 * t143 + t150 * t221) * t131;
t118 = (-t141 * t143 + t149 * t221) * t131;
t117 = (t136 * t143 + t146 * t221) * t131;
t1 = [-t137 * t143 * t131, t119, t118, t117, 0, 0; (-t133 * t121 - (-t127 + (t143 * t222 + t127) * t131) * t224) * t120 ((t159 * t179 - t183 * t196) * t121 - (t119 * t199 - t127 * t139 + t128 * t150) * t225) * t120 ((t157 * t213 - t179 * t189) * t121 - (t118 * t199 - t127 * t141 + t128 * t149) * t225) * t120 (t138 * t121 - (t117 * t199 + t127 * t136 + t128 * t146) * t225) * t120, 0, 0; ((t136 * t178 - t147 * t182) * t125 - (t136 * t182 + t147 * t178) * t223) * t124 ((t140 * t178 - t151 * t182) * t125 - (t140 * t182 + t151 * t178) * t223) * t124 ((t142 * t178 - t182 * t220) * t125 - (t142 * t182 + t178 * t220) * t223) * t124 (-t125 * t178 + t182 * t223) * t137 * t124, t200 * t124, 0;];
Ja_rot  = t1;
