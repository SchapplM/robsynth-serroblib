% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR13_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_jacobia_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:04
% EndTime: 2019-02-26 22:23:05
% DurationCPUTime: 0.58s
% Computational Cost: add. (2487->68), mult. (5494->156), div. (115->9), fcn. (7543->17), ass. (0->80)
t150 = cos(pkin(6));
t156 = cos(qJ(2));
t186 = sin(qJ(1));
t165 = t186 * t156;
t153 = sin(qJ(2));
t157 = cos(qJ(1));
t168 = t157 * t153;
t140 = t150 * t168 + t165;
t152 = sin(qJ(3));
t155 = cos(qJ(3));
t166 = t186 * t153;
t169 = t156 * t157;
t139 = -t150 * t169 + t166;
t147 = sin(pkin(7));
t149 = cos(pkin(7));
t148 = sin(pkin(6));
t173 = t148 * t157;
t162 = t139 * t149 + t147 * t173;
t126 = -t140 * t155 + t162 * t152;
t136 = -t139 * t147 + t149 * t173;
t146 = pkin(13) + qJ(5);
t144 = sin(t146);
t145 = cos(t146);
t187 = t126 * t144 - t136 * t145;
t114 = t126 * t145 + t136 * t144;
t172 = t149 * t152;
t174 = t147 * t150;
t135 = t152 * t174 + (t153 * t155 + t156 * t172) * t148;
t138 = -t147 * t148 * t156 + t149 * t150;
t121 = t135 * t144 - t138 * t145;
t108 = atan2(t187, t121);
t103 = sin(t108);
t104 = cos(t108);
t101 = t103 * t187 + t104 * t121;
t100 = 0.1e1 / t101 ^ 2;
t161 = t150 * t165 + t168;
t167 = t148 * t186;
t163 = t147 * t167;
t141 = -t150 * t166 + t169;
t176 = t141 * t155;
t128 = t176 + (-t161 * t149 + t163) * t152;
t158 = t161 * t147 + t149 * t167;
t115 = t128 * t144 - t158 * t145;
t185 = t100 * t115;
t184 = t100 * t115 ^ 2;
t116 = t128 * t145 + t158 * t144;
t154 = cos(qJ(6));
t160 = t161 * t155;
t127 = t141 * t152 + t149 * t160 - t155 * t163;
t151 = sin(qJ(6));
t181 = t127 * t151;
t110 = t116 * t154 + t181;
t107 = 0.1e1 / t110 ^ 2;
t180 = t127 * t154;
t109 = t116 * t151 - t180;
t183 = t107 * t109;
t120 = 0.1e1 / t121 ^ 2;
t182 = t187 * t120;
t175 = t145 * t147;
t171 = t152 * t153;
t170 = t155 * t156;
t164 = t107 * t109 ^ 2 + 0.1e1;
t159 = -t140 * t152 - t162 * t155;
t134 = t155 * t174 + (t149 * t170 - t171) * t148;
t131 = -t141 * t172 - t160;
t130 = t149 * t176 - t161 * t152;
t129 = ((-t149 * t171 + t170) * t144 - t153 * t175) * t148;
t122 = t135 * t145 + t138 * t144;
t119 = 0.1e1 / t121;
t118 = t141 * t144 * t147 + t131 * t145;
t117 = (-t139 * t155 - t140 * t172) * t144 - t140 * t175;
t106 = 0.1e1 / t110;
t105 = 0.1e1 / (t120 * t187 ^ 2 + 0.1e1);
t102 = 0.1e1 / t164;
t99 = 0.1e1 / t101;
t98 = 0.1e1 / (0.1e1 + t184);
t97 = (-t119 * t159 - t134 * t182) * t144 * t105;
t96 = (-t117 * t119 - t129 * t182) * t105;
t95 = (t114 * t119 - t122 * t182) * t105;
t1 = [-t115 * t119 * t105, t96, t97, 0, t95, 0; (t187 * t99 - (-t103 + (-t104 * t119 * t187 + t103) * t105) * t184) * t98 ((t131 * t144 - t141 * t175) * t99 - ((t187 * t96 + t129) * t104 + (-t121 * t96 - t117) * t103) * t185) * t98 (-t127 * t144 * t99 - ((t134 * t144 + t187 * t97) * t104 + (-t121 * t97 - t144 * t159) * t103) * t185) * t98, 0 (t116 * t99 - ((t187 * t95 + t122) * t104 + (-t121 * t95 + t114) * t103) * t185) * t98, 0; ((t114 * t151 - t154 * t159) * t106 - (t114 * t154 + t151 * t159) * t183) * t102 ((t118 * t151 - t130 * t154) * t106 - (t118 * t154 + t130 * t151) * t183) * t102 ((-t128 * t154 - t145 * t181) * t106 - (t128 * t151 - t145 * t180) * t183) * t102, 0 (-t106 * t151 + t154 * t183) * t115 * t102, t164 * t102;];
Ja_rot  = t1;
