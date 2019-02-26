% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRP11
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
% Datum: 2019-02-26 22:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRP11_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:45:44
% EndTime: 2019-02-26 22:45:44
% DurationCPUTime: 0.59s
% Computational Cost: add. (1916->67), mult. (5494->157), div. (115->9), fcn. (7543->17), ass. (0->78)
t141 = cos(pkin(6));
t149 = cos(qJ(2));
t178 = sin(qJ(1));
t157 = t178 * t149;
t145 = sin(qJ(2));
t150 = cos(qJ(1));
t162 = t150 * t145;
t134 = t141 * t162 + t157;
t144 = sin(qJ(3));
t148 = cos(qJ(3));
t158 = t178 * t145;
t161 = t150 * t149;
t133 = -t141 * t161 + t158;
t138 = sin(pkin(7));
t140 = cos(pkin(7));
t139 = sin(pkin(6));
t166 = t139 * t150;
t155 = t133 * t140 + t138 * t166;
t120 = -t134 * t148 + t155 * t144;
t130 = -t133 * t138 + t140 * t166;
t143 = sin(qJ(4));
t147 = cos(qJ(4));
t179 = t120 * t143 - t130 * t147;
t108 = t120 * t147 + t130 * t143;
t154 = t141 * t157 + t162;
t159 = t139 * t178;
t156 = t138 * t159;
t135 = -t141 * t158 + t161;
t169 = t135 * t148;
t122 = t169 + (-t154 * t140 + t156) * t144;
t151 = t154 * t138 + t140 * t159;
t110 = t122 * t147 + t151 * t143;
t153 = t154 * t148;
t121 = t135 * t144 + t140 * t153 - t148 * t156;
t142 = sin(qJ(5));
t146 = cos(qJ(5));
t100 = t110 * t146 + t121 * t142;
t98 = 0.1e1 / t100 ^ 2;
t99 = t110 * t142 - t121 * t146;
t177 = t98 * t99;
t109 = t122 * t143 - t151 * t147;
t165 = t140 * t144;
t168 = t138 * t141;
t129 = t144 * t168 + (t145 * t148 + t149 * t165) * t139;
t132 = -t138 * t139 * t149 + t140 * t141;
t115 = t129 * t143 - t132 * t147;
t104 = atan2(t179, t115);
t101 = sin(t104);
t102 = cos(t104);
t95 = t101 * t179 + t102 * t115;
t94 = 0.1e1 / t95 ^ 2;
t176 = t109 * t94;
t175 = t109 ^ 2 * t94;
t114 = 0.1e1 / t115 ^ 2;
t174 = t179 * t114;
t173 = t121 * t147;
t167 = t138 * t147;
t164 = t144 * t145;
t163 = t148 * t149;
t160 = t98 * t99 ^ 2 + 0.1e1;
t152 = -t134 * t144 - t155 * t148;
t128 = t148 * t168 + (t140 * t163 - t164) * t139;
t125 = ((-t140 * t164 + t163) * t143 - t145 * t167) * t139;
t124 = -t135 * t165 - t153;
t123 = t140 * t169 - t154 * t144;
t116 = t129 * t147 + t132 * t143;
t113 = 0.1e1 / t115;
t112 = t135 * t138 * t143 + t124 * t147;
t111 = (-t133 * t148 - t134 * t165) * t143 - t134 * t167;
t103 = 0.1e1 / (t114 * t179 ^ 2 + 0.1e1);
t97 = 0.1e1 / t100;
t96 = 0.1e1 / t160;
t93 = 0.1e1 / t95;
t92 = 0.1e1 / (0.1e1 + t175);
t91 = (-t113 * t152 - t128 * t174) * t143 * t103;
t90 = (-t111 * t113 - t125 * t174) * t103;
t89 = (t108 * t113 - t116 * t174) * t103;
t1 = [-t109 * t113 * t103, t90, t91, t89, 0, 0; (t179 * t93 - (-t101 + (-t102 * t113 * t179 + t101) * t103) * t175) * t92 ((t124 * t143 - t135 * t167) * t93 - ((t179 * t90 + t125) * t102 + (-t115 * t90 - t111) * t101) * t176) * t92 (-t121 * t143 * t93 - ((t128 * t143 + t179 * t91) * t102 + (-t115 * t91 - t143 * t152) * t101) * t176) * t92 (t110 * t93 - ((t179 * t89 + t116) * t102 + (-t115 * t89 + t108) * t101) * t176) * t92, 0, 0; ((t108 * t142 - t146 * t152) * t97 - (t108 * t146 + t142 * t152) * t177) * t96 ((t112 * t142 - t123 * t146) * t97 - (t112 * t146 + t123 * t142) * t177) * t96 ((-t122 * t146 - t142 * t173) * t97 - (t122 * t142 - t146 * t173) * t177) * t96 (-t142 * t97 + t146 * t177) * t96 * t109, t160 * t96, 0;];
Ja_rot  = t1;
