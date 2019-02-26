% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR15_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:24:15
% EndTime: 2019-02-26 22:24:16
% DurationCPUTime: 0.58s
% Computational Cost: add. (1916->67), mult. (5494->156), div. (115->9), fcn. (7543->17), ass. (0->78)
t142 = cos(pkin(6));
t151 = cos(qJ(2));
t152 = cos(qJ(1));
t158 = t152 * t151;
t146 = sin(qJ(2));
t147 = sin(qJ(1));
t161 = t147 * t146;
t133 = -t142 * t158 + t161;
t139 = sin(pkin(7));
t141 = cos(pkin(7));
t140 = sin(pkin(6));
t165 = t140 * t152;
t128 = -t133 * t139 + t141 * t165;
t144 = sin(qJ(5));
t149 = cos(qJ(5));
t159 = t152 * t146;
t160 = t147 * t151;
t134 = t142 * t159 + t160;
t145 = sin(qJ(3));
t150 = cos(qJ(3));
t154 = t133 * t141 + t139 * t165;
t153 = t134 * t145 + t150 * t154;
t108 = -t128 * t149 + t144 * t153;
t180 = t128 * t144 + t149 * t153;
t121 = -t134 * t150 + t145 * t154;
t135 = -t142 * t160 - t159;
t166 = t140 * t147;
t130 = -t135 * t139 + t141 * t166;
t156 = t139 * t166;
t164 = t141 * t150;
t136 = -t142 * t161 + t158;
t169 = t136 * t145;
t155 = -t135 * t164 - t150 * t156 + t169;
t110 = t130 * t149 + t144 * t155;
t148 = cos(qJ(6));
t122 = t136 * t150 + (t135 * t141 + t156) * t145;
t143 = sin(qJ(6));
t173 = t122 * t143;
t100 = t110 * t148 + t173;
t98 = 0.1e1 / t100 ^ 2;
t172 = t122 * t148;
t99 = t110 * t143 - t172;
t177 = t98 * t99;
t109 = t130 * t144 - t149 * t155;
t168 = t139 * t142;
t126 = -t150 * t168 + (t145 * t146 - t151 * t164) * t140;
t132 = -t140 * t151 * t139 + t142 * t141;
t117 = -t126 * t149 + t132 * t144;
t104 = atan2(t180, t117);
t101 = sin(t104);
t102 = cos(t104);
t95 = t101 * t180 + t102 * t117;
t94 = 0.1e1 / t95 ^ 2;
t176 = t109 * t94;
t175 = t109 ^ 2 * t94;
t116 = 0.1e1 / t117 ^ 2;
t174 = t180 * t116;
t167 = t139 * t144;
t163 = t145 * t151;
t162 = t146 * t150;
t157 = t99 ^ 2 * t98 + 0.1e1;
t127 = t145 * t168 + (t141 * t163 + t162) * t140;
t125 = (t146 * t167 - (t141 * t162 + t163) * t149) * t140;
t124 = t135 * t150 - t141 * t169;
t123 = t135 * t145 + t136 * t164;
t118 = t126 * t144 + t132 * t149;
t115 = 0.1e1 / t117;
t112 = t136 * t139 * t149 + t123 * t144;
t111 = -t134 * t167 + (-t133 * t145 + t134 * t164) * t149;
t103 = 0.1e1 / (t116 * t180 ^ 2 + 0.1e1);
t97 = 0.1e1 / t100;
t96 = 0.1e1 / t157;
t93 = 0.1e1 / t95;
t92 = 0.1e1 / (0.1e1 + t175);
t91 = (-t115 * t121 + t127 * t174) * t149 * t103;
t90 = (t111 * t115 - t125 * t174) * t103;
t89 = (-t108 * t115 - t118 * t174) * t103;
t1 = [-t109 * t115 * t103, t90, t91, 0, t89, 0; (t180 * t93 - (-t101 + (-t102 * t115 * t180 + t101) * t103) * t175) * t92 ((-t123 * t149 + t136 * t167) * t93 - ((t180 * t90 + t125) * t102 + (-t117 * t90 + t111) * t101) * t176) * t92 (-t122 * t149 * t93 - ((-t127 * t149 + t180 * t91) * t102 + (-t117 * t91 - t121 * t149) * t101) * t176) * t92, 0 (t110 * t93 - ((t180 * t89 + t118) * t102 + (-t117 * t89 - t108) * t101) * t176) * t92, 0; ((-t108 * t143 - t121 * t148) * t97 - (-t108 * t148 + t121 * t143) * t177) * t96 ((t112 * t143 - t124 * t148) * t97 - (t112 * t148 + t124 * t143) * t177) * t96 ((t144 * t173 + t148 * t155) * t97 - (-t143 * t155 + t144 * t172) * t177) * t96, 0 (-t143 * t97 + t148 * t177) * t96 * t109, t157 * t96;];
Ja_rot  = t1;
