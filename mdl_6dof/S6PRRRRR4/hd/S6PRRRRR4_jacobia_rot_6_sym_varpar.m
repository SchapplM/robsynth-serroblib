% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRR4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_jacobia_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:20:36
% EndTime: 2019-02-26 20:20:36
% DurationCPUTime: 0.38s
% Computational Cost: add. (2625->62), mult. (5778->144), div. (125->9), fcn. (7926->17), ass. (0->78)
t144 = sin(pkin(13));
t147 = cos(pkin(13));
t155 = cos(qJ(2));
t149 = cos(pkin(6));
t152 = sin(qJ(2));
t165 = t149 * t152;
t136 = t144 * t155 + t147 * t165;
t151 = sin(qJ(3));
t154 = cos(qJ(3));
t148 = cos(pkin(7));
t164 = t149 * t155;
t158 = -t144 * t152 + t147 * t164;
t145 = sin(pkin(7));
t146 = sin(pkin(6));
t169 = t146 * t145;
t156 = -t147 * t169 + t158 * t148;
t124 = t136 * t154 + t156 * t151;
t143 = qJ(4) + qJ(5);
t141 = sin(t143);
t142 = cos(t143);
t168 = t146 * t148;
t157 = -t158 * t145 - t147 * t168;
t112 = t124 * t141 - t157 * t142;
t167 = t148 * t151;
t170 = t145 * t149;
t133 = t151 * t170 + (t152 * t154 + t155 * t167) * t146;
t135 = t149 * t148 - t155 * t169;
t121 = t133 * t141 - t135 * t142;
t111 = atan2(-t112, t121);
t108 = sin(t111);
t109 = cos(t111);
t102 = -t108 * t112 + t109 * t121;
t101 = 0.1e1 / t102 ^ 2;
t137 = -t144 * t164 - t147 * t152;
t138 = -t144 * t165 + t147 * t155;
t161 = t144 * t169;
t126 = t138 * t154 + (t137 * t148 + t161) * t151;
t159 = -t137 * t145 + t144 * t168;
t115 = t126 * t141 - t159 * t142;
t176 = t101 * t115;
t116 = t126 * t142 + t159 * t141;
t153 = cos(qJ(6));
t166 = t148 * t154;
t125 = -t137 * t166 + t138 * t151 - t154 * t161;
t150 = sin(qJ(6));
t173 = t125 * t150;
t107 = t116 * t153 + t173;
t105 = 0.1e1 / t107 ^ 2;
t172 = t125 * t153;
t106 = t116 * t150 - t172;
t175 = t105 * t106;
t120 = 0.1e1 / t121 ^ 2;
t174 = t112 * t120;
t171 = t145 * t142;
t163 = t151 * t152;
t162 = t154 * t155;
t160 = t106 ^ 2 * t105 + 0.1e1;
t132 = t154 * t170 + (t148 * t162 - t163) * t146;
t129 = ((-t148 * t163 + t162) * t141 - t152 * t171) * t146;
t128 = t137 * t154 - t138 * t167;
t127 = t137 * t151 + t138 * t166;
t123 = -t136 * t151 + t156 * t154;
t122 = t133 * t142 + t135 * t141;
t119 = 0.1e1 / t121;
t118 = t138 * t145 * t141 + t128 * t142;
t117 = (-t136 * t167 + t158 * t154) * t141 - t136 * t171;
t114 = t124 * t142 + t157 * t141;
t110 = 0.1e1 / (t112 ^ 2 * t120 + 0.1e1);
t104 = 0.1e1 / t107;
t103 = 0.1e1 / t160;
t100 = 0.1e1 / t102;
t99 = 0.1e1 / (t115 ^ 2 * t101 + 0.1e1);
t98 = (-t119 * t123 + t132 * t174) * t141 * t110;
t97 = (-t117 * t119 + t129 * t174) * t110;
t96 = (-t114 * t119 + t122 * t174) * t110;
t95 = (-t150 * t104 + t153 * t175) * t115 * t103;
t94 = (t116 * t100 - ((-t112 * t96 + t122) * t109 + (-t121 * t96 - t114) * t108) * t176) * t99;
t1 = [0, t97, t98, t96, t96, 0; 0 ((t128 * t141 - t138 * t171) * t100 - ((-t112 * t97 + t129) * t109 + (-t121 * t97 - t117) * t108) * t176) * t99 (-t125 * t141 * t100 - ((-t112 * t98 + t132 * t141) * t109 + (-t121 * t98 - t123 * t141) * t108) * t176) * t99, t94, t94, 0; 0 ((t118 * t150 - t127 * t153) * t104 - (t118 * t153 + t127 * t150) * t175) * t103 ((-t126 * t153 - t142 * t173) * t104 - (t126 * t150 - t142 * t172) * t175) * t103, t95, t95, t160 * t103;];
Ja_rot  = t1;
