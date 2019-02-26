% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR12_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_jacobia_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:36:59
% EndTime: 2019-02-26 22:37:00
% DurationCPUTime: 0.59s
% Computational Cost: add. (2487->68), mult. (5494->156), div. (115->9), fcn. (7543->17), ass. (0->80)
t151 = cos(pkin(6));
t157 = cos(qJ(2));
t187 = sin(qJ(1));
t166 = t187 * t157;
t154 = sin(qJ(2));
t158 = cos(qJ(1));
t170 = t158 * t154;
t141 = t151 * t170 + t166;
t153 = sin(qJ(3));
t156 = cos(qJ(3));
t167 = t187 * t154;
t169 = t158 * t157;
t140 = -t151 * t169 + t167;
t148 = sin(pkin(7));
t150 = cos(pkin(7));
t149 = sin(pkin(6));
t174 = t149 * t158;
t163 = t140 * t150 + t148 * t174;
t127 = -t141 * t156 + t163 * t153;
t137 = -t140 * t148 + t150 * t174;
t147 = qJ(4) + pkin(13);
t145 = sin(t147);
t146 = cos(t147);
t188 = t127 * t145 - t137 * t146;
t115 = t127 * t146 + t137 * t145;
t173 = t150 * t153;
t175 = t148 * t151;
t136 = t153 * t175 + (t154 * t156 + t157 * t173) * t149;
t139 = -t149 * t157 * t148 + t151 * t150;
t122 = t136 * t145 - t139 * t146;
t109 = atan2(t188, t122);
t104 = sin(t109);
t105 = cos(t109);
t102 = t104 * t188 + t105 * t122;
t101 = 0.1e1 / t102 ^ 2;
t162 = t151 * t166 + t170;
t168 = t149 * t187;
t164 = t148 * t168;
t142 = -t151 * t167 + t169;
t177 = t142 * t156;
t129 = t177 + (-t162 * t150 + t164) * t153;
t159 = t162 * t148 + t150 * t168;
t116 = t129 * t145 - t159 * t146;
t186 = t101 * t116;
t117 = t129 * t146 + t159 * t145;
t155 = cos(qJ(6));
t161 = t162 * t156;
t128 = t142 * t153 + t150 * t161 - t156 * t164;
t152 = sin(qJ(6));
t182 = t128 * t152;
t111 = t117 * t155 + t182;
t108 = 0.1e1 / t111 ^ 2;
t181 = t128 * t155;
t110 = t117 * t152 - t181;
t185 = t108 * t110;
t121 = 0.1e1 / t122 ^ 2;
t184 = t188 * t121;
t183 = t116 ^ 2 * t101;
t176 = t148 * t146;
t172 = t153 * t154;
t171 = t156 * t157;
t165 = t110 ^ 2 * t108 + 0.1e1;
t160 = -t141 * t153 - t163 * t156;
t135 = t156 * t175 + (t150 * t171 - t172) * t149;
t132 = -t142 * t173 - t161;
t131 = t150 * t177 - t162 * t153;
t130 = ((-t150 * t172 + t171) * t145 - t154 * t176) * t149;
t123 = t136 * t146 + t139 * t145;
t120 = 0.1e1 / t122;
t119 = t142 * t148 * t145 + t132 * t146;
t118 = (-t140 * t156 - t141 * t173) * t145 - t141 * t176;
t107 = 0.1e1 / t111;
t106 = 0.1e1 / (t121 * t188 ^ 2 + 0.1e1);
t103 = 0.1e1 / t165;
t100 = 0.1e1 / t102;
t99 = 0.1e1 / (0.1e1 + t183);
t98 = (-t120 * t160 - t135 * t184) * t145 * t106;
t97 = (-t118 * t120 - t130 * t184) * t106;
t96 = (t115 * t120 - t123 * t184) * t106;
t1 = [-t116 * t120 * t106, t97, t98, t96, 0, 0; (t188 * t100 - (-t104 + (-t105 * t120 * t188 + t104) * t106) * t183) * t99 ((t132 * t145 - t142 * t176) * t100 - ((t188 * t97 + t130) * t105 + (-t122 * t97 - t118) * t104) * t186) * t99 (-t128 * t145 * t100 - ((t135 * t145 + t188 * t98) * t105 + (-t122 * t98 - t145 * t160) * t104) * t186) * t99 (t117 * t100 - ((t188 * t96 + t123) * t105 + (-t122 * t96 + t115) * t104) * t186) * t99, 0, 0; ((t115 * t152 - t155 * t160) * t107 - (t115 * t155 + t152 * t160) * t185) * t103 ((t119 * t152 - t131 * t155) * t107 - (t119 * t155 + t131 * t152) * t185) * t103 ((-t129 * t155 - t146 * t182) * t107 - (t129 * t152 - t146 * t181) * t185) * t103 (-t152 * t107 + t155 * t185) * t116 * t103, 0, t165 * t103;];
Ja_rot  = t1;
