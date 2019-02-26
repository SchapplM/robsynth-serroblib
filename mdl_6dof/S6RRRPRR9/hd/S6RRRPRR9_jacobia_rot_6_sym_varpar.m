% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR9
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
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR9_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobia_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:20:37
% EndTime: 2019-02-26 22:20:38
% DurationCPUTime: 0.60s
% Computational Cost: add. (2899->70), mult. (8111->162), div. (115->9), fcn. (11143->19), ass. (0->78)
t153 = sin(pkin(13));
t160 = sin(qJ(3));
t186 = cos(pkin(13));
t187 = cos(qJ(3));
t150 = -t187 * t153 - t160 * t186;
t154 = sin(pkin(7));
t141 = t150 * t154;
t156 = cos(pkin(7));
t143 = t150 * t156;
t157 = cos(pkin(6));
t165 = cos(qJ(2));
t166 = cos(qJ(1));
t172 = t166 * t165;
t161 = sin(qJ(2));
t162 = sin(qJ(1));
t175 = t162 * t161;
t145 = -t157 * t172 + t175;
t173 = t166 * t161;
t174 = t162 * t165;
t146 = t157 * t173 + t174;
t169 = -t160 * t153 + t187 * t186;
t155 = sin(pkin(6));
t176 = t155 * t166;
t128 = -t141 * t176 - t145 * t143 - t146 * t169;
t139 = -t145 * t154 + t156 * t176;
t159 = sin(qJ(5));
t164 = cos(qJ(5));
t188 = t128 * t159 - t139 * t164;
t116 = t128 * t164 + t139 * t159;
t133 = -t157 * t141 + (-t143 * t165 + t161 * t169) * t155;
t144 = -t155 * t165 * t154 + t157 * t156;
t121 = t133 * t159 - t144 * t164;
t112 = atan2(t188, t121);
t109 = sin(t112);
t110 = cos(t112);
t103 = t109 * t188 + t110 * t121;
t102 = 0.1e1 / t103 ^ 2;
t147 = -t157 * t174 - t173;
t148 = -t157 * t175 + t172;
t177 = t155 * t162;
t167 = -t141 * t177 - t147 * t143 + t148 * t169;
t170 = -t147 * t154 + t156 * t177;
t117 = t159 * t167 - t170 * t164;
t185 = t102 * t117;
t118 = t170 * t159 + t164 * t167;
t140 = t169 * t154;
t142 = t169 * t156;
t130 = t140 * t177 + t147 * t142 + t148 * t150;
t158 = sin(qJ(6));
t163 = cos(qJ(6));
t108 = t118 * t163 - t130 * t158;
t106 = 0.1e1 / t108 ^ 2;
t107 = t118 * t158 + t130 * t163;
t184 = t106 * t107;
t120 = 0.1e1 / t121 ^ 2;
t183 = t188 * t120;
t182 = t117 ^ 2 * t102;
t181 = t130 * t164;
t178 = t154 * t164;
t171 = t107 ^ 2 * t106 + 0.1e1;
t168 = -t140 * t176 - t145 * t142 + t146 * t150;
t136 = ((t143 * t161 + t165 * t169) * t159 - t161 * t178) * t155;
t135 = t148 * t143 + t147 * t169;
t134 = t148 * t142 - t147 * t150;
t132 = t157 * t140 + (t142 * t165 + t150 * t161) * t155;
t124 = t148 * t154 * t159 + t135 * t164;
t123 = (t146 * t143 - t145 * t169) * t159 - t146 * t178;
t122 = t133 * t164 + t144 * t159;
t119 = 0.1e1 / t121;
t111 = 0.1e1 / (t120 * t188 ^ 2 + 0.1e1);
t105 = 0.1e1 / t108;
t104 = 0.1e1 / t171;
t101 = 0.1e1 / t103;
t100 = 0.1e1 / (0.1e1 + t182);
t99 = (-t119 * t123 - t136 * t183) * t111;
t98 = (-t119 * t168 - t132 * t183) * t159 * t111;
t97 = (t116 * t119 - t122 * t183) * t111;
t1 = [-t117 * t119 * t111, t99, t98, 0, t97, 0; (t188 * t101 - (-t109 + (-t110 * t119 * t188 + t109) * t111) * t182) * t100 ((t135 * t159 - t148 * t178) * t101 - ((t188 * t99 + t136) * t110 + (-t121 * t99 - t123) * t109) * t185) * t100 (t130 * t159 * t101 - ((t132 * t159 + t188 * t98) * t110 + (-t121 * t98 - t159 * t168) * t109) * t185) * t100, 0 (t118 * t101 - ((t188 * t97 + t122) * t110 + (-t121 * t97 + t116) * t109) * t185) * t100, 0; ((t116 * t158 - t163 * t168) * t105 - (t116 * t163 + t158 * t168) * t184) * t104 ((t124 * t158 - t134 * t163) * t105 - (t124 * t163 + t134 * t158) * t184) * t104 ((t158 * t181 - t163 * t167) * t105 - (t158 * t167 + t163 * t181) * t184) * t104, 0 (-t158 * t105 + t163 * t184) * t117 * t104, t171 * t104;];
Ja_rot  = t1;
