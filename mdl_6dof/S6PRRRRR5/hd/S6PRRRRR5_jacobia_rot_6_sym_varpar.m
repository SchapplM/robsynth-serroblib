% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR5
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
% Datum: 2019-02-26 20:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRR5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_jacobia_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:21:10
% EndTime: 2019-02-26 20:21:11
% DurationCPUTime: 0.39s
% Computational Cost: add. (1691->62), mult. (4622->145), div. (100->9), fcn. (6331->17), ass. (0->76)
t138 = sin(pkin(13));
t141 = cos(pkin(13));
t146 = sin(qJ(2));
t143 = cos(pkin(6));
t149 = cos(qJ(2));
t158 = t143 * t149;
t131 = -t138 * t158 - t141 * t146;
t159 = t143 * t146;
t132 = -t138 * t159 + t141 * t149;
t142 = cos(pkin(7));
t145 = sin(qJ(3));
t148 = cos(qJ(3));
t139 = sin(pkin(7));
t140 = sin(pkin(6));
t163 = t140 * t139;
t155 = t138 * t163;
t118 = t132 * t148 + (t131 * t142 + t155) * t145;
t144 = sin(qJ(4));
t147 = cos(qJ(4));
t162 = t140 * t142;
t153 = -t131 * t139 + t138 * t162;
t110 = t118 * t147 + t153 * t144;
t160 = t142 * t148;
t117 = -t131 * t160 + t132 * t145 - t148 * t155;
t137 = qJ(5) + qJ(6);
t135 = sin(t137);
t136 = cos(t137);
t100 = t110 * t135 - t117 * t136;
t101 = t110 * t136 + t117 * t135;
t99 = 0.1e1 / t101 ^ 2;
t169 = t100 * t99;
t109 = t118 * t144 - t153 * t147;
t130 = t138 * t149 + t141 * t159;
t152 = -t138 * t146 + t141 * t158;
t150 = -t141 * t163 + t152 * t142;
t116 = t130 * t148 + t150 * t145;
t151 = -t152 * t139 - t141 * t162;
t106 = t116 * t144 - t151 * t147;
t161 = t142 * t145;
t165 = t139 * t143;
t127 = t145 * t165 + (t146 * t148 + t149 * t161) * t140;
t129 = t143 * t142 - t149 * t163;
t119 = t127 * t144 - t129 * t147;
t105 = atan2(-t106, t119);
t102 = sin(t105);
t103 = cos(t105);
t96 = -t102 * t106 + t103 * t119;
t95 = 0.1e1 / t96 ^ 2;
t168 = t109 * t95;
t114 = 0.1e1 / t119 ^ 2;
t167 = t106 * t114;
t166 = t117 * t147;
t164 = t139 * t147;
t157 = t145 * t146;
t156 = t148 * t149;
t154 = t100 ^ 2 * t99 + 0.1e1;
t126 = t148 * t165 + (t142 * t156 - t157) * t140;
t123 = ((-t142 * t157 + t156) * t144 - t146 * t164) * t140;
t122 = t131 * t148 - t132 * t161;
t121 = t131 * t145 + t132 * t160;
t120 = t127 * t147 + t129 * t144;
t115 = -t130 * t145 + t150 * t148;
t113 = 0.1e1 / t119;
t112 = t132 * t139 * t144 + t122 * t147;
t111 = (-t130 * t161 + t152 * t148) * t144 - t130 * t164;
t108 = t116 * t147 + t151 * t144;
t104 = 0.1e1 / (t106 ^ 2 * t114 + 0.1e1);
t98 = 0.1e1 / t101;
t97 = 0.1e1 / t154;
t94 = 0.1e1 / t96;
t93 = 0.1e1 / (t109 ^ 2 * t95 + 0.1e1);
t92 = (-t113 * t115 + t126 * t167) * t144 * t104;
t91 = (-t111 * t113 + t123 * t167) * t104;
t90 = (-t108 * t113 + t120 * t167) * t104;
t89 = t154 * t97;
t1 = [0, t91, t92, t90, 0, 0; 0 ((t122 * t144 - t132 * t164) * t94 - ((-t106 * t91 + t123) * t103 + (-t119 * t91 - t111) * t102) * t168) * t93 (-t117 * t144 * t94 - ((-t106 * t92 + t126 * t144) * t103 + (-t115 * t144 - t119 * t92) * t102) * t168) * t93 (t110 * t94 - ((-t106 * t90 + t120) * t103 + (-t119 * t90 - t108) * t102) * t168) * t93, 0, 0; 0 ((t112 * t135 - t121 * t136) * t98 - (t112 * t136 + t121 * t135) * t169) * t97 ((-t118 * t136 - t135 * t166) * t98 - (t118 * t135 - t136 * t166) * t169) * t97 (-t135 * t98 + t136 * t169) * t97 * t109, t89, t89;];
Ja_rot  = t1;
