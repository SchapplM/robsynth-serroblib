% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRP12_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:14:18
% EndTime: 2019-02-26 21:14:18
% DurationCPUTime: 0.64s
% Computational Cost: add. (2697->57), mult. (7742->135), div. (107->9), fcn. (10572->17), ass. (0->77)
t145 = cos(pkin(6));
t140 = sin(pkin(12));
t153 = cos(qJ(1));
t161 = t153 * t140;
t143 = cos(pkin(12));
t149 = sin(qJ(1));
t162 = t149 * t143;
t137 = t145 * t161 + t162;
t148 = sin(qJ(3));
t152 = cos(qJ(3));
t160 = t153 * t143;
t163 = t149 * t140;
t136 = -t145 * t160 + t163;
t141 = sin(pkin(7));
t144 = cos(pkin(7));
t142 = sin(pkin(6));
t166 = t142 * t153;
t157 = t136 * t144 + t141 * t166;
t126 = -t137 * t152 + t157 * t148;
t131 = -t136 * t141 + t144 * t166;
t147 = sin(qJ(4));
t151 = cos(qJ(4));
t116 = t126 * t151 + t131 * t147;
t146 = sin(qJ(5));
t150 = cos(qJ(5));
t154 = t137 * t148 + t157 * t152;
t184 = t116 * t146 + t154 * t150;
t183 = t116 * t150 - t154 * t146;
t114 = t126 * t147 - t131 * t151;
t156 = t145 * t162 + t161;
t167 = t142 * t149;
t178 = -t141 * t167 + t156 * t144;
t165 = t143 * t144;
t168 = t141 * t145;
t130 = t148 * t168 + (t140 * t152 + t148 * t165) * t142;
t135 = -t142 * t143 * t141 + t145 * t144;
t123 = t130 * t151 + t135 * t147;
t129 = -t152 * t168 + (t140 * t148 - t152 * t165) * t142;
t111 = t123 * t146 - t129 * t150;
t99 = atan2(t184, t111);
t97 = cos(t99);
t177 = t184 * t97;
t138 = -t145 * t163 + t160;
t128 = t138 * t152 - t178 * t148;
t133 = t156 * t141 + t144 * t167;
t118 = t128 * t151 + t133 * t147;
t127 = t138 * t148 + t178 * t152;
t171 = t127 * t150;
t106 = t118 * t146 - t171;
t96 = sin(t99);
t95 = t97 * t111 + t184 * t96;
t94 = 0.1e1 / t95 ^ 2;
t176 = t106 * t94;
t175 = t106 ^ 2 * t94;
t107 = t118 * t150 + t127 * t146;
t102 = 0.1e1 / t107 ^ 2;
t117 = -t128 * t147 + t133 * t151;
t174 = t102 * t117;
t110 = 0.1e1 / t111 ^ 2;
t173 = t184 * t110;
t172 = t117 ^ 2 * t102;
t164 = t146 * t151;
t158 = -t111 * t96 + t177;
t122 = -t130 * t147 + t135 * t151;
t119 = -t129 * t164 - t130 * t150;
t112 = t123 * t150 + t129 * t146;
t109 = 0.1e1 / t111;
t108 = t126 * t150 - t154 * t164;
t101 = 0.1e1 / t107;
t100 = 0.1e1 / (0.1e1 + t172);
t98 = 0.1e1 / (t110 * t184 ^ 2 + 0.1e1);
t93 = 0.1e1 / t95;
t92 = 0.1e1 / (0.1e1 + t175);
t91 = (-t109 * t114 - t122 * t173) * t98 * t146;
t90 = (-t108 * t109 - t119 * t173) * t98;
t89 = (t109 * t183 - t112 * t173) * t98;
t1 = [-t106 * t109 * t98, 0, t90, t91, t89, 0; (t184 * t93 - (-t96 + (-t109 * t177 + t96) * t98) * t175) * t92, 0 ((-t127 * t164 - t128 * t150) * t93 - (-t96 * t108 + t97 * t119 + t158 * t90) * t176) * t92 (t117 * t146 * t93 - (t158 * t91 + (-t114 * t96 + t122 * t97) * t146) * t176) * t92 (t107 * t93 - (t97 * t112 + t158 * t89 + t183 * t96) * t176) * t92, 0; (-t114 * t101 - t183 * t174) * t100, 0 (t127 * t147 * t101 - (t128 * t146 - t151 * t171) * t174) * t100 (-t101 * t118 - t150 * t172) * t100, t106 * t100 * t174, 0;];
Ja_rot  = t1;
