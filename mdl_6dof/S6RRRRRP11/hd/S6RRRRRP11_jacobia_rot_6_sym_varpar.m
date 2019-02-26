% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function Ja_rot = S6RRRRRP11_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:45:38
% EndTime: 2019-02-26 22:45:39
% DurationCPUTime: 0.62s
% Computational Cost: add. (1916->67), mult. (5494->157), div. (115->9), fcn. (7543->17), ass. (0->78)
t142 = cos(pkin(6));
t150 = cos(qJ(2));
t179 = sin(qJ(1));
t158 = t179 * t150;
t146 = sin(qJ(2));
t151 = cos(qJ(1));
t163 = t151 * t146;
t135 = t142 * t163 + t158;
t145 = sin(qJ(3));
t149 = cos(qJ(3));
t159 = t179 * t146;
t162 = t151 * t150;
t134 = -t142 * t162 + t159;
t139 = sin(pkin(7));
t141 = cos(pkin(7));
t140 = sin(pkin(6));
t167 = t140 * t151;
t156 = t134 * t141 + t139 * t167;
t121 = -t135 * t149 + t156 * t145;
t131 = -t134 * t139 + t141 * t167;
t144 = sin(qJ(4));
t148 = cos(qJ(4));
t180 = t121 * t144 - t131 * t148;
t109 = t121 * t148 + t131 * t144;
t155 = t142 * t158 + t163;
t160 = t140 * t179;
t157 = t139 * t160;
t136 = -t142 * t159 + t162;
t170 = t136 * t149;
t123 = t170 + (-t155 * t141 + t157) * t145;
t152 = t155 * t139 + t141 * t160;
t111 = t123 * t148 + t152 * t144;
t154 = t155 * t149;
t122 = t136 * t145 + t141 * t154 - t149 * t157;
t143 = sin(qJ(5));
t147 = cos(qJ(5));
t100 = t111 * t143 - t122 * t147;
t101 = t111 * t147 + t122 * t143;
t99 = 0.1e1 / t101 ^ 2;
t178 = t100 * t99;
t110 = t123 * t144 - t152 * t148;
t166 = t141 * t145;
t169 = t139 * t142;
t130 = t145 * t169 + (t146 * t149 + t150 * t166) * t140;
t133 = -t140 * t150 * t139 + t142 * t141;
t116 = t130 * t144 - t133 * t148;
t105 = atan2(t180, t116);
t102 = sin(t105);
t103 = cos(t105);
t96 = t102 * t180 + t103 * t116;
t95 = 0.1e1 / t96 ^ 2;
t177 = t110 * t95;
t176 = t110 ^ 2 * t95;
t115 = 0.1e1 / t116 ^ 2;
t175 = t180 * t115;
t174 = t122 * t148;
t168 = t139 * t148;
t165 = t145 * t146;
t164 = t149 * t150;
t161 = t100 ^ 2 * t99 + 0.1e1;
t153 = -t135 * t145 - t156 * t149;
t129 = t149 * t169 + (t141 * t164 - t165) * t140;
t126 = ((-t141 * t165 + t164) * t144 - t146 * t168) * t140;
t125 = -t136 * t166 - t154;
t124 = t141 * t170 - t155 * t145;
t117 = t130 * t148 + t133 * t144;
t114 = 0.1e1 / t116;
t113 = t136 * t139 * t144 + t125 * t148;
t112 = (-t134 * t149 - t135 * t166) * t144 - t135 * t168;
t104 = 0.1e1 / (t115 * t180 ^ 2 + 0.1e1);
t98 = 0.1e1 / t101;
t97 = 0.1e1 / t161;
t94 = 0.1e1 / t96;
t93 = 0.1e1 / (0.1e1 + t176);
t92 = (-t114 * t153 - t129 * t175) * t144 * t104;
t91 = (-t112 * t114 - t126 * t175) * t104;
t90 = (t109 * t114 - t117 * t175) * t104;
t1 = [-t110 * t114 * t104, t91, t92, t90, 0, 0; (t180 * t94 - (-t102 + (-t103 * t114 * t180 + t102) * t104) * t176) * t93 ((t125 * t144 - t136 * t168) * t94 - ((t180 * t91 + t126) * t103 + (-t116 * t91 - t112) * t102) * t177) * t93 (-t122 * t144 * t94 - ((t129 * t144 + t180 * t92) * t103 + (-t116 * t92 - t144 * t153) * t102) * t177) * t93 (t111 * t94 - ((t180 * t90 + t117) * t103 + (-t116 * t90 + t109) * t102) * t177) * t93, 0, 0; ((t109 * t143 - t147 * t153) * t98 - (t109 * t147 + t143 * t153) * t178) * t97 ((t113 * t143 - t124 * t147) * t98 - (t113 * t147 + t124 * t143) * t178) * t97 ((-t123 * t147 - t143 * t174) * t98 - (t123 * t143 - t147 * t174) * t178) * t97 (-t143 * t98 + t147 * t178) * t97 * t110, t161 * t97, 0;];
Ja_rot  = t1;
