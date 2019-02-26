% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR14
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
% Datum: 2019-02-26 22:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR14_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_jacobia_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:11
% EndTime: 2019-02-26 22:38:12
% DurationCPUTime: 0.60s
% Computational Cost: add. (2000->68), mult. (5494->157), div. (115->9), fcn. (7543->17), ass. (0->79)
t145 = cos(pkin(6));
t151 = cos(qJ(2));
t180 = sin(qJ(1));
t159 = t180 * t151;
t148 = sin(qJ(2));
t152 = cos(qJ(1));
t164 = t152 * t148;
t135 = t145 * t164 + t159;
t147 = sin(qJ(3));
t150 = cos(qJ(3));
t160 = t180 * t148;
t163 = t152 * t151;
t134 = -t145 * t163 + t160;
t142 = sin(pkin(7));
t144 = cos(pkin(7));
t143 = sin(pkin(6));
t168 = t143 * t152;
t157 = t134 * t144 + t142 * t168;
t121 = -t135 * t150 + t147 * t157;
t131 = -t134 * t142 + t144 * t168;
t146 = sin(qJ(4));
t149 = cos(qJ(4));
t181 = t121 * t146 - t131 * t149;
t109 = t121 * t149 + t131 * t146;
t156 = t145 * t159 + t164;
t161 = t143 * t180;
t158 = t142 * t161;
t136 = -t145 * t160 + t163;
t171 = t136 * t150;
t123 = t171 + (-t144 * t156 + t158) * t147;
t153 = t142 * t156 + t144 * t161;
t111 = t123 * t149 + t146 * t153;
t155 = t156 * t150;
t122 = t136 * t147 + t144 * t155 - t150 * t158;
t141 = pkin(13) + qJ(6);
t139 = sin(t141);
t140 = cos(t141);
t100 = t111 * t139 - t122 * t140;
t101 = t111 * t140 + t122 * t139;
t99 = 0.1e1 / t101 ^ 2;
t179 = t100 * t99;
t110 = t123 * t146 - t149 * t153;
t167 = t144 * t147;
t170 = t142 * t145;
t130 = t147 * t170 + (t148 * t150 + t151 * t167) * t143;
t133 = -t143 * t151 * t142 + t145 * t144;
t116 = t130 * t146 - t133 * t149;
t105 = atan2(t181, t116);
t102 = sin(t105);
t103 = cos(t105);
t96 = t102 * t181 + t103 * t116;
t95 = 0.1e1 / t96 ^ 2;
t178 = t110 * t95;
t177 = t110 ^ 2 * t95;
t115 = 0.1e1 / t116 ^ 2;
t176 = t181 * t115;
t175 = t122 * t149;
t169 = t142 * t149;
t166 = t147 * t148;
t165 = t150 * t151;
t162 = t100 ^ 2 * t99 + 0.1e1;
t154 = -t135 * t147 - t150 * t157;
t129 = t150 * t170 + (t144 * t165 - t166) * t143;
t126 = ((-t144 * t166 + t165) * t146 - t148 * t169) * t143;
t125 = -t136 * t167 - t155;
t124 = t144 * t171 - t147 * t156;
t117 = t130 * t149 + t133 * t146;
t114 = 0.1e1 / t116;
t113 = t136 * t142 * t146 + t125 * t149;
t112 = (-t134 * t150 - t135 * t167) * t146 - t135 * t169;
t104 = 0.1e1 / (t115 * t181 ^ 2 + 0.1e1);
t98 = 0.1e1 / t101;
t97 = 0.1e1 / t162;
t94 = 0.1e1 / t96;
t93 = 0.1e1 / (0.1e1 + t177);
t92 = (-t114 * t154 - t129 * t176) * t146 * t104;
t91 = (-t112 * t114 - t126 * t176) * t104;
t90 = (t109 * t114 - t117 * t176) * t104;
t1 = [-t110 * t114 * t104, t91, t92, t90, 0, 0; (t181 * t94 - (-t102 + (-t103 * t114 * t181 + t102) * t104) * t177) * t93 ((t125 * t146 - t136 * t169) * t94 - ((t181 * t91 + t126) * t103 + (-t116 * t91 - t112) * t102) * t178) * t93 (-t122 * t146 * t94 - ((t129 * t146 + t181 * t92) * t103 + (-t116 * t92 - t146 * t154) * t102) * t178) * t93 (t111 * t94 - ((t181 * t90 + t117) * t103 + (-t116 * t90 + t109) * t102) * t178) * t93, 0, 0; ((t109 * t139 - t140 * t154) * t98 - (t109 * t140 + t139 * t154) * t179) * t97 ((t113 * t139 - t124 * t140) * t98 - (t113 * t140 + t124 * t139) * t179) * t97 ((-t123 * t140 - t139 * t175) * t98 - (t123 * t139 - t140 * t175) * t179) * t97 (-t139 * t98 + t140 * t179) * t97 * t110, 0, t162 * t97;];
Ja_rot  = t1;
