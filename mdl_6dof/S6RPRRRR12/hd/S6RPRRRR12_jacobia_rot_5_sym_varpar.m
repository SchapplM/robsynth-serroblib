% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR12
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRR12_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobia_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:09
% EndTime: 2019-02-26 21:21:10
% DurationCPUTime: 0.68s
% Computational Cost: add. (2433->62), mult. (7147->138), div. (85->9), fcn. (9679->19), ass. (0->79)
t143 = cos(pkin(6));
t140 = cos(pkin(14));
t151 = cos(qJ(1));
t166 = t151 * t140;
t136 = sin(pkin(14));
t147 = sin(qJ(1));
t170 = t147 * t136;
t132 = -t143 * t166 + t170;
t167 = t151 * t136;
t168 = t147 * t140;
t133 = t143 * t167 + t168;
t146 = sin(qJ(3));
t150 = cos(qJ(3));
t138 = sin(pkin(7));
t139 = sin(pkin(6));
t173 = t139 * t151;
t164 = t138 * t173;
t142 = cos(pkin(7));
t171 = t142 * t146;
t124 = t132 * t171 - t133 * t150 + t146 * t164;
t145 = sin(qJ(4));
t149 = cos(qJ(4));
t123 = (t132 * t142 + t164) * t150 + t133 * t146;
t130 = -t132 * t138 + t142 * t173;
t137 = sin(pkin(8));
t141 = cos(pkin(8));
t163 = t123 * t141 + t130 * t137;
t108 = t124 * t149 + t163 * t145;
t105 = -t124 * t145 + t163 * t149;
t134 = -t143 * t170 + t166;
t160 = -t143 * t168 - t167;
t169 = t147 * t139;
t156 = t138 * t169 + t160 * t142;
t125 = t134 * t150 + t156 * t146;
t153 = t134 * t146 - t156 * t150;
t157 = -t160 * t138 + t142 * t169;
t154 = t157 * t137;
t110 = t125 * t149 + (-t153 * t141 + t154) * t145;
t118 = t153 * t137 + t157 * t141;
t144 = sin(qJ(5));
t148 = cos(qJ(5));
t100 = t110 * t148 + t118 * t144;
t98 = 0.1e1 / t100 ^ 2;
t99 = t110 * t144 - t118 * t148;
t182 = t98 * t99;
t152 = t153 * t149;
t177 = t125 * t145;
t109 = t141 * t152 - t149 * t154 + t177;
t174 = t138 * t143;
t129 = t146 * t174 + (t136 * t150 + t140 * t171) * t139;
t128 = t150 * t174 + (t140 * t142 * t150 - t136 * t146) * t139;
t162 = t128 * t141 + (-t139 * t140 * t138 + t143 * t142) * t137;
t115 = t129 * t145 - t162 * t149;
t104 = atan2(-t105, t115);
t101 = sin(t104);
t102 = cos(t104);
t95 = -t101 * t105 + t102 * t115;
t94 = 0.1e1 / t95 ^ 2;
t181 = t109 * t94;
t180 = t109 ^ 2 * t94;
t114 = 0.1e1 / t115 ^ 2;
t179 = t105 * t114;
t178 = t125 * t137;
t172 = t141 * t149;
t165 = t99 ^ 2 * t98 + 0.1e1;
t119 = t128 * t145 + t129 * t172;
t117 = -t123 * t137 + t130 * t141;
t116 = t129 * t149 + t162 * t145;
t113 = 0.1e1 / t115;
t112 = -t141 * t177 - t152;
t111 = -t123 * t145 - t124 * t172;
t103 = 0.1e1 / (t105 ^ 2 * t114 + 0.1e1);
t97 = 0.1e1 / t100;
t96 = 0.1e1 / t165;
t93 = 0.1e1 / t95;
t92 = 0.1e1 / (0.1e1 + t180);
t91 = (-t111 * t113 + t119 * t179) * t103;
t90 = (t108 * t113 + t116 * t179) * t103;
t1 = [-t109 * t113 * t103, 0, t91, t90, 0, 0; (-t105 * t93 - (-t101 + (t102 * t105 * t113 + t101) * t103) * t180) * t92, 0 ((t125 * t172 - t153 * t145) * t93 - ((-t105 * t91 + t119) * t102 + (-t115 * t91 - t111) * t101) * t181) * t92 (t110 * t93 - ((-t105 * t90 + t116) * t102 + (-t115 * t90 + t108) * t101) * t181) * t92, 0, 0; ((t108 * t144 - t117 * t148) * t97 - (t108 * t148 + t117 * t144) * t182) * t96, 0 ((t112 * t144 - t148 * t178) * t97 - (t112 * t148 + t144 * t178) * t182) * t96 (-t144 * t97 + t148 * t182) * t96 * t109, t165 * t96, 0;];
Ja_rot  = t1;
