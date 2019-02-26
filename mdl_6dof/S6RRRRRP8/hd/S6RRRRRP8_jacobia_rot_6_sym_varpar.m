% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP8
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRP8_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:43:48
% EndTime: 2019-02-26 22:43:49
% DurationCPUTime: 0.36s
% Computational Cost: add. (1836->47), mult. (3176->112), div. (137->9), fcn. (4481->13), ass. (0->64)
t140 = cos(pkin(6));
t142 = sin(qJ(2));
t146 = cos(qJ(1));
t151 = t146 * t142;
t143 = sin(qJ(1));
t145 = cos(qJ(2));
t152 = t143 * t145;
t132 = t140 * t151 + t152;
t138 = qJ(3) + qJ(4);
t136 = sin(t138);
t137 = cos(t138);
t139 = sin(pkin(6));
t155 = t139 * t146;
t124 = -t132 * t137 + t136 * t155;
t141 = sin(qJ(5));
t144 = cos(qJ(5));
t150 = t146 * t145;
t153 = t143 * t142;
t149 = -t140 * t150 + t153;
t167 = t124 * t141 + t149 * t144;
t148 = t149 * t141;
t166 = t124 * t144 - t148;
t157 = t139 * t142;
t129 = t140 * t136 + t137 * t157;
t120 = t139 * t145 * t144 + t129 * t141;
t108 = atan2(t167, t120);
t104 = sin(t108);
t105 = cos(t108);
t103 = t104 * t167 + t105 * t120;
t102 = 0.1e1 / t103 ^ 2;
t134 = -t140 * t153 + t150;
t156 = t139 * t143;
t126 = t134 * t137 + t136 * t156;
t133 = t140 * t152 + t151;
t158 = t133 * t144;
t114 = t126 * t141 - t158;
t164 = t102 * t114;
t159 = t133 * t141;
t115 = t126 * t144 + t159;
t110 = 0.1e1 / t115 ^ 2;
t125 = -t134 * t136 + t137 * t156;
t163 = t110 * t125;
t118 = 0.1e1 / t120 ^ 2;
t162 = t167 * t118;
t161 = t114 ^ 2 * t102;
t160 = t125 ^ 2 * t110;
t154 = t141 * t145;
t147 = t132 * t136 + t137 * t155;
t128 = -t136 * t157 + t140 * t137;
t127 = (t137 * t154 - t142 * t144) * t139;
t121 = t129 * t144 - t139 * t154;
t117 = 0.1e1 / t120;
t116 = -t132 * t144 - t137 * t148;
t109 = 0.1e1 / t115;
t107 = 0.1e1 / (0.1e1 + t160);
t106 = 0.1e1 / (t118 * t167 ^ 2 + 0.1e1);
t101 = 0.1e1 / t103;
t100 = 0.1e1 / (0.1e1 + t161);
t99 = (-t109 * t126 - t144 * t160) * t107;
t98 = (t117 * t147 - t128 * t162) * t141 * t106;
t97 = (-t116 * t117 - t127 * t162) * t106;
t96 = (t117 * t166 - t121 * t162) * t106;
t95 = (t125 * t141 * t101 - ((t128 * t141 + t167 * t98) * t105 + (-t120 * t98 + t141 * t147) * t104) * t164) * t100;
t1 = [-t114 * t117 * t106, t97, t98, t98, t96, 0; (t167 * t101 - (-t104 + (-t105 * t117 * t167 + t104) * t106) * t161) * t100 ((-t134 * t144 - t137 * t159) * t101 - ((t167 * t97 + t127) * t105 + (-t120 * t97 - t116) * t104) * t164) * t100, t95, t95 (t115 * t101 - ((t167 * t96 + t121) * t105 + (-t120 * t96 + t166) * t104) * t164) * t100, 0; (t147 * t109 - t166 * t163) * t107 (t133 * t136 * t109 - (t134 * t141 - t137 * t158) * t163) * t107, t99, t99, t114 * t107 * t163, 0;];
Ja_rot  = t1;
