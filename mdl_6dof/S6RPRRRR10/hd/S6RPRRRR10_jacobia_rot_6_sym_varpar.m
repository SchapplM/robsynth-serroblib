% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRR10_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_jacobia_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:19:50
% EndTime: 2019-02-26 21:19:50
% DurationCPUTime: 0.51s
% Computational Cost: add. (2501->51), mult. (5521->119), div. (115->9), fcn. (7584->17), ass. (0->71)
t145 = cos(pkin(6));
t143 = cos(pkin(13));
t175 = sin(qJ(1));
t158 = t175 * t143;
t140 = sin(pkin(13));
t150 = cos(qJ(1));
t162 = t150 * t140;
t134 = t145 * t162 + t158;
t147 = sin(qJ(3));
t149 = cos(qJ(3));
t159 = t175 * t140;
t161 = t150 * t143;
t133 = -t145 * t161 + t159;
t141 = sin(pkin(7));
t144 = cos(pkin(7));
t142 = sin(pkin(6));
t164 = t142 * t150;
t155 = t133 * t144 + t141 * t164;
t123 = -t134 * t149 + t155 * t147;
t130 = -t133 * t141 + t144 * t164;
t139 = qJ(4) + qJ(5);
t137 = sin(t139);
t138 = cos(t139);
t177 = t123 * t137 - t130 * t138;
t113 = t123 * t138 + t130 * t137;
t154 = t145 * t158 + t162;
t160 = t142 * t175;
t176 = -t141 * t160 + t154 * t144;
t135 = -t145 * t159 + t161;
t125 = t135 * t149 - t147 * t176;
t151 = t154 * t141 + t144 * t160;
t114 = t125 * t137 - t151 * t138;
t163 = t143 * t144;
t165 = t141 * t145;
t129 = t147 * t165 + (t140 * t149 + t147 * t163) * t142;
t132 = -t142 * t143 * t141 + t145 * t144;
t118 = t129 * t137 - t132 * t138;
t109 = atan2(t177, t118);
t102 = sin(t109);
t103 = cos(t109);
t100 = t102 * t177 + t103 * t118;
t99 = 0.1e1 / t100 ^ 2;
t174 = t114 * t99;
t173 = t114 ^ 2 * t99;
t115 = t125 * t138 + t151 * t137;
t148 = cos(qJ(6));
t124 = t135 * t147 + t149 * t176;
t146 = sin(qJ(6));
t170 = t124 * t146;
t108 = t115 * t148 + t170;
t105 = 0.1e1 / t108 ^ 2;
t169 = t124 * t148;
t107 = t115 * t146 - t169;
t172 = t105 * t107;
t117 = 0.1e1 / t118 ^ 2;
t171 = t177 * t117;
t157 = t107 ^ 2 * t105 + 0.1e1;
t152 = -t134 * t147 - t155 * t149;
t128 = t149 * t165 + (-t140 * t147 + t149 * t163) * t142;
t119 = t129 * t138 + t132 * t137;
t116 = 0.1e1 / t118;
t106 = 0.1e1 / (t117 * t177 ^ 2 + 0.1e1);
t104 = 0.1e1 / t108;
t101 = 0.1e1 / t157;
t98 = 0.1e1 / t100;
t97 = 0.1e1 / (0.1e1 + t173);
t96 = (-t116 * t152 - t128 * t171) * t137 * t106;
t95 = (t113 * t116 - t119 * t171) * t106;
t94 = (-t146 * t104 + t148 * t172) * t114 * t101;
t93 = (t115 * t98 - ((t177 * t95 + t119) * t103 + (-t118 * t95 + t113) * t102) * t174) * t97;
t1 = [-t114 * t116 * t106, 0, t96, t95, t95, 0; (t177 * t98 - (-t102 + (-t103 * t116 * t177 + t102) * t106) * t173) * t97, 0 (-t124 * t137 * t98 - ((t128 * t137 + t177 * t96) * t103 + (-t118 * t96 - t137 * t152) * t102) * t174) * t97, t93, t93, 0; ((t113 * t146 - t148 * t152) * t104 - (t113 * t148 + t146 * t152) * t172) * t101, 0 ((-t125 * t148 - t138 * t170) * t104 - (t125 * t146 - t138 * t169) * t172) * t101, t94, t94, t157 * t101;];
Ja_rot  = t1;
