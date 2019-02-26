% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRR3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_jacobia_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:14
% EndTime: 2019-02-26 20:05:14
% DurationCPUTime: 0.45s
% Computational Cost: add. (2312->62), mult. (6475->149), div. (95->9), fcn. (8887->19), ass. (0->77)
t166 = cos(qJ(3));
t141 = sin(qJ(5));
t145 = cos(qJ(5));
t135 = sin(pkin(7));
t133 = sin(pkin(13));
t142 = sin(qJ(3));
t162 = cos(pkin(13));
t150 = t166 * t133 + t142 * t162;
t122 = t150 * t135;
t138 = cos(pkin(7));
t124 = t150 * t138;
t134 = sin(pkin(12));
t137 = cos(pkin(12));
t143 = sin(qJ(2));
t139 = cos(pkin(6));
t146 = cos(qJ(2));
t155 = t139 * t146;
t127 = -t134 * t155 - t137 * t143;
t156 = t139 * t143;
t128 = -t134 * t156 + t137 * t146;
t149 = -t142 * t133 + t166 * t162;
t136 = sin(pkin(6));
t160 = t134 * t136;
t148 = t122 * t160 + t127 * t124 + t128 * t149;
t157 = t136 * t138;
t152 = -t127 * t135 + t134 * t157;
t102 = t152 * t141 + t145 * t148;
t121 = t149 * t135;
t123 = t149 * t138;
t110 = t121 * t160 + t127 * t123 - t128 * t150;
t140 = sin(qJ(6));
t144 = cos(qJ(6));
t93 = t102 * t144 - t110 * t140;
t91 = 0.1e1 / t93 ^ 2;
t92 = t102 * t140 + t110 * t144;
t165 = t91 * t92;
t101 = t141 * t148 - t152 * t145;
t115 = t139 * t122 + (t124 * t146 + t143 * t149) * t136;
t125 = -t136 * t146 * t135 + t139 * t138;
t112 = t115 * t141 - t125 * t145;
t126 = t134 * t146 + t137 * t156;
t151 = -t134 * t143 + t137 * t155;
t158 = t136 * t137;
t108 = -t122 * t158 + t151 * t124 + t126 * t149;
t147 = -t151 * t135 - t137 * t157;
t98 = t108 * t141 - t147 * t145;
t97 = atan2(-t98, t112);
t94 = sin(t97);
t95 = cos(t97);
t88 = t95 * t112 - t94 * t98;
t87 = 0.1e1 / t88 ^ 2;
t164 = t101 * t87;
t106 = 0.1e1 / t112 ^ 2;
t163 = t106 * t98;
t161 = t110 * t145;
t159 = t135 * t145;
t154 = t92 ^ 2 * t91 + 0.1e1;
t153 = -t112 * t94 - t95 * t98;
t118 = ((-t124 * t143 + t146 * t149) * t141 - t143 * t159) * t136;
t117 = -t128 * t124 + t127 * t149;
t116 = t128 * t123 + t127 * t150;
t114 = t139 * t121 + (t123 * t146 - t143 * t150) * t136;
t113 = t115 * t145 + t125 * t141;
t107 = -t121 * t158 + t151 * t123 - t126 * t150;
t105 = 0.1e1 / t112;
t104 = t128 * t135 * t141 + t117 * t145;
t103 = (-t126 * t124 + t149 * t151) * t141 - t126 * t159;
t100 = t108 * t145 + t147 * t141;
t96 = 0.1e1 / (t98 ^ 2 * t106 + 0.1e1);
t90 = 0.1e1 / t93;
t89 = 0.1e1 / t154;
t86 = 0.1e1 / t88;
t85 = 0.1e1 / (t101 ^ 2 * t87 + 0.1e1);
t84 = (-t103 * t105 + t118 * t163) * t96;
t83 = (-t105 * t107 + t114 * t163) * t96 * t141;
t82 = (-t100 * t105 + t113 * t163) * t96;
t1 = [0, t84, t83, 0, t82, 0; 0 ((t117 * t141 - t128 * t159) * t86 - (-t94 * t103 + t95 * t118 + t153 * t84) * t164) * t85 (t110 * t141 * t86 - (t153 * t83 + (-t107 * t94 + t114 * t95) * t141) * t164) * t85, 0 (t102 * t86 - (-t94 * t100 + t95 * t113 + t153 * t82) * t164) * t85, 0; 0 ((t104 * t140 - t116 * t144) * t90 - (t104 * t144 + t116 * t140) * t165) * t89 ((t140 * t161 - t144 * t148) * t90 - (t140 * t148 + t144 * t161) * t165) * t89, 0 (-t140 * t90 + t144 * t165) * t89 * t101, t154 * t89;];
Ja_rot  = t1;
