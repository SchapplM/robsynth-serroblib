% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Ja_rot = S6RRRRPR14_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_jacobia_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:11
% EndTime: 2019-02-26 22:38:12
% DurationCPUTime: 0.56s
% Computational Cost: add. (1829->65), mult. (5250->153), div. (110->9), fcn. (7214->17), ass. (0->79)
t138 = cos(pkin(6));
t144 = cos(qJ(2));
t174 = sin(qJ(1));
t153 = t174 * t144;
t141 = sin(qJ(2));
t145 = cos(qJ(1));
t157 = t145 * t141;
t129 = t138 * t157 + t153;
t140 = sin(qJ(3));
t143 = cos(qJ(3));
t154 = t174 * t141;
t156 = t145 * t144;
t128 = -t138 * t156 + t154;
t134 = sin(pkin(7));
t137 = cos(pkin(7));
t135 = sin(pkin(6));
t161 = t135 * t145;
t150 = t128 * t137 + t134 * t161;
t115 = -t129 * t143 + t140 * t150;
t125 = -t128 * t134 + t137 * t161;
t139 = sin(qJ(4));
t142 = cos(qJ(4));
t175 = t115 * t139 - t125 * t142;
t103 = t115 * t142 + t125 * t139;
t149 = t138 * t153 + t157;
t155 = t135 * t174;
t152 = t134 * t155;
t130 = -t138 * t154 + t156;
t164 = t130 * t143;
t117 = t164 + (-t137 * t149 + t152) * t140;
t146 = t134 * t149 + t137 * t155;
t105 = t117 * t142 + t139 * t146;
t148 = t149 * t143;
t116 = t130 * t140 + t137 * t148 - t143 * t152;
t133 = sin(pkin(13));
t136 = cos(pkin(13));
t95 = t105 * t136 + t116 * t133;
t93 = 0.1e1 / t95 ^ 2;
t94 = t105 * t133 - t116 * t136;
t173 = t93 * t94;
t160 = t137 * t140;
t163 = t134 * t138;
t124 = t140 * t163 + (t141 * t143 + t144 * t160) * t135;
t127 = -t135 * t144 * t134 + t138 * t137;
t110 = t124 * t139 - t127 * t142;
t99 = atan2(t175, t110);
t97 = cos(t99);
t172 = t175 * t97;
t104 = t117 * t139 - t142 * t146;
t96 = sin(t99);
t90 = t97 * t110 + t175 * t96;
t89 = 0.1e1 / t90 ^ 2;
t171 = t104 * t89;
t170 = t104 ^ 2 * t89;
t109 = 0.1e1 / t110 ^ 2;
t169 = t175 * t109;
t168 = t116 * t142;
t162 = t134 * t142;
t159 = t140 * t141;
t158 = t143 * t144;
t151 = -t110 * t96 + t172;
t147 = -t129 * t140 - t143 * t150;
t123 = t143 * t163 + (t137 * t158 - t159) * t135;
t120 = ((-t137 * t159 + t158) * t139 - t141 * t162) * t135;
t119 = -t130 * t160 - t148;
t118 = t137 * t164 - t140 * t149;
t111 = t124 * t142 + t127 * t139;
t108 = 0.1e1 / t110;
t107 = t130 * t134 * t139 + t119 * t142;
t106 = (-t128 * t143 - t129 * t160) * t139 - t129 * t162;
t98 = 0.1e1 / (t109 * t175 ^ 2 + 0.1e1);
t92 = 0.1e1 / t95;
t91 = 0.1e1 / (t94 ^ 2 * t93 + 0.1e1);
t88 = 0.1e1 / t90;
t87 = 0.1e1 / (0.1e1 + t170);
t86 = (-t108 * t147 - t123 * t169) * t98 * t139;
t85 = (-t106 * t108 - t120 * t169) * t98;
t84 = (t103 * t108 - t111 * t169) * t98;
t1 = [-t104 * t108 * t98, t85, t86, t84, 0, 0; (t175 * t88 - (-t96 + (-t108 * t172 + t96) * t98) * t170) * t87 ((t119 * t139 - t130 * t162) * t88 - (-t96 * t106 + t97 * t120 + t151 * t85) * t171) * t87 (-t116 * t139 * t88 - (t151 * t86 + (t123 * t97 - t147 * t96) * t139) * t171) * t87 (t105 * t88 - (t103 * t96 + t97 * t111 + t151 * t84) * t171) * t87, 0, 0; ((t103 * t133 - t136 * t147) * t92 - (t103 * t136 + t133 * t147) * t173) * t91 ((t107 * t133 - t118 * t136) * t92 - (t107 * t136 + t118 * t133) * t173) * t91 ((-t117 * t136 - t133 * t168) * t92 - (t117 * t133 - t136 * t168) * t173) * t91 (-t133 * t92 + t136 * t173) * t91 * t104, 0, 0;];
Ja_rot  = t1;
