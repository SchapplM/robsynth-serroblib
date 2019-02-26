% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR11_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:21:50
% EndTime: 2019-02-26 22:21:51
% DurationCPUTime: 0.52s
% Computational Cost: add. (1357->47), mult. (3691->106), div. (115->9), fcn. (5181->15), ass. (0->66)
t130 = cos(pkin(6));
t135 = sin(qJ(1));
t139 = cos(qJ(2));
t148 = t135 * t139;
t134 = sin(qJ(2));
t140 = cos(qJ(1));
t149 = t134 * t140;
t123 = t130 * t149 + t148;
t133 = sin(qJ(3));
t138 = cos(qJ(3));
t129 = sin(pkin(6));
t151 = t129 * t140;
t117 = -t123 * t138 + t133 * t151;
t132 = sin(qJ(5));
t137 = cos(qJ(5));
t141 = t123 * t133 + t138 * t151;
t102 = t117 * t137 - t141 * t132;
t147 = t139 * t140;
t150 = t134 * t135;
t125 = -t130 * t150 + t147;
t153 = t129 * t133;
t118 = t125 * t138 + t135 * t153;
t152 = t129 * t138;
t145 = t125 * t133 - t135 * t152;
t105 = t118 * t137 + t145 * t132;
t120 = -t130 * t138 + t134 * t153;
t121 = t130 * t133 + t134 * t152;
t112 = t120 * t132 + t121 * t137;
t111 = -t120 * t137 + t121 * t132;
t164 = t117 * t132 + t141 * t137;
t93 = atan2(t164, t111);
t91 = cos(t93);
t161 = t91 * t164;
t90 = sin(t93);
t144 = -t111 * t90 + t161;
t104 = t118 * t132 - t145 * t137;
t88 = t111 * t91 + t164 * t90;
t87 = 0.1e1 / t88 ^ 2;
t159 = t104 * t87;
t108 = 0.1e1 / t111;
t109 = 0.1e1 / t111 ^ 2;
t157 = t109 * t164;
t92 = 0.1e1 / (t109 * t164 ^ 2 + 0.1e1);
t169 = (t102 * t108 - t112 * t157) * t92;
t158 = t104 ^ 2 * t87;
t85 = 0.1e1 / (0.1e1 + t158);
t86 = 0.1e1 / t88;
t172 = (-t105 * t86 + (t102 * t90 + t91 * t112 + t144 * t169) * t159) * t85;
t131 = sin(qJ(6));
t136 = cos(qJ(6));
t124 = -t130 * t148 - t149;
t97 = t105 * t136 + t124 * t131;
t95 = 0.1e1 / t97 ^ 2;
t96 = t105 * t131 - t124 * t136;
t160 = t95 * t96;
t146 = t95 * t96 ^ 2 + 0.1e1;
t89 = 0.1e1 / t146;
t94 = 0.1e1 / t97;
t162 = (-t131 * t94 + t136 * t160) * t89 * t104;
t143 = t132 * t138 - t133 * t137;
t122 = -t130 * t147 + t150;
t119 = t143 * t139 * t129;
t107 = (t132 * t133 + t137 * t138) * t124;
t106 = t143 * t122;
t84 = (t106 * t108 - t119 * t157) * t92;
t1 = [-t104 * t108 * t92, t84, -t169, 0, t169, 0; (t164 * t86 - (-t90 + (-t108 * t161 + t90) * t92) * t158) * t85 (-(t90 * t106 + t91 * t119 + t144 * t84) * t159 + t143 * t86 * t124) * t85, t172, 0, -t172, 0; ((t102 * t131 - t122 * t136) * t94 - (t102 * t136 + t122 * t131) * t160) * t89 ((t107 * t131 + t125 * t136) * t94 - (t107 * t136 - t125 * t131) * t160) * t89, -t162, 0, t162, t146 * t89;];
Ja_rot  = t1;
