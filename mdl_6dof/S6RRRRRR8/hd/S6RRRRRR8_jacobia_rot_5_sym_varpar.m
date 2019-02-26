% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRR8_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_jacobia_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:51:32
% EndTime: 2019-02-26 22:51:33
% DurationCPUTime: 0.27s
% Computational Cost: add. (965->52), mult. (2645->114), div. (90->9), fcn. (3627->15), ass. (0->70)
t151 = sin(qJ(1));
t119 = sin(pkin(6));
t122 = sin(qJ(3));
t123 = sin(qJ(2));
t124 = cos(qJ(3));
t125 = cos(qJ(2));
t120 = cos(pkin(7));
t140 = t120 * t124;
t118 = sin(pkin(7));
t121 = cos(pkin(6));
t143 = t118 * t121;
t102 = -t124 * t143 + (t122 * t123 - t125 * t140) * t119;
t132 = t151 * t123;
t126 = cos(qJ(1));
t136 = t126 * t125;
t109 = -t121 * t136 + t132;
t142 = t119 * t126;
t134 = t118 * t142;
t131 = t151 * t125;
t137 = t126 * t123;
t110 = t121 * t137 + t131;
t145 = t110 * t122;
t92 = t109 * t140 + t124 * t134 + t145;
t91 = atan2(-t92, t102);
t88 = sin(t91);
t89 = cos(t91);
t82 = t89 * t102 - t88 * t92;
t81 = 0.1e1 / t82 ^ 2;
t111 = -t121 * t132 + t136;
t128 = t121 * t131 + t137;
t127 = t128 * t124;
t133 = t119 * t151;
t130 = t118 * t133;
t96 = t111 * t122 + t120 * t127 - t124 * t130;
t150 = t81 * t96;
t105 = t128 * t118 + t120 * t133;
t117 = qJ(4) + qJ(5);
t115 = sin(t117);
t116 = cos(t117);
t97 = t111 * t124 + (-t128 * t120 + t130) * t122;
t87 = t105 * t115 + t97 * t116;
t85 = 0.1e1 / t87 ^ 2;
t86 = -t105 * t116 + t97 * t115;
t149 = t85 * t86;
t148 = t89 * t92;
t147 = t96 ^ 2 * t81;
t101 = 0.1e1 / t102 ^ 2;
t146 = t101 * t92;
t144 = t111 * t118;
t141 = t120 * t122;
t139 = t122 * t125;
t138 = t123 * t124;
t135 = t86 ^ 2 * t85 + 0.1e1;
t129 = -t102 * t88 - t148;
t95 = t109 * t141 - t110 * t124 + t122 * t134;
t108 = (t120 * t138 + t139) * t119;
t104 = -t109 * t118 + t120 * t142;
t103 = t122 * t143 + (t120 * t139 + t138) * t119;
t100 = 0.1e1 / t102;
t99 = -t111 * t141 - t127;
t98 = -t109 * t122 + t110 * t140;
t90 = 0.1e1 / (t92 ^ 2 * t101 + 0.1e1);
t84 = 0.1e1 / t87;
t83 = 0.1e1 / t135;
t80 = 0.1e1 / t82;
t79 = 0.1e1 / (0.1e1 + t147);
t78 = (-t100 * t98 + t108 * t146) * t90;
t77 = (t100 * t95 + t103 * t146) * t90;
t76 = t135 * t83;
t1 = [-t96 * t100 * t90, t78, t77, 0, 0, 0; ((-t145 + (-t109 * t120 - t134) * t124) * t80 - (-t88 + (t100 * t148 + t88) * t90) * t147) * t79 ((t111 * t140 - t128 * t122) * t80 - (t89 * t108 + t129 * t78 - t88 * t98) * t150) * t79 (t97 * t80 - (t89 * t103 + t129 * t77 + t88 * t95) * t150) * t79, 0, 0, 0; ((-t104 * t116 + t95 * t115) * t84 - (t104 * t115 + t95 * t116) * t149) * t83 ((t99 * t115 - t116 * t144) * t84 - (t115 * t144 + t99 * t116) * t149) * t83 (-t115 * t84 + t116 * t149) * t96 * t83, t76, t76, 0;];
Ja_rot  = t1;
