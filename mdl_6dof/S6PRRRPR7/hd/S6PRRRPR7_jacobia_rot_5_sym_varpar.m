% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPR7_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_jacobia_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:04
% EndTime: 2019-02-26 20:14:05
% DurationCPUTime: 0.32s
% Computational Cost: add. (1437->59), mult. (4134->141), div. (90->9), fcn. (5673->17), ass. (0->75)
t120 = sin(pkin(6));
t126 = sin(qJ(3));
t127 = sin(qJ(2));
t129 = cos(qJ(3));
t130 = cos(qJ(2));
t123 = cos(pkin(7));
t142 = t123 * t126;
t119 = sin(pkin(7));
t124 = cos(pkin(6));
t146 = t119 * t124;
t109 = t126 * t146 + (t127 * t129 + t130 * t142) * t120;
t144 = t120 * t119;
t111 = t124 * t123 - t130 * t144;
t125 = sin(qJ(4));
t128 = cos(qJ(4));
t101 = t109 * t125 - t111 * t128;
t122 = cos(pkin(12));
t118 = sin(pkin(12));
t139 = t124 * t130;
t133 = -t118 * t127 + t122 * t139;
t143 = t120 * t123;
t132 = -t119 * t133 - t122 * t143;
t140 = t124 * t127;
t112 = t118 * t130 + t122 * t140;
t131 = -t122 * t144 + t123 * t133;
t98 = t112 * t129 + t126 * t131;
t88 = t98 * t125 - t128 * t132;
t87 = atan2(-t88, t101);
t84 = sin(t87);
t85 = cos(t87);
t78 = t85 * t101 - t84 * t88;
t77 = 0.1e1 / t78 ^ 2;
t113 = -t118 * t139 - t122 * t127;
t114 = -t118 * t140 + t122 * t130;
t136 = t118 * t144;
t100 = t114 * t129 + (t113 * t123 + t136) * t126;
t134 = -t113 * t119 + t118 * t143;
t91 = t100 * t125 - t128 * t134;
t151 = t77 * t91;
t121 = cos(pkin(13));
t117 = sin(pkin(13));
t141 = t123 * t129;
t99 = -t113 * t141 + t114 * t126 - t129 * t136;
t148 = t99 * t117;
t92 = t100 * t128 + t125 * t134;
t83 = t92 * t121 + t148;
t81 = 0.1e1 / t83 ^ 2;
t147 = t99 * t121;
t82 = t92 * t117 - t147;
t150 = t81 * t82;
t96 = 0.1e1 / t101 ^ 2;
t149 = t88 * t96;
t145 = t119 * t128;
t138 = t126 * t127;
t137 = t129 * t130;
t135 = -t101 * t84 - t85 * t88;
t108 = t129 * t146 + (t123 * t137 - t138) * t120;
t105 = ((-t123 * t138 + t137) * t125 - t127 * t145) * t120;
t104 = t113 * t129 - t114 * t142;
t103 = t113 * t126 + t114 * t141;
t102 = t109 * t128 + t111 * t125;
t97 = -t112 * t126 + t129 * t131;
t95 = 0.1e1 / t101;
t94 = t114 * t119 * t125 + t104 * t128;
t93 = (-t112 * t142 + t129 * t133) * t125 - t112 * t145;
t90 = t125 * t132 + t98 * t128;
t86 = 0.1e1 / (t88 ^ 2 * t96 + 0.1e1);
t80 = 0.1e1 / t83;
t79 = 0.1e1 / (t82 ^ 2 * t81 + 0.1e1);
t76 = 0.1e1 / t78;
t75 = 0.1e1 / (t91 ^ 2 * t77 + 0.1e1);
t74 = (t108 * t149 - t95 * t97) * t86 * t125;
t73 = (t105 * t149 - t93 * t95) * t86;
t72 = (t102 * t149 - t90 * t95) * t86;
t1 = [0, t73, t74, t72, 0, 0; 0 ((t104 * t125 - t114 * t145) * t76 - (t85 * t105 + t135 * t73 - t84 * t93) * t151) * t75 (-t99 * t125 * t76 - (t135 * t74 + (t108 * t85 - t84 * t97) * t125) * t151) * t75 (t92 * t76 - (t85 * t102 + t135 * t72 - t84 * t90) * t151) * t75, 0, 0; 0 ((-t103 * t121 + t94 * t117) * t80 - (t103 * t117 + t94 * t121) * t150) * t79 ((-t100 * t121 - t128 * t148) * t80 - (t100 * t117 - t128 * t147) * t150) * t79 (-t117 * t80 + t121 * t150) * t91 * t79, 0, 0;];
Ja_rot  = t1;
