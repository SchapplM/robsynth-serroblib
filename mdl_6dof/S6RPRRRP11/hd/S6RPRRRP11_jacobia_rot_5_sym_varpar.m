% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRP11_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:13:41
% EndTime: 2019-02-26 21:13:42
% DurationCPUTime: 0.38s
% Computational Cost: add. (1440->49), mult. (4121->118), div. (85->9), fcn. (5660->17), ass. (0->69)
t119 = cos(pkin(6));
t117 = cos(pkin(12));
t152 = sin(qJ(1));
t134 = t152 * t117;
t114 = sin(pkin(12));
t126 = cos(qJ(1));
t139 = t126 * t114;
t111 = t119 * t139 + t134;
t122 = sin(qJ(3));
t125 = cos(qJ(3));
t135 = t152 * t114;
t138 = t126 * t117;
t110 = -t119 * t138 + t135;
t115 = sin(pkin(7));
t118 = cos(pkin(7));
t116 = sin(pkin(6));
t141 = t116 * t126;
t131 = t110 * t118 + t115 * t141;
t100 = -t111 * t125 + t131 * t122;
t107 = -t110 * t115 + t118 * t141;
t121 = sin(qJ(4));
t124 = cos(qJ(4));
t154 = t100 * t121 - t107 * t124;
t90 = t100 * t124 + t107 * t121;
t130 = t119 * t134 + t139;
t136 = t116 * t152;
t153 = -t115 * t136 + t130 * t118;
t140 = t117 * t118;
t142 = t115 * t119;
t106 = t122 * t142 + (t114 * t125 + t122 * t140) * t116;
t109 = -t116 * t117 * t115 + t119 * t118;
t95 = t106 * t121 - t109 * t124;
t86 = atan2(t154, t95);
t83 = sin(t86);
t84 = cos(t86);
t77 = t154 * t83 + t84 * t95;
t76 = 0.1e1 / t77 ^ 2;
t112 = -t119 * t135 + t138;
t102 = t112 * t125 - t153 * t122;
t127 = t130 * t115 + t118 * t136;
t91 = t102 * t121 - t124 * t127;
t151 = t76 * t91;
t101 = t112 * t122 + t153 * t125;
t120 = sin(qJ(5));
t123 = cos(qJ(5));
t92 = t102 * t124 + t121 * t127;
t82 = t101 * t120 + t92 * t123;
t80 = 0.1e1 / t82 ^ 2;
t81 = -t101 * t123 + t92 * t120;
t150 = t80 * t81;
t149 = t84 * t154;
t94 = 0.1e1 / t95 ^ 2;
t148 = t154 * t94;
t147 = t91 ^ 2 * t76;
t146 = t101 * t124;
t137 = t81 ^ 2 * t80 + 0.1e1;
t132 = -t83 * t95 + t149;
t128 = -t111 * t122 - t131 * t125;
t105 = t125 * t142 + (-t114 * t122 + t125 * t140) * t116;
t96 = t106 * t124 + t109 * t121;
t93 = 0.1e1 / t95;
t85 = 0.1e1 / (t154 ^ 2 * t94 + 0.1e1);
t79 = 0.1e1 / t82;
t78 = 0.1e1 / t137;
t75 = 0.1e1 / t77;
t74 = 0.1e1 / (0.1e1 + t147);
t73 = (-t105 * t148 - t128 * t93) * t85 * t121;
t72 = (-t96 * t148 + t90 * t93) * t85;
t1 = [-t91 * t93 * t85, 0, t73, t72, 0, 0; (t154 * t75 - (-t83 + (-t93 * t149 + t83) * t85) * t147) * t74, 0 (-t101 * t121 * t75 - (t132 * t73 + (t105 * t84 - t128 * t83) * t121) * t151) * t74 (t92 * t75 - (t132 * t72 + t83 * t90 + t84 * t96) * t151) * t74, 0, 0; ((t90 * t120 - t123 * t128) * t79 - (t120 * t128 + t90 * t123) * t150) * t78, 0 ((-t102 * t123 - t120 * t146) * t79 - (t102 * t120 - t123 * t146) * t150) * t78 (-t120 * t79 + t123 * t150) * t91 * t78, t137 * t78, 0;];
Ja_rot  = t1;
