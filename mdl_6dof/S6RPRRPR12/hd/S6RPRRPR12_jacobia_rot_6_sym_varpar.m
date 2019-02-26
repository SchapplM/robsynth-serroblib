% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR12_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:07:16
% EndTime: 2019-02-26 21:07:16
% DurationCPUTime: 0.45s
% Computational Cost: add. (1440->49), mult. (4121->117), div. (85->9), fcn. (5660->17), ass. (0->70)
t127 = cos(pkin(6));
t122 = sin(pkin(12));
t135 = cos(qJ(1));
t144 = t135 * t122;
t125 = cos(pkin(12));
t131 = sin(qJ(1));
t145 = t131 * t125;
t118 = t127 * t144 + t145;
t130 = sin(qJ(3));
t134 = cos(qJ(3));
t143 = t135 * t125;
t146 = t131 * t122;
t117 = -t127 * t143 + t146;
t123 = sin(pkin(7));
t126 = cos(pkin(7));
t124 = sin(pkin(6));
t148 = t124 * t135;
t139 = t117 * t126 + t123 * t148;
t107 = -t118 * t134 + t130 * t139;
t112 = -t117 * t123 + t126 * t148;
t129 = sin(qJ(4));
t133 = cos(qJ(4));
t97 = t107 * t129 - t112 * t133;
t162 = t107 * t133 + t112 * t129;
t138 = t127 * t145 + t144;
t149 = t124 * t131;
t159 = -t123 * t149 + t126 * t138;
t147 = t125 * t126;
t150 = t123 * t127;
t111 = t130 * t150 + (t122 * t134 + t130 * t147) * t124;
t116 = -t123 * t124 * t125 + t126 * t127;
t103 = t111 * t133 + t116 * t129;
t93 = atan2(t162, t103);
t90 = sin(t93);
t91 = cos(t93);
t84 = t103 * t91 + t162 * t90;
t83 = 0.1e1 / t84 ^ 2;
t119 = -t127 * t146 + t143;
t109 = t119 * t134 - t130 * t159;
t114 = t123 * t138 + t126 * t149;
t99 = t109 * t133 + t114 * t129;
t158 = t83 * t99;
t128 = sin(qJ(6));
t108 = t119 * t130 + t134 * t159;
t132 = cos(qJ(6));
t152 = t108 * t132;
t98 = t109 * t129 - t114 * t133;
t89 = t128 * t98 + t152;
t87 = 0.1e1 / t89 ^ 2;
t153 = t108 * t128;
t88 = -t132 * t98 + t153;
t157 = t87 * t88;
t156 = t91 * t162;
t155 = t99 ^ 2 * t83;
t101 = 0.1e1 / t103 ^ 2;
t154 = t101 * t162;
t142 = t87 * t88 ^ 2 + 0.1e1;
t140 = -t103 * t90 + t156;
t136 = -t118 * t130 - t134 * t139;
t110 = t134 * t150 + (-t122 * t130 + t134 * t147) * t124;
t102 = -t111 * t129 + t116 * t133;
t100 = 0.1e1 / t103;
t92 = 0.1e1 / (t101 * t162 ^ 2 + 0.1e1);
t86 = 0.1e1 / t89;
t85 = 0.1e1 / t142;
t82 = 0.1e1 / t84;
t81 = 0.1e1 / (0.1e1 + t155);
t80 = (-t100 * t136 - t110 * t154) * t92 * t133;
t79 = (-t100 * t97 - t102 * t154) * t92;
t1 = [-t99 * t100 * t92, 0, t80, t79, 0, 0; (t162 * t82 - (-t90 + (-t100 * t156 + t90) * t92) * t155) * t81, 0 (-t108 * t133 * t82 - (t140 * t80 + (t110 * t91 - t136 * t90) * t133) * t158) * t81 (-t98 * t82 - (t91 * t102 + t140 * t79 - t90 * t97) * t158) * t81, 0, 0; ((t128 * t136 - t132 * t97) * t86 - (t128 * t97 + t132 * t136) * t157) * t85, 0 ((t109 * t128 + t129 * t152) * t86 - (t109 * t132 - t129 * t153) * t157) * t85 (-t128 * t157 - t132 * t86) * t99 * t85, 0, t142 * t85;];
Ja_rot  = t1;
