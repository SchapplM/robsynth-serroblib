% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PPRRRP2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:42:05
% EndTime: 2019-02-26 19:42:05
% DurationCPUTime: 0.40s
% Computational Cost: add. (2130->55), mult. (6124->131), div. (87->9), fcn. (8354->17), ass. (0->74)
t153 = cos(qJ(3));
t120 = sin(pkin(12));
t122 = sin(pkin(7));
t123 = sin(pkin(6));
t127 = cos(pkin(6));
t130 = sin(qJ(3));
t124 = cos(pkin(12));
t126 = cos(pkin(7));
t144 = t124 * t126;
t113 = t127 * t122 * t130 + (t153 * t120 + t130 * t144) * t123;
t146 = t123 * t122;
t117 = -t124 * t146 + t127 * t126;
t129 = sin(qJ(4));
t132 = cos(qJ(4));
t111 = t113 * t132 + t117 * t129;
t141 = t122 * t153;
t112 = -t127 * t141 + (t120 * t130 - t153 * t144) * t123;
t128 = sin(qJ(5));
t131 = cos(qJ(5));
t103 = t111 * t128 - t112 * t131;
t125 = cos(pkin(11));
t121 = sin(pkin(11));
t143 = t125 * t127;
t138 = -t121 * t120 + t124 * t143;
t135 = t138 * t126;
t136 = t120 * t143 + t121 * t124;
t107 = t136 * t153 + (-t125 * t146 + t135) * t130;
t145 = t123 * t126;
t114 = -t138 * t122 - t125 * t145;
t100 = t107 * t132 + t114 * t129;
t140 = t123 * t141;
t133 = t125 * t140 + t136 * t130 - t153 * t135;
t90 = t100 * t128 - t133 * t131;
t87 = atan2(-t90, t103);
t83 = sin(t87);
t84 = cos(t87);
t82 = t84 * t103 - t83 * t90;
t81 = 0.1e1 / t82 ^ 2;
t147 = t121 * t127;
t118 = -t120 * t147 + t125 * t124;
t137 = t125 * t120 + t124 * t147;
t134 = t137 * t126;
t109 = t118 * t153 + (t121 * t146 - t134) * t130;
t115 = t121 * t145 + t137 * t122;
t102 = t109 * t132 + t115 * t129;
t108 = t118 * t130 - t121 * t140 + t153 * t134;
t148 = t108 * t131;
t93 = t102 * t128 - t148;
t152 = t81 * t93;
t98 = 0.1e1 / t103 ^ 2;
t151 = t90 * t98;
t101 = -t109 * t129 + t115 * t132;
t94 = t102 * t131 + t108 * t128;
t89 = 0.1e1 / t94 ^ 2;
t150 = t101 ^ 2 * t89;
t149 = t101 * t89;
t142 = t128 * t132;
t139 = -t103 * t83 - t84 * t90;
t110 = -t113 * t129 + t117 * t132;
t105 = -t112 * t142 - t113 * t131;
t104 = t111 * t131 + t112 * t128;
t99 = -t107 * t129 + t114 * t132;
t97 = 0.1e1 / t103;
t95 = -t107 * t131 - t133 * t142;
t92 = t100 * t131 + t133 * t128;
t88 = 0.1e1 / t94;
t86 = 0.1e1 / (t90 ^ 2 * t98 + 0.1e1);
t85 = 0.1e1 / (0.1e1 + t150);
t80 = 0.1e1 / t82;
t79 = 0.1e1 / (t93 ^ 2 * t81 + 0.1e1);
t78 = (t110 * t151 - t97 * t99) * t86 * t128;
t77 = (t105 * t151 - t95 * t97) * t86;
t76 = (t104 * t151 - t92 * t97) * t86;
t1 = [0, 0, t77, t78, t76, 0; 0, 0 ((-t108 * t142 - t109 * t131) * t80 - (t84 * t105 + t139 * t77 - t83 * t95) * t152) * t79 (t101 * t128 * t80 - (t139 * t78 + (t110 * t84 - t83 * t99) * t128) * t152) * t79 (t94 * t80 - (t84 * t104 + t139 * t76 - t83 * t92) * t152) * t79, 0; 0, 0 (t108 * t129 * t88 - (t109 * t128 - t132 * t148) * t149) * t85 (-t102 * t88 - t131 * t150) * t85, t93 * t85 * t149, 0;];
Ja_rot  = t1;
