% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR11_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_jacobia_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:54:31
% EndTime: 2019-02-26 20:54:32
% DurationCPUTime: 0.47s
% Computational Cost: add. (1859->50), mult. (4121->117), div. (85->9), fcn. (5660->17), ass. (0->71)
t128 = cos(pkin(6));
t126 = cos(pkin(12));
t160 = sin(qJ(1));
t141 = t160 * t126;
t123 = sin(pkin(12));
t133 = cos(qJ(1));
t146 = t133 * t123;
t117 = t128 * t146 + t141;
t130 = sin(qJ(3));
t132 = cos(qJ(3));
t142 = t160 * t123;
t145 = t133 * t126;
t116 = -t128 * t145 + t142;
t124 = sin(pkin(7));
t127 = cos(pkin(7));
t125 = sin(pkin(6));
t148 = t125 * t133;
t138 = t116 * t127 + t124 * t148;
t106 = -t117 * t132 + t138 * t130;
t113 = -t116 * t124 + t127 * t148;
t122 = pkin(13) + qJ(5);
t120 = sin(t122);
t121 = cos(t122);
t162 = t106 * t120 - t113 * t121;
t96 = t106 * t121 + t113 * t120;
t137 = t128 * t141 + t146;
t143 = t125 * t160;
t161 = -t124 * t143 + t137 * t127;
t147 = t126 * t127;
t149 = t124 * t128;
t112 = t130 * t149 + (t123 * t132 + t130 * t147) * t125;
t115 = -t125 * t126 * t124 + t128 * t127;
t101 = t112 * t120 - t115 * t121;
t90 = atan2(t162, t101);
t85 = sin(t90);
t86 = cos(t90);
t83 = t86 * t101 + t162 * t85;
t82 = 0.1e1 / t83 ^ 2;
t118 = -t128 * t142 + t145;
t108 = t118 * t132 - t161 * t130;
t134 = t137 * t124 + t127 * t143;
t97 = t108 * t120 - t134 * t121;
t159 = t82 * t97;
t158 = t86 * t162;
t131 = cos(qJ(6));
t107 = t118 * t130 + t161 * t132;
t129 = sin(qJ(6));
t154 = t107 * t129;
t98 = t108 * t121 + t134 * t120;
t92 = t98 * t131 + t154;
t89 = 0.1e1 / t92 ^ 2;
t153 = t107 * t131;
t91 = t98 * t129 - t153;
t157 = t89 * t91;
t156 = t97 ^ 2 * t82;
t100 = 0.1e1 / t101 ^ 2;
t155 = t100 * t162;
t144 = t91 ^ 2 * t89 + 0.1e1;
t139 = -t101 * t85 + t158;
t135 = -t117 * t130 - t138 * t132;
t111 = t132 * t149 + (-t123 * t130 + t132 * t147) * t125;
t102 = t112 * t121 + t115 * t120;
t99 = 0.1e1 / t101;
t88 = 0.1e1 / t92;
t87 = 0.1e1 / (t100 * t162 ^ 2 + 0.1e1);
t84 = 0.1e1 / t144;
t81 = 0.1e1 / t83;
t80 = 0.1e1 / (0.1e1 + t156);
t79 = (-t111 * t155 - t135 * t99) * t87 * t120;
t78 = (-t102 * t155 + t96 * t99) * t87;
t1 = [-t97 * t99 * t87, 0, t79, 0, t78, 0; (t162 * t81 - (-t85 + (-t99 * t158 + t85) * t87) * t156) * t80, 0 (-t107 * t120 * t81 - (t139 * t79 + (t111 * t86 - t135 * t85) * t120) * t159) * t80, 0 (t98 * t81 - (t86 * t102 + t139 * t78 + t85 * t96) * t159) * t80, 0; ((t96 * t129 - t131 * t135) * t88 - (t129 * t135 + t96 * t131) * t157) * t84, 0 ((-t108 * t131 - t121 * t154) * t88 - (t108 * t129 - t121 * t153) * t157) * t84, 0 (-t129 * t88 + t131 * t157) * t97 * t84, t144 * t84;];
Ja_rot  = t1;
