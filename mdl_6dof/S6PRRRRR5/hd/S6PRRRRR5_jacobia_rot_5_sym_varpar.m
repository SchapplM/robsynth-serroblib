% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRR5_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_jacobia_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:21:10
% EndTime: 2019-02-26 20:21:11
% DurationCPUTime: 0.40s
% Computational Cost: add. (1524->59), mult. (4378->143), div. (95->9), fcn. (6002->17), ass. (0->75)
t122 = sin(pkin(6));
t128 = sin(qJ(3));
t129 = sin(qJ(2));
t132 = cos(qJ(3));
t133 = cos(qJ(2));
t124 = cos(pkin(7));
t146 = t124 * t128;
t121 = sin(pkin(7));
t125 = cos(pkin(6));
t150 = t121 * t125;
t112 = t128 * t150 + (t129 * t132 + t133 * t146) * t122;
t148 = t122 * t121;
t114 = t125 * t124 - t133 * t148;
t127 = sin(qJ(4));
t131 = cos(qJ(4));
t104 = t112 * t127 - t114 * t131;
t120 = sin(pkin(13));
t123 = cos(pkin(13));
t144 = t125 * t129;
t115 = t120 * t133 + t123 * t144;
t143 = t125 * t133;
t136 = -t120 * t129 + t123 * t143;
t134 = -t123 * t148 + t136 * t124;
t101 = t115 * t132 + t134 * t128;
t147 = t122 * t124;
t135 = -t136 * t121 - t123 * t147;
t91 = t101 * t127 - t135 * t131;
t90 = atan2(-t91, t104);
t87 = sin(t90);
t88 = cos(t90);
t81 = t88 * t104 - t87 * t91;
t80 = 0.1e1 / t81 ^ 2;
t116 = -t120 * t143 - t123 * t129;
t117 = -t120 * t144 + t123 * t133;
t139 = t120 * t148;
t103 = t117 * t132 + (t116 * t124 + t139) * t128;
t137 = -t116 * t121 + t120 * t147;
t94 = t103 * t127 - t137 * t131;
t154 = t80 * t94;
t145 = t124 * t132;
t102 = -t116 * t145 + t117 * t128 - t132 * t139;
t126 = sin(qJ(5));
t130 = cos(qJ(5));
t95 = t103 * t131 + t137 * t127;
t86 = t102 * t126 + t95 * t130;
t84 = 0.1e1 / t86 ^ 2;
t85 = -t102 * t130 + t95 * t126;
t153 = t84 * t85;
t99 = 0.1e1 / t104 ^ 2;
t152 = t91 * t99;
t151 = t102 * t131;
t149 = t121 * t131;
t142 = t128 * t129;
t141 = t132 * t133;
t140 = t85 ^ 2 * t84 + 0.1e1;
t138 = -t104 * t87 - t88 * t91;
t111 = t132 * t150 + (t124 * t141 - t142) * t122;
t108 = ((-t124 * t142 + t141) * t127 - t129 * t149) * t122;
t107 = t116 * t132 - t117 * t146;
t106 = t116 * t128 + t117 * t145;
t105 = t112 * t131 + t114 * t127;
t100 = -t115 * t128 + t134 * t132;
t98 = 0.1e1 / t104;
t97 = t117 * t121 * t127 + t107 * t131;
t96 = (-t115 * t146 + t136 * t132) * t127 - t115 * t149;
t93 = t101 * t131 + t135 * t127;
t89 = 0.1e1 / (t91 ^ 2 * t99 + 0.1e1);
t83 = 0.1e1 / t86;
t82 = 0.1e1 / t140;
t79 = 0.1e1 / t81;
t78 = 0.1e1 / (t94 ^ 2 * t80 + 0.1e1);
t77 = (-t100 * t98 + t111 * t152) * t89 * t127;
t76 = (t108 * t152 - t96 * t98) * t89;
t75 = (t105 * t152 - t93 * t98) * t89;
t1 = [0, t76, t77, t75, 0, 0; 0 ((t107 * t127 - t117 * t149) * t79 - (t88 * t108 + t138 * t76 - t87 * t96) * t154) * t78 (-t102 * t127 * t79 - (t138 * t77 + (-t100 * t87 + t111 * t88) * t127) * t154) * t78 (t95 * t79 - (t88 * t105 + t138 * t75 - t87 * t93) * t154) * t78, 0, 0; 0 ((-t106 * t130 + t97 * t126) * t83 - (t106 * t126 + t97 * t130) * t153) * t82 ((-t103 * t130 - t126 * t151) * t83 - (t103 * t126 - t130 * t151) * t153) * t82 (-t126 * t83 + t130 * t153) * t94 * t82, t140 * t82, 0;];
Ja_rot  = t1;
