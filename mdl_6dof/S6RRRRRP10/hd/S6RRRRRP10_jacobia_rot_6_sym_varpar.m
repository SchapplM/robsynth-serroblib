% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRP10_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:45:00
% EndTime: 2019-02-26 22:45:01
% DurationCPUTime: 0.34s
% Computational Cost: add. (1806->45), mult. (3178->109), div. (134->9), fcn. (4488->13), ass. (0->65)
t124 = cos(pkin(6));
t126 = sin(qJ(2));
t130 = cos(qJ(1));
t139 = t130 * t126;
t127 = sin(qJ(1));
t129 = cos(qJ(2));
t140 = t127 * t129;
t117 = t124 * t139 + t140;
t125 = sin(qJ(3));
t128 = cos(qJ(3));
t123 = sin(pkin(6));
t142 = t123 * t130;
t109 = -t117 * t128 + t125 * t142;
t122 = qJ(4) + qJ(5);
t120 = sin(t122);
t121 = cos(t122);
t138 = t130 * t129;
t141 = t127 * t126;
t136 = -t124 * t138 + t141;
t154 = t109 * t120 + t121 * t136;
t134 = t136 * t120;
t153 = t109 * t121 - t134;
t144 = t123 * t128;
t116 = t124 * t125 + t126 * t144;
t143 = t123 * t129;
t104 = t116 * t120 + t121 * t143;
t92 = atan2(t154, t104);
t89 = sin(t92);
t90 = cos(t92);
t88 = t90 * t104 + t154 * t89;
t87 = 0.1e1 / t88 ^ 2;
t118 = -t124 * t141 + t138;
t145 = t123 * t125;
t111 = t118 * t128 + t127 * t145;
t132 = t124 * t140 + t139;
t99 = t111 * t120 - t121 * t132;
t151 = t87 * t99;
t150 = t90 * t154;
t149 = t99 ^ 2 * t87;
t103 = 0.1e1 / t104 ^ 2;
t148 = t103 * t154;
t110 = -t118 * t125 + t127 * t144;
t100 = t111 * t121 + t120 * t132;
t95 = 0.1e1 / t100 ^ 2;
t147 = t110 ^ 2 * t95;
t146 = t110 * t95;
t93 = 0.1e1 / (0.1e1 + t147);
t137 = t99 * t93 * t146;
t135 = -t104 * t89 + t150;
t133 = t117 * t125 + t128 * t142;
t131 = t128 * t132;
t115 = t124 * t128 - t126 * t145;
t112 = (t120 * t128 * t129 - t121 * t126) * t123;
t105 = t116 * t121 - t120 * t143;
t102 = 0.1e1 / t104;
t101 = -t117 * t121 - t128 * t134;
t94 = 0.1e1 / t100;
t91 = 0.1e1 / (t103 * t154 ^ 2 + 0.1e1);
t86 = 0.1e1 / t88;
t85 = 0.1e1 / (0.1e1 + t149);
t84 = (t102 * t133 - t115 * t148) * t91 * t120;
t83 = (-t101 * t102 - t112 * t148) * t91;
t82 = (t102 * t153 - t105 * t148) * t91;
t81 = (t100 * t86 - (t90 * t105 + t135 * t82 + t153 * t89) * t151) * t85;
t1 = [-t99 * t102 * t91, t83, t84, t82, t82, 0; (t154 * t86 - (-t89 + (-t102 * t150 + t89) * t91) * t149) * t85 ((-t118 * t121 - t120 * t131) * t86 - (-t89 * t101 + t90 * t112 + t135 * t83) * t151) * t85 (t110 * t120 * t86 - (t135 * t84 + (t115 * t90 + t133 * t89) * t120) * t151) * t85, t81, t81, 0; (t133 * t94 - t153 * t146) * t93 (t132 * t125 * t94 - (t118 * t120 - t121 * t131) * t146) * t93 (-t111 * t94 - t121 * t147) * t93, t137, t137, 0;];
Ja_rot  = t1;
