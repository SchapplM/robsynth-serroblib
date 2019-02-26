% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PPRRRR1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_jacobia_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:42:43
% EndTime: 2019-02-26 19:42:43
% DurationCPUTime: 0.34s
% Computational Cost: add. (1997->44), mult. (4405->105), div. (95->9), fcn. (6043->17), ass. (0->69)
t120 = sin(pkin(12));
t125 = cos(pkin(7));
t119 = sin(pkin(13));
t123 = cos(pkin(13));
t124 = cos(pkin(12));
t126 = cos(pkin(6));
t145 = t120 * t126;
t135 = t124 * t119 + t123 * t145;
t121 = sin(pkin(7));
t122 = sin(pkin(6));
t143 = t122 * t121;
t151 = -t120 * t143 + t135 * t125;
t128 = sin(qJ(3));
t130 = cos(qJ(3));
t141 = t123 * t125;
t144 = t121 * t126;
t110 = t128 * t144 + (t119 * t130 + t128 * t141) * t122;
t112 = -t123 * t143 + t126 * t125;
t118 = qJ(4) + qJ(5);
t116 = sin(t118);
t117 = cos(t118);
t101 = t110 * t116 - t112 * t117;
t140 = t124 * t126;
t113 = t119 * t140 + t120 * t123;
t136 = -t120 * t119 + t123 * t140;
t132 = -t124 * t143 + t125 * t136;
t104 = t113 * t130 + t128 * t132;
t142 = t122 * t125;
t133 = -t121 * t136 - t124 * t142;
t94 = t104 * t116 - t117 * t133;
t93 = atan2(-t94, t101);
t90 = sin(t93);
t91 = cos(t93);
t84 = t91 * t101 - t90 * t94;
t83 = 0.1e1 / t84 ^ 2;
t114 = -t119 * t145 + t124 * t123;
t106 = t114 * t130 - t151 * t128;
t131 = t120 * t142 + t121 * t135;
t97 = t106 * t116 - t117 * t131;
t150 = t83 * t97;
t129 = cos(qJ(6));
t105 = t114 * t128 + t151 * t130;
t127 = sin(qJ(6));
t147 = t105 * t127;
t98 = t106 * t117 + t116 * t131;
t89 = t98 * t129 + t147;
t87 = 0.1e1 / t89 ^ 2;
t146 = t105 * t129;
t88 = t98 * t127 - t146;
t149 = t87 * t88;
t100 = 0.1e1 / t101 ^ 2;
t148 = t100 * t94;
t139 = t88 ^ 2 * t87 + 0.1e1;
t137 = -t101 * t90 - t91 * t94;
t109 = t130 * t144 + (-t119 * t128 + t130 * t141) * t122;
t103 = -t113 * t128 + t130 * t132;
t102 = t110 * t117 + t112 * t116;
t99 = 0.1e1 / t101;
t96 = t104 * t117 + t116 * t133;
t92 = 0.1e1 / (t94 ^ 2 * t100 + 0.1e1);
t86 = 0.1e1 / t89;
t85 = 0.1e1 / t139;
t82 = 0.1e1 / t84;
t81 = 0.1e1 / (t97 ^ 2 * t83 + 0.1e1);
t80 = (-t103 * t99 + t109 * t148) * t92 * t116;
t79 = (t102 * t148 - t96 * t99) * t92;
t78 = (-t127 * t86 + t129 * t149) * t97 * t85;
t77 = (t98 * t82 - (t91 * t102 + t137 * t79 - t90 * t96) * t150) * t81;
t1 = [0, 0, t80, t79, t79, 0; 0, 0 (-t105 * t116 * t82 - (t137 * t80 + (-t103 * t90 + t109 * t91) * t116) * t150) * t81, t77, t77, 0; 0, 0 ((-t106 * t129 - t117 * t147) * t86 - (t106 * t127 - t117 * t146) * t149) * t85, t78, t78, t139 * t85;];
Ja_rot  = t1;
