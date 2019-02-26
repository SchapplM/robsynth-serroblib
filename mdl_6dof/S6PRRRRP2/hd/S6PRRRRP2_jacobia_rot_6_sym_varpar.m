% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRP2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:15:42
% EndTime: 2019-02-26 20:15:43
% DurationCPUTime: 0.25s
% Computational Cost: add. (1533->43), mult. (2650->105), div. (117->9), fcn. (3731->13), ass. (0->60)
t130 = sin(pkin(11));
t132 = cos(pkin(11));
t137 = cos(qJ(2));
t133 = cos(pkin(6));
t135 = sin(qJ(2));
t141 = t133 * t135;
t126 = -t130 * t141 + t132 * t137;
t129 = qJ(3) + qJ(4);
t127 = sin(t129);
t128 = cos(t129);
t131 = sin(pkin(6));
t144 = t130 * t131;
t115 = t126 * t128 + t127 * t144;
t134 = sin(qJ(5));
t140 = t133 * t137;
t125 = t130 * t140 + t132 * t135;
t136 = cos(qJ(5));
t146 = t125 * t136;
t108 = t115 * t134 - t146;
t124 = t130 * t137 + t132 * t141;
t143 = t131 * t132;
t113 = t124 * t128 - t127 * t143;
t138 = -t130 * t135 + t132 * t140;
t105 = t113 * t134 + t136 * t138;
t142 = t131 * t135;
t122 = t127 * t133 + t128 * t142;
t118 = t131 * t136 * t137 + t122 * t134;
t102 = atan2(-t105, t118);
t100 = cos(t102);
t99 = sin(t102);
t97 = t100 * t118 - t105 * t99;
t96 = 0.1e1 / t97 ^ 2;
t150 = t108 * t96;
t109 = t115 * t136 + t125 * t134;
t104 = 0.1e1 / t109 ^ 2;
t114 = -t126 * t127 + t128 * t144;
t149 = t104 * t114;
t117 = 0.1e1 / t118 ^ 2;
t148 = t105 * t117;
t147 = t114 ^ 2 * t104;
t145 = t128 * t134;
t139 = t134 * t137;
t121 = -t127 * t142 + t128 * t133;
t120 = (t128 * t139 - t135 * t136) * t131;
t119 = t122 * t136 - t131 * t139;
t116 = 0.1e1 / t118;
t112 = -t124 * t127 - t128 * t143;
t110 = -t124 * t136 + t138 * t145;
t107 = t113 * t136 - t134 * t138;
t103 = 0.1e1 / t109;
t101 = 0.1e1 / (t105 ^ 2 * t117 + 0.1e1);
t98 = 0.1e1 / (0.1e1 + t147);
t95 = 0.1e1 / t97;
t94 = 0.1e1 / (t108 ^ 2 * t96 + 0.1e1);
t93 = (-t112 * t116 + t121 * t148) * t134 * t101;
t92 = (-t103 * t115 - t136 * t147) * t98;
t91 = (-t110 * t116 + t120 * t148) * t101;
t90 = (-t107 * t116 + t119 * t148) * t101;
t89 = (t114 * t134 * t95 - ((-t112 * t134 - t118 * t93) * t99 + (-t105 * t93 + t121 * t134) * t100) * t150) * t94;
t1 = [0, t91, t93, t93, t90, 0; 0 ((-t125 * t145 - t126 * t136) * t95 - ((-t118 * t91 - t110) * t99 + (-t105 * t91 + t120) * t100) * t150) * t94, t89, t89 (t109 * t95 - ((-t118 * t90 - t107) * t99 + (-t105 * t90 + t119) * t100) * t150) * t94, 0; 0 (t125 * t127 * t103 - (t126 * t134 - t128 * t146) * t149) * t98, t92, t92, t108 * t98 * t149, 0;];
Ja_rot  = t1;
