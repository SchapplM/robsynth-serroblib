% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPP1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPP1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:02:47
% EndTime: 2019-02-26 22:02:48
% DurationCPUTime: 0.29s
% Computational Cost: add. (665->43), mult. (2053->105), div. (80->9), fcn. (2819->13), ass. (0->60)
t118 = sin(pkin(6));
t120 = cos(pkin(6));
t125 = cos(qJ(2));
t121 = sin(qJ(3));
t122 = sin(qJ(2));
t140 = t121 * t122;
t127 = t118 * t125 + t120 * t140;
t117 = sin(pkin(10));
t119 = cos(pkin(10));
t124 = cos(qJ(3));
t138 = t122 * t124;
t151 = t127 * t117 - t119 * t138;
t123 = sin(qJ(1));
t126 = cos(qJ(1));
t134 = t126 * t121;
t114 = t123 * t124 - t125 * t134;
t137 = t122 * t126;
t128 = t114 * t120 + t118 * t137;
t133 = t126 * t124;
t136 = t123 * t121;
t115 = t125 * t133 + t136;
t144 = t115 * t119;
t101 = t128 * t117 + t144;
t135 = t124 * t125;
t113 = t123 * t135 - t134;
t112 = t125 * t136 + t133;
t139 = t122 * t123;
t129 = t112 * t120 - t118 * t139;
t98 = t113 * t119 - t129 * t117;
t96 = atan2(-t98, -t151);
t93 = sin(t96);
t94 = cos(t96);
t92 = -t151 * t94 - t93 * t98;
t91 = 0.1e1 / t92 ^ 2;
t147 = t101 ^ 2 * t91;
t89 = 0.1e1 / (0.1e1 + t147);
t90 = 0.1e1 / t92;
t150 = t89 * t90;
t149 = t94 * t98;
t148 = t101 * t91;
t105 = 0.1e1 / t151 ^ 2;
t146 = t105 * t98;
t100 = -t115 * t117 + t128 * t119;
t110 = -t114 * t118 + t120 * t137;
t109 = 0.1e1 / t110 ^ 2;
t145 = t100 * t109;
t143 = t117 * t120;
t141 = t120 * t125;
t130 = t151 * t93 - t149;
t111 = (-t119 * t121 - t124 * t143) * t122;
t108 = 0.1e1 / t110;
t107 = t119 * t135 + (t118 * t122 - t121 * t141) * t117;
t104 = 0.1e1 / t151;
t103 = t151 * t123;
t102 = -t112 * t119 - t113 * t143;
t97 = 0.1e1 / (t100 ^ 2 * t109 + 0.1e1);
t95 = 0.1e1 / (t98 ^ 2 * t105 + 0.1e1);
t88 = (t102 * t104 + t111 * t146) * t95;
t87 = (t103 * t104 + t107 * t146) * t95;
t1 = [t101 * t104 * t95, t87, t88, 0, 0, 0; -t98 * t150 - (-t93 + (-t104 * t149 + t93) * t95) * t89 * t147 -(-t93 * t103 + t94 * t107 + t130 * t87) * t89 * t148 + t151 * t126 * t150 ((t114 * t119 - t115 * t143) * t90 - (-t93 * t102 + t94 * t111 + t130 * t88) * t148) * t89, 0, 0, 0; ((t113 * t117 + t129 * t119) * t108 - (-t112 * t118 - t120 * t139) * t145) * t97 ((t117 * t138 + t127 * t119) * t108 - (-t118 * t140 + t141) * t145) * t97 * t126 ((-t114 * t117 - t120 * t144) * t108 - t115 * t118 * t145) * t97, 0, 0, 0;];
Ja_rot  = t1;
