% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRP5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:17:28
% EndTime: 2019-02-26 20:17:29
% DurationCPUTime: 0.34s
% Computational Cost: add. (1524->59), mult. (4378->143), div. (95->9), fcn. (6002->17), ass. (0->75)
t129 = sin(pkin(12));
t132 = cos(pkin(12));
t138 = sin(qJ(2));
t134 = cos(pkin(6));
t142 = cos(qJ(2));
t152 = t134 * t142;
t125 = -t129 * t152 - t132 * t138;
t153 = t134 * t138;
t126 = -t129 * t153 + t132 * t142;
t133 = cos(pkin(7));
t137 = sin(qJ(3));
t141 = cos(qJ(3));
t130 = sin(pkin(7));
t131 = sin(pkin(6));
t157 = t131 * t130;
t148 = t129 * t157;
t112 = t126 * t141 + (t125 * t133 + t148) * t137;
t136 = sin(qJ(4));
t140 = cos(qJ(4));
t156 = t131 * t133;
t146 = -t125 * t130 + t129 * t156;
t104 = t112 * t140 + t146 * t136;
t154 = t133 * t141;
t111 = -t125 * t154 + t126 * t137 - t141 * t148;
t135 = sin(qJ(5));
t139 = cos(qJ(5));
t95 = t104 * t139 + t111 * t135;
t93 = 0.1e1 / t95 ^ 2;
t94 = t104 * t135 - t111 * t139;
t163 = t93 * t94;
t103 = t112 * t136 - t146 * t140;
t124 = t129 * t142 + t132 * t153;
t145 = -t129 * t138 + t132 * t152;
t143 = -t132 * t157 + t145 * t133;
t110 = t124 * t141 + t143 * t137;
t144 = -t145 * t130 - t132 * t156;
t100 = t110 * t136 - t144 * t140;
t155 = t133 * t137;
t159 = t130 * t134;
t121 = t137 * t159 + (t138 * t141 + t142 * t155) * t131;
t123 = t134 * t133 - t142 * t157;
t113 = t121 * t136 - t123 * t140;
t99 = atan2(-t100, t113);
t96 = sin(t99);
t97 = cos(t99);
t90 = -t100 * t96 + t113 * t97;
t89 = 0.1e1 / t90 ^ 2;
t162 = t103 * t89;
t108 = 0.1e1 / t113 ^ 2;
t161 = t100 * t108;
t160 = t111 * t140;
t158 = t130 * t140;
t151 = t137 * t138;
t150 = t141 * t142;
t149 = t93 * t94 ^ 2 + 0.1e1;
t147 = -t100 * t97 - t113 * t96;
t120 = t141 * t159 + (t133 * t150 - t151) * t131;
t117 = ((-t133 * t151 + t150) * t136 - t138 * t158) * t131;
t116 = t125 * t141 - t126 * t155;
t115 = t125 * t137 + t126 * t154;
t114 = t121 * t140 + t123 * t136;
t109 = -t124 * t137 + t143 * t141;
t107 = 0.1e1 / t113;
t106 = t126 * t130 * t136 + t116 * t140;
t105 = (-t124 * t155 + t145 * t141) * t136 - t124 * t158;
t102 = t110 * t140 + t144 * t136;
t98 = 0.1e1 / (t100 ^ 2 * t108 + 0.1e1);
t92 = 0.1e1 / t95;
t91 = 0.1e1 / t149;
t88 = 0.1e1 / t90;
t87 = 0.1e1 / (t103 ^ 2 * t89 + 0.1e1);
t86 = (-t107 * t109 + t120 * t161) * t98 * t136;
t85 = (-t105 * t107 + t117 * t161) * t98;
t84 = (-t102 * t107 + t114 * t161) * t98;
t1 = [0, t85, t86, t84, 0, 0; 0 ((t116 * t136 - t126 * t158) * t88 - (-t105 * t96 + t117 * t97 + t147 * t85) * t162) * t87 (-t111 * t136 * t88 - (t147 * t86 + (-t109 * t96 + t120 * t97) * t136) * t162) * t87 (t104 * t88 - (-t102 * t96 + t114 * t97 + t147 * t84) * t162) * t87, 0, 0; 0 ((t106 * t135 - t115 * t139) * t92 - (t106 * t139 + t115 * t135) * t163) * t91 ((-t112 * t139 - t135 * t160) * t92 - (t112 * t135 - t139 * t160) * t163) * t91 (-t135 * t92 + t139 * t163) * t91 * t103, t149 * t91, 0;];
Ja_rot  = t1;
