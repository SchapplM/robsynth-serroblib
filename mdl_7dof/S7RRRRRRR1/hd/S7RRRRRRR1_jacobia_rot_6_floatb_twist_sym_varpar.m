% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% Ja_rot [3x7]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_rot = S7RRRRRRR1_jacobia_rot_6_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobia_rot_6_floatb_twist_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobia_rot_6_floatb_twist_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-26 21:20:53
% EndTime: 2018-11-26 21:20:54
% DurationCPUTime: 0.51s
% Computational Cost: add. (1373->65), mult. (3817->158), div. (145->9), fcn. (5357->15), ass. (0->78)
t141 = cos(qJ(3));
t135 = sin(qJ(3));
t143 = cos(qJ(1));
t152 = t143 * t135;
t137 = sin(qJ(1));
t142 = cos(qJ(2));
t154 = t137 * t142;
t125 = t141 * t154 + t152;
t134 = sin(qJ(4));
t140 = cos(qJ(4));
t136 = sin(qJ(2));
t157 = t136 * t137;
t112 = t125 * t140 + t134 * t157;
t133 = sin(qJ(5));
t139 = cos(qJ(5));
t151 = t143 * t141;
t146 = -t135 * t154 + t151;
t97 = t112 * t133 - t146 * t139;
t144 = t146 * t133;
t99 = t112 * t139 + t144;
t127 = -t137 * t135 + t142 * t151;
t155 = t136 * t143;
t117 = t127 * t140 + t134 * t155;
t126 = -t137 * t141 - t142 * t152;
t102 = t117 * t139 + t126 * t133;
t116 = t127 * t134 - t140 * t155;
t132 = sin(qJ(6));
t138 = cos(qJ(6));
t92 = t102 * t138 + t116 * t132;
t90 = 0.1e1 / t92 ^ 2;
t91 = t102 * t132 - t116 * t138;
t167 = t90 * t91;
t156 = t136 * t141;
t124 = -t142 * t134 + t140 * t156;
t158 = t135 * t139;
t148 = t136 * t158;
t109 = t124 * t133 + t148;
t96 = atan2(-t97, t109);
t94 = cos(t96);
t166 = t94 * t97;
t160 = t126 * t139;
t101 = t117 * t133 - t160;
t93 = sin(t96);
t87 = t94 * t109 - t93 * t97;
t86 = 0.1e1 / t87 ^ 2;
t165 = t101 * t86;
t164 = t101 ^ 2 * t86;
t108 = 0.1e1 / t109 ^ 2;
t163 = t108 * t97;
t162 = t116 * t139;
t161 = t126 * t134;
t159 = t133 * t140;
t153 = t142 * t140;
t150 = t91 ^ 2 * t90 + 0.1e1;
t149 = t136 * t135 * t133;
t147 = -t125 * t134 + t140 * t157;
t145 = -t109 * t93 - t166;
t123 = -t134 * t156 - t153;
t120 = t124 * t143;
t119 = t123 * t143;
t118 = (-t135 * t159 + t139 * t141) * t136;
t115 = (t136 * t134 + t141 * t153) * t133 + t142 * t158;
t110 = t124 * t139 - t149;
t107 = 0.1e1 / t109;
t106 = -t120 * t139 + t143 * t149;
t105 = t109 * t137;
t104 = -t127 * t133 + t140 * t160;
t103 = t125 * t139 + t140 * t144;
t95 = 0.1e1 / (t97 ^ 2 * t108 + 0.1e1);
t89 = 0.1e1 / t92;
t88 = 0.1e1 / t150;
t85 = 0.1e1 / t87;
t84 = 0.1e1 / (0.1e1 + t164);
t83 = (-t107 * t147 + t123 * t163) * t95 * t133;
t82 = (-t103 * t107 + t118 * t163) * t95;
t81 = (t105 * t107 + t115 * t163) * t95;
t80 = (-t107 * t99 + t110 * t163) * t95;
t1 = [-t101 * t107 * t95, t81, t82, t83, t80, 0, 0; (-t97 * t85 - (-t93 + (t107 * t166 + t93) * t95) * t164) * t84 ((-t120 * t133 - t143 * t148) * t85 - (t93 * t105 + t94 * t115 + t145 * t81) * t165) * t84 ((t126 * t159 + t127 * t139) * t85 - (-t93 * t103 + t94 * t118 + t145 * t82) * t165) * t84 (-t116 * t133 * t85 - (t145 * t83 + (t123 * t94 - t147 * t93) * t133) * t165) * t84 (t102 * t85 - (t94 * t110 + t145 * t80 - t93 * t99) * t165) * t84, 0, 0; ((-t132 * t99 - t138 * t147) * t89 - (t132 * t147 - t138 * t99) * t167) * t88 ((t106 * t132 - t119 * t138) * t89 - (t106 * t138 + t119 * t132) * t167) * t88 ((t104 * t132 - t138 * t161) * t89 - (t104 * t138 + t132 * t161) * t167) * t88 ((-t117 * t138 - t132 * t162) * t89 - (t117 * t132 - t138 * t162) * t167) * t88 (-t132 * t89 + t138 * t167) * t88 * t101, t150 * t88, 0;];
Ja_rot  = t1;
