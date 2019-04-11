% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR10V2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRR10V2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_jacobia_rot_6_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:56:32
% EndTime: 2019-04-11 14:56:32
% DurationCPUTime: 0.28s
% Computational Cost: add. (1532->48), mult. (2174->116), div. (145->9), fcn. (3100->13), ass. (0->60)
t128 = qJ(2) + qJ(3);
t127 = cos(t128);
t135 = cos(qJ(4));
t136 = cos(qJ(1));
t139 = t136 * t135;
t131 = sin(qJ(4));
t132 = sin(qJ(1));
t142 = t132 * t131;
t121 = t127 * t139 + t142;
t130 = sin(qJ(5));
t134 = cos(qJ(5));
t126 = sin(t128);
t144 = t126 * t136;
t110 = t121 * t134 + t130 * t144;
t140 = t136 * t131;
t141 = t132 * t135;
t120 = t127 * t140 - t141;
t129 = sin(qJ(6));
t133 = cos(qJ(6));
t100 = t110 * t133 + t120 * t129;
t98 = 0.1e1 / t100 ^ 2;
t99 = t110 * t129 - t120 * t133;
t151 = t98 * t99;
t109 = t121 * t130 - t134 * t144;
t119 = t127 * t141 - t140;
t145 = t126 * t134;
t105 = t119 * t130 - t132 * t145;
t143 = t130 * t135;
t115 = t126 * t143 + t127 * t134;
t104 = atan2(-t105, t115);
t101 = sin(t104);
t102 = cos(t104);
t95 = -t101 * t105 + t102 * t115;
t94 = 0.1e1 / t95 ^ 2;
t150 = t109 * t94;
t149 = t109 ^ 2 * t94;
t114 = 0.1e1 / t115 ^ 2;
t148 = t105 * t114;
t147 = t120 * t134;
t146 = t126 * t131;
t138 = t99 ^ 2 * t98 + 0.1e1;
t137 = t126 * t140;
t107 = t132 * t126 * t130 + t119 * t134;
t116 = -t127 * t130 + t135 * t145;
t118 = -t127 * t142 - t139;
t117 = t127 * t143 - t145;
t113 = 0.1e1 / t115;
t112 = t116 * t136;
t111 = t115 * t132;
t103 = 0.1e1 / (t105 ^ 2 * t114 + 0.1e1);
t97 = 0.1e1 / t100;
t96 = 0.1e1 / t138;
t93 = 0.1e1 / t95;
t92 = 0.1e1 / (0.1e1 + t149);
t91 = (-t113 * t118 - t146 * t148) * t130 * t103;
t90 = (t111 * t113 + t117 * t148) * t103;
t89 = (-t107 * t113 + t116 * t148) * t103;
t88 = ((-t112 * t129 + t133 * t137) * t97 - (-t112 * t133 - t129 * t137) * t151) * t96;
t87 = (-((-t105 * t90 + t117) * t102 + (-t115 * t90 + t111) * t101) * t150 - t115 * t93 * t136) * t92;
t1 = [-t109 * t113 * t103, t90, t90, t91, t89, 0; (-t105 * t93 - (-t101 + (t102 * t105 * t113 + t101) * t103) * t149) * t92, t87, t87 (-t120 * t130 * t93 - ((-t105 * t91 - t130 * t146) * t102 + (-t115 * t91 - t118 * t130) * t101) * t150) * t92 (t110 * t93 - ((-t105 * t89 + t116) * t102 + (-t115 * t89 - t107) * t101) * t150) * t92, 0; ((-t107 * t129 - t118 * t133) * t97 - (-t107 * t133 + t118 * t129) * t151) * t96, t88, t88 ((-t121 * t133 - t129 * t147) * t97 - (t121 * t129 - t133 * t147) * t151) * t96 (-t129 * t97 + t133 * t151) * t96 * t109, t138 * t96;];
Ja_rot  = t1;
