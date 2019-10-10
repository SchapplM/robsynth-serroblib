% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PPPRRR1
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
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6PPPRRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 08:49
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PPPRRR1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPPRRR1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_jacobia_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (87->0), mult. (251->0), div. (5->0), fcn. (336->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:22
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (968->40), mult. (2836->90), div. (35->9), fcn. (3836->19), ass. (0->56)
	t95 = cos(pkin(12));
	t98 = cos(pkin(6));
	t111 = t95 * t98;
	t88 = sin(pkin(13));
	t89 = sin(pkin(12));
	t94 = cos(pkin(13));
	t107 = t94 * t111 - t89 * t88;
	t92 = sin(pkin(6));
	t112 = t95 * t92;
	t91 = sin(pkin(7));
	t97 = cos(pkin(7));
	t105 = t107 * t97 - t91 * t112;
	t106 = t88 * t111 + t89 * t94;
	t87 = sin(pkin(14));
	t90 = sin(pkin(8));
	t93 = cos(pkin(14));
	t96 = cos(pkin(8));
	t120 = (t105 * t93 - t106 * t87) * t96 + (-t107 * t91 - t97 * t112) * t90;
	t114 = t92 * t91;
	t116 = t89 * t98;
	t85 = -t94 * t116 - t95 * t88;
	t108 = t89 * t114 + t85 * t97;
	t86 = -t88 * t116 + t95 * t94;
	t77 = t108 * t93 - t86 * t87;
	t83 = t89 * t92 * t97 - t85 * t91;
	t119 = t77 * t96 + t83 * t90;
	t115 = t91 * t98;
	t113 = t94 * t97;
	t101 = cos(qJ(5));
	t100 = sin(qJ(4));
	t102 = cos(qJ(4));
	t78 = t108 * t87 + t86 * t93;
	t69 = t119 * t100 + t78 * t102;
	t73 = -t77 * t90 + t83 * t96;
	t99 = sin(qJ(5));
	t60 = t69 * t101 + t73 * t99;
	t58 = 0.1e1 / t60 ^ 2;
	t59 = -t73 * t101 + t69 * t99;
	t110 = t59 ^ 2 * t58 + 0.1e1;
	t109 = (t93 * t115 + (t93 * t113 - t87 * t88) * t92) * t96 + (-t94 * t114 + t98 * t97) * t90;
	t82 = t92 * t88 * t93 + (t92 * t113 + t115) * t87;
	t76 = t105 * t87 + t106 * t93;
	t72 = t109 * t100 + t82 * t102;
	t71 = t82 * t100 - t109 * t102;
	t70 = 0.1e1 / t71 ^ 2;
	t68 = t78 * t100 - t119 * t102;
	t67 = t120 * t100 + t76 * t102;
	t65 = t76 * t100 - t120 * t102;
	t64 = atan2(-t65, t71);
	t62 = cos(t64);
	t61 = sin(t64);
	t57 = 0.1e1 / t110;
	t56 = -t61 * t65 + t62 * t71;
	t55 = 0.1e1 / t56 ^ 2;
	t53 = (-t67 / t71 + t72 * t65 * t70) / (t65 ^ 2 * t70 + 0.1e1);
	t1 = [0, 0, 0, t53, 0, 0; 0, 0, 0, (t69 / t56 - (-t61 * t67 + t62 * t72 + (-t61 * t71 - t62 * t65) * t53) * t68 * t55) / (t68 ^ 2 * t55 + 0.1e1), 0, 0; 0, 0, 0, (-t99 / t60 + t101 * t59 * t58) * t68 * t57, t110 * t57, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:22
	% DurationCPUTime: 0.73s
	% Computational Cost: add. (2922->55), mult. (8425->130), div. (65->9), fcn. (11413->21), ass. (0->81)
	t128 = sin(pkin(8));
	t134 = cos(pkin(8));
	t126 = sin(pkin(13));
	t132 = cos(pkin(13));
	t133 = cos(pkin(12));
	t127 = sin(pkin(12));
	t136 = cos(pkin(6));
	t162 = t127 * t136;
	t124 = -t126 * t162 + t133 * t132;
	t125 = sin(pkin(14));
	t131 = cos(pkin(14));
	t123 = -t133 * t126 - t132 * t162;
	t135 = cos(pkin(7));
	t129 = sin(pkin(7));
	t130 = sin(pkin(6));
	t161 = t129 * t130;
	t150 = t123 * t135 + t127 * t161;
	t147 = t124 * t125 - t150 * t131;
	t159 = t130 * t135;
	t151 = -t123 * t129 + t127 * t159;
	t167 = -t151 * t128 + t147 * t134;
	t158 = t132 * t135;
	t160 = t129 * t136;
	t119 = t130 * t126 * t131 + (t130 * t158 + t160) * t125;
	t139 = sin(qJ(4));
	t142 = cos(qJ(4));
	t118 = t131 * t160 + (-t125 * t126 + t131 * t158) * t130;
	t120 = -t132 * t161 + t136 * t135;
	t154 = t118 * t134 + t120 * t128;
	t112 = t119 * t142 + t154 * t139;
	t116 = -t118 * t128 + t120 * t134;
	t138 = sin(qJ(5));
	t141 = cos(qJ(5));
	t103 = t112 * t138 - t116 * t141;
	t157 = t133 * t136;
	t122 = t126 * t157 + t127 * t132;
	t121 = -t127 * t126 + t132 * t157;
	t152 = t121 * t135 - t133 * t161;
	t114 = t122 * t131 + t152 * t125;
	t148 = -t122 * t125 + t152 * t131;
	t153 = -t121 * t129 - t133 * t159;
	t144 = t153 * t128 + t148 * t134;
	t106 = t114 * t142 + t144 * t139;
	t145 = -t148 * t128 + t153 * t134;
	t96 = t106 * t138 - t145 * t141;
	t95 = atan2(-t96, t103);
	t92 = sin(t95);
	t93 = cos(t95);
	t86 = t93 * t103 - t92 * t96;
	t85 = 0.1e1 / t86 ^ 2;
	t115 = t124 * t131 + t150 * t125;
	t108 = t115 * t142 - t167 * t139;
	t143 = t147 * t128 + t151 * t134;
	t99 = t108 * t138 - t143 * t141;
	t166 = t85 * t99;
	t100 = t108 * t141 + t143 * t138;
	t107 = t115 * t139 + t167 * t142;
	t137 = sin(qJ(6));
	t140 = cos(qJ(6));
	t91 = t100 * t140 + t107 * t137;
	t89 = 0.1e1 / t91 ^ 2;
	t90 = t100 * t137 - t107 * t140;
	t165 = t89 * t90;
	t102 = 0.1e1 / t103 ^ 2;
	t164 = t102 * t96;
	t163 = t107 * t141;
	t156 = t90 ^ 2 * t89 + 0.1e1;
	t155 = -t103 * t92 - t93 * t96;
	t111 = -t119 * t139 + t154 * t142;
	t105 = -t114 * t139 + t144 * t142;
	t104 = t112 * t141 + t116 * t138;
	t101 = 0.1e1 / t103;
	t98 = t106 * t141 + t145 * t138;
	t94 = 0.1e1 / (t96 ^ 2 * t102 + 0.1e1);
	t88 = 0.1e1 / t91;
	t87 = 0.1e1 / t156;
	t84 = 0.1e1 / t86;
	t83 = 0.1e1 / (t99 ^ 2 * t85 + 0.1e1);
	t82 = (-t101 * t105 + t111 * t164) * t94 * t138;
	t81 = (-t101 * t98 + t104 * t164) * t94;
	t1 = [0, 0, 0, t82, t81, 0; 0, 0, 0, (-t107 * t138 * t84 - (t155 * t82 + (-t105 * t92 + t111 * t93) * t138) * t166) * t83, (t100 * t84 - (t93 * t104 + t155 * t81 - t92 * t98) * t166) * t83, 0; 0, 0, 0, ((-t108 * t140 - t137 * t163) * t88 - (t108 * t137 - t140 * t163) * t165) * t87, (-t137 * t88 + t140 * t165) * t99 * t87, t156 * t87;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end