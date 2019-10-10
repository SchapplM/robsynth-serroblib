% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
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
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6RRRRRR10V2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:38
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRR10V2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR10V2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_jacobia_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (358->21), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
	t67 = qJ(2) + qJ(3);
	t66 = cos(t67);
	t65 = sin(t67);
	t69 = sin(qJ(1));
	t77 = t69 * t65;
	t60 = atan2(-t77, -t66);
	t58 = sin(t60);
	t59 = cos(t60);
	t51 = -t58 * t77 - t59 * t66;
	t50 = 0.1e1 / t51 ^ 2;
	t71 = cos(qJ(1));
	t83 = t50 * t71 ^ 2;
	t70 = cos(qJ(4));
	t73 = t71 * t70;
	t68 = sin(qJ(4));
	t76 = t69 * t68;
	t57 = t66 * t73 + t76;
	t55 = 0.1e1 / t57 ^ 2;
	t74 = t71 * t68;
	t75 = t69 * t70;
	t56 = t66 * t74 - t75;
	t82 = t55 * t56;
	t81 = t58 * t66;
	t62 = t65 ^ 2;
	t80 = t62 / t66 ^ 2;
	t79 = t65 * t71;
	t61 = 0.1e1 / (t69 ^ 2 * t80 + 0.1e1);
	t78 = t69 * t61;
	t72 = t56 ^ 2 * t55 + 0.1e1;
	t63 = 0.1e1 / t66;
	t54 = 0.1e1 / t57;
	t53 = 0.1e1 / t72;
	t52 = (0.1e1 + t80) * t78;
	t49 = 0.1e1 / t51;
	t48 = 0.1e1 / (t62 * t83 + 0.1e1);
	t47 = (-t54 * t68 + t70 * t82) * t53 * t79;
	t46 = (t66 * t49 - (-t69 * t81 + t59 * t65 + (-t59 * t77 + t81) * t52) * t65 * t50) * t71 * t48;
	t1 = [t63 * t61 * t79, t52, t52, 0, 0, 0; (-t49 * t77 - (-t59 * t62 * t63 * t78 + (t61 - 0.1e1) * t65 * t58) * t65 * t83) * t48, t46, t46, 0, 0, 0; ((-t66 * t76 - t73) * t54 - (-t66 * t75 + t74) * t82) * t53, t47, t47, t72 * t53, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (600->33), mult. (873->84), div. (144->11), fcn. (1297->11), ass. (0->49)
	t85 = qJ(2) + qJ(3);
	t81 = sin(t85);
	t87 = sin(qJ(4));
	t100 = t81 * t87;
	t82 = cos(t85);
	t90 = cos(qJ(4));
	t91 = cos(qJ(1));
	t94 = t90 * t91;
	t88 = sin(qJ(1));
	t96 = t88 * t87;
	t72 = t82 * t96 + t94;
	t71 = atan2(-t72, t100);
	t68 = sin(t71);
	t69 = cos(t71);
	t62 = t69 * t100 - t68 * t72;
	t61 = 0.1e1 / t62 ^ 2;
	t93 = t91 * t87;
	t95 = t88 * t90;
	t75 = t82 * t93 - t95;
	t105 = t61 * t75;
	t104 = t61 * t75 ^ 2;
	t76 = t82 * t94 + t96;
	t86 = sin(qJ(5));
	t89 = cos(qJ(5));
	t97 = t81 * t91;
	t67 = t76 * t89 + t86 * t97;
	t65 = 0.1e1 / t67 ^ 2;
	t66 = t76 * t86 - t89 * t97;
	t103 = t65 * t66;
	t102 = t69 * t72;
	t79 = 0.1e1 / t81;
	t83 = 0.1e1 / t87;
	t101 = t79 * t83;
	t99 = t81 * t88;
	t98 = t81 * t90;
	t92 = t65 * t66 ^ 2 + 0.1e1;
	t84 = 0.1e1 / t87 ^ 2;
	t80 = 0.1e1 / t81 ^ 2;
	t74 = t82 * t95 - t93;
	t70 = 0.1e1 / (t72 ^ 2 * t80 * t84 + 0.1e1);
	t64 = 0.1e1 / t67;
	t63 = 0.1e1 / t92;
	t60 = 0.1e1 / t62;
	t59 = (t72 * t80 * t82 * t83 + t88) * t70;
	t58 = 0.1e1 / (0.1e1 + t104);
	t57 = (t72 * t84 * t90 - t74 * t83) * t79 * t70;
	t56 = ((-t82 * t89 - t86 * t98) * t64 - (t82 * t86 - t89 * t98) * t103) * t63 * t91;
	t55 = (t59 * t102 * t105 + (-t60 * t97 - (t69 * t82 + (-t59 * t81 + t99) * t68) * t105) * t87) * t58;
	t1 = [-t75 * t70 * t101, t59, t59, t57, 0, 0; (-t72 * t60 - (-t68 + (t101 * t102 + t68) * t70) * t104) * t58, t55, t55, (t76 * t60 - (t69 * t98 - t68 * t74 + (-t68 * t100 - t102) * t57) * t105) * t58, 0, 0; ((-t74 * t86 + t89 * t99) * t64 - (-t74 * t89 - t86 * t99) * t103) * t63, t56, t56, (t89 * t103 - t64 * t86) * t75 * t63, t92 * t63, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.40s
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
	t138 = t98 * t99 ^ 2 + 0.1e1;
	t137 = t126 * t140;
	t107 = t126 * t130 * t132 + t119 * t134;
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
	t1 = [-t109 * t113 * t103, t90, t90, t91, t89, 0; (-t105 * t93 - (-t101 + (t102 * t105 * t113 + t101) * t103) * t149) * t92, t87, t87, (-t120 * t130 * t93 - ((-t105 * t91 - t130 * t146) * t102 + (-t115 * t91 - t118 * t130) * t101) * t150) * t92, (t110 * t93 - ((-t105 * t89 + t116) * t102 + (-t115 * t89 - t107) * t101) * t150) * t92, 0; ((-t107 * t129 - t118 * t133) * t97 - (-t107 * t133 + t118 * t129) * t151) * t96, t88, t88, ((-t121 * t133 - t129 * t147) * t97 - (t121 * t129 - t133 * t147) * t151) * t96, (-t129 * t97 + t133 * t151) * t96 * t109, t138 * t96;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end