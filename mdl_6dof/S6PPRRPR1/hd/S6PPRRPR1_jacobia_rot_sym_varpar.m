% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PPRRPR1
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
%   Wie in S6PPRRPR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:10
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PPRRPR1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRPR1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_jacobia_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (31->0), mult. (90->0), div. (5->0), fcn. (119->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (333->29), mult. (995->67), div. (35->9), fcn. (1360->15), ass. (0->42)
	t64 = sin(pkin(12));
	t65 = sin(pkin(11));
	t68 = cos(pkin(12));
	t69 = cos(pkin(11));
	t70 = cos(pkin(7));
	t71 = cos(pkin(6));
	t81 = t69 * t71;
	t66 = sin(pkin(7));
	t67 = sin(pkin(6));
	t82 = t67 * t66;
	t85 = (-t65 * t64 + t68 * t81) * t70 - t69 * t82;
	t84 = t65 * t71;
	t83 = t66 * t71;
	t75 = cos(qJ(3));
	t80 = t70 * t75;
	t79 = t65 * t82;
	t61 = -t69 * t64 - t68 * t84;
	t62 = -t64 * t84 + t69 * t68;
	t73 = sin(qJ(3));
	t53 = t62 * t75 + (t61 * t70 + t79) * t73;
	t57 = t65 * t67 * t70 - t61 * t66;
	t72 = sin(qJ(4));
	t74 = cos(qJ(4));
	t44 = t53 * t74 + t57 * t72;
	t42 = 0.1e1 / t44 ^ 2;
	t43 = t53 * t72 - t57 * t74;
	t77 = t43 ^ 2 * t42 + 0.1e1;
	t60 = t64 * t81 + t65 * t68;
	t56 = t73 * t83 + (t68 * t70 * t73 + t64 * t75) * t67;
	t55 = -t75 * t83 + (t64 * t73 - t68 * t80) * t67;
	t54 = 0.1e1 / t55 ^ 2;
	t52 = -t61 * t80 + t62 * t73 - t75 * t79;
	t51 = t60 * t75 + t85 * t73;
	t49 = t60 * t73 - t85 * t75;
	t48 = atan2(-t49, t55);
	t46 = cos(t48);
	t45 = sin(t48);
	t41 = 0.1e1 / t77;
	t40 = -t45 * t49 + t46 * t55;
	t39 = 0.1e1 / t40 ^ 2;
	t37 = (-t51 / t55 + t56 * t49 * t54) / (t49 ^ 2 * t54 + 0.1e1);
	t1 = [0, 0, t37, 0, 0, 0; 0, 0, (t53 / t40 - (-t45 * t51 + t46 * t56 + (-t45 * t55 - t46 * t49) * t37) * t52 * t39) / (t52 ^ 2 * t39 + 0.1e1), 0, 0, 0; 0, 0, (-t72 / t44 + t74 * t43 * t42) * t52 * t41, t77 * t41, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.37s
	% Computational Cost: add. (961->43), mult. (2761->105), div. (60->9), fcn. (3790->17), ass. (0->64)
	t89 = sin(pkin(11));
	t96 = cos(pkin(6));
	t114 = t89 * t96;
	t88 = sin(pkin(12));
	t93 = cos(pkin(12));
	t94 = cos(pkin(11));
	t105 = t93 * t114 + t94 * t88;
	t90 = sin(pkin(7));
	t91 = sin(pkin(6));
	t112 = t91 * t90;
	t95 = cos(pkin(7));
	t119 = t105 * t95 - t89 * t112;
	t109 = t94 * t96;
	t106 = t93 * t109 - t89 * t88;
	t111 = t91 * t95;
	t103 = -t106 * t90 - t94 * t111;
	t100 = cos(qJ(3));
	t102 = t106 * t95 - t94 * t112;
	t84 = t88 * t109 + t89 * t93;
	t98 = sin(qJ(3));
	t71 = t84 * t100 + t102 * t98;
	t97 = sin(qJ(4));
	t99 = cos(qJ(4));
	t65 = -t103 * t99 + t71 * t97;
	t110 = t93 * t95;
	t113 = t90 * t96;
	t81 = t98 * t113 + (t100 * t88 + t98 * t110) * t91;
	t83 = -t93 * t112 + t96 * t95;
	t76 = t81 * t97 - t83 * t99;
	t64 = atan2(-t65, t76);
	t61 = sin(t64);
	t62 = cos(t64);
	t55 = -t61 * t65 + t62 * t76;
	t54 = 0.1e1 / t55 ^ 2;
	t101 = t105 * t90 + t89 * t111;
	t85 = -t88 * t114 + t94 * t93;
	t73 = t85 * t100 - t119 * t98;
	t68 = -t101 * t99 + t73 * t97;
	t118 = t54 * t68;
	t69 = t101 * t97 + t73 * t99;
	t72 = t119 * t100 + t85 * t98;
	t87 = sin(pkin(13));
	t92 = cos(pkin(13));
	t60 = t69 * t92 + t72 * t87;
	t58 = 0.1e1 / t60 ^ 2;
	t59 = t69 * t87 - t72 * t92;
	t117 = t58 * t59;
	t75 = 0.1e1 / t76 ^ 2;
	t116 = t65 * t75;
	t115 = t72 * t99;
	t107 = -t61 * t76 - t62 * t65;
	t80 = -t91 * t88 * t98 + (t91 * t110 + t113) * t100;
	t77 = t81 * t99 + t83 * t97;
	t74 = 0.1e1 / t76;
	t70 = t102 * t100 - t84 * t98;
	t67 = t103 * t97 + t71 * t99;
	t63 = 0.1e1 / (t65 ^ 2 * t75 + 0.1e1);
	t57 = 0.1e1 / t60;
	t56 = 0.1e1 / (t59 ^ 2 * t58 + 0.1e1);
	t53 = 0.1e1 / t55;
	t52 = 0.1e1 / (t68 ^ 2 * t54 + 0.1e1);
	t51 = (t80 * t116 - t70 * t74) * t97 * t63;
	t50 = (t77 * t116 - t67 * t74) * t63;
	t1 = [0, 0, t51, t50, 0, 0; 0, 0, (-t72 * t97 * t53 - ((-t61 * t70 + t62 * t80) * t97 + t107 * t51) * t118) * t52, (t69 * t53 - (t107 * t50 - t61 * t67 + t62 * t77) * t118) * t52, 0, 0; 0, 0, ((-t87 * t115 - t73 * t92) * t57 - (-t92 * t115 + t73 * t87) * t117) * t56, (t92 * t117 - t57 * t87) * t68 * t56, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (1096->44), mult. (3005->105), div. (65->9), fcn. (4119->17), ass. (0->67)
	t108 = sin(pkin(11));
	t113 = cos(pkin(7));
	t107 = sin(pkin(12));
	t111 = cos(pkin(12));
	t112 = cos(pkin(11));
	t114 = cos(pkin(6));
	t133 = t108 * t114;
	t123 = t112 * t107 + t111 * t133;
	t109 = sin(pkin(7));
	t110 = sin(pkin(6));
	t131 = t110 * t109;
	t139 = -t108 * t131 + t123 * t113;
	t115 = sin(qJ(4));
	t117 = cos(qJ(4));
	t128 = t112 * t114;
	t124 = -t108 * t107 + t111 * t128;
	t130 = t110 * t113;
	t121 = -t124 * t109 - t112 * t130;
	t101 = t107 * t128 + t108 * t111;
	t116 = sin(qJ(3));
	t118 = cos(qJ(3));
	t120 = -t112 * t131 + t124 * t113;
	t88 = t101 * t118 + t120 * t116;
	t82 = t88 * t115 - t121 * t117;
	t100 = -t111 * t131 + t114 * t113;
	t129 = t111 * t113;
	t132 = t109 * t114;
	t98 = t116 * t132 + (t107 * t118 + t116 * t129) * t110;
	t93 = -t100 * t117 + t98 * t115;
	t81 = atan2(-t82, t93);
	t78 = sin(t81);
	t79 = cos(t81);
	t72 = -t78 * t82 + t79 * t93;
	t71 = 0.1e1 / t72 ^ 2;
	t119 = t108 * t130 + t123 * t109;
	t102 = -t107 * t133 + t112 * t111;
	t90 = t102 * t118 - t139 * t116;
	t85 = t90 * t115 - t119 * t117;
	t138 = t71 * t85;
	t106 = pkin(13) + qJ(6);
	t105 = cos(t106);
	t104 = sin(t106);
	t89 = t102 * t116 + t139 * t118;
	t135 = t89 * t104;
	t86 = t119 * t115 + t90 * t117;
	t77 = t86 * t105 + t135;
	t75 = 0.1e1 / t77 ^ 2;
	t134 = t89 * t105;
	t76 = t86 * t104 - t134;
	t137 = t75 * t76;
	t92 = 0.1e1 / t93 ^ 2;
	t136 = t82 * t92;
	t127 = t76 ^ 2 * t75 + 0.1e1;
	t125 = -t78 * t93 - t79 * t82;
	t97 = t118 * t132 + (-t107 * t116 + t118 * t129) * t110;
	t94 = t100 * t115 + t98 * t117;
	t91 = 0.1e1 / t93;
	t87 = -t101 * t116 + t120 * t118;
	t84 = t121 * t115 + t88 * t117;
	t80 = 0.1e1 / (t82 ^ 2 * t92 + 0.1e1);
	t74 = 0.1e1 / t77;
	t73 = 0.1e1 / t127;
	t70 = 0.1e1 / t72;
	t69 = 0.1e1 / (t85 ^ 2 * t71 + 0.1e1);
	t68 = (t97 * t136 - t87 * t91) * t80 * t115;
	t67 = (t94 * t136 - t84 * t91) * t80;
	t1 = [0, 0, t68, t67, 0, 0; 0, 0, (-t89 * t115 * t70 - (t125 * t68 + (-t78 * t87 + t79 * t97) * t115) * t138) * t69, (t86 * t70 - (t125 * t67 - t78 * t84 + t79 * t94) * t138) * t69, 0, 0; 0, 0, ((-t90 * t105 - t117 * t135) * t74 - (t90 * t104 - t117 * t134) * t137) * t73, (-t104 * t74 + t105 * t137) * t85 * t73, 0, t127 * t73;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end