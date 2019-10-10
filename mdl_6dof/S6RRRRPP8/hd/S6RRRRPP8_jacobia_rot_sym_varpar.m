% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRPP8
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
%   Wie in S6RRRRPP8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:31
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPP8_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP8_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:31:15
	% EndTime: 2019-10-10 12:31:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:31:15
	% EndTime: 2019-10-10 12:31:15
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:31:15
	% EndTime: 2019-10-10 12:31:16
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (35->13), mult. (116->34), div. (25->9), fcn. (175->9), ass. (0->25)
	t30 = cos(pkin(6));
	t29 = sin(pkin(6));
	t34 = cos(qJ(1));
	t38 = t34 * t29;
	t26 = atan2(t38, t30);
	t23 = sin(t26);
	t24 = cos(t26);
	t18 = t23 * t38 + t24 * t30;
	t32 = sin(qJ(1));
	t42 = 0.1e1 / t18 ^ 2 * t32 ^ 2;
	t27 = t29 ^ 2;
	t25 = 0.1e1 / (0.1e1 + t34 ^ 2 * t27 / t30 ^ 2);
	t41 = t25 / t30;
	t31 = sin(qJ(2));
	t40 = t32 * t31;
	t33 = cos(qJ(2));
	t39 = t32 * t33;
	t37 = t34 * t31;
	t36 = t34 * t33;
	t22 = -t30 * t40 + t36;
	t20 = 0.1e1 / t22 ^ 2;
	t21 = t30 * t39 + t37;
	t35 = t21 ^ 2 * t20 + 0.1e1;
	t19 = 0.1e1 / t35;
	t1 = [-t32 * t29 * t41, 0, 0, 0, 0, 0; (0.1e1 / t18 * t38 - (-t24 * t27 * t34 * t41 + (t25 - 0.1e1) * t29 * t23) * t29 * t42) / (t27 * t42 + 0.1e1), 0, 0, 0, 0, 0; ((t30 * t36 - t40) / t22 - (-t30 * t37 - t39) * t21 * t20) * t19, t35 * t19, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:31:15
	% EndTime: 2019-10-10 12:31:16
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (173->24), mult. (451->64), div. (72->11), fcn. (673->11), ass. (0->40)
	t53 = sin(qJ(2));
	t56 = cos(qJ(2));
	t57 = cos(qJ(1));
	t54 = sin(qJ(1));
	t61 = cos(pkin(6));
	t59 = t54 * t61;
	t46 = -t53 * t59 + t57 * t56;
	t52 = sin(qJ(3));
	t55 = cos(qJ(3));
	t51 = sin(pkin(6));
	t64 = t51 * t54;
	t37 = t46 * t55 + t52 * t64;
	t35 = 0.1e1 / t37 ^ 2;
	t36 = t46 * t52 - t55 * t64;
	t68 = t35 * t36;
	t58 = t57 * t61;
	t42 = t54 * t53 - t56 * t58;
	t63 = t51 * t56;
	t40 = atan2(-t42, -t63);
	t39 = cos(t40);
	t67 = t39 * t42;
	t38 = sin(t40);
	t32 = -t38 * t42 - t39 * t63;
	t31 = 0.1e1 / t32 ^ 2;
	t45 = t57 * t53 + t56 * t59;
	t66 = t45 ^ 2 * t31;
	t48 = 0.1e1 / t51;
	t49 = 0.1e1 / t56;
	t65 = t48 * t49;
	t62 = t51 * t57;
	t60 = t36 ^ 2 * t35 + 0.1e1;
	t50 = 0.1e1 / t56 ^ 2;
	t44 = t53 * t58 + t54 * t56;
	t41 = 0.1e1 / (0.1e1 + t42 ^ 2 / t51 ^ 2 * t50);
	t34 = 0.1e1 / t37;
	t33 = 0.1e1 / t60;
	t30 = 0.1e1 / t32;
	t29 = 0.1e1 / (0.1e1 + t66);
	t28 = (t42 * t50 * t53 + t44 * t49) * t48 * t41;
	t1 = [t45 * t41 * t65, t28, 0, 0, 0, 0; (-t42 * t30 - (-t38 + (-t65 * t67 + t38) * t41) * t66) * t29, (t46 * t30 - (t39 * t51 * t53 - t38 * t44 + (t38 * t63 - t67) * t28) * t45 * t31) * t29, 0, 0, 0, 0; ((-t44 * t52 - t55 * t62) * t34 - (-t44 * t55 + t52 * t62) * t68) * t33, (-t52 * t34 + t55 * t68) * t45 * t33, t60 * t33, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:31:16
	% EndTime: 2019-10-10 12:31:16
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (459->37), mult. (1300->91), div. (85->9), fcn. (1858->13), ass. (0->54)
	t76 = cos(pkin(6));
	t79 = sin(qJ(2));
	t84 = cos(qJ(1));
	t88 = t84 * t79;
	t80 = sin(qJ(1));
	t83 = cos(qJ(2));
	t89 = t80 * t83;
	t70 = t76 * t88 + t89;
	t78 = sin(qJ(3));
	t82 = cos(qJ(3));
	t75 = sin(pkin(6));
	t91 = t75 * t84;
	t59 = t70 * t78 + t82 * t91;
	t94 = t75 * t78;
	t67 = -t76 * t82 + t79 * t94;
	t58 = atan2(-t59, t67);
	t55 = sin(t58);
	t56 = cos(t58);
	t49 = -t55 * t59 + t56 * t67;
	t48 = 0.1e1 / t49 ^ 2;
	t87 = t84 * t83;
	t90 = t80 * t79;
	t72 = -t76 * t90 + t87;
	t93 = t75 * t82;
	t63 = t72 * t78 - t80 * t93;
	t100 = t48 * t63;
	t64 = t72 * t82 + t80 * t94;
	t71 = t76 * t89 + t88;
	t77 = sin(qJ(4));
	t81 = cos(qJ(4));
	t54 = t64 * t81 + t71 * t77;
	t52 = 0.1e1 / t54 ^ 2;
	t53 = t64 * t77 - t71 * t81;
	t99 = t52 * t53;
	t98 = t56 * t59;
	t66 = 0.1e1 / t67 ^ 2;
	t97 = t59 * t66;
	t96 = t63 ^ 2 * t48;
	t95 = t71 * t82;
	t92 = t75 * t83;
	t86 = t53 ^ 2 * t52 + 0.1e1;
	t61 = t70 * t82 - t78 * t91;
	t85 = -t55 * t67 - t98;
	t69 = t76 * t87 - t90;
	t68 = t76 * t78 + t79 * t93;
	t65 = 0.1e1 / t67;
	t57 = 0.1e1 / (t59 ^ 2 * t66 + 0.1e1);
	t51 = 0.1e1 / t54;
	t50 = 0.1e1 / t86;
	t47 = 0.1e1 / t49;
	t46 = 0.1e1 / (0.1e1 + t96);
	t45 = (-t65 * t69 + t92 * t97) * t78 * t57;
	t44 = (-t61 * t65 + t68 * t97) * t57;
	t1 = [-t63 * t65 * t57, t45, t44, 0, 0, 0; (-t59 * t47 - (-t55 + (t65 * t98 + t55) * t57) * t96) * t46, (-t71 * t78 * t47 - ((-t55 * t69 + t56 * t92) * t78 + t85 * t45) * t100) * t46, (t64 * t47 - (t85 * t44 - t55 * t61 + t56 * t68) * t100) * t46, 0, 0, 0; ((-t61 * t77 - t69 * t81) * t51 - (-t61 * t81 + t69 * t77) * t99) * t50, ((-t72 * t81 - t77 * t95) * t51 - (t72 * t77 - t81 * t95) * t99) * t50, (-t77 * t51 + t81 * t99) * t63 * t50, t86 * t50, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:31:16
	% EndTime: 2019-10-10 12:31:16
	% DurationCPUTime: 0.45s
	% Computational Cost: add. (905->44), mult. (2519->109), div. (107->9), fcn. (3557->13), ass. (0->63)
	t103 = sin(qJ(4));
	t107 = cos(qJ(4));
	t102 = cos(pkin(6));
	t109 = cos(qJ(2));
	t110 = cos(qJ(1));
	t115 = t110 * t109;
	t105 = sin(qJ(2));
	t106 = sin(qJ(1));
	t118 = t106 * t105;
	t114 = -t102 * t115 + t118;
	t104 = sin(qJ(3));
	t108 = cos(qJ(3));
	t101 = sin(pkin(6));
	t120 = t101 * t110;
	t116 = t110 * t105;
	t117 = t106 * t109;
	t97 = t102 * t116 + t117;
	t89 = t104 * t120 - t97 * t108;
	t133 = t89 * t103 + t114 * t107;
	t112 = t114 * t103;
	t132 = t89 * t107 - t112;
	t121 = t101 * t108;
	t96 = t102 * t104 + t105 * t121;
	t85 = t101 * t109 * t107 + t96 * t103;
	t73 = atan2(t133, t85);
	t69 = sin(t73);
	t70 = cos(t73);
	t68 = t133 * t69 + t70 * t85;
	t67 = 0.1e1 / t68 ^ 2;
	t98 = t102 * t117 + t116;
	t123 = t98 * t107;
	t122 = t101 * t104;
	t99 = -t102 * t118 + t115;
	t91 = t106 * t122 + t99 * t108;
	t79 = t91 * t103 - t123;
	t130 = t67 * t79;
	t129 = t70 * t133;
	t124 = t98 * t103;
	t80 = t91 * t107 + t124;
	t75 = 0.1e1 / t80 ^ 2;
	t90 = -t99 * t104 + t106 * t121;
	t128 = t75 * t90;
	t83 = 0.1e1 / t85 ^ 2;
	t127 = t133 * t83;
	t126 = t79 ^ 2 * t67;
	t125 = t90 ^ 2 * t75;
	t119 = t103 * t109;
	t113 = -t69 * t85 + t129;
	t111 = t97 * t104 + t108 * t120;
	t95 = t102 * t108 - t105 * t122;
	t92 = (-t105 * t107 + t108 * t119) * t101;
	t86 = -t101 * t119 + t96 * t107;
	t82 = 0.1e1 / t85;
	t81 = -t97 * t107 - t108 * t112;
	t74 = 0.1e1 / t80;
	t72 = 0.1e1 / (0.1e1 + t125);
	t71 = 0.1e1 / (t133 ^ 2 * t83 + 0.1e1);
	t66 = 0.1e1 / t68;
	t65 = 0.1e1 / (0.1e1 + t126);
	t64 = (t111 * t82 - t95 * t127) * t71 * t103;
	t63 = (-t92 * t127 - t81 * t82) * t71;
	t62 = (-t86 * t127 + t132 * t82) * t71;
	t1 = [-t79 * t82 * t71, t63, t64, t62, 0, 0; (t133 * t66 - (-t69 + (-t82 * t129 + t69) * t71) * t126) * t65, ((-t99 * t107 - t108 * t124) * t66 - (t113 * t63 - t69 * t81 + t70 * t92) * t130) * t65, (t90 * t103 * t66 - (t113 * t64 + (t111 * t69 + t70 * t95) * t103) * t130) * t65, (t80 * t66 - (t113 * t62 + t132 * t69 + t70 * t86) * t130) * t65, 0, 0; (t111 * t74 - t132 * t128) * t72, (t98 * t104 * t74 - (t99 * t103 - t108 * t123) * t128) * t72, (-t107 * t125 - t74 * t91) * t72, t79 * t72 * t128, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:31:16
	% EndTime: 2019-10-10 12:31:16
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (463->36), mult. (1311->91), div. (86->9), fcn. (1874->13), ass. (0->54)
	t119 = sin(qJ(1));
	t100 = cos(qJ(3));
	t102 = cos(qJ(1));
	t94 = sin(pkin(6));
	t109 = t102 * t94;
	t101 = cos(qJ(2));
	t104 = t119 * t101;
	t98 = sin(qJ(2));
	t108 = t102 * t98;
	t95 = cos(pkin(6));
	t89 = t95 * t108 + t104;
	t97 = sin(qJ(3));
	t78 = t100 * t109 + t89 * t97;
	t112 = t94 * t98;
	t86 = t100 * t95 - t97 * t112;
	t77 = atan2(t78, t86);
	t74 = sin(t77);
	t75 = cos(t77);
	t67 = t74 * t78 + t75 * t86;
	t66 = 0.1e1 / t67 ^ 2;
	t106 = t94 * t119;
	t105 = t119 * t98;
	t107 = t102 * t101;
	t91 = -t95 * t105 + t107;
	t81 = -t100 * t106 + t91 * t97;
	t118 = t66 * t81;
	t117 = t66 * t81 ^ 2;
	t83 = t91 * t100 + t97 * t106;
	t90 = t95 * t104 + t108;
	t96 = sin(qJ(4));
	t99 = cos(qJ(4));
	t73 = t83 * t99 + t90 * t96;
	t71 = 0.1e1 / t73 ^ 2;
	t72 = -t83 * t96 + t90 * t99;
	t116 = t72 ^ 2 * t71;
	t115 = t71 * t72;
	t114 = t75 * t78;
	t85 = 0.1e1 / t86 ^ 2;
	t113 = t78 * t85;
	t111 = t100 * t90;
	t110 = t101 * t94;
	t79 = t100 * t89 - t97 * t109;
	t103 = -t74 * t86 + t114;
	t88 = t95 * t107 - t105;
	t87 = -t100 * t112 - t95 * t97;
	t84 = 0.1e1 / t86;
	t76 = 0.1e1 / (t78 ^ 2 * t85 + 0.1e1);
	t70 = 0.1e1 / t73;
	t68 = 0.1e1 / (0.1e1 + t116);
	t65 = 0.1e1 / t67;
	t64 = 0.1e1 / (0.1e1 + t117);
	t63 = (t110 * t113 + t84 * t88) * t97 * t76;
	t62 = (-t87 * t113 + t79 * t84) * t76;
	t1 = [t81 * t84 * t76, t63, t62, 0, 0, 0; (t78 * t65 + (t74 + (t84 * t114 - t74) * t76) * t117) * t64, (t90 * t97 * t65 + ((-t75 * t110 + t74 * t88) * t97 + t103 * t63) * t118) * t64, (-t83 * t65 + (t103 * t62 + t74 * t79 + t75 * t87) * t118) * t64, 0, 0, 0; ((t79 * t96 + t88 * t99) * t70 - (-t79 * t99 + t88 * t96) * t115) * t68, ((t96 * t111 + t91 * t99) * t70 - (-t99 * t111 + t91 * t96) * t115) * t68, (t99 * t115 + t96 * t70) * t81 * t68, (-t70 * t73 - t116) * t68, 0, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end