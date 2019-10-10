% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRP12
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
%   Wie in S6RRRPRP12_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:53
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRP12_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP12_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:53:24
	% EndTime: 2019-10-10 11:53:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:53:24
	% EndTime: 2019-10-10 11:53:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:53:24
	% EndTime: 2019-10-10 11:53:24
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
	% StartTime: 2019-10-10 11:53:24
	% EndTime: 2019-10-10 11:53:25
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
	% StartTime: 2019-10-10 11:53:24
	% EndTime: 2019-10-10 11:53:25
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (355->31), mult. (1022->76), div. (77->9), fcn. (1479->11), ass. (0->48)
	t63 = cos(pkin(6));
	t65 = sin(qJ(2));
	t69 = cos(qJ(1));
	t72 = t69 * t65;
	t66 = sin(qJ(1));
	t68 = cos(qJ(2));
	t73 = t66 * t68;
	t57 = t63 * t72 + t73;
	t64 = sin(qJ(3));
	t67 = cos(qJ(3));
	t62 = sin(pkin(6));
	t75 = t62 * t69;
	t45 = t57 * t64 + t67 * t75;
	t78 = t62 * t64;
	t54 = -t63 * t67 + t65 * t78;
	t44 = atan2(-t45, t54);
	t40 = sin(t44);
	t41 = cos(t44);
	t39 = -t40 * t45 + t41 * t54;
	t38 = 0.1e1 / t39 ^ 2;
	t71 = t69 * t68;
	t74 = t66 * t65;
	t59 = -t63 * t74 + t71;
	t77 = t62 * t67;
	t48 = t59 * t64 - t66 * t77;
	t83 = t38 * t48;
	t82 = t41 * t45;
	t51 = 0.1e1 / t54 ^ 2;
	t81 = t45 * t51;
	t80 = t48 ^ 2 * t38;
	t49 = t59 * t67 + t66 * t78;
	t58 = t63 * t73 + t72;
	t53 = 0.1e1 / t58 ^ 2;
	t79 = t49 * t53;
	t76 = t62 * t68;
	t47 = t57 * t67 - t64 * t75;
	t70 = -t40 * t54 - t82;
	t56 = t63 * t71 - t74;
	t55 = t63 * t64 + t65 * t77;
	t52 = 0.1e1 / t58;
	t50 = 0.1e1 / t54;
	t43 = 0.1e1 / (t49 ^ 2 * t53 + 0.1e1);
	t42 = 0.1e1 / (t45 ^ 2 * t51 + 0.1e1);
	t37 = 0.1e1 / t39;
	t36 = 0.1e1 / (0.1e1 + t80);
	t35 = (-t50 * t56 + t76 * t81) * t64 * t42;
	t34 = (-t47 * t50 + t55 * t81) * t42;
	t1 = [-t48 * t50 * t42, t35, t34, 0, 0, 0; (-t45 * t37 - (-t40 + (t50 * t82 + t40) * t42) * t80) * t36, (-t58 * t64 * t37 - ((-t40 * t56 + t41 * t76) * t64 + t70 * t35) * t83) * t36, (t49 * t37 - (t70 * t34 - t40 * t47 + t41 * t55) * t83) * t36, 0, 0, 0; (-t47 * t52 - t56 * t79) * t43, (-t52 * t58 * t67 - t59 * t79) * t43, -t48 * t52 * t43, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:53:24
	% EndTime: 2019-10-10 11:53:25
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (459->37), mult. (1300->90), div. (85->9), fcn. (1858->13), ass. (0->55)
	t72 = cos(pkin(6));
	t75 = sin(qJ(2));
	t80 = cos(qJ(1));
	t84 = t80 * t75;
	t76 = sin(qJ(1));
	t79 = cos(qJ(2));
	t85 = t76 * t79;
	t67 = t72 * t84 + t85;
	t74 = sin(qJ(3));
	t78 = cos(qJ(3));
	t71 = sin(pkin(6));
	t87 = t71 * t80;
	t57 = t67 * t78 - t74 * t87;
	t89 = t71 * t78;
	t65 = t72 * t74 + t75 * t89;
	t55 = atan2(-t57, t65);
	t52 = sin(t55);
	t53 = cos(t55);
	t46 = -t52 * t57 + t53 * t65;
	t45 = 0.1e1 / t46 ^ 2;
	t83 = t80 * t79;
	t86 = t76 * t75;
	t69 = -t72 * t86 + t83;
	t90 = t71 * t74;
	t61 = t69 * t78 + t76 * t90;
	t97 = t45 * t61;
	t60 = t69 * t74 - t76 * t89;
	t73 = sin(qJ(5));
	t68 = t72 * t85 + t84;
	t77 = cos(qJ(5));
	t91 = t68 * t77;
	t51 = t60 * t73 + t91;
	t49 = 0.1e1 / t51 ^ 2;
	t92 = t68 * t73;
	t50 = -t60 * t77 + t92;
	t96 = t49 * t50;
	t95 = t53 * t57;
	t63 = 0.1e1 / t65 ^ 2;
	t94 = t57 * t63;
	t93 = t61 ^ 2 * t45;
	t88 = t71 * t79;
	t82 = t50 ^ 2 * t49 + 0.1e1;
	t81 = -t52 * t65 - t95;
	t56 = t67 * t74 + t78 * t87;
	t66 = t72 * t83 - t86;
	t64 = t72 * t78 - t75 * t90;
	t62 = 0.1e1 / t65;
	t54 = 0.1e1 / (t57 ^ 2 * t63 + 0.1e1);
	t48 = 0.1e1 / t51;
	t47 = 0.1e1 / t82;
	t44 = 0.1e1 / t46;
	t43 = 0.1e1 / (0.1e1 + t93);
	t42 = (-t62 * t66 + t88 * t94) * t78 * t54;
	t41 = (t56 * t62 + t64 * t94) * t54;
	t1 = [-t61 * t62 * t54, t42, t41, 0, 0, 0; (-t57 * t44 - (-t52 + (t62 * t95 + t52) * t54) * t93) * t43, (-t68 * t78 * t44 - ((-t52 * t66 + t53 * t88) * t78 + t81 * t42) * t97) * t43, (-t60 * t44 - (t81 * t41 + t52 * t56 + t53 * t64) * t97) * t43, 0, 0, 0; ((t56 * t77 + t66 * t73) * t48 - (-t56 * t73 + t66 * t77) * t96) * t47, ((t69 * t73 + t74 * t91) * t48 - (t69 * t77 - t74 * t92) * t96) * t47, (-t77 * t48 - t73 * t96) * t61 * t47, 0, t82 * t47, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:53:25
	% EndTime: 2019-10-10 11:53:25
	% DurationCPUTime: 0.47s
	% Computational Cost: add. (907->45), mult. (2525->110), div. (107->9), fcn. (3566->13), ass. (0->62)
	t104 = sin(qJ(5));
	t108 = cos(qJ(5));
	t105 = sin(qJ(3));
	t109 = cos(qJ(3));
	t102 = sin(pkin(6));
	t111 = cos(qJ(1));
	t121 = t102 * t111;
	t103 = cos(pkin(6));
	t106 = sin(qJ(2));
	t116 = t111 * t106;
	t107 = sin(qJ(1));
	t110 = cos(qJ(2));
	t118 = t107 * t110;
	t98 = t103 * t116 + t118;
	t113 = t98 * t105 + t109 * t121;
	t115 = t111 * t110;
	t119 = t107 * t106;
	t97 = -t103 * t115 + t119;
	t80 = t113 * t104 + t97 * t108;
	t133 = -t97 * t104 + t113 * t108;
	t123 = t102 * t105;
	t95 = -t103 * t109 + t106 * t123;
	t88 = -t102 * t110 * t104 - t95 * t108;
	t75 = atan2(t133, t88);
	t71 = sin(t75);
	t72 = cos(t75);
	t70 = t133 * t71 + t72 * t88;
	t69 = 0.1e1 / t70 ^ 2;
	t100 = -t103 * t119 + t115;
	t122 = t102 * t109;
	t112 = t100 * t105 - t107 * t122;
	t99 = t103 * t118 + t116;
	t124 = t99 * t104;
	t81 = -t112 * t108 + t124;
	t132 = t69 * t81;
	t131 = t72 * t133;
	t82 = t112 * t104 + t99 * t108;
	t77 = 0.1e1 / t82 ^ 2;
	t92 = t100 * t109 + t107 * t123;
	t130 = t77 * t92;
	t87 = 0.1e1 / t88 ^ 2;
	t129 = t133 * t87;
	t128 = t81 ^ 2 * t69;
	t127 = t92 ^ 2 * t77;
	t120 = t105 * t108;
	t117 = t108 * t110;
	t114 = -t71 * t88 + t131;
	t90 = -t105 * t121 + t98 * t109;
	t96 = t103 * t105 + t106 * t122;
	t94 = (t104 * t106 - t105 * t117) * t102;
	t89 = -t102 * t117 + t95 * t104;
	t86 = 0.1e1 / t88;
	t83 = t98 * t104 + t97 * t120;
	t76 = 0.1e1 / t82;
	t74 = 0.1e1 / (0.1e1 + t127);
	t73 = 0.1e1 / (t133 ^ 2 * t87 + 0.1e1);
	t68 = 0.1e1 / t70;
	t67 = 0.1e1 / (0.1e1 + t128);
	t66 = (t96 * t129 + t86 * t90) * t73 * t108;
	t65 = (-t94 * t129 - t83 * t86) * t73;
	t64 = (-t89 * t129 - t80 * t86) * t73;
	t1 = [-t81 * t86 * t73, t65, t66, 0, t64, 0; (t133 * t68 - (-t71 + (-t86 * t131 + t71) * t73) * t128) * t67, ((t100 * t104 + t99 * t120) * t68 - (t114 * t65 - t71 * t83 + t72 * t94) * t132) * t67, (-t92 * t108 * t68 - (t114 * t66 + (t71 * t90 - t72 * t96) * t108) * t132) * t67, 0, (t82 * t68 - (t114 * t64 - t71 * t80 + t72 * t89) * t132) * t67, 0; (-t80 * t130 + t90 * t76) * t74, (t99 * t109 * t76 + (t100 * t108 - t105 * t124) * t130) * t74, (t104 * t127 + t112 * t76) * t74, 0, -t81 * t74 * t130, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end