% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR9
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
%   Wie in S6RRRPPR9_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:31
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPR9_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR9_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:24
	% EndTime: 2019-10-10 11:31:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:24
	% EndTime: 2019-10-10 11:31:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:24
	% EndTime: 2019-10-10 11:31:24
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
	% StartTime: 2019-10-10 11:31:24
	% EndTime: 2019-10-10 11:31:24
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
	% StartTime: 2019-10-10 11:31:24
	% EndTime: 2019-10-10 11:31:25
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (428->37), mult. (1217->90), div. (80->9), fcn. (1746->13), ass. (0->53)
	t72 = cos(pkin(6));
	t74 = sin(qJ(2));
	t78 = cos(qJ(1));
	t81 = t78 * t74;
	t75 = sin(qJ(1));
	t77 = cos(qJ(2));
	t82 = t75 * t77;
	t64 = t72 * t81 + t82;
	t73 = sin(qJ(3));
	t76 = cos(qJ(3));
	t70 = sin(pkin(6));
	t84 = t70 * t78;
	t53 = t64 * t73 + t76 * t84;
	t87 = t70 * t73;
	t61 = -t72 * t76 + t74 * t87;
	t52 = atan2(-t53, t61);
	t49 = sin(t52);
	t50 = cos(t52);
	t43 = -t49 * t53 + t50 * t61;
	t42 = 0.1e1 / t43 ^ 2;
	t80 = t78 * t77;
	t83 = t75 * t74;
	t66 = -t72 * t83 + t80;
	t86 = t70 * t76;
	t57 = t66 * t73 - t75 * t86;
	t93 = t42 * t57;
	t58 = t66 * t76 + t75 * t87;
	t65 = t72 * t82 + t81;
	t69 = sin(pkin(11));
	t71 = cos(pkin(11));
	t48 = t58 * t71 + t65 * t69;
	t46 = 0.1e1 / t48 ^ 2;
	t47 = t58 * t69 - t65 * t71;
	t92 = t46 * t47;
	t91 = t50 * t53;
	t60 = 0.1e1 / t61 ^ 2;
	t90 = t53 * t60;
	t89 = t57 ^ 2 * t42;
	t88 = t65 * t76;
	t85 = t70 * t77;
	t55 = t64 * t76 - t73 * t84;
	t79 = -t49 * t61 - t91;
	t63 = t72 * t80 - t83;
	t62 = t72 * t73 + t74 * t86;
	t59 = 0.1e1 / t61;
	t51 = 0.1e1 / (t53 ^ 2 * t60 + 0.1e1);
	t45 = 0.1e1 / t48;
	t44 = 0.1e1 / (t46 * t47 ^ 2 + 0.1e1);
	t41 = 0.1e1 / t43;
	t40 = 0.1e1 / (0.1e1 + t89);
	t39 = (-t59 * t63 + t85 * t90) * t73 * t51;
	t38 = (-t55 * t59 + t62 * t90) * t51;
	t1 = [-t57 * t59 * t51, t39, t38, 0, 0, 0; (-t53 * t41 - (-t49 + (t59 * t91 + t49) * t51) * t89) * t40, (-t65 * t73 * t41 - ((-t49 * t63 + t50 * t85) * t73 + t79 * t39) * t93) * t40, (t58 * t41 - (t38 * t79 - t49 * t55 + t50 * t62) * t93) * t40, 0, 0, 0; ((-t55 * t69 - t63 * t71) * t45 - (-t55 * t71 + t63 * t69) * t92) * t44, ((-t66 * t71 - t69 * t88) * t45 - (t66 * t69 - t71 * t88) * t92) * t44, (-t45 * t69 + t71 * t92) * t57 * t44, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:24
	% EndTime: 2019-10-10 11:31:25
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (666->39), mult. (1860->96), div. (80->9), fcn. (2626->13), ass. (0->58)
	t96 = cos(qJ(2));
	t97 = cos(qJ(1));
	t102 = t97 * t96;
	t93 = sin(qJ(2));
	t94 = sin(qJ(1));
	t105 = t94 * t93;
	t91 = cos(pkin(6));
	t101 = -t91 * t102 + t105;
	t89 = sin(pkin(6));
	t106 = t89 * t97;
	t103 = t97 * t93;
	t104 = t94 * t96;
	t84 = t91 * t103 + t104;
	t92 = sin(qJ(3));
	t95 = cos(qJ(3));
	t77 = t92 * t106 - t84 * t95;
	t88 = sin(pkin(11));
	t90 = cos(pkin(11));
	t118 = t101 * t90 + t77 * t88;
	t107 = t89 * t95;
	t74 = (t93 * t107 + t91 * t92) * t88 + t89 * t96 * t90;
	t63 = atan2(t118, t74);
	t60 = sin(t63);
	t61 = cos(t63);
	t59 = t118 * t60 + t61 * t74;
	t58 = 0.1e1 / t59 ^ 2;
	t85 = t91 * t104 + t103;
	t110 = t85 * t90;
	t108 = t89 * t92;
	t86 = -t91 * t105 + t102;
	t79 = t94 * t108 + t86 * t95;
	t69 = t79 * t88 - t110;
	t117 = t58 * t69;
	t116 = t61 * t118;
	t70 = t79 * t90 + t85 * t88;
	t66 = 0.1e1 / t70 ^ 2;
	t78 = t94 * t107 - t86 * t92;
	t115 = t66 * t78;
	t73 = 0.1e1 / t74 ^ 2;
	t114 = t118 * t73;
	t113 = t69 ^ 2 * t58;
	t112 = t78 ^ 2 * t66;
	t109 = t88 * t95;
	t100 = t101 * t88;
	t99 = -t60 * t74 + t116;
	t98 = t95 * t106 + t84 * t92;
	t83 = -t93 * t108 + t91 * t95;
	t80 = (t96 * t109 - t90 * t93) * t89;
	t72 = 0.1e1 / t74;
	t71 = -t95 * t100 - t84 * t90;
	t65 = 0.1e1 / t70;
	t64 = 0.1e1 / (0.1e1 + t112);
	t62 = 0.1e1 / (t118 ^ 2 * t73 + 0.1e1);
	t57 = 0.1e1 / t59;
	t56 = 0.1e1 / (0.1e1 + t113);
	t55 = (-t83 * t114 + t72 * t98) * t88 * t62;
	t54 = (-t80 * t114 - t71 * t72) * t62;
	t1 = [-t69 * t72 * t62, t54, t55, 0, 0, 0; (t118 * t57 - (-t60 + (-t72 * t116 + t60) * t62) * t113) * t56, ((-t85 * t109 - t86 * t90) * t57 - (t99 * t54 - t60 * t71 + t61 * t80) * t117) * t56, (t78 * t88 * t57 - ((t60 * t98 + t61 * t83) * t88 + t99 * t55) * t117) * t56, 0, 0, 0; (t98 * t65 - (t77 * t90 - t100) * t115) * t64, (t85 * t92 * t65 - (-t95 * t110 + t86 * t88) * t115) * t64, (-t90 * t112 - t65 * t79) * t64, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:24
	% EndTime: 2019-10-10 11:31:25
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (623->43), mult. (1730->106), div. (85->9), fcn. (2452->15), ass. (0->63)
	t129 = sin(qJ(1));
	t105 = sin(qJ(3));
	t108 = cos(qJ(3));
	t110 = cos(qJ(1));
	t121 = sin(pkin(6));
	t113 = t110 * t121;
	t103 = cos(pkin(6));
	t109 = cos(qJ(2));
	t116 = t129 * t109;
	t106 = sin(qJ(2));
	t120 = t110 * t106;
	t96 = t103 * t120 + t116;
	t85 = t96 * t105 + t108 * t113;
	t115 = t106 * t121;
	t93 = t103 * t108 - t105 * t115;
	t82 = atan2(t85, t93);
	t79 = sin(t82);
	t80 = cos(t82);
	t70 = t79 * t85 + t80 * t93;
	t69 = 0.1e1 / t70 ^ 2;
	t112 = t121 * t129;
	t117 = t129 * t106;
	t119 = t110 * t109;
	t98 = -t103 * t117 + t119;
	t88 = t98 * t105 - t108 * t112;
	t128 = t69 * t88;
	t104 = sin(qJ(6));
	t107 = cos(qJ(6));
	t101 = sin(pkin(11));
	t102 = cos(pkin(11));
	t97 = t103 * t116 + t120;
	t122 = t97 * t102;
	t90 = t105 * t112 + t98 * t108;
	t77 = t90 * t101 - t122;
	t123 = t97 * t101;
	t78 = t90 * t102 + t123;
	t74 = t77 * t104 + t78 * t107;
	t72 = 0.1e1 / t74 ^ 2;
	t73 = t78 * t104 - t77 * t107;
	t127 = t72 * t73;
	t126 = t80 * t85;
	t92 = 0.1e1 / t93 ^ 2;
	t125 = t85 * t92;
	t124 = t88 ^ 2 * t69;
	t118 = t73 ^ 2 * t72 + 0.1e1;
	t114 = t109 * t121;
	t86 = -t105 * t113 + t96 * t108;
	t111 = -t79 * t93 + t126;
	t95 = t103 * t119 - t117;
	t94 = -t103 * t105 - t108 * t115;
	t91 = 0.1e1 / t93;
	t84 = t98 * t101 - t108 * t122;
	t83 = -t98 * t102 - t108 * t123;
	t81 = 0.1e1 / (t85 ^ 2 * t92 + 0.1e1);
	t76 = t95 * t101 - t102 * t86;
	t75 = -t101 * t86 - t95 * t102;
	t71 = 0.1e1 / t74;
	t68 = 0.1e1 / t70;
	t67 = 0.1e1 / (0.1e1 + t124);
	t66 = (t114 * t125 + t91 * t95) * t81 * t105;
	t65 = (-t94 * t125 + t86 * t91) * t81;
	t64 = 0.1e1 / t118;
	t1 = [t88 * t91 * t81, t66, t65, 0, 0, 0; (t85 * t68 + (t79 + (t91 * t126 - t79) * t81) * t124) * t67, (t97 * t105 * t68 + (t111 * t66 + (-t80 * t114 + t79 * t95) * t105) * t128) * t67, (-t90 * t68 + (t111 * t65 + t79 * t86 + t80 * t94) * t128) * t67, 0, 0, 0; ((t76 * t104 - t75 * t107) * t71 - (t75 * t104 + t76 * t107) * t127) * t64, ((t84 * t104 - t83 * t107) * t71 - (t83 * t104 + t84 * t107) * t127) * t64, ((t101 * t107 - t102 * t104) * t71 - (-t101 * t104 - t102 * t107) * t127) * t64 * t88, 0, 0, t118 * t64;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end