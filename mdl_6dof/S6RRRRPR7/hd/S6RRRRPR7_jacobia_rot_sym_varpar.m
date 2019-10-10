% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR7
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
%   Wie in S6RRRRPR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:42
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR7_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR7_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:42:14
	% EndTime: 2019-10-10 12:42:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:42:14
	% EndTime: 2019-10-10 12:42:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:42:14
	% EndTime: 2019-10-10 12:42:14
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
	% StartTime: 2019-10-10 12:42:14
	% EndTime: 2019-10-10 12:42:14
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
	% StartTime: 2019-10-10 12:42:14
	% EndTime: 2019-10-10 12:42:14
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (252->25), mult. (499->64), div. (77->11), fcn. (736->11), ass. (0->42)
	t69 = sin(qJ(2));
	t71 = cos(qJ(2));
	t72 = cos(qJ(1));
	t70 = sin(qJ(1));
	t76 = cos(pkin(6));
	t74 = t70 * t76;
	t60 = -t69 * t74 + t72 * t71;
	t67 = qJ(3) + qJ(4);
	t62 = sin(t67);
	t63 = cos(t67);
	t68 = sin(pkin(6));
	t79 = t68 * t70;
	t51 = t60 * t63 + t62 * t79;
	t49 = 0.1e1 / t51 ^ 2;
	t50 = t60 * t62 - t63 * t79;
	t83 = t49 * t50;
	t73 = t72 * t76;
	t56 = t70 * t69 - t71 * t73;
	t78 = t68 * t71;
	t54 = atan2(-t56, -t78);
	t53 = cos(t54);
	t82 = t53 * t56;
	t52 = sin(t54);
	t46 = -t52 * t56 - t53 * t78;
	t45 = 0.1e1 / t46 ^ 2;
	t59 = t72 * t69 + t71 * t74;
	t81 = t59 ^ 2 * t45;
	t64 = 0.1e1 / t68;
	t65 = 0.1e1 / t71;
	t80 = t64 * t65;
	t77 = t68 * t72;
	t75 = t50 ^ 2 * t49 + 0.1e1;
	t66 = 0.1e1 / t71 ^ 2;
	t58 = t69 * t73 + t70 * t71;
	t55 = 0.1e1 / (0.1e1 + t56 ^ 2 / t68 ^ 2 * t66);
	t48 = 0.1e1 / t51;
	t47 = 0.1e1 / t75;
	t44 = 0.1e1 / t46;
	t43 = 0.1e1 / (0.1e1 + t81);
	t42 = (t56 * t66 * t69 + t58 * t65) * t64 * t55;
	t41 = t75 * t47;
	t1 = [t59 * t55 * t80, t42, 0, 0, 0, 0; (-t56 * t44 - (-t52 + (-t80 * t82 + t52) * t55) * t81) * t43, (t60 * t44 - (t53 * t68 * t69 - t52 * t58 + (t52 * t78 - t82) * t42) * t59 * t45) * t43, 0, 0, 0, 0; ((-t58 * t62 - t63 * t77) * t48 - (-t58 * t63 + t62 * t77) * t83) * t47, (-t62 * t48 + t63 * t83) * t59 * t47, t41, t41, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:42:14
	% EndTime: 2019-10-10 12:42:14
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (314->25), mult. (499->64), div. (77->11), fcn. (736->11), ass. (0->42)
	t74 = sin(qJ(2));
	t76 = cos(qJ(2));
	t77 = cos(qJ(1));
	t75 = sin(qJ(1));
	t81 = cos(pkin(6));
	t79 = t75 * t81;
	t65 = -t74 * t79 + t77 * t76;
	t69 = qJ(3) + qJ(4) + pkin(12);
	t67 = sin(t69);
	t68 = cos(t69);
	t73 = sin(pkin(6));
	t84 = t73 * t75;
	t56 = t65 * t68 + t67 * t84;
	t54 = 0.1e1 / t56 ^ 2;
	t55 = t65 * t67 - t68 * t84;
	t88 = t54 * t55;
	t78 = t77 * t81;
	t61 = t75 * t74 - t76 * t78;
	t83 = t73 * t76;
	t59 = atan2(-t61, -t83);
	t58 = cos(t59);
	t87 = t58 * t61;
	t57 = sin(t59);
	t51 = -t57 * t61 - t58 * t83;
	t50 = 0.1e1 / t51 ^ 2;
	t64 = t77 * t74 + t76 * t79;
	t86 = t64 ^ 2 * t50;
	t70 = 0.1e1 / t73;
	t71 = 0.1e1 / t76;
	t85 = t70 * t71;
	t82 = t73 * t77;
	t80 = t55 ^ 2 * t54 + 0.1e1;
	t72 = 0.1e1 / t76 ^ 2;
	t63 = t74 * t78 + t75 * t76;
	t60 = 0.1e1 / (0.1e1 + t61 ^ 2 / t73 ^ 2 * t72);
	t53 = 0.1e1 / t56;
	t52 = 0.1e1 / t80;
	t49 = 0.1e1 / t51;
	t48 = 0.1e1 / (0.1e1 + t86);
	t47 = (t61 * t72 * t74 + t63 * t71) * t70 * t60;
	t46 = t80 * t52;
	t1 = [t64 * t60 * t85, t47, 0, 0, 0, 0; (-t61 * t49 - (-t57 + (-t85 * t87 + t57) * t60) * t86) * t48, (t65 * t49 - (t58 * t73 * t74 - t57 * t63 + (t57 * t83 - t87) * t47) * t64 * t50) * t48, 0, 0, 0, 0; ((-t63 * t67 - t68 * t82) * t53 - (-t63 * t68 + t67 * t82) * t88) * t52, (-t67 * t53 + t68 * t88) * t64 * t52, t46, t46, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:42:14
	% EndTime: 2019-10-10 12:42:15
	% DurationCPUTime: 0.37s
	% Computational Cost: add. (1756->38), mult. (1738->90), div. (115->9), fcn. (2485->13), ass. (0->58)
	t104 = qJ(3) + qJ(4) + pkin(12);
	t102 = sin(t104);
	t103 = cos(t104);
	t105 = sin(pkin(6));
	t112 = cos(qJ(1));
	t119 = t105 * t112;
	t106 = cos(pkin(6));
	t108 = sin(qJ(2));
	t116 = t112 * t108;
	t109 = sin(qJ(1));
	t111 = cos(qJ(2));
	t117 = t109 * t111;
	t97 = t106 * t116 + t117;
	t86 = t97 * t102 + t103 * t119;
	t122 = t105 * t108;
	t94 = t102 * t122 - t106 * t103;
	t81 = atan2(-t86, t94);
	t78 = sin(t81);
	t79 = cos(t81);
	t76 = -t78 * t86 + t79 * t94;
	t75 = 0.1e1 / t76 ^ 2;
	t121 = t105 * t109;
	t115 = t112 * t111;
	t118 = t109 * t108;
	t99 = -t106 * t118 + t115;
	t90 = t99 * t102 - t103 * t121;
	t129 = t75 * t90;
	t128 = t79 * t86;
	t110 = cos(qJ(6));
	t107 = sin(qJ(6));
	t98 = t106 * t117 + t116;
	t124 = t98 * t107;
	t91 = t102 * t121 + t99 * t103;
	t85 = t91 * t110 + t124;
	t83 = 0.1e1 / t85 ^ 2;
	t123 = t98 * t110;
	t84 = t91 * t107 - t123;
	t127 = t83 * t84;
	t93 = 0.1e1 / t94 ^ 2;
	t126 = t86 * t93;
	t125 = t90 ^ 2 * t75;
	t120 = t105 * t111;
	t114 = t84 ^ 2 * t83 + 0.1e1;
	t88 = -t102 * t119 + t97 * t103;
	t113 = -t78 * t94 - t128;
	t96 = t106 * t115 - t118;
	t95 = t106 * t102 + t103 * t122;
	t92 = 0.1e1 / t94;
	t82 = 0.1e1 / t85;
	t80 = 0.1e1 / (t86 ^ 2 * t93 + 0.1e1);
	t77 = 0.1e1 / t114;
	t74 = 0.1e1 / t76;
	t73 = 0.1e1 / (0.1e1 + t125);
	t72 = (t120 * t126 - t92 * t96) * t80 * t102;
	t71 = (t95 * t126 - t88 * t92) * t80;
	t70 = (-t107 * t82 + t110 * t127) * t90 * t77;
	t69 = (t91 * t74 - (t113 * t71 - t78 * t88 + t79 * t95) * t129) * t73;
	t1 = [-t90 * t92 * t80, t72, t71, t71, 0, 0; (-t86 * t74 - (-t78 + (t92 * t128 + t78) * t80) * t125) * t73, (-t98 * t102 * t74 - (t113 * t72 + (t79 * t120 - t78 * t96) * t102) * t129) * t73, t69, t69, 0, 0; ((-t107 * t88 - t96 * t110) * t82 - (t96 * t107 - t110 * t88) * t127) * t77, ((-t103 * t124 - t99 * t110) * t82 - (-t103 * t123 + t99 * t107) * t127) * t77, t70, t70, 0, t114 * t77;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end