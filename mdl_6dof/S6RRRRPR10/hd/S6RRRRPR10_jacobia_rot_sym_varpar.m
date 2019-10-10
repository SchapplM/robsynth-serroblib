% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR10
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
%   Wie in S6RRRRPR10_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:46
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR10_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR10_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:05
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
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:05
	% DurationCPUTime: 0.18s
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
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:05
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
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:05
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (992->33), mult. (1383->75), div. (104->9), fcn. (2002->11), ass. (0->51)
	t88 = cos(pkin(6));
	t89 = sin(qJ(2));
	t92 = cos(qJ(1));
	t95 = t92 * t89;
	t90 = sin(qJ(1));
	t91 = cos(qJ(2));
	t96 = t90 * t91;
	t79 = t88 * t95 + t96;
	t86 = qJ(3) + qJ(4);
	t84 = sin(t86);
	t85 = cos(t86);
	t87 = sin(pkin(6));
	t98 = t87 * t92;
	t67 = t79 * t84 + t85 * t98;
	t101 = t87 * t89;
	t74 = t84 * t101 - t88 * t85;
	t65 = atan2(-t67, t74);
	t62 = sin(t65);
	t63 = cos(t65);
	t60 = -t62 * t67 + t63 * t74;
	t59 = 0.1e1 / t60 ^ 2;
	t100 = t87 * t90;
	t94 = t92 * t91;
	t97 = t90 * t89;
	t81 = -t88 * t97 + t94;
	t70 = -t85 * t100 + t81 * t84;
	t106 = t59 * t70;
	t105 = t63 * t67;
	t73 = 0.1e1 / t74 ^ 2;
	t104 = t67 * t73;
	t103 = t70 ^ 2 * t59;
	t71 = t84 * t100 + t81 * t85;
	t80 = t88 * t96 + t95;
	t77 = 0.1e1 / t80 ^ 2;
	t102 = t71 * t77;
	t99 = t87 * t91;
	t69 = t79 * t85 - t84 * t98;
	t93 = -t62 * t74 - t105;
	t78 = t88 * t94 - t97;
	t76 = 0.1e1 / t80;
	t75 = t85 * t101 + t88 * t84;
	t72 = 0.1e1 / t74;
	t66 = 0.1e1 / (t71 ^ 2 * t77 + 0.1e1);
	t64 = 0.1e1 / (t67 ^ 2 * t73 + 0.1e1);
	t61 = t70 * t76 * t66;
	t58 = 0.1e1 / t60;
	t57 = 0.1e1 / (0.1e1 + t103);
	t56 = (t99 * t104 - t72 * t78) * t84 * t64;
	t55 = (t75 * t104 - t69 * t72) * t64;
	t54 = (t71 * t58 - (t93 * t55 - t62 * t69 + t63 * t75) * t106) * t57;
	t1 = [-t70 * t72 * t64, t56, t55, t55, 0, 0; (-t67 * t58 - (-t62 + (t72 * t105 + t62) * t64) * t103) * t57, (-t80 * t84 * t58 - ((-t62 * t78 + t63 * t99) * t84 + t93 * t56) * t106) * t57, t54, t54, 0, 0; (-t78 * t102 - t69 * t76) * t66, (-t76 * t80 * t85 - t81 * t102) * t66, -t61, -t61, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:06
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (1185->38), mult. (1738->90), div. (115->9), fcn. (2485->13), ass. (0->58)
	t106 = cos(qJ(1));
	t99 = sin(pkin(6));
	t115 = t106 * t99;
	t100 = cos(pkin(6));
	t102 = sin(qJ(2));
	t110 = t106 * t102;
	t103 = sin(qJ(1));
	t105 = cos(qJ(2));
	t111 = t103 * t105;
	t92 = t100 * t110 + t111;
	t98 = qJ(3) + qJ(4);
	t96 = sin(t98);
	t97 = cos(t98);
	t82 = -t96 * t115 + t92 * t97;
	t118 = t102 * t99;
	t90 = t100 * t96 + t97 * t118;
	t80 = atan2(-t82, t90);
	t75 = sin(t80);
	t76 = cos(t80);
	t71 = -t75 * t82 + t76 * t90;
	t70 = 0.1e1 / t71 ^ 2;
	t117 = t103 * t99;
	t109 = t106 * t105;
	t112 = t103 * t102;
	t94 = -t100 * t112 + t109;
	t86 = t96 * t117 + t94 * t97;
	t123 = t70 * t86;
	t101 = sin(qJ(6));
	t104 = cos(qJ(6));
	t93 = t100 * t111 + t110;
	t113 = t93 * t104;
	t85 = -t97 * t117 + t94 * t96;
	t78 = t85 * t101 + t113;
	t74 = 0.1e1 / t78 ^ 2;
	t114 = t93 * t101;
	t77 = -t85 * t104 + t114;
	t122 = t74 * t77;
	t121 = t76 * t82;
	t88 = 0.1e1 / t90 ^ 2;
	t120 = t82 * t88;
	t119 = t86 ^ 2 * t70;
	t116 = t105 * t99;
	t108 = t77 ^ 2 * t74 + 0.1e1;
	t107 = -t75 * t90 - t121;
	t81 = t97 * t115 + t92 * t96;
	t91 = t100 * t109 - t112;
	t89 = t100 * t97 - t96 * t118;
	t87 = 0.1e1 / t90;
	t79 = 0.1e1 / (t82 ^ 2 * t88 + 0.1e1);
	t73 = 0.1e1 / t78;
	t72 = 0.1e1 / t108;
	t69 = 0.1e1 / t71;
	t68 = 0.1e1 / (0.1e1 + t119);
	t67 = (t116 * t120 - t87 * t91) * t97 * t79;
	t66 = (t89 * t120 + t81 * t87) * t79;
	t65 = (-t101 * t122 - t104 * t73) * t86 * t72;
	t64 = (-t85 * t69 - (t107 * t66 + t75 * t81 + t76 * t89) * t123) * t68;
	t1 = [-t86 * t87 * t79, t67, t66, t66, 0, 0; (-t82 * t69 - (-t75 + (t87 * t121 + t75) * t79) * t119) * t68, (-t93 * t97 * t69 - ((t76 * t116 - t75 * t91) * t97 + t107 * t67) * t123) * t68, t64, t64, 0, 0; ((t91 * t101 + t104 * t81) * t73 - (-t101 * t81 + t91 * t104) * t122) * t72, ((t94 * t101 + t96 * t113) * t73 - (t94 * t104 - t96 * t114) * t122) * t72, t65, t65, 0, t108 * t72;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end