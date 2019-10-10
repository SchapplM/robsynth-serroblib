% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR5
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
%   Wie in S6RRRPPR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:24
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPR5_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR5_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:24:03
	% EndTime: 2019-10-10 11:24:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:24:03
	% EndTime: 2019-10-10 11:24:03
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:24:03
	% EndTime: 2019-10-10 11:24:03
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
	% StartTime: 2019-10-10 11:24:03
	% EndTime: 2019-10-10 11:24:03
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
	% StartTime: 2019-10-10 11:24:03
	% EndTime: 2019-10-10 11:24:03
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (221->25), mult. (451->64), div. (72->11), fcn. (673->11), ass. (0->41)
	t59 = sin(qJ(2));
	t61 = cos(qJ(2));
	t62 = cos(qJ(1));
	t60 = sin(qJ(1));
	t66 = cos(pkin(6));
	t64 = t60 * t66;
	t50 = -t59 * t64 + t62 * t61;
	t55 = qJ(3) + pkin(11);
	t52 = sin(t55);
	t53 = cos(t55);
	t58 = sin(pkin(6));
	t69 = t58 * t60;
	t41 = t50 * t53 + t52 * t69;
	t39 = 0.1e1 / t41 ^ 2;
	t40 = t50 * t52 - t53 * t69;
	t73 = t39 * t40;
	t63 = t62 * t66;
	t46 = t60 * t59 - t61 * t63;
	t68 = t58 * t61;
	t44 = atan2(-t46, -t68);
	t43 = cos(t44);
	t72 = t43 * t46;
	t42 = sin(t44);
	t36 = -t42 * t46 - t43 * t68;
	t35 = 0.1e1 / t36 ^ 2;
	t49 = t62 * t59 + t61 * t64;
	t71 = t49 ^ 2 * t35;
	t54 = 0.1e1 / t58;
	t56 = 0.1e1 / t61;
	t70 = t54 * t56;
	t67 = t58 * t62;
	t65 = t40 ^ 2 * t39 + 0.1e1;
	t57 = 0.1e1 / t61 ^ 2;
	t48 = t59 * t63 + t60 * t61;
	t45 = 0.1e1 / (0.1e1 + t46 ^ 2 / t58 ^ 2 * t57);
	t38 = 0.1e1 / t41;
	t37 = 0.1e1 / t65;
	t34 = 0.1e1 / t36;
	t33 = 0.1e1 / (0.1e1 + t71);
	t32 = (t46 * t57 * t59 + t48 * t56) * t54 * t45;
	t1 = [t49 * t45 * t70, t32, 0, 0, 0, 0; (-t46 * t34 - (-t42 + (-t70 * t72 + t42) * t45) * t71) * t33, (t50 * t34 - (t43 * t58 * t59 - t42 * t48 + (t42 * t68 - t72) * t32) * t49 * t35) * t33, 0, 0, 0, 0; ((-t48 * t52 - t53 * t67) * t38 - (-t48 * t53 + t52 * t67) * t73) * t37, (-t52 * t38 + t53 * t73) * t49 * t37, t65 * t37, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:24:03
	% EndTime: 2019-10-10 11:24:03
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (833->38), mult. (1217->89), div. (80->9), fcn. (1746->13), ass. (0->55)
	t81 = cos(pkin(6));
	t82 = sin(qJ(2));
	t85 = cos(qJ(1));
	t88 = t85 * t82;
	t83 = sin(qJ(1));
	t84 = cos(qJ(2));
	t89 = t83 * t84;
	t70 = t81 * t88 + t89;
	t77 = qJ(3) + pkin(11);
	t75 = sin(t77);
	t76 = cos(t77);
	t79 = sin(pkin(6));
	t91 = t79 * t85;
	t59 = t70 * t75 + t76 * t91;
	t94 = t79 * t82;
	t67 = t75 * t94 - t81 * t76;
	t58 = atan2(-t59, t67);
	t53 = sin(t58);
	t54 = cos(t58);
	t49 = -t53 * t59 + t54 * t67;
	t48 = 0.1e1 / t49 ^ 2;
	t87 = t85 * t84;
	t90 = t83 * t82;
	t72 = -t81 * t90 + t87;
	t93 = t79 * t83;
	t63 = t72 * t75 - t76 * t93;
	t101 = t48 * t63;
	t64 = t72 * t76 + t75 * t93;
	t80 = cos(pkin(12));
	t71 = t81 * t89 + t88;
	t78 = sin(pkin(12));
	t96 = t71 * t78;
	t56 = t64 * t80 + t96;
	t52 = 0.1e1 / t56 ^ 2;
	t95 = t71 * t80;
	t55 = t64 * t78 - t95;
	t100 = t52 * t55;
	t99 = t54 * t59;
	t66 = 0.1e1 / t67 ^ 2;
	t98 = t59 * t66;
	t97 = t63 ^ 2 * t48;
	t92 = t79 * t84;
	t61 = t70 * t76 - t75 * t91;
	t86 = -t53 * t67 - t99;
	t69 = t81 * t87 - t90;
	t68 = t81 * t75 + t76 * t94;
	t65 = 0.1e1 / t67;
	t57 = 0.1e1 / (t59 ^ 2 * t66 + 0.1e1);
	t51 = 0.1e1 / t56;
	t50 = 0.1e1 / (t55 ^ 2 * t52 + 0.1e1);
	t47 = 0.1e1 / t49;
	t46 = 0.1e1 / (0.1e1 + t97);
	t45 = (-t65 * t69 + t92 * t98) * t75 * t57;
	t44 = (-t61 * t65 + t68 * t98) * t57;
	t1 = [-t63 * t65 * t57, t45, t44, 0, 0, 0; (-t59 * t47 - (-t53 + (t65 * t99 + t53) * t57) * t97) * t46, (-t71 * t75 * t47 - ((-t53 * t69 + t54 * t92) * t75 + t86 * t45) * t101) * t46, (t64 * t47 - (t86 * t44 - t53 * t61 + t54 * t68) * t101) * t46, 0, 0, 0; ((-t61 * t78 - t69 * t80) * t51 - (-t61 * t80 + t69 * t78) * t100) * t50, ((-t72 * t80 - t76 * t96) * t51 - (t72 * t78 - t76 * t95) * t100) * t50, (t80 * t100 - t78 * t51) * t63 * t50, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:24:03
	% EndTime: 2019-10-10 11:24:03
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (944->39), mult. (1300->91), div. (85->9), fcn. (1858->13), ass. (0->56)
	t88 = sin(pkin(6));
	t93 = cos(qJ(1));
	t100 = t88 * t93;
	t89 = cos(pkin(6));
	t90 = sin(qJ(2));
	t97 = t93 * t90;
	t91 = sin(qJ(1));
	t92 = cos(qJ(2));
	t98 = t91 * t92;
	t77 = t89 * t97 + t98;
	t87 = qJ(3) + pkin(11);
	t83 = sin(t87);
	t85 = cos(t87);
	t66 = t85 * t100 + t77 * t83;
	t103 = t88 * t90;
	t74 = t83 * t103 - t89 * t85;
	t65 = atan2(-t66, t74);
	t62 = sin(t65);
	t63 = cos(t65);
	t56 = -t62 * t66 + t63 * t74;
	t55 = 0.1e1 / t56 ^ 2;
	t102 = t88 * t91;
	t96 = t93 * t92;
	t99 = t91 * t90;
	t79 = -t89 * t99 + t96;
	t70 = -t85 * t102 + t79 * t83;
	t109 = t55 * t70;
	t71 = t83 * t102 + t79 * t85;
	t78 = t89 * t98 + t97;
	t86 = pkin(12) + qJ(6);
	t82 = sin(t86);
	t84 = cos(t86);
	t61 = t71 * t84 + t78 * t82;
	t59 = 0.1e1 / t61 ^ 2;
	t60 = t71 * t82 - t78 * t84;
	t108 = t59 * t60;
	t107 = t63 * t66;
	t73 = 0.1e1 / t74 ^ 2;
	t106 = t66 * t73;
	t105 = t70 ^ 2 * t55;
	t104 = t78 * t85;
	t101 = t88 * t92;
	t95 = t60 ^ 2 * t59 + 0.1e1;
	t68 = -t83 * t100 + t77 * t85;
	t94 = -t62 * t74 - t107;
	t76 = t89 * t96 - t99;
	t75 = t85 * t103 + t89 * t83;
	t72 = 0.1e1 / t74;
	t64 = 0.1e1 / (t66 ^ 2 * t73 + 0.1e1);
	t58 = 0.1e1 / t61;
	t57 = 0.1e1 / t95;
	t54 = 0.1e1 / t56;
	t53 = 0.1e1 / (0.1e1 + t105);
	t52 = (t101 * t106 - t72 * t76) * t83 * t64;
	t51 = (t75 * t106 - t68 * t72) * t64;
	t1 = [-t70 * t72 * t64, t52, t51, 0, 0, 0; (-t66 * t54 - (-t62 + (t72 * t107 + t62) * t64) * t105) * t53, (-t78 * t83 * t54 - ((t63 * t101 - t62 * t76) * t83 + t94 * t52) * t109) * t53, (t71 * t54 - (t94 * t51 - t62 * t68 + t63 * t75) * t109) * t53, 0, 0, 0; ((-t68 * t82 - t76 * t84) * t58 - (-t68 * t84 + t76 * t82) * t108) * t57, ((-t82 * t104 - t79 * t84) * t58 - (-t84 * t104 + t79 * t82) * t108) * t57, (t84 * t108 - t82 * t58) * t70 * t57, 0, 0, t95 * t57;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end