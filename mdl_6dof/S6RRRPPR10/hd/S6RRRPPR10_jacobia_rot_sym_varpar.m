% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR10
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
%   Wie in S6RRRPPR10_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:33
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPR10_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR10_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:33:19
	% EndTime: 2019-10-10 11:33:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:33:19
	% EndTime: 2019-10-10 11:33:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:33:19
	% EndTime: 2019-10-10 11:33:19
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
	% StartTime: 2019-10-10 11:33:19
	% EndTime: 2019-10-10 11:33:20
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
	% StartTime: 2019-10-10 11:33:19
	% EndTime: 2019-10-10 11:33:20
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
	% StartTime: 2019-10-10 11:33:19
	% EndTime: 2019-10-10 11:33:20
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (428->37), mult. (1217->90), div. (80->9), fcn. (1746->13), ass. (0->53)
	t69 = cos(pkin(6));
	t71 = sin(qJ(2));
	t75 = cos(qJ(1));
	t78 = t75 * t71;
	t72 = sin(qJ(1));
	t74 = cos(qJ(2));
	t79 = t72 * t74;
	t62 = t69 * t78 + t79;
	t70 = sin(qJ(3));
	t73 = cos(qJ(3));
	t67 = sin(pkin(6));
	t81 = t67 * t75;
	t52 = t62 * t73 - t70 * t81;
	t83 = t67 * t73;
	t60 = t69 * t70 + t71 * t83;
	t50 = atan2(-t52, t60);
	t47 = sin(t50);
	t48 = cos(t50);
	t41 = -t47 * t52 + t48 * t60;
	t40 = 0.1e1 / t41 ^ 2;
	t77 = t75 * t74;
	t80 = t72 * t71;
	t64 = -t69 * t80 + t77;
	t84 = t67 * t70;
	t56 = t64 * t73 + t72 * t84;
	t90 = t40 * t56;
	t55 = t64 * t70 - t72 * t83;
	t63 = t69 * t79 + t78;
	t66 = sin(pkin(11));
	t68 = cos(pkin(11));
	t46 = t55 * t66 + t63 * t68;
	t44 = 0.1e1 / t46 ^ 2;
	t45 = -t55 * t68 + t63 * t66;
	t89 = t44 * t45;
	t88 = t48 * t52;
	t58 = 0.1e1 / t60 ^ 2;
	t87 = t52 * t58;
	t86 = t56 ^ 2 * t40;
	t85 = t63 * t70;
	t82 = t67 * t74;
	t76 = -t47 * t60 - t88;
	t51 = t62 * t70 + t73 * t81;
	t61 = t69 * t77 - t80;
	t59 = t69 * t73 - t71 * t84;
	t57 = 0.1e1 / t60;
	t49 = 0.1e1 / (t52 ^ 2 * t58 + 0.1e1);
	t43 = 0.1e1 / t46;
	t42 = 0.1e1 / (t44 * t45 ^ 2 + 0.1e1);
	t39 = 0.1e1 / t41;
	t38 = 0.1e1 / (0.1e1 + t86);
	t37 = (-t57 * t61 + t82 * t87) * t73 * t49;
	t36 = (t51 * t57 + t59 * t87) * t49;
	t1 = [-t56 * t57 * t49, t37, t36, 0, 0, 0; (-t52 * t39 - (-t47 + (t57 * t88 + t47) * t49) * t86) * t38, (-t63 * t73 * t39 - ((-t47 * t61 + t48 * t82) * t73 + t76 * t37) * t90) * t38, (-t55 * t39 - (t36 * t76 + t47 * t51 + t48 * t59) * t90) * t38, 0, 0, 0; ((t51 * t68 + t61 * t66) * t43 - (-t51 * t66 + t61 * t68) * t89) * t42, ((t64 * t66 + t68 * t85) * t43 - (t64 * t68 - t66 * t85) * t89) * t42, (-t43 * t68 - t66 * t89) * t56 * t42, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:33:20
	% EndTime: 2019-10-10 11:33:20
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (525->38), mult. (1300->91), div. (85->9), fcn. (1858->13), ass. (0->55)
	t79 = cos(pkin(6));
	t81 = sin(qJ(2));
	t85 = cos(qJ(1));
	t89 = t85 * t81;
	t82 = sin(qJ(1));
	t84 = cos(qJ(2));
	t90 = t82 * t84;
	t71 = t79 * t89 + t90;
	t80 = sin(qJ(3));
	t83 = cos(qJ(3));
	t78 = sin(pkin(6));
	t92 = t78 * t85;
	t61 = t71 * t83 - t80 * t92;
	t94 = t78 * t83;
	t69 = t79 * t80 + t81 * t94;
	t59 = atan2(-t61, t69);
	t56 = sin(t59);
	t57 = cos(t59);
	t50 = -t56 * t61 + t57 * t69;
	t49 = 0.1e1 / t50 ^ 2;
	t88 = t85 * t84;
	t91 = t82 * t81;
	t73 = -t79 * t91 + t88;
	t95 = t78 * t80;
	t65 = t73 * t83 + t82 * t95;
	t101 = t49 * t65;
	t64 = t73 * t80 - t82 * t94;
	t72 = t79 * t90 + t89;
	t77 = pkin(11) + qJ(6);
	t75 = sin(t77);
	t76 = cos(t77);
	t55 = t64 * t75 + t72 * t76;
	t53 = 0.1e1 / t55 ^ 2;
	t54 = -t64 * t76 + t72 * t75;
	t100 = t53 * t54;
	t99 = t57 * t61;
	t67 = 0.1e1 / t69 ^ 2;
	t98 = t61 * t67;
	t97 = t65 ^ 2 * t49;
	t96 = t72 * t80;
	t93 = t78 * t84;
	t87 = t54 ^ 2 * t53 + 0.1e1;
	t86 = -t56 * t69 - t99;
	t60 = t71 * t80 + t83 * t92;
	t70 = t79 * t88 - t91;
	t68 = t79 * t83 - t81 * t95;
	t66 = 0.1e1 / t69;
	t58 = 0.1e1 / (t61 ^ 2 * t67 + 0.1e1);
	t52 = 0.1e1 / t55;
	t51 = 0.1e1 / t87;
	t48 = 0.1e1 / t50;
	t47 = 0.1e1 / (0.1e1 + t97);
	t46 = (-t66 * t70 + t93 * t98) * t83 * t58;
	t45 = (t60 * t66 + t68 * t98) * t58;
	t1 = [-t65 * t66 * t58, t46, t45, 0, 0, 0; (-t61 * t48 - (-t56 + (t66 * t99 + t56) * t58) * t97) * t47, (-t72 * t83 * t48 - ((-t56 * t70 + t57 * t93) * t83 + t86 * t46) * t101) * t47, (-t64 * t48 - (t86 * t45 + t56 * t60 + t57 * t68) * t101) * t47, 0, 0, 0; ((t60 * t76 + t70 * t75) * t52 - (-t60 * t75 + t70 * t76) * t100) * t51, ((t73 * t75 + t76 * t96) * t52 - (t73 * t76 - t75 * t96) * t100) * t51, (-t75 * t100 - t76 * t52) * t65 * t51, 0, 0, t87 * t51;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end