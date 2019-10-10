% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR8
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
%   Wie in S6RRRPPR8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:29
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPR8_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR8_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
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
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
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
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (359->30), mult. (1033->74), div. (77->9), fcn. (1491->11), ass. (0->49)
	t67 = cos(pkin(6));
	t69 = sin(qJ(2));
	t73 = cos(qJ(1));
	t76 = t73 * t69;
	t70 = sin(qJ(1));
	t72 = cos(qJ(2));
	t77 = t70 * t72;
	t60 = t67 * t76 + t77;
	t68 = sin(qJ(3));
	t71 = cos(qJ(3));
	t66 = sin(pkin(6));
	t79 = t66 * t73;
	t49 = t60 * t68 + t71 * t79;
	t82 = t66 * t68;
	t57 = -t67 * t71 + t69 * t82;
	t46 = atan2(-t49, t57);
	t42 = sin(t46);
	t43 = cos(t46);
	t41 = -t42 * t49 + t43 * t57;
	t40 = 0.1e1 / t41 ^ 2;
	t75 = t73 * t72;
	t78 = t70 * t69;
	t62 = -t67 * t78 + t75;
	t81 = t66 * t71;
	t52 = t62 * t68 - t70 * t81;
	t88 = t40 * t52;
	t87 = t43 * t49;
	t53 = t62 * t71 + t70 * t82;
	t48 = 0.1e1 / t53 ^ 2;
	t61 = -t67 * t77 - t76;
	t86 = t48 * t61;
	t55 = 0.1e1 / t57 ^ 2;
	t85 = t49 * t55;
	t84 = t52 ^ 2 * t40;
	t83 = t61 ^ 2 * t48;
	t80 = t66 * t72;
	t51 = t60 * t71 - t68 * t79;
	t74 = -t42 * t57 - t87;
	t59 = -t67 * t75 + t78;
	t58 = t67 * t68 + t69 * t81;
	t54 = 0.1e1 / t57;
	t47 = 0.1e1 / t53;
	t45 = 0.1e1 / (0.1e1 + t83);
	t44 = 0.1e1 / (t49 ^ 2 * t55 + 0.1e1);
	t39 = 0.1e1 / t41;
	t38 = 0.1e1 / (0.1e1 + t84);
	t37 = (t54 * t59 + t80 * t85) * t68 * t44;
	t36 = (-t51 * t54 + t58 * t85) * t44;
	t1 = [-t52 * t54 * t44, t37, t36, 0, 0, 0; (-t49 * t39 - (-t42 + (t54 * t87 + t42) * t44) * t84) * t38, (t61 * t68 * t39 - ((t42 * t59 + t43 * t80) * t68 + t74 * t37) * t88) * t38, (t53 * t39 - (t36 * t74 - t42 * t51 + t43 * t58) * t88) * t38, 0, 0, 0; (t59 * t47 + t51 * t86) * t45, (-t47 * t62 - t71 * t83) * t45, t52 * t45 * t86, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (151->23), mult. (457->62), div. (73->11), fcn. (682->11), ass. (0->42)
	t61 = cos(pkin(6));
	t66 = cos(qJ(2));
	t67 = cos(qJ(1));
	t68 = t67 * t66;
	t63 = sin(qJ(2));
	t64 = sin(qJ(1));
	t71 = t64 * t63;
	t55 = -t61 * t71 + t68;
	t62 = sin(qJ(3));
	t65 = cos(qJ(3));
	t60 = sin(pkin(6));
	t74 = t60 * t64;
	t45 = t55 * t62 - t65 * t74;
	t43 = 0.1e1 / t45 ^ 2;
	t46 = t55 * t65 + t62 * t74;
	t79 = t43 * t46;
	t78 = t46 ^ 2 * t43;
	t51 = -t61 * t68 + t71;
	t73 = t60 * t66;
	t50 = atan2(t51, t73);
	t48 = cos(t50);
	t77 = t48 * t51;
	t47 = sin(t50);
	t40 = t47 * t51 + t48 * t73;
	t39 = 0.1e1 / t40 ^ 2;
	t69 = t67 * t63;
	t70 = t64 * t66;
	t53 = t61 * t70 + t69;
	t76 = t53 ^ 2 * t39;
	t57 = 0.1e1 / t60;
	t58 = 0.1e1 / t66;
	t75 = t57 * t58;
	t72 = t60 * t67;
	t59 = 0.1e1 / t66 ^ 2;
	t52 = t61 * t69 + t70;
	t49 = 0.1e1 / (0.1e1 + t51 ^ 2 / t60 ^ 2 * t59);
	t42 = 0.1e1 / t45;
	t41 = 0.1e1 / (0.1e1 + t78);
	t38 = 0.1e1 / t40;
	t37 = 0.1e1 / (0.1e1 + t76);
	t36 = (t51 * t59 * t63 + t52 * t58) * t57 * t49;
	t1 = [t53 * t49 * t75, t36, 0, 0, 0, 0; (t51 * t38 + (t47 + (t75 * t77 - t47) * t49) * t76) * t37, (-t55 * t38 + (-t48 * t60 * t63 + t47 * t52 + (-t47 * t73 + t77) * t36) * t53 * t39) * t37, 0, 0, 0, 0; ((-t52 * t65 + t62 * t72) * t42 - (-t52 * t62 - t65 * t72) * t79) * t41, (-t65 * t42 + t62 * t79) * t53 * t41, (-t42 * t45 - t78) * t41, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (459->37), mult. (1300->90), div. (85->9), fcn. (1858->13), ass. (0->55)
	t76 = cos(pkin(6));
	t79 = sin(qJ(2));
	t84 = cos(qJ(1));
	t88 = t84 * t79;
	t80 = sin(qJ(1));
	t83 = cos(qJ(2));
	t89 = t80 * t83;
	t69 = t76 * t88 + t89;
	t78 = sin(qJ(3));
	t82 = cos(qJ(3));
	t75 = sin(pkin(6));
	t91 = t75 * t84;
	t59 = t69 * t82 - t78 * t91;
	t93 = t75 * t82;
	t67 = t76 * t78 + t79 * t93;
	t57 = atan2(-t59, t67);
	t54 = sin(t57);
	t55 = cos(t57);
	t48 = -t54 * t59 + t55 * t67;
	t47 = 0.1e1 / t48 ^ 2;
	t87 = t84 * t83;
	t90 = t80 * t79;
	t71 = -t76 * t90 + t87;
	t94 = t75 * t78;
	t63 = t71 * t82 + t80 * t94;
	t101 = t47 * t63;
	t62 = t71 * t78 - t80 * t93;
	t81 = cos(qJ(6));
	t70 = -t76 * t89 - t88;
	t77 = sin(qJ(6));
	t96 = t70 * t77;
	t53 = t62 * t81 + t96;
	t51 = 0.1e1 / t53 ^ 2;
	t95 = t70 * t81;
	t52 = t62 * t77 - t95;
	t100 = t51 * t52;
	t99 = t55 * t59;
	t65 = 0.1e1 / t67 ^ 2;
	t98 = t59 * t65;
	t97 = t63 ^ 2 * t47;
	t92 = t75 * t83;
	t86 = t52 ^ 2 * t51 + 0.1e1;
	t85 = -t54 * t67 - t99;
	t58 = t69 * t78 + t82 * t91;
	t68 = -t76 * t87 + t90;
	t66 = t76 * t82 - t79 * t94;
	t64 = 0.1e1 / t67;
	t56 = 0.1e1 / (t59 ^ 2 * t65 + 0.1e1);
	t50 = 0.1e1 / t53;
	t49 = 0.1e1 / t86;
	t46 = 0.1e1 / t48;
	t45 = 0.1e1 / (0.1e1 + t97);
	t44 = (t64 * t68 + t92 * t98) * t82 * t56;
	t43 = (t58 * t64 + t66 * t98) * t56;
	t1 = [-t63 * t64 * t56, t44, t43, 0, 0, 0; (-t59 * t46 - (-t54 + (t64 * t99 + t54) * t56) * t97) * t45, (t70 * t82 * t46 - ((t54 * t68 + t55 * t92) * t82 + t85 * t44) * t101) * t45, (-t62 * t46 - (t85 * t43 + t54 * t58 + t55 * t66) * t101) * t45, 0, 0, 0; ((-t58 * t77 - t68 * t81) * t50 - (-t58 * t81 + t68 * t77) * t100) * t49, ((t71 * t81 + t78 * t96) * t50 - (-t71 * t77 + t78 * t95) * t100) * t49, (-t81 * t100 + t77 * t50) * t63 * t49, 0, 0, t86 * t49;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end