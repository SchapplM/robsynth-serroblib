% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5RRRPR12_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRRPR12_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR12_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:42:20
	% EndTime: 2019-12-31 21:42:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:42:20
	% EndTime: 2019-12-31 21:42:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:42:20
	% EndTime: 2019-12-31 21:42:20
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (35->13), mult. (116->34), div. (25->9), fcn. (175->9), ass. (0->25)
	t30 = cos(pkin(5));
	t29 = sin(pkin(5));
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
	t1 = [-t32 * t29 * t41, 0, 0, 0, 0; (0.1e1 / t18 * t38 - (-t24 * t27 * t34 * t41 + (t25 - 0.1e1) * t29 * t23) * t29 * t42) / (t27 * t42 + 0.1e1), 0, 0, 0, 0; ((t30 * t36 - t40) / t22 - (-t30 * t37 - t39) * t21 * t20) * t19, t35 * t19, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:42:20
	% EndTime: 2019-12-31 21:42:20
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (173->24), mult. (451->64), div. (72->11), fcn. (673->11), ass. (0->40)
	t53 = sin(qJ(2));
	t56 = cos(qJ(2));
	t57 = cos(qJ(1));
	t54 = sin(qJ(1));
	t61 = cos(pkin(5));
	t59 = t54 * t61;
	t46 = -t53 * t59 + t57 * t56;
	t52 = sin(qJ(3));
	t55 = cos(qJ(3));
	t51 = sin(pkin(5));
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
	t1 = [t45 * t41 * t65, t28, 0, 0, 0; (-t42 * t30 - (-t38 + (-t65 * t67 + t38) * t41) * t66) * t29, (t46 * t30 - (t39 * t51 * t53 - t38 * t44 + (t38 * t63 - t67) * t28) * t45 * t31) * t29, 0, 0, 0; ((-t44 * t52 - t55 * t62) * t34 - (-t44 * t55 + t52 * t62) * t68) * t33, (-t52 * t34 + t55 * t68) * t45 * t33, t60 * t33, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:42:20
	% EndTime: 2019-12-31 21:42:21
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (428->37), mult. (1217->90), div. (80->9), fcn. (1746->13), ass. (0->53)
	t72 = cos(pkin(5));
	t74 = sin(qJ(2));
	t78 = cos(qJ(1));
	t81 = t78 * t74;
	t75 = sin(qJ(1));
	t77 = cos(qJ(2));
	t82 = t75 * t77;
	t64 = t72 * t81 + t82;
	t73 = sin(qJ(3));
	t76 = cos(qJ(3));
	t70 = sin(pkin(5));
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
	t69 = sin(pkin(10));
	t71 = cos(pkin(10));
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
	t44 = 0.1e1 / (t47 ^ 2 * t46 + 0.1e1);
	t41 = 0.1e1 / t43;
	t40 = 0.1e1 / (0.1e1 + t89);
	t39 = (-t59 * t63 + t85 * t90) * t73 * t51;
	t38 = (-t55 * t59 + t62 * t90) * t51;
	t1 = [-t57 * t59 * t51, t39, t38, 0, 0; (-t53 * t41 - (-t49 + (t59 * t91 + t49) * t51) * t89) * t40, (-t65 * t73 * t41 - ((-t49 * t63 + t50 * t85) * t73 + t79 * t39) * t93) * t40, (t58 * t41 - (t79 * t38 - t49 * t55 + t50 * t62) * t93) * t40, 0, 0; ((-t55 * t69 - t63 * t71) * t45 - (-t55 * t71 + t63 * t69) * t92) * t44, ((-t66 * t71 - t69 * t88) * t45 - (t66 * t69 - t71 * t88) * t92) * t44, (-t69 * t45 + t71 * t92) * t57 * t44, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:42:20
	% EndTime: 2019-12-31 21:42:21
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (525->38), mult. (1300->91), div. (85->9), fcn. (1858->13), ass. (0->55)
	t80 = cos(pkin(5));
	t82 = sin(qJ(2));
	t86 = cos(qJ(1));
	t90 = t86 * t82;
	t83 = sin(qJ(1));
	t85 = cos(qJ(2));
	t91 = t83 * t85;
	t71 = t80 * t90 + t91;
	t81 = sin(qJ(3));
	t84 = cos(qJ(3));
	t79 = sin(pkin(5));
	t93 = t79 * t86;
	t60 = t71 * t81 + t84 * t93;
	t96 = t79 * t81;
	t68 = -t80 * t84 + t82 * t96;
	t59 = atan2(-t60, t68);
	t56 = sin(t59);
	t57 = cos(t59);
	t50 = -t56 * t60 + t57 * t68;
	t49 = 0.1e1 / t50 ^ 2;
	t89 = t86 * t85;
	t92 = t83 * t82;
	t73 = -t80 * t92 + t89;
	t95 = t79 * t84;
	t64 = t73 * t81 - t83 * t95;
	t102 = t49 * t64;
	t65 = t73 * t84 + t83 * t96;
	t72 = t80 * t91 + t90;
	t78 = pkin(10) + qJ(5);
	t76 = sin(t78);
	t77 = cos(t78);
	t55 = t65 * t77 + t72 * t76;
	t53 = 0.1e1 / t55 ^ 2;
	t54 = t65 * t76 - t72 * t77;
	t101 = t53 * t54;
	t100 = t57 * t60;
	t67 = 0.1e1 / t68 ^ 2;
	t99 = t60 * t67;
	t98 = t64 ^ 2 * t49;
	t97 = t72 * t84;
	t94 = t79 * t85;
	t88 = t54 ^ 2 * t53 + 0.1e1;
	t62 = t71 * t84 - t81 * t93;
	t87 = -t56 * t68 - t100;
	t70 = t80 * t89 - t92;
	t69 = t80 * t81 + t82 * t95;
	t66 = 0.1e1 / t68;
	t58 = 0.1e1 / (t60 ^ 2 * t67 + 0.1e1);
	t52 = 0.1e1 / t55;
	t51 = 0.1e1 / t88;
	t48 = 0.1e1 / t50;
	t47 = 0.1e1 / (0.1e1 + t98);
	t46 = (-t66 * t70 + t94 * t99) * t81 * t58;
	t45 = (-t62 * t66 + t69 * t99) * t58;
	t1 = [-t64 * t66 * t58, t46, t45, 0, 0; (-t60 * t48 - (-t56 + (t66 * t100 + t56) * t58) * t98) * t47, (-t72 * t81 * t48 - ((-t56 * t70 + t57 * t94) * t81 + t87 * t46) * t102) * t47, (t65 * t48 - (t87 * t45 - t56 * t62 + t57 * t69) * t102) * t47, 0, 0; ((-t62 * t76 - t70 * t77) * t52 - (-t62 * t77 + t70 * t76) * t101) * t51, ((-t73 * t77 - t76 * t97) * t52 - (t73 * t76 - t77 * t97) * t101) * t51, (t77 * t101 - t76 * t52) * t64 * t51, 0, t88 * t51;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end