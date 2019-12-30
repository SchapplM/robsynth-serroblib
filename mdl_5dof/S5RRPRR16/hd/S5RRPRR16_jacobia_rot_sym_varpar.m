% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRPRR16
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
%   Wie in S5RRPRR16_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 19:30
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRPRR16_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR16_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:30:31
	% EndTime: 2019-12-29 19:30:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:30:29
	% EndTime: 2019-12-29 19:30:29
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:30:35
	% EndTime: 2019-12-29 19:30:35
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-12-29 19:30:29
	% EndTime: 2019-12-29 19:30:30
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (126->20), mult. (316->53), div. (70->11), fcn. (497->9), ass. (0->33)
	t43 = sin(qJ(2));
	t44 = sin(qJ(1));
	t45 = cos(qJ(2));
	t46 = cos(qJ(1));
	t49 = cos(pkin(5));
	t47 = t46 * t49;
	t30 = t44 * t43 - t45 * t47;
	t42 = sin(pkin(5));
	t50 = t42 * t45;
	t27 = atan2(-t30, -t50);
	t26 = cos(t27);
	t54 = t26 * t30;
	t48 = t44 * t49;
	t34 = -t43 * t48 + t46 * t45;
	t36 = 0.1e1 / t42;
	t37 = 0.1e1 / t42 ^ 2;
	t39 = 0.1e1 / t44 ^ 2;
	t53 = 0.1e1 / (t34 ^ 2 * t39 * t37 + 0.1e1) * t36;
	t25 = sin(t27);
	t24 = -t25 * t30 - t26 * t50;
	t23 = 0.1e1 / t24 ^ 2;
	t33 = t46 * t43 + t45 * t48;
	t52 = t33 ^ 2 * t23;
	t40 = 0.1e1 / t45;
	t51 = t36 * t40;
	t41 = 0.1e1 / t45 ^ 2;
	t38 = 0.1e1 / t44;
	t32 = t43 * t47 + t44 * t45;
	t28 = 0.1e1 / (t30 ^ 2 * t37 * t41 + 0.1e1);
	t22 = 0.1e1 / t24;
	t21 = 0.1e1 / (0.1e1 + t52);
	t20 = (t30 * t41 * t43 + t32 * t40) * t36 * t28;
	t1 = [t33 * t28 * t51, t20, 0, 0, 0; (-t30 * t22 - (-t25 + (-t51 * t54 + t25) * t28) * t52) * t21, (t34 * t22 - (t26 * t42 * t43 - t25 * t32 + (t25 * t50 - t54) * t20) * t33 * t23) * t21, 0, 0, 0; (-t34 * t39 * t46 - t32 * t38) * t53, -t33 * t38 * t53, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:30:30
	% EndTime: 2019-12-29 19:30:30
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (149->24), mult. (451->64), div. (72->11), fcn. (673->11), ass. (0->40)
	t54 = sin(qJ(2));
	t57 = cos(qJ(2));
	t58 = cos(qJ(1));
	t55 = sin(qJ(1));
	t62 = cos(pkin(5));
	t60 = t55 * t62;
	t45 = t58 * t54 + t57 * t60;
	t53 = sin(qJ(4));
	t56 = cos(qJ(4));
	t52 = sin(pkin(5));
	t64 = t52 * t55;
	t37 = t45 * t53 + t56 * t64;
	t35 = 0.1e1 / t37 ^ 2;
	t36 = -t45 * t56 + t53 * t64;
	t69 = t35 * t36;
	t59 = t58 * t62;
	t43 = t54 * t59 + t55 * t57;
	t65 = t52 * t54;
	t41 = atan2(-t43, t65);
	t39 = cos(t41);
	t68 = t39 * t43;
	t38 = sin(t41);
	t32 = -t38 * t43 + t39 * t65;
	t31 = 0.1e1 / t32 ^ 2;
	t46 = -t54 * t60 + t58 * t57;
	t67 = t46 ^ 2 * t31;
	t49 = 0.1e1 / t52;
	t50 = 0.1e1 / t54;
	t66 = t49 * t50;
	t63 = t52 * t58;
	t61 = t35 * t36 ^ 2 + 0.1e1;
	t51 = 0.1e1 / t54 ^ 2;
	t42 = t55 * t54 - t57 * t59;
	t40 = 0.1e1 / (0.1e1 + t43 ^ 2 / t52 ^ 2 * t51);
	t34 = 0.1e1 / t37;
	t33 = 0.1e1 / t61;
	t30 = 0.1e1 / t32;
	t29 = 0.1e1 / (0.1e1 + t67);
	t28 = (t43 * t51 * t57 + t42 * t50) * t49 * t40;
	t1 = [-t46 * t40 * t66, t28, 0, 0, 0; (-t43 * t30 - (-t38 + (t66 * t68 + t38) * t40) * t67) * t29, (-t45 * t30 - (t39 * t52 * t57 + t38 * t42 + (-t38 * t65 - t68) * t28) * t46 * t31) * t29, 0, 0, 0; ((t42 * t56 + t53 * t63) * t34 - (-t42 * t53 + t56 * t63) * t69) * t33, (-t34 * t56 - t53 * t69) * t46 * t33, 0, t61 * t33, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:30:30
	% EndTime: 2019-12-29 19:30:30
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (459->36), mult. (1300->90), div. (85->9), fcn. (1858->13), ass. (0->55)
	t72 = cos(pkin(5));
	t79 = cos(qJ(2));
	t80 = cos(qJ(1));
	t85 = t80 * t79;
	t75 = sin(qJ(2));
	t76 = sin(qJ(1));
	t88 = t76 * t75;
	t67 = -t72 * t85 + t88;
	t74 = sin(qJ(4));
	t78 = cos(qJ(4));
	t71 = sin(pkin(5));
	t89 = t71 * t80;
	t59 = t67 * t78 + t74 * t89;
	t90 = t71 * t79;
	t65 = t72 * t74 + t78 * t90;
	t56 = atan2(t59, t65);
	t53 = sin(t56);
	t54 = cos(t56);
	t47 = t53 * t59 + t54 * t65;
	t46 = 0.1e1 / t47 ^ 2;
	t86 = t80 * t75;
	t87 = t76 * t79;
	t81 = t72 * t87 + t86;
	t91 = t71 * t76;
	t57 = t74 * t91 - t81 * t78;
	t99 = t46 * t57;
	t58 = t81 * t74 + t78 * t91;
	t77 = cos(qJ(5));
	t69 = -t72 * t88 + t85;
	t73 = sin(qJ(5));
	t94 = t69 * t73;
	t52 = t58 * t77 + t94;
	t50 = 0.1e1 / t52 ^ 2;
	t93 = t69 * t77;
	t51 = t58 * t73 - t93;
	t98 = t50 * t51;
	t97 = t54 * t59;
	t96 = t57 ^ 2 * t46;
	t64 = 0.1e1 / t65 ^ 2;
	t95 = t59 * t64;
	t92 = t71 * t75;
	t84 = t51 ^ 2 * t50 + 0.1e1;
	t83 = -t53 * t65 + t97;
	t82 = -t67 * t74 + t78 * t89;
	t68 = t72 * t86 + t87;
	t66 = t72 * t78 - t74 * t90;
	t63 = 0.1e1 / t65;
	t55 = 0.1e1 / (t59 ^ 2 * t64 + 0.1e1);
	t49 = 0.1e1 / t52;
	t48 = 0.1e1 / t84;
	t45 = 0.1e1 / t47;
	t44 = 0.1e1 / (0.1e1 + t96);
	t43 = (t63 * t68 + t92 * t95) * t78 * t55;
	t42 = (t63 * t82 - t66 * t95) * t55;
	t1 = [-t57 * t63 * t55, t43, 0, t42, 0; (t59 * t45 - (-t53 + (-t63 * t97 + t53) * t55) * t96) * t44, (-t69 * t78 * t45 - ((t53 * t68 - t54 * t92) * t78 + t83 * t43) * t99) * t44, 0, (t58 * t45 - (t83 * t42 + t53 * t82 + t54 * t66) * t99) * t44, 0; ((t68 * t77 + t73 * t82) * t49 - (-t68 * t73 + t77 * t82) * t98) * t48, ((t74 * t94 + t81 * t77) * t49 - (-t81 * t73 + t74 * t93) * t98) * t48, 0, (-t73 * t49 + t77 * t98) * t57 * t48, t84 * t48;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end