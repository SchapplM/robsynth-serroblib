% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPPR3
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
%   Wie in S6PRRPPR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:10
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPPR3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (101->18), mult. (268->43), div. (47->11), fcn. (395->11), ass. (0->29)
	t41 = sin(pkin(10));
	t42 = sin(pkin(6));
	t53 = t41 * t42;
	t48 = cos(qJ(2));
	t52 = t42 * t48;
	t44 = cos(pkin(6));
	t46 = sin(qJ(2));
	t51 = t44 * t46;
	t50 = t44 * t48;
	t43 = cos(pkin(10));
	t37 = -t41 * t51 + t43 * t48;
	t45 = sin(qJ(3));
	t47 = cos(qJ(3));
	t28 = t37 * t47 + t45 * t53;
	t26 = 0.1e1 / t28 ^ 2;
	t27 = t37 * t45 - t47 * t53;
	t49 = t27 ^ 2 * t26 + 0.1e1;
	t40 = 0.1e1 / t48 ^ 2;
	t36 = t41 * t50 + t43 * t46;
	t35 = t41 * t48 + t43 * t51;
	t33 = t41 * t46 - t43 * t50;
	t31 = atan2(-t33, -t52);
	t30 = cos(t31);
	t29 = sin(t31);
	t25 = 0.1e1 / t49;
	t24 = -t29 * t33 - t30 * t52;
	t23 = 0.1e1 / t24 ^ 2;
	t21 = (t35 / t48 + t46 * t33 * t40) / t42 / (0.1e1 + t33 ^ 2 / t42 ^ 2 * t40);
	t1 = [0, t21, 0, 0, 0, 0; 0, (t37 / t24 - (t30 * t42 * t46 - t29 * t35 + (t29 * t52 - t30 * t33) * t21) * t36 * t23) / (t36 ^ 2 * t23 + 0.1e1), 0, 0, 0, 0; 0, (-t45 / t28 + t47 * t27 * t26) * t36 * t25, t49 * t25, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (256->25), mult. (739->65), div. (57->9), fcn. (1061->11), ass. (0->42)
	t55 = sin(pkin(10));
	t57 = cos(pkin(10));
	t62 = cos(qJ(2));
	t58 = cos(pkin(6));
	t60 = sin(qJ(2));
	t65 = t58 * t60;
	t50 = t55 * t62 + t57 * t65;
	t59 = sin(qJ(3));
	t56 = sin(pkin(6));
	t61 = cos(qJ(3));
	t67 = t56 * t61;
	t41 = t50 * t59 + t57 * t67;
	t68 = t56 * t59;
	t53 = -t58 * t61 + t60 * t68;
	t39 = atan2(-t41, t53);
	t35 = sin(t39);
	t36 = cos(t39);
	t34 = -t35 * t41 + t36 * t53;
	t33 = 0.1e1 / t34 ^ 2;
	t52 = -t55 * t65 + t57 * t62;
	t44 = t52 * t59 - t55 * t67;
	t71 = t33 * t44;
	t48 = 0.1e1 / t53 ^ 2;
	t70 = t41 * t48;
	t45 = t52 * t61 + t55 * t68;
	t40 = 0.1e1 / t45 ^ 2;
	t64 = t58 * t62;
	t51 = -t55 * t64 - t57 * t60;
	t69 = t51 ^ 2 * t40;
	t66 = t56 * t62;
	t63 = -t35 * t53 - t36 * t41;
	t54 = t58 * t59 + t60 * t67;
	t49 = -t55 * t60 + t57 * t64;
	t47 = 0.1e1 / t53;
	t43 = t50 * t61 - t57 * t68;
	t38 = 0.1e1 / (t41 ^ 2 * t48 + 0.1e1);
	t37 = 0.1e1 / (0.1e1 + t69);
	t32 = 0.1e1 / t34;
	t31 = 0.1e1 / (t44 ^ 2 * t33 + 0.1e1);
	t30 = (-t47 * t49 + t66 * t70) * t59 * t38;
	t29 = (-t43 * t47 + t54 * t70) * t38;
	t1 = [0, t30, t29, 0, 0, 0; 0, (t51 * t59 * t32 - ((-t35 * t49 + t36 * t66) * t59 + t63 * t30) * t71) * t31, (t45 * t32 - (t63 * t29 - t35 * t43 + t36 * t54) * t71) * t31, 0, 0, 0; 0, (-t52 / t45 - t61 * t69) * t37, t44 * t51 * t40 * t37, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (91->17), mult. (274->44), div. (48->11), fcn. (404->11), ass. (0->30)
	t51 = sin(pkin(10));
	t53 = cos(pkin(10));
	t58 = cos(qJ(2));
	t54 = cos(pkin(6));
	t56 = sin(qJ(2));
	t60 = t54 * t56;
	t47 = -t51 * t60 + t53 * t58;
	t55 = sin(qJ(3));
	t57 = cos(qJ(3));
	t52 = sin(pkin(6));
	t62 = t51 * t52;
	t38 = t47 * t55 - t57 * t62;
	t36 = 0.1e1 / t38 ^ 2;
	t39 = t47 * t57 + t55 * t62;
	t63 = t39 ^ 2 * t36;
	t61 = t52 * t58;
	t59 = t54 * t58;
	t50 = 0.1e1 / t58 ^ 2;
	t46 = -t51 * t59 - t53 * t56;
	t45 = t51 * t58 + t53 * t60;
	t44 = t51 * t56 - t53 * t59;
	t43 = atan2(t44, t61);
	t41 = cos(t43);
	t40 = sin(t43);
	t35 = 0.1e1 / t38;
	t34 = 0.1e1 / (0.1e1 + t63);
	t33 = t40 * t44 + t41 * t61;
	t32 = 0.1e1 / t33 ^ 2;
	t30 = (t45 / t58 + t56 * t44 * t50) / t52 / (0.1e1 + t44 ^ 2 / t52 ^ 2 * t50);
	t1 = [0, t30, 0, 0, 0, 0; 0, (-t47 / t33 - (-t41 * t52 * t56 + t40 * t45 + (-t40 * t61 + t41 * t44) * t30) * t46 * t32) / (t46 ^ 2 * t32 + 0.1e1), 0, 0, 0, 0; 0, (-t36 * t39 * t55 + t35 * t57) * t46 * t34, (-t35 * t38 - t63) * t34, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (334->30), mult. (949->78), div. (65->9), fcn. (1349->13), ass. (0->50)
	t63 = sin(pkin(10));
	t65 = cos(pkin(10));
	t72 = cos(qJ(2));
	t66 = cos(pkin(6));
	t69 = sin(qJ(2));
	t76 = t66 * t69;
	t56 = t63 * t72 + t65 * t76;
	t71 = cos(qJ(3));
	t64 = sin(pkin(6));
	t68 = sin(qJ(3));
	t79 = t64 * t68;
	t49 = t56 * t71 - t65 * t79;
	t78 = t64 * t71;
	t60 = t66 * t68 + t69 * t78;
	t47 = atan2(-t49, t60);
	t44 = sin(t47);
	t45 = cos(t47);
	t38 = -t44 * t49 + t45 * t60;
	t37 = 0.1e1 / t38 ^ 2;
	t58 = -t63 * t76 + t65 * t72;
	t52 = t58 * t71 + t63 * t79;
	t84 = t37 * t52;
	t51 = t58 * t68 - t63 * t78;
	t70 = cos(qJ(6));
	t75 = t66 * t72;
	t57 = -t63 * t75 - t65 * t69;
	t67 = sin(qJ(6));
	t81 = t57 * t67;
	t43 = t51 * t70 + t81;
	t41 = 0.1e1 / t43 ^ 2;
	t80 = t57 * t70;
	t42 = t51 * t67 - t80;
	t83 = t41 * t42;
	t54 = 0.1e1 / t60 ^ 2;
	t82 = t49 * t54;
	t77 = t64 * t72;
	t74 = t42 ^ 2 * t41 + 0.1e1;
	t73 = -t44 * t60 - t45 * t49;
	t59 = t66 * t71 - t69 * t79;
	t55 = -t63 * t69 + t65 * t75;
	t53 = 0.1e1 / t60;
	t48 = t56 * t68 + t65 * t78;
	t46 = 0.1e1 / (t49 ^ 2 * t54 + 0.1e1);
	t40 = 0.1e1 / t43;
	t39 = 0.1e1 / t74;
	t36 = 0.1e1 / t38;
	t35 = 0.1e1 / (t52 ^ 2 * t37 + 0.1e1);
	t34 = (-t53 * t55 + t77 * t82) * t71 * t46;
	t33 = (t48 * t53 + t59 * t82) * t46;
	t1 = [0, t34, t33, 0, 0, 0; 0, (t57 * t71 * t36 - ((-t44 * t55 + t45 * t77) * t71 + t73 * t34) * t84) * t35, (-t51 * t36 - (t73 * t33 + t44 * t48 + t45 * t59) * t84) * t35, 0, 0, 0; 0, ((t58 * t70 + t68 * t81) * t40 - (-t58 * t67 + t68 * t80) * t83) * t39, (t40 * t67 - t70 * t83) * t52 * t39, 0, 0, t74 * t39;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end