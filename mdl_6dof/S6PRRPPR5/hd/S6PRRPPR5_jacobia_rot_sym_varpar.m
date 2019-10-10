% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPPR5
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
%   Wie in S6PRRPPR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:14
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPPR5_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR5_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:14:35
	% EndTime: 2019-10-09 22:14:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:14:35
	% EndTime: 2019-10-09 22:14:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:14:35
	% EndTime: 2019-10-09 22:14:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:14:35
	% EndTime: 2019-10-09 22:14:35
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
	% StartTime: 2019-10-09 22:14:35
	% EndTime: 2019-10-09 22:14:35
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (253->26), mult. (731->68), div. (57->9), fcn. (1053->11), ass. (0->42)
	t53 = sin(pkin(10));
	t55 = cos(pkin(10));
	t60 = cos(qJ(2));
	t56 = cos(pkin(6));
	t58 = sin(qJ(2));
	t63 = t56 * t58;
	t48 = t53 * t60 + t55 * t63;
	t57 = sin(qJ(3));
	t54 = sin(pkin(6));
	t59 = cos(qJ(3));
	t65 = t54 * t59;
	t38 = t48 * t57 + t55 * t65;
	t66 = t54 * t57;
	t51 = -t56 * t59 + t58 * t66;
	t37 = atan2(-t38, t51);
	t33 = sin(t37);
	t34 = cos(t37);
	t32 = -t33 * t38 + t34 * t51;
	t31 = 0.1e1 / t32 ^ 2;
	t50 = -t53 * t63 + t55 * t60;
	t41 = t50 * t57 - t53 * t65;
	t68 = t31 * t41;
	t46 = 0.1e1 / t51 ^ 2;
	t67 = t38 * t46;
	t64 = t54 * t60;
	t62 = t56 * t60;
	t61 = -t33 * t51 - t34 * t38;
	t52 = t56 * t57 + t58 * t65;
	t49 = t53 * t62 + t55 * t58;
	t47 = -t53 * t58 + t55 * t62;
	t45 = 0.1e1 / t51;
	t44 = 0.1e1 / t49 ^ 2;
	t43 = 0.1e1 / t49;
	t42 = t50 * t59 + t53 * t66;
	t40 = t48 * t59 - t55 * t66;
	t36 = 0.1e1 / (t38 ^ 2 * t46 + 0.1e1);
	t35 = 0.1e1 / (t42 ^ 2 * t44 + 0.1e1);
	t30 = 0.1e1 / t32;
	t29 = 0.1e1 / (t41 ^ 2 * t31 + 0.1e1);
	t28 = (-t45 * t47 + t64 * t67) * t57 * t36;
	t27 = (-t40 * t45 + t52 * t67) * t36;
	t1 = [0, t28, t27, 0, 0, 0; 0, (-t49 * t57 * t30 - ((-t33 * t47 + t34 * t64) * t57 + t61 * t28) * t68) * t29, (t42 * t30 - (t61 * t27 - t33 * t40 + t34 * t52) * t68) * t29, 0, 0, 0; 0, (-t42 * t44 * t50 - t43 * t49 * t59) * t35, -t41 * t43 * t35, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:14:35
	% EndTime: 2019-10-09 22:14:35
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (303->30), mult. (866->78), div. (60->9), fcn. (1237->13), ass. (0->48)
	t61 = sin(pkin(10));
	t64 = cos(pkin(10));
	t69 = cos(qJ(2));
	t65 = cos(pkin(6));
	t67 = sin(qJ(2));
	t72 = t65 * t67;
	t54 = t61 * t69 + t64 * t72;
	t68 = cos(qJ(3));
	t62 = sin(pkin(6));
	t66 = sin(qJ(3));
	t75 = t62 * t66;
	t47 = t54 * t68 - t64 * t75;
	t74 = t62 * t68;
	t58 = t65 * t66 + t67 * t74;
	t45 = atan2(-t47, t58);
	t42 = sin(t45);
	t43 = cos(t45);
	t36 = -t42 * t47 + t43 * t58;
	t35 = 0.1e1 / t36 ^ 2;
	t56 = -t61 * t72 + t64 * t69;
	t50 = t56 * t68 + t61 * t75;
	t79 = t35 * t50;
	t49 = t56 * t66 - t61 * t74;
	t71 = t65 * t69;
	t55 = t61 * t71 + t64 * t67;
	t60 = sin(pkin(11));
	t63 = cos(pkin(11));
	t41 = t49 * t60 + t55 * t63;
	t39 = 0.1e1 / t41 ^ 2;
	t40 = -t49 * t63 + t55 * t60;
	t78 = t39 * t40;
	t52 = 0.1e1 / t58 ^ 2;
	t77 = t47 * t52;
	t76 = t55 * t66;
	t73 = t62 * t69;
	t70 = -t42 * t58 - t43 * t47;
	t57 = t65 * t68 - t67 * t75;
	t53 = -t61 * t67 + t64 * t71;
	t51 = 0.1e1 / t58;
	t46 = t54 * t66 + t64 * t74;
	t44 = 0.1e1 / (t47 ^ 2 * t52 + 0.1e1);
	t38 = 0.1e1 / t41;
	t37 = 0.1e1 / (t40 ^ 2 * t39 + 0.1e1);
	t34 = 0.1e1 / t36;
	t33 = 0.1e1 / (t50 ^ 2 * t35 + 0.1e1);
	t32 = (-t51 * t53 + t73 * t77) * t68 * t44;
	t31 = (t46 * t51 + t57 * t77) * t44;
	t1 = [0, t32, t31, 0, 0, 0; 0, (-t55 * t68 * t34 - ((-t42 * t53 + t43 * t73) * t68 + t70 * t32) * t79) * t33, (-t49 * t34 - (t70 * t31 + t42 * t46 + t43 * t57) * t79) * t33, 0, 0, 0; 0, ((t56 * t60 + t63 * t76) * t38 - (t56 * t63 - t60 * t76) * t78) * t37, (-t38 * t63 - t60 * t78) * t50 * t37, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:14:35
	% EndTime: 2019-10-09 22:14:35
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (382->31), mult. (949->79), div. (65->9), fcn. (1349->13), ass. (0->50)
	t71 = sin(pkin(10));
	t73 = cos(pkin(10));
	t78 = cos(qJ(2));
	t74 = cos(pkin(6));
	t76 = sin(qJ(2));
	t82 = t74 * t76;
	t62 = t71 * t78 + t73 * t82;
	t77 = cos(qJ(3));
	t72 = sin(pkin(6));
	t75 = sin(qJ(3));
	t85 = t72 * t75;
	t55 = t62 * t77 - t73 * t85;
	t84 = t72 * t77;
	t66 = t74 * t75 + t76 * t84;
	t53 = atan2(-t55, t66);
	t50 = sin(t53);
	t51 = cos(t53);
	t44 = -t50 * t55 + t51 * t66;
	t43 = 0.1e1 / t44 ^ 2;
	t64 = -t71 * t82 + t73 * t78;
	t58 = t64 * t77 + t71 * t85;
	t89 = t43 * t58;
	t57 = t64 * t75 - t71 * t84;
	t81 = t74 * t78;
	t63 = t71 * t81 + t73 * t76;
	t70 = pkin(11) + qJ(6);
	t68 = sin(t70);
	t69 = cos(t70);
	t49 = t57 * t68 + t63 * t69;
	t47 = 0.1e1 / t49 ^ 2;
	t48 = -t57 * t69 + t63 * t68;
	t88 = t47 * t48;
	t60 = 0.1e1 / t66 ^ 2;
	t87 = t55 * t60;
	t86 = t63 * t75;
	t83 = t72 * t78;
	t80 = t48 ^ 2 * t47 + 0.1e1;
	t79 = -t50 * t66 - t51 * t55;
	t65 = t74 * t77 - t76 * t85;
	t61 = -t71 * t76 + t73 * t81;
	t59 = 0.1e1 / t66;
	t54 = t62 * t75 + t73 * t84;
	t52 = 0.1e1 / (t55 ^ 2 * t60 + 0.1e1);
	t46 = 0.1e1 / t49;
	t45 = 0.1e1 / t80;
	t42 = 0.1e1 / t44;
	t41 = 0.1e1 / (t58 ^ 2 * t43 + 0.1e1);
	t40 = (-t59 * t61 + t83 * t87) * t77 * t52;
	t39 = (t54 * t59 + t65 * t87) * t52;
	t1 = [0, t40, t39, 0, 0, 0; 0, (-t63 * t77 * t42 - ((-t50 * t61 + t51 * t83) * t77 + t79 * t40) * t89) * t41, (-t57 * t42 - (t79 * t39 + t50 * t54 + t51 * t65) * t89) * t41, 0, 0, 0; 0, ((t64 * t68 + t69 * t86) * t46 - (t64 * t69 - t68 * t86) * t88) * t45, (-t46 * t69 - t68 * t88) * t58 * t45, 0, 0, t80 * t45;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end