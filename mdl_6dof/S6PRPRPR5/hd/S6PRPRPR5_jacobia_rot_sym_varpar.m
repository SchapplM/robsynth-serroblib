% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR5
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
%   Wie in S6PRPRPR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:37
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRPR5_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR5_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (84->18), mult. (220->41), div. (42->11), fcn. (332->11), ass. (0->27)
	t39 = sin(pkin(10));
	t40 = sin(pkin(6));
	t49 = t39 * t40;
	t45 = cos(qJ(2));
	t48 = t40 * t45;
	t43 = cos(pkin(6));
	t44 = sin(qJ(2));
	t47 = t43 * t44;
	t46 = t43 * t45;
	t42 = cos(pkin(10));
	t41 = cos(pkin(11));
	t38 = sin(pkin(11));
	t37 = 0.1e1 / t45 ^ 2;
	t34 = -t39 * t47 + t42 * t45;
	t33 = t39 * t46 + t42 * t44;
	t32 = t39 * t45 + t42 * t47;
	t30 = t39 * t44 - t42 * t46;
	t28 = atan2(-t30, -t48);
	t27 = cos(t28);
	t26 = sin(t28);
	t25 = t34 * t41 + t38 * t49;
	t24 = t34 * t38 - t41 * t49;
	t23 = 0.1e1 / t25 ^ 2;
	t21 = -t26 * t30 - t27 * t48;
	t20 = 0.1e1 / t21 ^ 2;
	t18 = (t32 / t45 + t44 * t30 * t37) / t40 / (0.1e1 + t30 ^ 2 / t40 ^ 2 * t37);
	t1 = [0, t18, 0, 0, 0, 0; 0, (t34 / t21 - (t27 * t40 * t44 - t26 * t32 + (t26 * t48 - t27 * t30) * t18) * t33 * t20) / (t33 ^ 2 * t20 + 0.1e1), 0, 0, 0, 0; 0, (-t38 / t25 + t41 * t24 * t23) * t33 / (t24 ^ 2 * t23 + 0.1e1), 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (131->19), mult. (268->43), div. (47->11), fcn. (395->11), ass. (0->30)
	t47 = sin(pkin(10));
	t48 = sin(pkin(6));
	t57 = t47 * t48;
	t52 = cos(qJ(2));
	t56 = t48 * t52;
	t50 = cos(pkin(6));
	t51 = sin(qJ(2));
	t55 = t50 * t51;
	t54 = t50 * t52;
	t49 = cos(pkin(10));
	t40 = -t47 * t55 + t49 * t52;
	t45 = pkin(11) + qJ(4);
	t42 = sin(t45);
	t43 = cos(t45);
	t31 = t40 * t43 + t42 * t57;
	t29 = 0.1e1 / t31 ^ 2;
	t30 = t40 * t42 - t43 * t57;
	t53 = t30 ^ 2 * t29 + 0.1e1;
	t46 = 0.1e1 / t52 ^ 2;
	t39 = t47 * t54 + t49 * t51;
	t38 = t47 * t52 + t49 * t55;
	t36 = t47 * t51 - t49 * t54;
	t34 = atan2(-t36, -t56);
	t33 = cos(t34);
	t32 = sin(t34);
	t28 = 0.1e1 / t53;
	t27 = -t32 * t36 - t33 * t56;
	t26 = 0.1e1 / t27 ^ 2;
	t24 = (t38 / t52 + t51 * t36 * t46) / t48 / (0.1e1 + t36 ^ 2 / t48 ^ 2 * t46);
	t1 = [0, t24, 0, 0, 0, 0; 0, (t40 / t27 - (t33 * t48 * t51 - t32 * t38 + (t32 * t56 - t33 * t36) * t24) * t39 * t26) / (t39 ^ 2 * t26 + 0.1e1), 0, 0, 0, 0; 0, (-t42 / t31 + t43 * t30 * t29) * t39 * t28, 0, t53 * t28, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (523->27), mult. (731->69), div. (57->9), fcn. (1053->11), ass. (0->44)
	t62 = sin(pkin(10));
	t64 = cos(pkin(10));
	t67 = cos(qJ(2));
	t65 = cos(pkin(6));
	t66 = sin(qJ(2));
	t70 = t65 * t66;
	t56 = t62 * t67 + t64 * t70;
	t61 = pkin(11) + qJ(4);
	t59 = sin(t61);
	t60 = cos(t61);
	t63 = sin(pkin(6));
	t73 = t63 * t64;
	t44 = t56 * t59 + t60 * t73;
	t72 = t63 * t66;
	t51 = t59 * t72 - t65 * t60;
	t42 = atan2(-t44, t51);
	t39 = sin(t42);
	t40 = cos(t42);
	t38 = -t39 * t44 + t40 * t51;
	t37 = 0.1e1 / t38 ^ 2;
	t58 = -t62 * t70 + t64 * t67;
	t74 = t62 * t63;
	t47 = t58 * t59 - t60 * t74;
	t76 = t37 * t47;
	t50 = 0.1e1 / t51 ^ 2;
	t75 = t44 * t50;
	t71 = t63 * t67;
	t69 = t65 * t67;
	t68 = -t39 * t51 - t40 * t44;
	t57 = t62 * t69 + t64 * t66;
	t55 = -t62 * t66 + t64 * t69;
	t54 = 0.1e1 / t57 ^ 2;
	t53 = 0.1e1 / t57;
	t52 = t65 * t59 + t60 * t72;
	t49 = 0.1e1 / t51;
	t48 = t58 * t60 + t59 * t74;
	t46 = t56 * t60 - t59 * t73;
	t43 = 0.1e1 / (t48 ^ 2 * t54 + 0.1e1);
	t41 = 0.1e1 / (t44 ^ 2 * t50 + 0.1e1);
	t36 = 0.1e1 / t38;
	t35 = 0.1e1 / (t47 ^ 2 * t37 + 0.1e1);
	t34 = (-t49 * t55 + t71 * t75) * t59 * t41;
	t33 = (-t46 * t49 + t52 * t75) * t41;
	t1 = [0, t34, 0, t33, 0, 0; 0, (-t57 * t59 * t36 - ((-t39 * t55 + t40 * t71) * t59 + t68 * t34) * t76) * t35, 0, (t48 * t36 - (t68 * t33 - t39 * t46 + t40 * t52) * t76) * t35, 0, 0; 0, (-t48 * t54 * t58 - t53 * t57 * t60) * t43, 0, -t47 * t53 * t43, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (641->31), mult. (949->79), div. (65->9), fcn. (1349->13), ass. (0->52)
	t70 = sin(pkin(10));
	t72 = cos(pkin(10));
	t77 = cos(qJ(2));
	t73 = cos(pkin(6));
	t75 = sin(qJ(2));
	t81 = t73 * t75;
	t63 = t70 * t77 + t72 * t81;
	t69 = pkin(11) + qJ(4);
	t67 = sin(t69);
	t68 = cos(t69);
	t71 = sin(pkin(6));
	t84 = t71 * t72;
	t54 = t63 * t68 - t67 * t84;
	t83 = t71 * t75;
	t61 = t73 * t67 + t68 * t83;
	t52 = atan2(-t54, t61);
	t47 = sin(t52);
	t48 = cos(t52);
	t43 = -t47 * t54 + t48 * t61;
	t42 = 0.1e1 / t43 ^ 2;
	t65 = -t70 * t81 + t72 * t77;
	t85 = t70 * t71;
	t57 = t65 * t68 + t67 * t85;
	t90 = t42 * t57;
	t56 = t65 * t67 - t68 * t85;
	t74 = sin(qJ(6));
	t80 = t73 * t77;
	t64 = t70 * t80 + t72 * t75;
	t76 = cos(qJ(6));
	t86 = t64 * t76;
	t50 = t56 * t74 + t86;
	t46 = 0.1e1 / t50 ^ 2;
	t87 = t64 * t74;
	t49 = -t56 * t76 + t87;
	t89 = t46 * t49;
	t59 = 0.1e1 / t61 ^ 2;
	t88 = t54 * t59;
	t82 = t71 * t77;
	t79 = t49 ^ 2 * t46 + 0.1e1;
	t78 = -t47 * t61 - t48 * t54;
	t62 = -t70 * t75 + t72 * t80;
	t60 = -t67 * t83 + t73 * t68;
	t58 = 0.1e1 / t61;
	t53 = t63 * t67 + t68 * t84;
	t51 = 0.1e1 / (t54 ^ 2 * t59 + 0.1e1);
	t45 = 0.1e1 / t50;
	t44 = 0.1e1 / t79;
	t41 = 0.1e1 / t43;
	t40 = 0.1e1 / (t57 ^ 2 * t42 + 0.1e1);
	t39 = (-t58 * t62 + t82 * t88) * t68 * t51;
	t38 = (t53 * t58 + t60 * t88) * t51;
	t1 = [0, t39, 0, t38, 0, 0; 0, (-t64 * t68 * t41 - ((-t47 * t62 + t48 * t82) * t68 + t78 * t39) * t90) * t40, 0, (-t56 * t41 - (t78 * t38 + t47 * t53 + t48 * t60) * t90) * t40, 0, 0; 0, ((t65 * t74 + t67 * t86) * t45 - (t65 * t76 - t67 * t87) * t89) * t44, 0, (-t45 * t76 - t74 * t89) * t57 * t44, 0, t79 * t44;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end