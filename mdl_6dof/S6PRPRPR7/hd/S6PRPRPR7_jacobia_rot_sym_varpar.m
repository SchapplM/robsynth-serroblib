% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR7
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
%   Wie in S6PRPRPR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:40
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRPR7_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR7_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (68->16), mult. (173->36), div. (41->11), fcn. (271->9), ass. (0->23)
	t34 = sin(pkin(6));
	t38 = cos(qJ(2));
	t41 = t34 * t38;
	t36 = cos(pkin(6));
	t37 = sin(qJ(2));
	t40 = t36 * t37;
	t39 = t36 * t38;
	t35 = cos(pkin(10));
	t33 = sin(pkin(10));
	t32 = 0.1e1 / t38 ^ 2;
	t31 = 0.1e1 / t34 ^ 2;
	t30 = 0.1e1 / t34;
	t28 = -t33 * t40 + t35 * t38;
	t27 = t33 * t39 + t35 * t37;
	t26 = t33 * t38 + t35 * t40;
	t24 = t33 * t37 - t35 * t39;
	t22 = atan2(-t24, -t41);
	t21 = cos(t22);
	t20 = sin(t22);
	t19 = -t20 * t24 - t21 * t41;
	t18 = 0.1e1 / t19 ^ 2;
	t16 = (t26 / t38 + t37 * t24 * t32) * t30 / (t24 ^ 2 * t31 * t32 + 0.1e1);
	t1 = [0, t16, 0, 0, 0, 0; 0, (t28 / t19 - (t21 * t34 * t37 - t20 * t26 + (t20 * t41 - t21 * t24) * t16) * t27 * t18) / (t27 ^ 2 * t18 + 0.1e1), 0, 0, 0, 0; 0, -t27 / t33 * t30 / (0.1e1 + t28 ^ 2 / t33 ^ 2 * t31), 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (89->17), mult. (268->43), div. (47->11), fcn. (395->11), ass. (0->29)
	t41 = sin(pkin(10));
	t42 = sin(pkin(6));
	t53 = t41 * t42;
	t46 = sin(qJ(2));
	t52 = t42 * t46;
	t44 = cos(pkin(6));
	t51 = t44 * t46;
	t48 = cos(qJ(2));
	t50 = t44 * t48;
	t43 = cos(pkin(10));
	t37 = t41 * t50 + t43 * t46;
	t45 = sin(qJ(4));
	t47 = cos(qJ(4));
	t29 = t37 * t45 + t47 * t53;
	t27 = 0.1e1 / t29 ^ 2;
	t28 = -t37 * t47 + t45 * t53;
	t49 = t28 ^ 2 * t27 + 0.1e1;
	t40 = 0.1e1 / t46 ^ 2;
	t38 = -t41 * t51 + t43 * t48;
	t35 = t41 * t48 + t43 * t51;
	t34 = t41 * t46 - t43 * t50;
	t33 = atan2(-t35, t52);
	t31 = cos(t33);
	t30 = sin(t33);
	t26 = 0.1e1 / t49;
	t25 = -t30 * t35 + t31 * t52;
	t24 = 0.1e1 / t25 ^ 2;
	t22 = (t34 / t46 + t48 * t35 * t40) / t42 / (0.1e1 + t35 ^ 2 / t42 ^ 2 * t40);
	t1 = [0, t22, 0, 0, 0, 0; 0, (-t37 / t25 - (t31 * t42 * t48 + t30 * t34 + (-t30 * t52 - t31 * t35) * t22) * t38 * t24) / (t38 ^ 2 * t24 + 0.1e1), 0, 0, 0, 0; 0, (-t47 / t29 - t45 * t28 * t27) * t38 * t26, 0, t49 * t26, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (251->25), mult. (724->66), div. (56->9), fcn. (1043->11), ass. (0->42)
	t53 = sin(pkin(10));
	t55 = cos(pkin(10));
	t58 = sin(qJ(2));
	t56 = cos(pkin(6));
	t60 = cos(qJ(2));
	t62 = t56 * t60;
	t47 = t53 * t58 - t55 * t62;
	t59 = cos(qJ(4));
	t54 = sin(pkin(6));
	t57 = sin(qJ(4));
	t67 = t54 * t57;
	t42 = t47 * t59 + t55 * t67;
	t64 = t54 * t60;
	t51 = t56 * t57 + t59 * t64;
	t39 = atan2(t42, t51);
	t35 = sin(t39);
	t36 = cos(t39);
	t34 = t35 * t42 + t36 * t51;
	t33 = 0.1e1 / t34 ^ 2;
	t49 = t53 * t62 + t55 * t58;
	t40 = -t49 * t59 + t53 * t67;
	t69 = t33 * t40;
	t46 = 0.1e1 / t51 ^ 2;
	t68 = t42 * t46;
	t66 = t54 * t58;
	t65 = t54 * t59;
	t63 = t56 * t58;
	t61 = -t35 * t51 + t36 * t42;
	t52 = t56 * t59 - t57 * t64;
	t50 = -t53 * t63 + t55 * t60;
	t48 = t53 * t60 + t55 * t63;
	t45 = 0.1e1 / t51;
	t44 = 0.1e1 / t50 ^ 2;
	t43 = -t47 * t57 + t55 * t65;
	t41 = t49 * t57 + t53 * t65;
	t38 = 0.1e1 / (t42 ^ 2 * t46 + 0.1e1);
	t37 = 0.1e1 / (t41 ^ 2 * t44 + 0.1e1);
	t32 = 0.1e1 / t34;
	t31 = 0.1e1 / (t40 ^ 2 * t33 + 0.1e1);
	t30 = (t45 * t48 + t66 * t68) * t59 * t38;
	t29 = (t43 * t45 - t52 * t68) * t38;
	t1 = [0, t30, 0, t29, 0, 0; 0, (-t50 * t59 * t32 - ((t35 * t48 - t36 * t66) * t59 + t61 * t30) * t69) * t31, 0, (t41 * t32 - (t61 * t29 + t35 * t43 + t36 * t52) * t69) * t31, 0, 0; 0, (t41 * t44 * t49 + t57) * t37, 0, -t40 / t50 * t37, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (334->29), mult. (949->80), div. (65->9), fcn. (1349->13), ass. (0->50)
	t63 = cos(pkin(10));
	t66 = sin(qJ(4));
	t61 = sin(pkin(10));
	t67 = sin(qJ(2));
	t64 = cos(pkin(6));
	t70 = cos(qJ(2));
	t74 = t64 * t70;
	t71 = t61 * t67 - t63 * t74;
	t62 = sin(pkin(6));
	t69 = cos(qJ(4));
	t77 = t62 * t69;
	t52 = t63 * t77 - t71 * t66;
	t76 = t62 * t70;
	t60 = t64 * t69 - t66 * t76;
	t48 = atan2(t52, t60);
	t45 = sin(t48);
	t46 = cos(t48);
	t39 = t45 * t52 + t46 * t60;
	t38 = 0.1e1 / t39 ^ 2;
	t57 = t61 * t74 + t63 * t67;
	t50 = t57 * t66 + t61 * t77;
	t83 = t38 * t50;
	t79 = t62 * t66;
	t49 = -t57 * t69 + t61 * t79;
	t75 = t64 * t67;
	t58 = -t61 * t75 + t63 * t70;
	t65 = sin(qJ(6));
	t68 = cos(qJ(6));
	t44 = t49 * t65 + t58 * t68;
	t42 = 0.1e1 / t44 ^ 2;
	t43 = -t49 * t68 + t58 * t65;
	t82 = t42 * t43;
	t55 = 0.1e1 / t60 ^ 2;
	t81 = t52 * t55;
	t80 = t58 * t69;
	t78 = t62 * t67;
	t73 = t43 ^ 2 * t42 + 0.1e1;
	t72 = -t45 * t60 + t46 * t52;
	t59 = -t64 * t66 - t69 * t76;
	t56 = t61 * t70 + t63 * t75;
	t54 = 0.1e1 / t60;
	t51 = t63 * t79 + t71 * t69;
	t47 = 0.1e1 / (t52 ^ 2 * t55 + 0.1e1);
	t41 = 0.1e1 / t44;
	t40 = 0.1e1 / t73;
	t37 = 0.1e1 / t39;
	t36 = 0.1e1 / (t50 ^ 2 * t38 + 0.1e1);
	t35 = (-t54 * t56 - t78 * t81) * t66 * t47;
	t34 = (-t51 * t54 - t59 * t81) * t47;
	t1 = [0, t35, 0, t34, 0, 0; 0, (t58 * t66 * t37 - ((-t45 * t56 + t46 * t78) * t66 + t72 * t35) * t83) * t36, 0, (-t49 * t37 - (t72 * t34 - t45 * t51 + t46 * t59) * t83) * t36, 0, 0; 0, ((-t57 * t65 + t68 * t80) * t41 - (-t57 * t68 - t65 * t80) * t82) * t40, 0, (-t41 * t68 - t65 * t82) * t50 * t40, 0, t73 * t40;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end