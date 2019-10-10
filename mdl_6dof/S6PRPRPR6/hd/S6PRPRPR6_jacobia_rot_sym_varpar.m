% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR6
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
%   Wie in S6PRPRPR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRPR6_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR6_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.09s
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
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.14s
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
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (303->29), mult. (866->79), div. (60->9), fcn. (1237->13), ass. (0->49)
	t61 = sin(pkin(10));
	t64 = cos(pkin(10));
	t67 = sin(qJ(2));
	t65 = cos(pkin(6));
	t69 = cos(qJ(2));
	t71 = t65 * t69;
	t54 = t61 * t67 - t64 * t71;
	t68 = cos(qJ(4));
	t62 = sin(pkin(6));
	t66 = sin(qJ(4));
	t76 = t62 * t66;
	t49 = t54 * t68 + t64 * t76;
	t73 = t62 * t69;
	t58 = t65 * t66 + t68 * t73;
	t46 = atan2(t49, t58);
	t43 = sin(t46);
	t44 = cos(t46);
	t37 = t43 * t49 + t44 * t58;
	t36 = 0.1e1 / t37 ^ 2;
	t56 = t61 * t71 + t64 * t67;
	t47 = -t56 * t68 + t61 * t76;
	t80 = t36 * t47;
	t74 = t62 * t68;
	t48 = t56 * t66 + t61 * t74;
	t72 = t65 * t67;
	t57 = -t61 * t72 + t64 * t69;
	t60 = sin(pkin(11));
	t63 = cos(pkin(11));
	t42 = t48 * t63 + t57 * t60;
	t40 = 0.1e1 / t42 ^ 2;
	t41 = t48 * t60 - t57 * t63;
	t79 = t40 * t41;
	t53 = 0.1e1 / t58 ^ 2;
	t78 = t49 * t53;
	t77 = t57 * t66;
	t75 = t62 * t67;
	t70 = -t43 * t58 + t44 * t49;
	t59 = t65 * t68 - t66 * t73;
	t55 = t61 * t69 + t64 * t72;
	t52 = 0.1e1 / t58;
	t50 = -t54 * t66 + t64 * t74;
	t45 = 0.1e1 / (t49 ^ 2 * t53 + 0.1e1);
	t39 = 0.1e1 / t42;
	t38 = 0.1e1 / (t41 ^ 2 * t40 + 0.1e1);
	t35 = 0.1e1 / t37;
	t34 = 0.1e1 / (t47 ^ 2 * t36 + 0.1e1);
	t33 = (t52 * t55 + t75 * t78) * t68 * t45;
	t32 = (t50 * t52 - t59 * t78) * t45;
	t1 = [0, t33, 0, t32, 0, 0; 0, (-t57 * t68 * t35 - ((t43 * t55 - t44 * t75) * t68 + t70 * t33) * t80) * t34, 0, (t48 * t35 - (t70 * t32 + t43 * t50 + t44 * t59) * t80) * t34, 0, 0; 0, ((t56 * t63 + t60 * t77) * t39 - (-t56 * t60 + t63 * t77) * t79) * t38, 0, (-t39 * t60 + t63 * t79) * t47 * t38, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (382->30), mult. (949->80), div. (65->9), fcn. (1349->13), ass. (0->51)
	t69 = sin(pkin(10));
	t71 = cos(pkin(10));
	t74 = sin(qJ(2));
	t72 = cos(pkin(6));
	t76 = cos(qJ(2));
	t79 = t72 * t76;
	t60 = t69 * t74 - t71 * t79;
	t75 = cos(qJ(4));
	t70 = sin(pkin(6));
	t73 = sin(qJ(4));
	t84 = t70 * t73;
	t55 = t60 * t75 + t71 * t84;
	t81 = t70 * t76;
	t64 = t72 * t73 + t75 * t81;
	t52 = atan2(t55, t64);
	t49 = sin(t52);
	t50 = cos(t52);
	t43 = t49 * t55 + t50 * t64;
	t42 = 0.1e1 / t43 ^ 2;
	t62 = t69 * t79 + t71 * t74;
	t53 = -t62 * t75 + t69 * t84;
	t88 = t42 * t53;
	t82 = t70 * t75;
	t54 = t62 * t73 + t69 * t82;
	t80 = t72 * t74;
	t63 = -t69 * t80 + t71 * t76;
	t68 = pkin(11) + qJ(6);
	t66 = sin(t68);
	t67 = cos(t68);
	t48 = t54 * t67 + t63 * t66;
	t46 = 0.1e1 / t48 ^ 2;
	t47 = t54 * t66 - t63 * t67;
	t87 = t46 * t47;
	t59 = 0.1e1 / t64 ^ 2;
	t86 = t55 * t59;
	t85 = t63 * t73;
	t83 = t70 * t74;
	t78 = t47 ^ 2 * t46 + 0.1e1;
	t77 = -t49 * t64 + t50 * t55;
	t65 = t72 * t75 - t73 * t81;
	t61 = t69 * t76 + t71 * t80;
	t58 = 0.1e1 / t64;
	t56 = -t60 * t73 + t71 * t82;
	t51 = 0.1e1 / (t55 ^ 2 * t59 + 0.1e1);
	t45 = 0.1e1 / t48;
	t44 = 0.1e1 / t78;
	t41 = 0.1e1 / t43;
	t40 = 0.1e1 / (t53 ^ 2 * t42 + 0.1e1);
	t39 = (t58 * t61 + t83 * t86) * t75 * t51;
	t38 = (t56 * t58 - t65 * t86) * t51;
	t1 = [0, t39, 0, t38, 0, 0; 0, (-t63 * t75 * t41 - ((t49 * t61 - t50 * t83) * t75 + t77 * t39) * t88) * t40, 0, (t54 * t41 - (t77 * t38 + t49 * t56 + t50 * t65) * t88) * t40, 0, 0; 0, ((t62 * t67 + t66 * t85) * t45 - (-t62 * t66 + t67 * t85) * t87) * t44, 0, (-t45 * t66 + t67 * t87) * t53 * t44, 0, t78 * t44;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end