% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRRP5
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
%   Wie in S6PRPRRP5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:50
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRP5_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRP5_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:50:08
	% EndTime: 2019-10-09 21:50:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:50:08
	% EndTime: 2019-10-09 21:50:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:50:08
	% EndTime: 2019-10-09 21:50:08
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:50:08
	% EndTime: 2019-10-09 21:50:08
	% DurationCPUTime: 0.07s
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
	% StartTime: 2019-10-09 21:50:08
	% EndTime: 2019-10-09 21:50:08
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
	% StartTime: 2019-10-09 21:50:08
	% EndTime: 2019-10-09 21:50:08
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (334->29), mult. (949->79), div. (65->9), fcn. (1349->13), ass. (0->51)
	t61 = sin(pkin(10));
	t63 = cos(pkin(10));
	t67 = sin(qJ(2));
	t64 = cos(pkin(6));
	t70 = cos(qJ(2));
	t73 = t64 * t70;
	t55 = t61 * t67 - t63 * t73;
	t69 = cos(qJ(4));
	t62 = sin(pkin(6));
	t66 = sin(qJ(4));
	t78 = t62 * t66;
	t50 = t55 * t69 + t63 * t78;
	t75 = t62 * t70;
	t59 = t64 * t66 + t69 * t75;
	t47 = atan2(t50, t59);
	t44 = sin(t47);
	t45 = cos(t47);
	t38 = t44 * t50 + t45 * t59;
	t37 = 0.1e1 / t38 ^ 2;
	t57 = t61 * t73 + t63 * t67;
	t48 = -t57 * t69 + t61 * t78;
	t83 = t37 * t48;
	t76 = t62 * t69;
	t49 = t57 * t66 + t61 * t76;
	t68 = cos(qJ(5));
	t74 = t64 * t67;
	t58 = -t61 * t74 + t63 * t70;
	t65 = sin(qJ(5));
	t80 = t58 * t65;
	t43 = t49 * t68 + t80;
	t41 = 0.1e1 / t43 ^ 2;
	t79 = t58 * t68;
	t42 = t49 * t65 - t79;
	t82 = t41 * t42;
	t54 = 0.1e1 / t59 ^ 2;
	t81 = t50 * t54;
	t77 = t62 * t67;
	t72 = t42 ^ 2 * t41 + 0.1e1;
	t71 = -t44 * t59 + t45 * t50;
	t60 = t64 * t69 - t66 * t75;
	t56 = t61 * t70 + t63 * t74;
	t53 = 0.1e1 / t59;
	t51 = -t55 * t66 + t63 * t76;
	t46 = 0.1e1 / (t50 ^ 2 * t54 + 0.1e1);
	t40 = 0.1e1 / t43;
	t39 = 0.1e1 / t72;
	t36 = 0.1e1 / t38;
	t35 = 0.1e1 / (t48 ^ 2 * t37 + 0.1e1);
	t34 = (t53 * t56 + t77 * t81) * t69 * t46;
	t33 = (t51 * t53 - t60 * t81) * t46;
	t1 = [0, t34, 0, t33, 0, 0; 0, (-t58 * t69 * t36 - ((t44 * t56 - t45 * t77) * t69 + t71 * t34) * t83) * t35, 0, (t49 * t36 - (t71 * t33 + t44 * t51 + t45 * t60) * t83) * t35, 0, 0; 0, ((t57 * t68 + t66 * t80) * t40 - (-t57 * t65 + t66 * t79) * t82) * t39, 0, (-t40 * t65 + t68 * t82) * t48 * t39, t72 * t39, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:50:08
	% EndTime: 2019-10-09 21:50:08
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (334->29), mult. (949->79), div. (65->9), fcn. (1349->13), ass. (0->51)
	t70 = sin(pkin(10));
	t72 = cos(pkin(10));
	t76 = sin(qJ(2));
	t73 = cos(pkin(6));
	t79 = cos(qJ(2));
	t82 = t73 * t79;
	t64 = t70 * t76 - t72 * t82;
	t78 = cos(qJ(4));
	t71 = sin(pkin(6));
	t75 = sin(qJ(4));
	t87 = t71 * t75;
	t59 = t64 * t78 + t72 * t87;
	t84 = t71 * t79;
	t68 = t73 * t75 + t78 * t84;
	t56 = atan2(t59, t68);
	t53 = sin(t56);
	t54 = cos(t56);
	t47 = t53 * t59 + t54 * t68;
	t46 = 0.1e1 / t47 ^ 2;
	t66 = t70 * t82 + t72 * t76;
	t57 = -t66 * t78 + t70 * t87;
	t92 = t46 * t57;
	t85 = t71 * t78;
	t58 = t66 * t75 + t70 * t85;
	t77 = cos(qJ(5));
	t83 = t73 * t76;
	t67 = -t70 * t83 + t72 * t79;
	t74 = sin(qJ(5));
	t89 = t67 * t74;
	t52 = t58 * t77 + t89;
	t50 = 0.1e1 / t52 ^ 2;
	t88 = t67 * t77;
	t51 = t58 * t74 - t88;
	t91 = t50 * t51;
	t63 = 0.1e1 / t68 ^ 2;
	t90 = t59 * t63;
	t86 = t71 * t76;
	t81 = t51 ^ 2 * t50 + 0.1e1;
	t80 = -t53 * t68 + t54 * t59;
	t69 = t73 * t78 - t75 * t84;
	t65 = t70 * t79 + t72 * t83;
	t62 = 0.1e1 / t68;
	t60 = -t64 * t75 + t72 * t85;
	t55 = 0.1e1 / (t59 ^ 2 * t63 + 0.1e1);
	t49 = 0.1e1 / t52;
	t48 = 0.1e1 / t81;
	t45 = 0.1e1 / t47;
	t44 = 0.1e1 / (t57 ^ 2 * t46 + 0.1e1);
	t43 = (t62 * t65 + t86 * t90) * t78 * t55;
	t42 = (t60 * t62 - t69 * t90) * t55;
	t1 = [0, t43, 0, t42, 0, 0; 0, (-t67 * t78 * t45 - ((t53 * t65 - t54 * t86) * t78 + t80 * t43) * t92) * t44, 0, (t58 * t45 - (t80 * t42 + t53 * t60 + t54 * t69) * t92) * t44, 0, 0; 0, ((t66 * t77 + t75 * t89) * t49 - (-t66 * t74 + t75 * t88) * t91) * t48, 0, (-t49 * t74 + t77 * t91) * t57 * t48, t81 * t48, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end