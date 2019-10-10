% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRRR6
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
%   Wie in S6PRPRRR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:03
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRR6_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR6_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:06
	% EndTime: 2019-10-09 22:03:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:06
	% EndTime: 2019-10-09 22:03:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:06
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:07
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (68->16), mult. (173->36), div. (41->11), fcn. (271->9), ass. (0->23)
	t34 = sin(pkin(6));
	t38 = cos(qJ(2));
	t41 = t34 * t38;
	t36 = cos(pkin(6));
	t37 = sin(qJ(2));
	t40 = t36 * t37;
	t39 = t36 * t38;
	t35 = cos(pkin(11));
	t33 = sin(pkin(11));
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
	% StartTime: 2019-10-09 22:03:07
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (89->17), mult. (268->43), div. (47->11), fcn. (395->11), ass. (0->29)
	t41 = sin(pkin(11));
	t42 = sin(pkin(6));
	t53 = t41 * t42;
	t46 = sin(qJ(2));
	t52 = t42 * t46;
	t44 = cos(pkin(6));
	t51 = t44 * t46;
	t48 = cos(qJ(2));
	t50 = t44 * t48;
	t43 = cos(pkin(11));
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
	% StartTime: 2019-10-09 22:03:07
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (334->29), mult. (949->79), div. (65->9), fcn. (1349->13), ass. (0->51)
	t61 = sin(pkin(11));
	t63 = cos(pkin(11));
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
	% StartTime: 2019-10-09 22:03:07
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (427->30), mult. (1032->80), div. (70->9), fcn. (1461->13), ass. (0->52)
	t79 = sin(pkin(11));
	t81 = cos(pkin(11));
	t84 = sin(qJ(2));
	t82 = cos(pkin(6));
	t86 = cos(qJ(2));
	t89 = t82 * t86;
	t70 = t79 * t84 - t81 * t89;
	t85 = cos(qJ(4));
	t80 = sin(pkin(6));
	t83 = sin(qJ(4));
	t94 = t80 * t83;
	t65 = t70 * t85 + t81 * t94;
	t91 = t80 * t86;
	t74 = t82 * t83 + t85 * t91;
	t62 = atan2(t65, t74);
	t59 = sin(t62);
	t60 = cos(t62);
	t53 = t59 * t65 + t60 * t74;
	t52 = 0.1e1 / t53 ^ 2;
	t72 = t79 * t89 + t81 * t84;
	t63 = -t72 * t85 + t79 * t94;
	t98 = t52 * t63;
	t92 = t80 * t85;
	t64 = t72 * t83 + t79 * t92;
	t90 = t82 * t84;
	t73 = -t79 * t90 + t81 * t86;
	t78 = qJ(5) + qJ(6);
	t76 = sin(t78);
	t77 = cos(t78);
	t58 = t64 * t77 + t73 * t76;
	t56 = 0.1e1 / t58 ^ 2;
	t57 = t64 * t76 - t73 * t77;
	t97 = t56 * t57;
	t69 = 0.1e1 / t74 ^ 2;
	t96 = t65 * t69;
	t95 = t73 * t83;
	t93 = t80 * t84;
	t88 = t57 ^ 2 * t56 + 0.1e1;
	t87 = -t59 * t74 + t60 * t65;
	t75 = t82 * t85 - t83 * t91;
	t71 = t79 * t86 + t81 * t90;
	t68 = 0.1e1 / t74;
	t66 = -t70 * t83 + t81 * t92;
	t61 = 0.1e1 / (t65 ^ 2 * t69 + 0.1e1);
	t55 = 0.1e1 / t58;
	t54 = 0.1e1 / t88;
	t51 = 0.1e1 / t53;
	t50 = 0.1e1 / (t63 ^ 2 * t52 + 0.1e1);
	t49 = (t68 * t71 + t93 * t96) * t85 * t61;
	t48 = (t66 * t68 - t75 * t96) * t61;
	t47 = t88 * t54;
	t1 = [0, t49, 0, t48, 0, 0; 0, (-t73 * t85 * t51 - ((t59 * t71 - t60 * t93) * t85 + t87 * t49) * t98) * t50, 0, (t64 * t51 - (t87 * t48 + t59 * t66 + t60 * t75) * t98) * t50, 0, 0; 0, ((t72 * t77 + t76 * t95) * t55 - (-t72 * t76 + t77 * t95) * t97) * t54, 0, (-t55 * t76 + t77 * t97) * t63 * t54, t47, t47;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end