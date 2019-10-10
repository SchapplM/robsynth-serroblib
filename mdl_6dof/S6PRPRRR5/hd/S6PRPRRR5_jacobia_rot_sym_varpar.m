% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRRR5
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
%   Wie in S6PRPRRR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:01
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRR5_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR5_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
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
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
	% DurationCPUTime: 0.19s
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
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (150->18), mult. (316->43), div. (52->11), fcn. (458->11), ass. (0->31)
	t58 = sin(pkin(11));
	t59 = sin(pkin(6));
	t68 = t58 * t59;
	t62 = sin(qJ(2));
	t67 = t59 * t62;
	t61 = cos(pkin(6));
	t66 = t61 * t62;
	t63 = cos(qJ(2));
	t65 = t61 * t63;
	t60 = cos(pkin(11));
	t51 = t58 * t65 + t60 * t62;
	t57 = qJ(4) + qJ(5);
	t53 = sin(t57);
	t54 = cos(t57);
	t43 = t51 * t53 + t54 * t68;
	t41 = 0.1e1 / t43 ^ 2;
	t42 = -t51 * t54 + t53 * t68;
	t64 = t42 ^ 2 * t41 + 0.1e1;
	t56 = 0.1e1 / t62 ^ 2;
	t52 = -t58 * t66 + t60 * t63;
	t49 = t58 * t63 + t60 * t66;
	t48 = t58 * t62 - t60 * t65;
	t47 = atan2(-t49, t67);
	t45 = cos(t47);
	t44 = sin(t47);
	t40 = 0.1e1 / t64;
	t39 = -t44 * t49 + t45 * t67;
	t38 = 0.1e1 / t39 ^ 2;
	t36 = (t48 / t62 + t63 * t49 * t56) / t59 / (0.1e1 + t49 ^ 2 / t59 ^ 2 * t56);
	t35 = t64 * t40;
	t1 = [0, t36, 0, 0, 0, 0; 0, (-t51 / t39 - (t45 * t59 * t63 + t44 * t48 + (-t44 * t67 - t45 * t49) * t36) * t52 * t38) / (t52 ^ 2 * t38 + 0.1e1), 0, 0, 0, 0; 0, (-t54 / t43 - t53 * t42 * t41) * t52 * t40, 0, t35, t35, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:17
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (948->30), mult. (1387->79), div. (95->9), fcn. (1976->13), ass. (0->54)
	t87 = sin(pkin(6));
	t88 = cos(pkin(11));
	t100 = t87 * t88;
	t86 = sin(pkin(11));
	t91 = sin(qJ(2));
	t89 = cos(pkin(6));
	t93 = cos(qJ(2));
	t96 = t89 * t93;
	t79 = t86 * t91 - t88 * t96;
	t85 = qJ(4) + qJ(5);
	t83 = sin(t85);
	t84 = cos(t85);
	t72 = t83 * t100 + t79 * t84;
	t98 = t87 * t93;
	t77 = t89 * t83 + t84 * t98;
	t69 = atan2(t72, t77);
	t66 = sin(t69);
	t67 = cos(t69);
	t60 = t66 * t72 + t67 * t77;
	t59 = 0.1e1 / t60 ^ 2;
	t101 = t86 * t87;
	t81 = t86 * t96 + t88 * t91;
	t70 = t83 * t101 - t81 * t84;
	t106 = t59 * t70;
	t97 = t89 * t91;
	t82 = -t86 * t97 + t88 * t93;
	t90 = sin(qJ(6));
	t103 = t82 * t90;
	t71 = t84 * t101 + t81 * t83;
	t92 = cos(qJ(6));
	t65 = t71 * t92 + t103;
	t63 = 0.1e1 / t65 ^ 2;
	t102 = t82 * t92;
	t64 = t71 * t90 - t102;
	t105 = t63 * t64;
	t76 = 0.1e1 / t77 ^ 2;
	t104 = t72 * t76;
	t99 = t87 * t91;
	t95 = t64 ^ 2 * t63 + 0.1e1;
	t94 = -t66 * t77 + t67 * t72;
	t80 = t86 * t93 + t88 * t97;
	t78 = -t83 * t98 + t89 * t84;
	t75 = 0.1e1 / t77;
	t73 = t84 * t100 - t79 * t83;
	t68 = 0.1e1 / (t72 ^ 2 * t76 + 0.1e1);
	t62 = 0.1e1 / t65;
	t61 = 0.1e1 / t95;
	t58 = 0.1e1 / t60;
	t57 = 0.1e1 / (t70 ^ 2 * t59 + 0.1e1);
	t56 = (t99 * t104 + t75 * t80) * t84 * t68;
	t55 = (-t78 * t104 + t73 * t75) * t68;
	t54 = (t92 * t105 - t62 * t90) * t70 * t61;
	t53 = (t71 * t58 - (t94 * t55 + t66 * t73 + t67 * t78) * t106) * t57;
	t1 = [0, t56, 0, t55, t55, 0; 0, (-t82 * t84 * t58 - ((t66 * t80 - t67 * t99) * t84 + t94 * t56) * t106) * t57, 0, t53, t53, 0; 0, ((t83 * t103 + t81 * t92) * t62 - (t83 * t102 - t81 * t90) * t105) * t61, 0, t54, t54, t95 * t61;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end