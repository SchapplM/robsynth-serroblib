% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR13
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
%   Wie in S6RRPRPR13_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:26
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR13_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR13_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:26:24
	% EndTime: 2019-10-10 10:26:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:26:24
	% EndTime: 2019-10-10 10:26:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:26:24
	% EndTime: 2019-10-10 10:26:24
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (35->13), mult. (116->34), div. (25->9), fcn. (175->9), ass. (0->25)
	t30 = cos(pkin(6));
	t29 = sin(pkin(6));
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
	t1 = [-t32 * t29 * t41, 0, 0, 0, 0, 0; (0.1e1 / t18 * t38 - (-t24 * t27 * t34 * t41 + (t25 - 0.1e1) * t29 * t23) * t29 * t42) / (t27 * t42 + 0.1e1), 0, 0, 0, 0, 0; ((t30 * t36 - t40) / t22 - (-t30 * t37 - t39) * t21 * t20) * t19, t35 * t19, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:26:24
	% EndTime: 2019-10-10 10:26:24
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (126->20), mult. (316->53), div. (70->11), fcn. (497->9), ass. (0->33)
	t43 = sin(qJ(2));
	t44 = sin(qJ(1));
	t45 = cos(qJ(2));
	t46 = cos(qJ(1));
	t49 = cos(pkin(6));
	t47 = t46 * t49;
	t30 = t44 * t43 - t45 * t47;
	t42 = sin(pkin(6));
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
	t1 = [t33 * t28 * t51, t20, 0, 0, 0, 0; (-t30 * t22 - (-t25 + (-t51 * t54 + t25) * t28) * t52) * t21, (t34 * t22 - (t26 * t42 * t43 - t25 * t32 + (t25 * t50 - t54) * t20) * t33 * t23) * t21, 0, 0, 0, 0; (-t34 * t39 * t46 - t32 * t38) * t53, -t33 * t38 * t53, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:26:24
	% EndTime: 2019-10-10 10:26:24
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (149->24), mult. (451->64), div. (72->11), fcn. (673->11), ass. (0->40)
	t54 = sin(qJ(2));
	t57 = cos(qJ(2));
	t58 = cos(qJ(1));
	t55 = sin(qJ(1));
	t62 = cos(pkin(6));
	t60 = t55 * t62;
	t45 = t58 * t54 + t57 * t60;
	t53 = sin(qJ(4));
	t56 = cos(qJ(4));
	t52 = sin(pkin(6));
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
	t61 = t36 ^ 2 * t35 + 0.1e1;
	t51 = 0.1e1 / t54 ^ 2;
	t42 = t55 * t54 - t57 * t59;
	t40 = 0.1e1 / (0.1e1 + t43 ^ 2 / t52 ^ 2 * t51);
	t34 = 0.1e1 / t37;
	t33 = 0.1e1 / t61;
	t30 = 0.1e1 / t32;
	t29 = 0.1e1 / (0.1e1 + t67);
	t28 = (t43 * t51 * t57 + t42 * t50) * t49 * t40;
	t1 = [-t46 * t40 * t66, t28, 0, 0, 0, 0; (-t43 * t30 - (-t38 + (t66 * t68 + t38) * t40) * t67) * t29, (-t45 * t30 - (t39 * t52 * t57 + t38 * t42 + (-t38 * t65 - t68) * t28) * t46 * t31) * t29, 0, 0, 0, 0; ((t42 * t56 + t53 * t63) * t34 - (-t42 * t53 + t56 * t63) * t69) * t33, (-t56 * t34 - t53 * t69) * t46 * t33, 0, t61 * t33, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:26:24
	% EndTime: 2019-10-10 10:26:24
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (428->36), mult. (1217->90), div. (80->9), fcn. (1746->13), ass. (0->53)
	t69 = cos(pkin(6));
	t74 = cos(qJ(2));
	t75 = cos(qJ(1));
	t79 = t75 * t74;
	t71 = sin(qJ(2));
	t72 = sin(qJ(1));
	t82 = t72 * t71;
	t62 = -t69 * t79 + t82;
	t70 = sin(qJ(4));
	t73 = cos(qJ(4));
	t67 = sin(pkin(6));
	t83 = t67 * t75;
	t54 = t62 * t73 + t70 * t83;
	t84 = t67 * t74;
	t60 = t69 * t70 + t73 * t84;
	t51 = atan2(t54, t60);
	t48 = sin(t51);
	t49 = cos(t51);
	t42 = t48 * t54 + t49 * t60;
	t41 = 0.1e1 / t42 ^ 2;
	t80 = t75 * t71;
	t81 = t72 * t74;
	t76 = t69 * t81 + t80;
	t85 = t67 * t72;
	t52 = t70 * t85 - t73 * t76;
	t92 = t41 * t52;
	t53 = t70 * t76 + t73 * t85;
	t64 = -t69 * t82 + t79;
	t66 = sin(pkin(11));
	t68 = cos(pkin(11));
	t47 = t53 * t68 + t64 * t66;
	t45 = 0.1e1 / t47 ^ 2;
	t46 = t53 * t66 - t64 * t68;
	t91 = t45 * t46;
	t90 = t49 * t54;
	t89 = t52 ^ 2 * t41;
	t59 = 0.1e1 / t60 ^ 2;
	t88 = t54 * t59;
	t87 = t64 * t70;
	t86 = t67 * t71;
	t78 = -t48 * t60 + t90;
	t77 = -t62 * t70 + t73 * t83;
	t63 = t69 * t80 + t81;
	t61 = t69 * t73 - t70 * t84;
	t58 = 0.1e1 / t60;
	t50 = 0.1e1 / (t54 ^ 2 * t59 + 0.1e1);
	t44 = 0.1e1 / t47;
	t43 = 0.1e1 / (t45 * t46 ^ 2 + 0.1e1);
	t40 = 0.1e1 / t42;
	t39 = 0.1e1 / (0.1e1 + t89);
	t38 = (t58 * t63 + t86 * t88) * t73 * t50;
	t37 = (t58 * t77 - t61 * t88) * t50;
	t1 = [-t52 * t58 * t50, t38, 0, t37, 0, 0; (t54 * t40 - (-t48 + (-t58 * t90 + t48) * t50) * t89) * t39, (-t64 * t73 * t40 - ((t48 * t63 - t49 * t86) * t73 + t78 * t38) * t92) * t39, 0, (t53 * t40 - (t37 * t78 + t48 * t77 + t49 * t61) * t92) * t39, 0, 0; ((t63 * t68 + t66 * t77) * t44 - (-t63 * t66 + t68 * t77) * t91) * t43, ((t66 * t87 + t68 * t76) * t44 - (-t66 * t76 + t68 * t87) * t91) * t43, 0, (-t44 * t66 + t68 * t91) * t52 * t43, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:26:24
	% EndTime: 2019-10-10 10:26:24
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (525->37), mult. (1300->91), div. (85->9), fcn. (1858->13), ass. (0->55)
	t78 = cos(pkin(6));
	t83 = cos(qJ(2));
	t84 = cos(qJ(1));
	t89 = t84 * t83;
	t80 = sin(qJ(2));
	t81 = sin(qJ(1));
	t92 = t81 * t80;
	t70 = -t78 * t89 + t92;
	t79 = sin(qJ(4));
	t82 = cos(qJ(4));
	t77 = sin(pkin(6));
	t93 = t77 * t84;
	t62 = t70 * t82 + t79 * t93;
	t94 = t77 * t83;
	t68 = t78 * t79 + t82 * t94;
	t59 = atan2(t62, t68);
	t56 = sin(t59);
	t57 = cos(t59);
	t50 = t56 * t62 + t57 * t68;
	t49 = 0.1e1 / t50 ^ 2;
	t90 = t84 * t80;
	t91 = t81 * t83;
	t85 = t78 * t91 + t90;
	t95 = t77 * t81;
	t60 = t79 * t95 - t85 * t82;
	t102 = t49 * t60;
	t61 = t85 * t79 + t82 * t95;
	t72 = -t78 * t92 + t89;
	t76 = pkin(11) + qJ(6);
	t74 = sin(t76);
	t75 = cos(t76);
	t55 = t61 * t75 + t72 * t74;
	t53 = 0.1e1 / t55 ^ 2;
	t54 = t61 * t74 - t72 * t75;
	t101 = t53 * t54;
	t100 = t57 * t62;
	t99 = t60 ^ 2 * t49;
	t67 = 0.1e1 / t68 ^ 2;
	t98 = t62 * t67;
	t97 = t72 * t79;
	t96 = t77 * t80;
	t88 = t54 ^ 2 * t53 + 0.1e1;
	t87 = -t56 * t68 + t100;
	t86 = -t70 * t79 + t82 * t93;
	t71 = t78 * t90 + t91;
	t69 = t78 * t82 - t79 * t94;
	t66 = 0.1e1 / t68;
	t58 = 0.1e1 / (t62 ^ 2 * t67 + 0.1e1);
	t52 = 0.1e1 / t55;
	t51 = 0.1e1 / t88;
	t48 = 0.1e1 / t50;
	t47 = 0.1e1 / (0.1e1 + t99);
	t46 = (t66 * t71 + t96 * t98) * t82 * t58;
	t45 = (t66 * t86 - t69 * t98) * t58;
	t1 = [-t60 * t66 * t58, t46, 0, t45, 0, 0; (t62 * t48 - (-t56 + (-t66 * t100 + t56) * t58) * t99) * t47, (-t72 * t82 * t48 - ((t56 * t71 - t57 * t96) * t82 + t87 * t46) * t102) * t47, 0, (t61 * t48 - (t87 * t45 + t56 * t86 + t57 * t69) * t102) * t47, 0, 0; ((t71 * t75 + t74 * t86) * t52 - (-t71 * t74 + t75 * t86) * t101) * t51, ((t74 * t97 + t85 * t75) * t52 - (-t85 * t74 + t75 * t97) * t101) * t51, 0, (t75 * t101 - t74 * t52) * t60 * t51, 0, t88 * t51;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end