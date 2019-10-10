% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR12
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
%   Wie in S6RRPRPR12_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:24
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR12_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR12_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
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
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
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
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
	% DurationCPUTime: 0.13s
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
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (197->25), mult. (451->64), div. (72->11), fcn. (673->11), ass. (0->41)
	t62 = sin(qJ(2));
	t64 = cos(qJ(2));
	t65 = cos(qJ(1));
	t63 = sin(qJ(1));
	t69 = cos(pkin(6));
	t67 = t63 * t69;
	t51 = t65 * t62 + t64 * t67;
	t58 = qJ(4) + pkin(11);
	t55 = sin(t58);
	t56 = cos(t58);
	t61 = sin(pkin(6));
	t71 = t61 * t63;
	t43 = t51 * t55 + t56 * t71;
	t41 = 0.1e1 / t43 ^ 2;
	t42 = -t51 * t56 + t55 * t71;
	t76 = t41 * t42;
	t66 = t65 * t69;
	t49 = t62 * t66 + t63 * t64;
	t72 = t61 * t62;
	t47 = atan2(-t49, t72);
	t45 = cos(t47);
	t75 = t45 * t49;
	t44 = sin(t47);
	t38 = -t44 * t49 + t45 * t72;
	t37 = 0.1e1 / t38 ^ 2;
	t52 = -t62 * t67 + t65 * t64;
	t74 = t52 ^ 2 * t37;
	t57 = 0.1e1 / t61;
	t59 = 0.1e1 / t62;
	t73 = t57 * t59;
	t70 = t61 * t65;
	t68 = t42 ^ 2 * t41 + 0.1e1;
	t60 = 0.1e1 / t62 ^ 2;
	t48 = t63 * t62 - t64 * t66;
	t46 = 0.1e1 / (0.1e1 + t49 ^ 2 / t61 ^ 2 * t60);
	t40 = 0.1e1 / t43;
	t39 = 0.1e1 / t68;
	t36 = 0.1e1 / t38;
	t35 = 0.1e1 / (0.1e1 + t74);
	t34 = (t49 * t60 * t64 + t48 * t59) * t57 * t46;
	t1 = [-t52 * t46 * t73, t34, 0, 0, 0, 0; (-t49 * t36 - (-t44 + (t73 * t75 + t44) * t46) * t74) * t35, (-t51 * t36 - (t45 * t61 * t64 + t44 * t48 + (-t44 * t72 - t75) * t34) * t52 * t37) * t35, 0, 0, 0, 0; ((t48 * t56 + t55 * t70) * t40 - (-t48 * t55 + t56 * t70) * t76) * t39, (-t56 * t40 - t55 * t76) * t52 * t39, 0, t68 * t39, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (878->37), mult. (1300->90), div. (85->9), fcn. (1858->13), ass. (0->56)
	t81 = cos(pkin(6));
	t86 = cos(qJ(2));
	t87 = cos(qJ(1));
	t92 = t87 * t86;
	t83 = sin(qJ(2));
	t84 = sin(qJ(1));
	t95 = t84 * t83;
	t73 = -t81 * t92 + t95;
	t79 = qJ(4) + pkin(11);
	t77 = sin(t79);
	t78 = cos(t79);
	t80 = sin(pkin(6));
	t96 = t80 * t87;
	t65 = t73 * t78 + t77 * t96;
	t97 = t80 * t86;
	t71 = t81 * t77 + t78 * t97;
	t62 = atan2(t65, t71);
	t55 = sin(t62);
	t56 = cos(t62);
	t53 = t55 * t65 + t56 * t71;
	t52 = 0.1e1 / t53 ^ 2;
	t93 = t87 * t83;
	t94 = t84 * t86;
	t88 = t81 * t94 + t93;
	t98 = t80 * t84;
	t63 = t77 * t98 - t78 * t88;
	t106 = t52 * t63;
	t105 = t56 * t65;
	t75 = -t81 * t95 + t92;
	t82 = sin(qJ(6));
	t101 = t75 * t82;
	t64 = t77 * t88 + t78 * t98;
	t85 = cos(qJ(6));
	t61 = t64 * t85 + t101;
	t58 = 0.1e1 / t61 ^ 2;
	t100 = t75 * t85;
	t60 = t64 * t82 - t100;
	t104 = t58 * t60;
	t103 = t63 ^ 2 * t52;
	t70 = 0.1e1 / t71 ^ 2;
	t102 = t65 * t70;
	t99 = t80 * t83;
	t91 = t60 ^ 2 * t58 + 0.1e1;
	t90 = -t55 * t71 + t105;
	t89 = -t73 * t77 + t78 * t96;
	t74 = t81 * t93 + t94;
	t72 = -t77 * t97 + t81 * t78;
	t69 = 0.1e1 / t71;
	t59 = 0.1e1 / (t65 ^ 2 * t70 + 0.1e1);
	t57 = 0.1e1 / t61;
	t54 = 0.1e1 / t91;
	t51 = 0.1e1 / t53;
	t50 = 0.1e1 / (0.1e1 + t103);
	t49 = (t102 * t99 + t69 * t74) * t78 * t59;
	t48 = (-t102 * t72 + t69 * t89) * t59;
	t1 = [-t63 * t69 * t59, t49, 0, t48, 0, 0; (t65 * t51 - (-t55 + (-t105 * t69 + t55) * t59) * t103) * t50, (-t75 * t78 * t51 - ((t55 * t74 - t56 * t99) * t78 + t90 * t49) * t106) * t50, 0, (t64 * t51 - (t48 * t90 + t55 * t89 + t56 * t72) * t106) * t50, 0, 0; ((t74 * t85 + t82 * t89) * t57 - (-t74 * t82 + t85 * t89) * t104) * t54, ((t101 * t77 + t85 * t88) * t57 - (t100 * t77 - t82 * t88) * t104) * t54, 0, (t104 * t85 - t82 * t57) * t63 * t54, 0, t91 * t54;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end