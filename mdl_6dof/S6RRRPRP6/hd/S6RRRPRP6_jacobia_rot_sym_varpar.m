% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRP6
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
%   Wie in S6RRRPRP6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:44
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRP6_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP6_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:43:59
	% EndTime: 2019-10-10 11:43:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:43:59
	% EndTime: 2019-10-10 11:43:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:43:59
	% EndTime: 2019-10-10 11:43:59
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
	% StartTime: 2019-10-10 11:43:59
	% EndTime: 2019-10-10 11:43:59
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (173->24), mult. (451->64), div. (72->11), fcn. (673->11), ass. (0->40)
	t53 = sin(qJ(2));
	t56 = cos(qJ(2));
	t57 = cos(qJ(1));
	t54 = sin(qJ(1));
	t61 = cos(pkin(6));
	t59 = t54 * t61;
	t46 = -t53 * t59 + t57 * t56;
	t52 = sin(qJ(3));
	t55 = cos(qJ(3));
	t51 = sin(pkin(6));
	t64 = t51 * t54;
	t37 = t46 * t55 + t52 * t64;
	t35 = 0.1e1 / t37 ^ 2;
	t36 = t46 * t52 - t55 * t64;
	t68 = t35 * t36;
	t58 = t57 * t61;
	t42 = t54 * t53 - t56 * t58;
	t63 = t51 * t56;
	t40 = atan2(-t42, -t63);
	t39 = cos(t40);
	t67 = t39 * t42;
	t38 = sin(t40);
	t32 = -t38 * t42 - t39 * t63;
	t31 = 0.1e1 / t32 ^ 2;
	t45 = t57 * t53 + t56 * t59;
	t66 = t45 ^ 2 * t31;
	t48 = 0.1e1 / t51;
	t49 = 0.1e1 / t56;
	t65 = t48 * t49;
	t62 = t51 * t57;
	t60 = t36 ^ 2 * t35 + 0.1e1;
	t50 = 0.1e1 / t56 ^ 2;
	t44 = t53 * t58 + t54 * t56;
	t41 = 0.1e1 / (0.1e1 + t42 ^ 2 / t51 ^ 2 * t50);
	t34 = 0.1e1 / t37;
	t33 = 0.1e1 / t60;
	t30 = 0.1e1 / t32;
	t29 = 0.1e1 / (0.1e1 + t66);
	t28 = (t42 * t50 * t53 + t44 * t49) * t48 * t41;
	t1 = [t45 * t41 * t65, t28, 0, 0, 0, 0; (-t42 * t30 - (-t38 + (-t65 * t67 + t38) * t41) * t66) * t29, (t46 * t30 - (t39 * t51 * t53 - t38 * t44 + (t38 * t63 - t67) * t28) * t45 * t31) * t29, 0, 0, 0, 0; ((-t44 * t52 - t55 * t62) * t34 - (-t44 * t55 + t52 * t62) * t68) * t33, (-t52 * t34 + t55 * t68) * t45 * t33, t60 * t33, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:43:59
	% EndTime: 2019-10-10 11:43:59
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (221->25), mult. (451->64), div. (72->11), fcn. (673->11), ass. (0->41)
	t59 = sin(qJ(2));
	t61 = cos(qJ(2));
	t62 = cos(qJ(1));
	t60 = sin(qJ(1));
	t66 = cos(pkin(6));
	t64 = t60 * t66;
	t50 = -t59 * t64 + t62 * t61;
	t55 = qJ(3) + pkin(11);
	t52 = sin(t55);
	t53 = cos(t55);
	t58 = sin(pkin(6));
	t69 = t58 * t60;
	t41 = t50 * t53 + t52 * t69;
	t39 = 0.1e1 / t41 ^ 2;
	t40 = t50 * t52 - t53 * t69;
	t73 = t39 * t40;
	t63 = t62 * t66;
	t46 = t60 * t59 - t61 * t63;
	t68 = t58 * t61;
	t44 = atan2(-t46, -t68);
	t43 = cos(t44);
	t72 = t43 * t46;
	t42 = sin(t44);
	t36 = -t42 * t46 - t43 * t68;
	t35 = 0.1e1 / t36 ^ 2;
	t49 = t62 * t59 + t61 * t64;
	t71 = t49 ^ 2 * t35;
	t54 = 0.1e1 / t58;
	t56 = 0.1e1 / t61;
	t70 = t54 * t56;
	t67 = t58 * t62;
	t65 = t40 ^ 2 * t39 + 0.1e1;
	t57 = 0.1e1 / t61 ^ 2;
	t48 = t59 * t63 + t60 * t61;
	t45 = 0.1e1 / (0.1e1 + t46 ^ 2 / t58 ^ 2 * t57);
	t38 = 0.1e1 / t41;
	t37 = 0.1e1 / t65;
	t34 = 0.1e1 / t36;
	t33 = 0.1e1 / (0.1e1 + t71);
	t32 = (t46 * t57 * t59 + t48 * t56) * t54 * t45;
	t1 = [t49 * t45 * t70, t32, 0, 0, 0, 0; (-t46 * t34 - (-t42 + (-t70 * t72 + t42) * t45) * t71) * t33, (t50 * t34 - (t43 * t58 * t59 - t42 * t48 + (t42 * t68 - t72) * t32) * t49 * t35) * t33, 0, 0, 0, 0; ((-t48 * t52 - t53 * t67) * t38 - (-t48 * t53 + t52 * t67) * t73) * t37, (-t52 * t38 + t53 * t73) * t49 * t37, t65 * t37, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:43:59
	% EndTime: 2019-10-10 11:44:00
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (878->38), mult. (1300->90), div. (85->9), fcn. (1858->13), ass. (0->56)
	t85 = cos(pkin(6));
	t87 = sin(qJ(2));
	t91 = cos(qJ(1));
	t95 = t91 * t87;
	t88 = sin(qJ(1));
	t90 = cos(qJ(2));
	t96 = t88 * t90;
	t76 = t85 * t95 + t96;
	t83 = qJ(3) + pkin(11);
	t81 = sin(t83);
	t82 = cos(t83);
	t84 = sin(pkin(6));
	t98 = t84 * t91;
	t65 = t76 * t81 + t82 * t98;
	t101 = t84 * t87;
	t73 = t81 * t101 - t85 * t82;
	t64 = atan2(-t65, t73);
	t57 = sin(t64);
	t58 = cos(t64);
	t55 = -t57 * t65 + t58 * t73;
	t54 = 0.1e1 / t55 ^ 2;
	t100 = t84 * t88;
	t94 = t91 * t90;
	t97 = t88 * t87;
	t78 = -t85 * t97 + t94;
	t69 = -t82 * t100 + t78 * t81;
	t108 = t54 * t69;
	t107 = t58 * t65;
	t77 = t85 * t96 + t95;
	t86 = sin(qJ(5));
	t103 = t77 * t86;
	t70 = t81 * t100 + t78 * t82;
	t89 = cos(qJ(5));
	t63 = t70 * t89 + t103;
	t60 = 0.1e1 / t63 ^ 2;
	t102 = t77 * t89;
	t62 = t70 * t86 - t102;
	t106 = t60 * t62;
	t72 = 0.1e1 / t73 ^ 2;
	t105 = t65 * t72;
	t104 = t69 ^ 2 * t54;
	t99 = t84 * t90;
	t93 = t62 ^ 2 * t60 + 0.1e1;
	t67 = t76 * t82 - t81 * t98;
	t92 = -t57 * t73 - t107;
	t75 = t85 * t94 - t97;
	t74 = t82 * t101 + t85 * t81;
	t71 = 0.1e1 / t73;
	t61 = 0.1e1 / (t65 ^ 2 * t72 + 0.1e1);
	t59 = 0.1e1 / t63;
	t56 = 0.1e1 / t93;
	t53 = 0.1e1 / t55;
	t52 = 0.1e1 / (0.1e1 + t104);
	t51 = (t99 * t105 - t71 * t75) * t81 * t61;
	t50 = (t74 * t105 - t67 * t71) * t61;
	t1 = [-t69 * t71 * t61, t51, t50, 0, 0, 0; (-t65 * t53 - (-t57 + (t71 * t107 + t57) * t61) * t104) * t52, (-t77 * t81 * t53 - ((-t57 * t75 + t58 * t99) * t81 + t92 * t51) * t108) * t52, (t70 * t53 - (t92 * t50 - t57 * t67 + t58 * t74) * t108) * t52, 0, 0, 0; ((-t67 * t86 - t75 * t89) * t59 - (-t67 * t89 + t75 * t86) * t106) * t56, ((-t82 * t103 - t78 * t89) * t59 - (-t82 * t102 + t78 * t86) * t106) * t56, (t89 * t106 - t86 * t59) * t69 * t56, 0, t93 * t56, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:43:59
	% EndTime: 2019-10-10 11:44:00
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (878->38), mult. (1300->90), div. (85->9), fcn. (1858->13), ass. (0->56)
	t86 = sin(pkin(6));
	t93 = cos(qJ(1));
	t100 = t86 * t93;
	t87 = cos(pkin(6));
	t90 = sin(qJ(1));
	t92 = cos(qJ(2));
	t97 = t90 * t92;
	t89 = sin(qJ(2));
	t99 = t89 * t93;
	t78 = t87 * t99 + t97;
	t85 = qJ(3) + pkin(11);
	t83 = sin(t85);
	t84 = cos(t85);
	t67 = t84 * t100 + t78 * t83;
	t103 = t86 * t89;
	t75 = t83 * t103 - t84 * t87;
	t66 = atan2(-t67, t75);
	t59 = sin(t66);
	t60 = cos(t66);
	t57 = -t59 * t67 + t60 * t75;
	t56 = 0.1e1 / t57 ^ 2;
	t102 = t86 * t90;
	t96 = t92 * t93;
	t98 = t90 * t89;
	t80 = -t87 * t98 + t96;
	t71 = -t84 * t102 + t80 * t83;
	t110 = t56 * t71;
	t109 = t56 * t71 ^ 2;
	t108 = t60 * t67;
	t79 = t87 * t97 + t99;
	t88 = sin(qJ(5));
	t105 = t79 * t88;
	t72 = t83 * t102 + t80 * t84;
	t91 = cos(qJ(5));
	t65 = t72 * t91 + t105;
	t62 = 0.1e1 / t65 ^ 2;
	t104 = t79 * t91;
	t64 = t72 * t88 - t104;
	t107 = t62 * t64;
	t74 = 0.1e1 / t75 ^ 2;
	t106 = t67 * t74;
	t101 = t86 * t92;
	t95 = t62 * t64 ^ 2 + 0.1e1;
	t69 = -t83 * t100 + t78 * t84;
	t94 = -t59 * t75 - t108;
	t77 = t87 * t96 - t98;
	t76 = t84 * t103 + t83 * t87;
	t73 = 0.1e1 / t75;
	t63 = 0.1e1 / (t67 ^ 2 * t74 + 0.1e1);
	t61 = 0.1e1 / t65;
	t58 = 0.1e1 / t95;
	t55 = 0.1e1 / t57;
	t54 = 0.1e1 / (0.1e1 + t109);
	t53 = (t101 * t106 - t73 * t77) * t83 * t63;
	t52 = (t76 * t106 - t69 * t73) * t63;
	t1 = [-t71 * t73 * t63, t53, t52, 0, 0, 0; (-t67 * t55 - (-t59 + (t73 * t108 + t59) * t63) * t109) * t54, (-t79 * t83 * t55 - ((t60 * t101 - t59 * t77) * t83 + t94 * t53) * t110) * t54, (t72 * t55 - (t94 * t52 - t59 * t69 + t60 * t76) * t110) * t54, 0, 0, 0; ((-t69 * t88 - t77 * t91) * t61 - (-t69 * t91 + t77 * t88) * t107) * t58, ((-t84 * t105 - t80 * t91) * t61 - (-t84 * t104 + t80 * t88) * t107) * t58, (t91 * t107 - t88 * t61) * t71 * t58, 0, t95 * t58, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end