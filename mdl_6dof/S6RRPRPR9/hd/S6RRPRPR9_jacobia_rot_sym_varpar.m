% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR9
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
%   Wie in S6RRPRPR9_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:19
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR9_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR9_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:18:59
	% EndTime: 2019-10-10 10:18:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:18:59
	% EndTime: 2019-10-10 10:18:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:18:59
	% EndTime: 2019-10-10 10:18:59
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
	% StartTime: 2019-10-10 10:18:59
	% EndTime: 2019-10-10 10:18:59
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (156->24), mult. (403->63), div. (67->11), fcn. (610->11), ass. (0->39)
	t46 = sin(qJ(2));
	t48 = cos(qJ(2));
	t49 = cos(qJ(1));
	t47 = sin(qJ(1));
	t52 = cos(pkin(6));
	t51 = t47 * t52;
	t38 = -t46 * t51 + t49 * t48;
	t43 = sin(pkin(11));
	t45 = cos(pkin(11));
	t44 = sin(pkin(6));
	t55 = t44 * t47;
	t29 = t38 * t45 + t43 * t55;
	t27 = 0.1e1 / t29 ^ 2;
	t28 = t38 * t43 - t45 * t55;
	t59 = t27 * t28;
	t50 = t49 * t52;
	t34 = t47 * t46 - t48 * t50;
	t54 = t44 * t48;
	t32 = atan2(-t34, -t54);
	t31 = cos(t32);
	t58 = t31 * t34;
	t30 = sin(t32);
	t24 = -t30 * t34 - t31 * t54;
	t23 = 0.1e1 / t24 ^ 2;
	t37 = t49 * t46 + t48 * t51;
	t57 = t37 ^ 2 * t23;
	t40 = 0.1e1 / t44;
	t41 = 0.1e1 / t48;
	t56 = t40 * t41;
	t53 = t44 * t49;
	t42 = 0.1e1 / t48 ^ 2;
	t36 = t46 * t50 + t47 * t48;
	t33 = 0.1e1 / (0.1e1 + t34 ^ 2 / t44 ^ 2 * t42);
	t26 = 0.1e1 / t29;
	t25 = 0.1e1 / (t28 ^ 2 * t27 + 0.1e1);
	t22 = 0.1e1 / t24;
	t21 = 0.1e1 / (0.1e1 + t57);
	t20 = (t34 * t42 * t46 + t36 * t41) * t40 * t33;
	t1 = [t37 * t33 * t56, t20, 0, 0, 0, 0; (-t34 * t22 - (-t30 + (-t56 * t58 + t30) * t33) * t57) * t21, (t38 * t22 - (t31 * t44 * t46 - t30 * t36 + (t30 * t54 - t58) * t20) * t37 * t23) * t21, 0, 0, 0, 0; ((-t36 * t43 - t45 * t53) * t26 - (-t36 * t45 + t43 * t53) * t59) * t25, (-t43 * t26 + t45 * t59) * t37 * t25, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:18:59
	% EndTime: 2019-10-10 10:18:59
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (221->25), mult. (451->64), div. (72->11), fcn. (673->11), ass. (0->41)
	t58 = sin(qJ(2));
	t60 = cos(qJ(2));
	t61 = cos(qJ(1));
	t59 = sin(qJ(1));
	t65 = cos(pkin(6));
	t63 = t59 * t65;
	t49 = -t58 * t63 + t61 * t60;
	t54 = pkin(11) + qJ(4);
	t51 = sin(t54);
	t52 = cos(t54);
	t57 = sin(pkin(6));
	t68 = t57 * t59;
	t40 = t49 * t52 + t51 * t68;
	t38 = 0.1e1 / t40 ^ 2;
	t39 = t49 * t51 - t52 * t68;
	t72 = t38 * t39;
	t62 = t61 * t65;
	t45 = t59 * t58 - t60 * t62;
	t67 = t57 * t60;
	t43 = atan2(-t45, -t67);
	t42 = cos(t43);
	t71 = t42 * t45;
	t41 = sin(t43);
	t35 = -t41 * t45 - t42 * t67;
	t34 = 0.1e1 / t35 ^ 2;
	t48 = t61 * t58 + t60 * t63;
	t70 = t48 ^ 2 * t34;
	t53 = 0.1e1 / t57;
	t55 = 0.1e1 / t60;
	t69 = t53 * t55;
	t66 = t57 * t61;
	t64 = t39 ^ 2 * t38 + 0.1e1;
	t56 = 0.1e1 / t60 ^ 2;
	t47 = t58 * t62 + t59 * t60;
	t44 = 0.1e1 / (0.1e1 + t45 ^ 2 / t57 ^ 2 * t56);
	t37 = 0.1e1 / t40;
	t36 = 0.1e1 / t64;
	t33 = 0.1e1 / t35;
	t32 = 0.1e1 / (0.1e1 + t70);
	t31 = (t45 * t56 * t58 + t47 * t55) * t53 * t44;
	t1 = [t48 * t44 * t69, t31, 0, 0, 0, 0; (-t45 * t33 - (-t41 + (-t69 * t71 + t41) * t44) * t70) * t32, (t49 * t33 - (t42 * t57 * t58 - t41 * t47 + (t41 * t67 - t71) * t31) * t48 * t34) * t32, 0, 0, 0, 0; ((-t47 * t51 - t52 * t66) * t37 - (-t47 * t52 + t51 * t66) * t72) * t36, (-t51 * t37 + t52 * t72) * t48 * t36, 0, t64 * t36, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:18:59
	% EndTime: 2019-10-10 10:18:59
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (833->38), mult. (1217->89), div. (80->9), fcn. (1746->13), ass. (0->55)
	t80 = cos(pkin(6));
	t81 = sin(qJ(2));
	t84 = cos(qJ(1));
	t87 = t84 * t81;
	t82 = sin(qJ(1));
	t83 = cos(qJ(2));
	t88 = t82 * t83;
	t69 = t80 * t87 + t88;
	t76 = pkin(11) + qJ(4);
	t74 = sin(t76);
	t75 = cos(t76);
	t78 = sin(pkin(6));
	t90 = t78 * t84;
	t58 = t69 * t74 + t75 * t90;
	t93 = t78 * t81;
	t66 = t74 * t93 - t80 * t75;
	t57 = atan2(-t58, t66);
	t52 = sin(t57);
	t53 = cos(t57);
	t48 = -t52 * t58 + t53 * t66;
	t47 = 0.1e1 / t48 ^ 2;
	t86 = t84 * t83;
	t89 = t82 * t81;
	t71 = -t80 * t89 + t86;
	t92 = t78 * t82;
	t62 = t71 * t74 - t75 * t92;
	t100 = t47 * t62;
	t63 = t71 * t75 + t74 * t92;
	t79 = cos(pkin(12));
	t70 = t80 * t88 + t87;
	t77 = sin(pkin(12));
	t95 = t70 * t77;
	t55 = t63 * t79 + t95;
	t51 = 0.1e1 / t55 ^ 2;
	t94 = t70 * t79;
	t54 = t63 * t77 - t94;
	t99 = t51 * t54;
	t98 = t53 * t58;
	t65 = 0.1e1 / t66 ^ 2;
	t97 = t58 * t65;
	t96 = t62 ^ 2 * t47;
	t91 = t78 * t83;
	t60 = t69 * t75 - t74 * t90;
	t85 = -t52 * t66 - t98;
	t68 = t80 * t86 - t89;
	t67 = t80 * t74 + t75 * t93;
	t64 = 0.1e1 / t66;
	t56 = 0.1e1 / (t58 ^ 2 * t65 + 0.1e1);
	t50 = 0.1e1 / t55;
	t49 = 0.1e1 / (t54 ^ 2 * t51 + 0.1e1);
	t46 = 0.1e1 / t48;
	t45 = 0.1e1 / (0.1e1 + t96);
	t44 = (-t64 * t68 + t91 * t97) * t74 * t56;
	t43 = (-t60 * t64 + t67 * t97) * t56;
	t1 = [-t62 * t64 * t56, t44, 0, t43, 0, 0; (-t58 * t46 - (-t52 + (t64 * t98 + t52) * t56) * t96) * t45, (-t70 * t74 * t46 - ((-t52 * t68 + t53 * t91) * t74 + t85 * t44) * t100) * t45, 0, (t63 * t46 - (t85 * t43 - t52 * t60 + t53 * t67) * t100) * t45, 0, 0; ((-t60 * t77 - t68 * t79) * t50 - (-t60 * t79 + t68 * t77) * t99) * t49, ((-t71 * t79 - t75 * t95) * t50 - (t71 * t77 - t75 * t94) * t99) * t49, 0, (-t77 * t50 + t79 * t99) * t62 * t49, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:18:59
	% EndTime: 2019-10-10 10:18:59
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (944->39), mult. (1300->91), div. (85->9), fcn. (1858->13), ass. (0->56)
	t88 = cos(pkin(6));
	t89 = sin(qJ(2));
	t92 = cos(qJ(1));
	t96 = t92 * t89;
	t90 = sin(qJ(1));
	t91 = cos(qJ(2));
	t97 = t90 * t91;
	t76 = t88 * t96 + t97;
	t86 = pkin(11) + qJ(4);
	t82 = sin(t86);
	t84 = cos(t86);
	t87 = sin(pkin(6));
	t99 = t87 * t92;
	t65 = t76 * t82 + t84 * t99;
	t102 = t87 * t89;
	t73 = t82 * t102 - t88 * t84;
	t64 = atan2(-t65, t73);
	t61 = sin(t64);
	t62 = cos(t64);
	t55 = -t61 * t65 + t62 * t73;
	t54 = 0.1e1 / t55 ^ 2;
	t101 = t87 * t90;
	t95 = t92 * t91;
	t98 = t90 * t89;
	t78 = -t88 * t98 + t95;
	t69 = -t84 * t101 + t78 * t82;
	t108 = t54 * t69;
	t70 = t82 * t101 + t78 * t84;
	t77 = t88 * t97 + t96;
	t85 = pkin(12) + qJ(6);
	t81 = sin(t85);
	t83 = cos(t85);
	t60 = t70 * t83 + t77 * t81;
	t58 = 0.1e1 / t60 ^ 2;
	t59 = t70 * t81 - t77 * t83;
	t107 = t58 * t59;
	t106 = t62 * t65;
	t72 = 0.1e1 / t73 ^ 2;
	t105 = t65 * t72;
	t104 = t69 ^ 2 * t54;
	t103 = t77 * t84;
	t100 = t87 * t91;
	t94 = t59 ^ 2 * t58 + 0.1e1;
	t67 = t76 * t84 - t82 * t99;
	t93 = -t61 * t73 - t106;
	t75 = t88 * t95 - t98;
	t74 = t84 * t102 + t88 * t82;
	t71 = 0.1e1 / t73;
	t63 = 0.1e1 / (t65 ^ 2 * t72 + 0.1e1);
	t57 = 0.1e1 / t60;
	t56 = 0.1e1 / t94;
	t53 = 0.1e1 / t55;
	t52 = 0.1e1 / (0.1e1 + t104);
	t51 = (t100 * t105 - t71 * t75) * t82 * t63;
	t50 = (t74 * t105 - t67 * t71) * t63;
	t1 = [-t69 * t71 * t63, t51, 0, t50, 0, 0; (-t65 * t53 - (-t61 + (t71 * t106 + t61) * t63) * t104) * t52, (-t77 * t82 * t53 - ((t62 * t100 - t61 * t75) * t82 + t93 * t51) * t108) * t52, 0, (t70 * t53 - (t93 * t50 - t61 * t67 + t62 * t74) * t108) * t52, 0, 0; ((-t67 * t81 - t75 * t83) * t57 - (-t67 * t83 + t75 * t81) * t107) * t56, ((-t81 * t103 - t78 * t83) * t57 - (-t83 * t103 + t78 * t81) * t107) * t56, 0, (t83 * t107 - t81 * t57) * t69 * t56, 0, t94 * t56;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end