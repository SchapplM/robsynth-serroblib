% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRRR6
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:23
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRRR6_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR6_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR6_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_jacobia_transl_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:35
	% EndTime: 2019-10-09 23:23:35
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:35
	% EndTime: 2019-10-09 23:23:35
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:35
	% EndTime: 2019-10-09 23:23:35
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (7->7), mult. (20->17), div. (0->0), fcn. (24->6), ass. (0->8)
	t4 = cos(pkin(6));
	t5 = sin(qJ(2));
	t8 = t4 * t5;
	t6 = cos(qJ(2));
	t7 = t4 * t6;
	t3 = cos(pkin(14));
	t1 = sin(pkin(14));
	t2 = [0, (-t1 * t7 - t3 * t5) * r_i_i_C(1) + (t1 * t8 - t3 * t6) * r_i_i_C(2), 0, 0, 0, 0; 0, (-t1 * t5 + t3 * t7) * r_i_i_C(1) + (-t1 * t6 - t3 * t8) * r_i_i_C(2), 0, 0, 0, 0; 1, (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * sin(pkin(6)), 0, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:35
	% EndTime: 2019-10-09 23:23:35
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (55->24), mult. (156->50), div. (0->0), fcn. (196->10), ass. (0->24)
	t7 = sin(pkin(7));
	t8 = sin(pkin(6));
	t24 = t7 * t8;
	t10 = cos(pkin(7));
	t15 = cos(qJ(2));
	t23 = t10 * t15;
	t11 = cos(pkin(6));
	t13 = sin(qJ(2));
	t22 = t11 * t13;
	t21 = t11 * t15;
	t12 = sin(qJ(3));
	t14 = cos(qJ(3));
	t20 = r_i_i_C(1) * t14 - r_i_i_C(2) * t12;
	t6 = sin(pkin(14));
	t9 = cos(pkin(14));
	t1 = -t6 * t13 + t9 * t21;
	t19 = -t1 * t10 + t9 * t24;
	t3 = -t9 * t13 - t6 * t21;
	t18 = t10 * t3 + t6 * t24;
	t17 = pkin(2) + t20;
	t16 = (pkin(10) + r_i_i_C(3)) * t7 + (-t12 * r_i_i_C(1) - t14 * r_i_i_C(2)) * t10;
	t4 = t15 * t9 - t6 * t22;
	t2 = t15 * t6 + t9 * t22;
	t5 = [0, t16 * t4 + t17 * t3, (-t12 * t4 + t18 * t14) * r_i_i_C(1) + (-t18 * t12 - t14 * t4) * r_i_i_C(2), 0, 0, 0; 0, t17 * t1 + t16 * t2, (-t12 * t2 - t19 * t14) * r_i_i_C(1) + (t19 * t12 - t14 * t2) * r_i_i_C(2), 0, 0, 0; 1, (t16 * t13 + t17 * t15) * t8, t20 * t7 * t11 + ((-t12 * t13 + t14 * t23) * r_i_i_C(1) + (-t12 * t23 - t13 * t14) * r_i_i_C(2)) * t8, 0, 0, 0;];
	Ja_transl = t5;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:35
	% EndTime: 2019-10-09 23:23:36
	% DurationCPUTime: 0.45s
	% Computational Cost: add. (250->71), mult. (730->136), div. (0->0), fcn. (944->14), ass. (0->57)
	t26 = cos(pkin(8));
	t28 = sin(qJ(4));
	t31 = cos(qJ(4));
	t45 = t28 * r_i_i_C(1) + t31 * r_i_i_C(2);
	t22 = sin(pkin(8));
	t65 = pkin(11) + r_i_i_C(3);
	t49 = t65 * t22;
	t34 = -t45 * t26 + t49;
	t23 = sin(pkin(7));
	t64 = t23 * pkin(10);
	t25 = cos(pkin(14));
	t30 = sin(qJ(2));
	t33 = cos(qJ(2));
	t51 = sin(pkin(14));
	t52 = cos(pkin(6));
	t41 = t52 * t51;
	t19 = t25 * t33 - t30 * t41;
	t29 = sin(qJ(3));
	t32 = cos(qJ(3));
	t18 = -t25 * t30 - t33 * t41;
	t27 = cos(pkin(7));
	t24 = sin(pkin(6));
	t48 = t24 * t51;
	t35 = t18 * t27 + t23 * t48;
	t4 = t19 * t32 + t35 * t29;
	t63 = t4 * t28;
	t62 = t4 * t31;
	t61 = t22 * t23;
	t60 = t23 * t26;
	t59 = t24 * t25;
	t58 = t27 * t29;
	t57 = t27 * t32;
	t56 = t29 * t30;
	t55 = t29 * t33;
	t54 = t30 * t32;
	t53 = t32 * t33;
	t50 = t23 * t59;
	t47 = t25 * t52;
	t46 = t52 * t23;
	t16 = -t51 * t30 + t33 * t47;
	t17 = t30 * t47 + t51 * t33;
	t1 = -t17 * t29 + (t16 * t27 - t50) * t32;
	t44 = t1 * t26 + (-t16 * t23 - t27 * t59) * t22;
	t3 = -t19 * t29 + t35 * t32;
	t43 = (-t18 * t23 + t27 * t48) * t22 + t26 * t3;
	t9 = t32 * t46 + (t27 * t53 - t56) * t24;
	t42 = (-t24 * t33 * t23 + t52 * t27) * t22 + t26 * t9;
	t40 = t31 * r_i_i_C(1) - t28 * r_i_i_C(2) + pkin(3);
	t5 = -t16 * t29 - t17 * t57;
	t38 = t17 * t61 + t26 * t5;
	t7 = -t18 * t29 - t19 * t57;
	t36 = t19 * t61 + t26 * t7;
	t10 = t29 * t46 + (t27 * t55 + t54) * t24;
	t8 = t18 * t32 - t19 * t58;
	t6 = t16 * t32 - t17 * t58;
	t2 = t16 * t58 + t17 * t32 - t29 * t50;
	t11 = [0, (t36 * t28 + t8 * t31) * r_i_i_C(1) + (-t8 * t28 + t36 * t31) * r_i_i_C(2) + t8 * pkin(3) + t18 * pkin(2) + t19 * t64 + t65 * (t19 * t60 - t7 * t22), (-t26 * t63 + t3 * t31) * r_i_i_C(1) + (-t26 * t62 - t3 * t28) * r_i_i_C(2) + t3 * pkin(3) + t4 * t49, (t43 * t31 - t63) * r_i_i_C(1) + (-t43 * t28 - t62) * r_i_i_C(2), 0, 0; 0, (t38 * t28 + t6 * t31) * r_i_i_C(1) + (-t6 * t28 + t38 * t31) * r_i_i_C(2) + t6 * pkin(3) + t16 * pkin(2) + t17 * t64 + t65 * (t17 * t60 - t5 * t22), t40 * t1 + t34 * t2, (-t2 * t28 + t44 * t31) * r_i_i_C(1) + (-t2 * t31 - t44 * t28) * r_i_i_C(2), 0, 0; 1, (t40 * (-t27 * t56 + t53) - t34 * (-t27 * t54 - t55) + t33 * pkin(2) + (t45 * t22 + t65 * t26 + pkin(10)) * t30 * t23) * t24, t34 * t10 + t40 * t9, (-t10 * t28 + t42 * t31) * r_i_i_C(1) + (-t10 * t31 - t42 * t28) * r_i_i_C(2), 0, 0;];
	Ja_transl = t11;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:35
	% EndTime: 2019-10-09 23:23:36
	% DurationCPUTime: 0.54s
	% Computational Cost: add. (640->118), mult. (1862->221), div. (0->0), fcn. (2442->16), ass. (0->83)
	t93 = r_i_i_C(3) + pkin(12);
	t46 = sin(pkin(8));
	t92 = t46 * pkin(11);
	t47 = sin(pkin(7));
	t91 = t47 * pkin(10);
	t48 = sin(pkin(6));
	t51 = cos(pkin(7));
	t55 = sin(qJ(2));
	t58 = cos(qJ(3));
	t77 = t55 * t58;
	t54 = sin(qJ(3));
	t59 = cos(qJ(2));
	t78 = t54 * t59;
	t37 = (-t51 * t77 - t78) * t48;
	t90 = t37 * t46;
	t89 = t46 * t47;
	t52 = sin(qJ(5));
	t88 = t46 * t52;
	t56 = cos(qJ(5));
	t87 = t46 * t56;
	t50 = cos(pkin(8));
	t86 = t47 * t50;
	t85 = t47 * t55;
	t49 = cos(pkin(14));
	t84 = t48 * t49;
	t53 = sin(qJ(4));
	t83 = t50 * t53;
	t57 = cos(qJ(4));
	t82 = t50 * t57;
	t81 = t51 * t54;
	t80 = t51 * t58;
	t79 = t54 * t55;
	t76 = t58 * t59;
	t75 = cos(pkin(6));
	t74 = sin(pkin(14));
	t73 = t47 * t84;
	t72 = t48 * t85;
	t71 = t48 * t74;
	t70 = t49 * t75;
	t69 = t75 * t47;
	t40 = -t55 * t74 + t59 * t70;
	t41 = t55 * t70 + t59 * t74;
	t24 = -t41 * t54 + (t40 * t51 - t73) * t58;
	t35 = -t40 * t47 - t51 * t84;
	t68 = t24 * t50 + t35 * t46;
	t65 = t75 * t74;
	t43 = t49 * t59 - t55 * t65;
	t42 = -t49 * t55 - t59 * t65;
	t60 = t42 * t51 + t47 * t71;
	t26 = -t43 * t54 + t58 * t60;
	t36 = -t42 * t47 + t51 * t71;
	t67 = t26 * t50 + t36 * t46;
	t33 = t58 * t69 + (t51 * t76 - t79) * t48;
	t39 = -t48 * t59 * t47 + t51 * t75;
	t66 = t33 * t50 + t39 * t46;
	t64 = r_i_i_C(1) * t56 - r_i_i_C(2) * t52 + pkin(4);
	t28 = -t40 * t54 - t41 * t80;
	t19 = -t28 * t46 + t41 * t86;
	t63 = t28 * t50 + t41 * t89;
	t30 = -t42 * t54 - t43 * t80;
	t20 = -t30 * t46 + t43 * t86;
	t62 = t30 * t50 + t43 * t89;
	t61 = t37 * t50 + t46 * t72;
	t38 = (-t51 * t79 + t76) * t48;
	t34 = t54 * t69 + (t51 * t78 + t77) * t48;
	t32 = t50 * t72 - t90;
	t31 = t42 * t58 - t43 * t81;
	t29 = t40 * t58 - t41 * t81;
	t27 = t43 * t58 + t54 * t60;
	t25 = t40 * t81 + t41 * t58 - t54 * t73;
	t23 = -t33 * t46 + t39 * t50;
	t22 = t38 * t57 + t53 * t61;
	t18 = t33 * t57 - t34 * t83;
	t16 = -t26 * t46 + t36 * t50;
	t15 = -t24 * t46 + t35 * t50;
	t14 = t34 * t57 + t53 * t66;
	t12 = t26 * t57 - t27 * t83;
	t10 = t24 * t57 - t25 * t83;
	t8 = t31 * t57 + t53 * t62;
	t6 = t29 * t57 + t53 * t63;
	t4 = t27 * t57 + t53 * t67;
	t2 = t25 * t57 + t53 * t68;
	t1 = [0, (t20 * t52 + t56 * t8) * r_i_i_C(1) + (t20 * t56 - t52 * t8) * r_i_i_C(2) + t8 * pkin(4) + t31 * pkin(3) + t42 * pkin(2) + t43 * t91 + t93 * (t31 * t53 - t57 * t62) + t20 * pkin(11), (t12 * t56 + t27 * t88) * r_i_i_C(1) + (-t12 * t52 + t27 * t87) * r_i_i_C(2) + t12 * pkin(4) + t26 * pkin(3) + t27 * t92 + t93 * (t26 * t53 + t27 * t82), t93 * t4 + t64 * (-t27 * t53 + t57 * t67), (t16 * t56 - t4 * t52) * r_i_i_C(1) + (-t16 * t52 - t4 * t56) * r_i_i_C(2), 0; 0, (t19 * t52 + t56 * t6) * r_i_i_C(1) + (t19 * t56 - t52 * t6) * r_i_i_C(2) + t6 * pkin(4) + t29 * pkin(3) + t40 * pkin(2) + t41 * t91 + t93 * (t29 * t53 - t57 * t63) + t19 * pkin(11), (t10 * t56 + t25 * t88) * r_i_i_C(1) + (-t10 * t52 + t25 * t87) * r_i_i_C(2) + t10 * pkin(4) + t24 * pkin(3) + t25 * t92 + t93 * (t24 * t53 + t25 * t82), t93 * t2 + t64 * (-t25 * t53 + t57 * t68), (t15 * t56 - t2 * t52) * r_i_i_C(1) + (-t15 * t52 - t2 * t56) * r_i_i_C(2), 0; 1, (t22 * t56 + t32 * t52) * r_i_i_C(1) + (-t22 * t52 + t32 * t56) * r_i_i_C(2) + t22 * pkin(4) + t38 * pkin(3) - pkin(11) * t90 + t93 * (t38 * t53 - t57 * t61) + (t59 * pkin(2) + (pkin(11) * t50 + pkin(10)) * t85) * t48, (t18 * t56 + t34 * t88) * r_i_i_C(1) + (-t18 * t52 + t34 * t87) * r_i_i_C(2) + t18 * pkin(4) + t33 * pkin(3) + t34 * t92 + t93 * (t33 * t53 + t34 * t82), t93 * t14 + t64 * (-t34 * t53 + t57 * t66), (-t14 * t52 + t23 * t56) * r_i_i_C(1) + (-t14 * t56 - t23 * t52) * r_i_i_C(2), 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:35
	% EndTime: 2019-10-09 23:23:36
	% DurationCPUTime: 0.96s
	% Computational Cost: add. (1404->161), mult. (4083->287), div. (0->0), fcn. (5397->18), ass. (0->106)
	t122 = r_i_i_C(3) + pkin(13);
	t64 = sin(pkin(8));
	t121 = pkin(11) * t64;
	t65 = sin(pkin(7));
	t120 = t65 * pkin(10);
	t119 = cos(qJ(4));
	t73 = sin(qJ(2));
	t76 = cos(qJ(3));
	t107 = t73 * t76;
	t72 = sin(qJ(3));
	t77 = cos(qJ(2));
	t108 = t72 * t77;
	t66 = sin(pkin(6));
	t68 = cos(pkin(7));
	t56 = (-t68 * t107 - t108) * t66;
	t118 = t56 * t64;
	t70 = sin(qJ(5));
	t117 = t64 * t70;
	t75 = cos(qJ(5));
	t116 = t64 * t75;
	t115 = t65 * t64;
	t67 = cos(pkin(8));
	t114 = t65 * t67;
	t113 = t65 * t73;
	t71 = sin(qJ(4));
	t112 = t67 * t71;
	t111 = t68 * t72;
	t110 = t68 * t76;
	t109 = t72 * t73;
	t106 = t76 * t77;
	t105 = cos(pkin(6));
	t104 = cos(pkin(14));
	t103 = sin(pkin(14));
	t102 = t66 * t113;
	t101 = t67 * t119;
	t100 = t66 * t104;
	t99 = t66 * t103;
	t98 = t105 * t65;
	t97 = t119 * t115;
	t96 = t65 * t100;
	t95 = t105 * t104;
	t94 = t105 * t103;
	t69 = sin(qJ(6));
	t74 = cos(qJ(6));
	t93 = t74 * r_i_i_C(1) - t69 * r_i_i_C(2) + pkin(5);
	t92 = t69 * r_i_i_C(1) + t74 * r_i_i_C(2) + pkin(12);
	t58 = -t103 * t73 + t77 * t95;
	t59 = t103 * t77 + t73 * t95;
	t46 = -t59 * t110 - t58 * t72;
	t37 = t59 * t114 - t46 * t64;
	t60 = -t104 * t73 - t77 * t94;
	t61 = t104 * t77 - t73 * t94;
	t48 = -t61 * t110 - t60 * t72;
	t38 = t61 * t114 - t48 * t64;
	t91 = -t68 * t100 - t58 * t65;
	t90 = -t60 * t65 + t68 * t99;
	t89 = -t66 * t77 * t65 + t105 * t68;
	t88 = t60 * t68 + t65 * t99;
	t87 = t91 * t64;
	t86 = t90 * t64;
	t85 = t89 * t64;
	t84 = -t122 * t70 - t93 * t75 - pkin(4);
	t83 = t59 * t72 - (t58 * t68 - t96) * t76;
	t82 = t61 * t72 - t88 * t76;
	t81 = t76 * t98 + (t68 * t106 - t109) * t66;
	t80 = t83 * t119;
	t79 = t82 * t119;
	t78 = t81 * t119;
	t57 = (-t68 * t109 + t106) * t66;
	t54 = t72 * t98 + (t68 * t108 + t107) * t66;
	t50 = t67 * t102 - t118;
	t49 = -t61 * t111 + t60 * t76;
	t47 = -t59 * t111 + t58 * t76;
	t45 = t61 * t76 + t88 * t72;
	t44 = t58 * t111 + t59 * t76 - t72 * t96;
	t43 = -t81 * t64 + t89 * t67;
	t40 = t57 * t119 + (t64 * t102 + t56 * t67) * t71;
	t39 = -t66 * t73 * t97 - t56 * t101 + t57 * t71;
	t36 = -t54 * t112 + t78;
	t35 = t54 * t101 + t81 * t71;
	t34 = t82 * t64 + t90 * t67;
	t33 = t83 * t64 + t91 * t67;
	t32 = t54 * t119 + (t81 * t67 + t85) * t71;
	t31 = -t119 * t85 + t54 * t71 - t67 * t78;
	t30 = t54 * t117 + t36 * t75;
	t28 = t40 * t75 + t50 * t70;
	t26 = -t45 * t112 - t79;
	t25 = t45 * t101 - t82 * t71;
	t24 = -t44 * t112 - t80;
	t23 = t44 * t101 - t83 * t71;
	t22 = t49 * t119 + (t61 * t115 + t48 * t67) * t71;
	t21 = -t48 * t101 + t49 * t71 - t61 * t97;
	t20 = t47 * t119 + (t59 * t115 + t46 * t67) * t71;
	t19 = -t46 * t101 + t47 * t71 - t59 * t97;
	t18 = t45 * t119 + (-t82 * t67 + t86) * t71;
	t17 = -t119 * t86 + t45 * t71 + t67 * t79;
	t16 = t44 * t119 + (-t83 * t67 + t87) * t71;
	t15 = -t119 * t87 + t44 * t71 + t67 * t80;
	t14 = t32 * t75 + t43 * t70;
	t12 = t45 * t117 + t26 * t75;
	t10 = t44 * t117 + t24 * t75;
	t8 = t22 * t75 + t38 * t70;
	t6 = t20 * t75 + t37 * t70;
	t4 = t18 * t75 + t34 * t70;
	t2 = t16 * t75 + t33 * t70;
	t1 = [0, (t21 * t69 + t8 * t74) * r_i_i_C(1) + (t21 * t74 - t8 * t69) * r_i_i_C(2) + t8 * pkin(5) + t22 * pkin(4) + t21 * pkin(12) + t49 * pkin(3) + t60 * pkin(2) + t61 * t120 + t122 * (t22 * t70 - t38 * t75) + t38 * pkin(11), (t12 * t74 + t25 * t69) * r_i_i_C(1) + (-t12 * t69 + t25 * t74) * r_i_i_C(2) + t12 * pkin(5) + t26 * pkin(4) + t25 * pkin(12) - t82 * pkin(3) + t45 * t121 + t122 * (-t45 * t116 + t26 * t70), t84 * t17 + t92 * t18, t122 * t4 + t93 * (-t18 * t70 + t34 * t75), (t17 * t74 - t4 * t69) * r_i_i_C(1) + (-t17 * t69 - t4 * t74) * r_i_i_C(2); 0, (t19 * t69 + t6 * t74) * r_i_i_C(1) + (t19 * t74 - t6 * t69) * r_i_i_C(2) + t6 * pkin(5) + t20 * pkin(4) + t19 * pkin(12) + t47 * pkin(3) + t58 * pkin(2) + t59 * t120 + t122 * (t20 * t70 - t37 * t75) + t37 * pkin(11), (t10 * t74 + t23 * t69) * r_i_i_C(1) + (-t10 * t69 + t23 * t74) * r_i_i_C(2) + t10 * pkin(5) + t24 * pkin(4) + t23 * pkin(12) - t83 * pkin(3) + t44 * t121 + t122 * (-t44 * t116 + t24 * t70), t84 * t15 + t92 * t16, t122 * t2 + t93 * (-t16 * t70 + t33 * t75), (t15 * t74 - t2 * t69) * r_i_i_C(1) + (-t15 * t69 - t2 * t74) * r_i_i_C(2); 1, (t28 * t74 + t39 * t69) * r_i_i_C(1) + (-t28 * t69 + t39 * t74) * r_i_i_C(2) + t28 * pkin(5) + t40 * pkin(4) + t39 * pkin(12) + t57 * pkin(3) - pkin(11) * t118 + t122 * (t40 * t70 - t50 * t75) + (t77 * pkin(2) + (pkin(11) * t67 + pkin(10)) * t113) * t66, (t30 * t74 + t35 * t69) * r_i_i_C(1) + (-t30 * t69 + t35 * t74) * r_i_i_C(2) + t30 * pkin(5) + t36 * pkin(4) + t35 * pkin(12) + t81 * pkin(3) + t54 * t121 + t122 * (-t54 * t116 + t36 * t70), t84 * t31 + t92 * t32, t122 * t14 + t93 * (-t32 * t70 + t43 * t75), (-t14 * t69 + t31 * t74) * r_i_i_C(1) + (-t14 * t74 - t31 * t69) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end