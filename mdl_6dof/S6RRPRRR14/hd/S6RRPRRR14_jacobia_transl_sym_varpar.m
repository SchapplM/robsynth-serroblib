% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:11
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR14_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR14_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobia_transl_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:58
	% EndTime: 2019-10-10 11:10:58
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:58
	% EndTime: 2019-10-10 11:10:58
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:58
	% EndTime: 2019-10-10 11:10:58
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (19->14), mult. (46->24), div. (0->0), fcn. (54->6), ass. (0->16)
	t7 = sin(qJ(2));
	t8 = sin(qJ(1));
	t15 = t8 * t7;
	t9 = cos(qJ(2));
	t14 = t8 * t9;
	t10 = cos(qJ(1));
	t13 = t10 * t7;
	t12 = t10 * t9;
	t5 = sin(pkin(6));
	t11 = (pkin(10) + r_i_i_C(3)) * t5;
	t6 = cos(pkin(6));
	t4 = -t6 * t15 + t12;
	t3 = -t6 * t14 - t13;
	t2 = -t6 * t13 - t14;
	t1 = -t6 * t12 + t15;
	t16 = [-t8 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t10 * t11, t3 * r_i_i_C(1) - t4 * r_i_i_C(2), 0, 0, 0, 0; t10 * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t11 * t8, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0, 0, 0, 0; 0, (r_i_i_C(1) * t9 - r_i_i_C(2) * t7) * t5, 0, 0, 0, 0;];
	Ja_transl = t16;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:58
	% EndTime: 2019-10-10 11:10:58
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (75->32), mult. (199->56), div. (0->0), fcn. (251->10), ass. (0->26)
	t26 = r_i_i_C(3) + qJ(3);
	t11 = cos(pkin(14));
	t12 = cos(pkin(7));
	t8 = sin(pkin(14));
	t9 = sin(pkin(7));
	t29 = (r_i_i_C(1) * t8 + r_i_i_C(2) * t11) * t12 - t26 * t9;
	t10 = sin(pkin(6));
	t14 = sin(qJ(1));
	t28 = t10 * t14;
	t16 = cos(qJ(1));
	t27 = t10 * t16;
	t25 = cos(pkin(6));
	t24 = t14 * t25;
	t23 = t16 * t25;
	t21 = r_i_i_C(1) * t11 - r_i_i_C(2) * t8 + pkin(2);
	t13 = sin(qJ(2));
	t15 = cos(qJ(2));
	t4 = -t16 * t13 - t15 * t24;
	t20 = t12 * t4 + t9 * t28;
	t1 = t12 * t28 - t4 * t9;
	t2 = t13 * t14 - t15 * t23;
	t19 = t12 * t2 + t9 * t27;
	t18 = t12 * t27 - t2 * t9;
	t5 = -t13 * t24 + t15 * t16;
	t3 = -t13 * t23 - t14 * t15;
	t6 = [(t3 * t11 + t19 * t8) * r_i_i_C(1) + (t19 * t11 - t3 * t8) * r_i_i_C(2) + t3 * pkin(2) - t14 * pkin(1) + pkin(10) * t27 + t26 * t18, t21 * t4 - t29 * t5, t1, 0, 0, 0; (t11 * t5 + t20 * t8) * r_i_i_C(1) + (t20 * t11 - t5 * t8) * r_i_i_C(2) + t5 * pkin(2) + t16 * pkin(1) + pkin(10) * t28 + t26 * t1, -t21 * t2 + t29 * t3, -t18, 0, 0, 0; 0, (-t13 * t29 + t21 * t15) * t10, -t10 * t15 * t9 + t25 * t12, 0, 0, 0;];
	Ja_transl = t6;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:58
	% EndTime: 2019-10-10 11:10:58
	% DurationCPUTime: 0.58s
	% Computational Cost: add. (273->80), mult. (775->145), div. (0->0), fcn. (1007->14), ass. (0->56)
	t34 = sin(qJ(4));
	t37 = cos(qJ(4));
	t35 = sin(qJ(2));
	t36 = sin(qJ(1));
	t38 = cos(qJ(2));
	t39 = cos(qJ(1));
	t56 = cos(pkin(6));
	t54 = t39 * t56;
	t21 = t36 * t35 - t38 * t54;
	t29 = sin(pkin(7));
	t33 = cos(pkin(7));
	t30 = sin(pkin(6));
	t61 = t30 * t39;
	t15 = -t21 * t29 + t33 * t61;
	t28 = sin(pkin(8));
	t32 = cos(pkin(8));
	t22 = t35 * t54 + t36 * t38;
	t27 = sin(pkin(14));
	t31 = cos(pkin(14));
	t41 = t21 * t33 + t29 * t61;
	t5 = t22 * t27 + t31 * t41;
	t50 = t15 * t28 + t32 * t5;
	t6 = -t22 * t31 + t27 * t41;
	t73 = t50 * t34 + t6 * t37;
	t72 = -t6 * t34 + t50 * t37;
	t68 = pkin(11) + r_i_i_C(3);
	t65 = t27 * t33;
	t64 = t28 * t29;
	t63 = t29 * t32;
	t62 = t30 * t36;
	t60 = t31 * t33;
	t59 = t33 * t35;
	t58 = t33 * t38;
	t57 = t29 * qJ(3);
	t55 = t36 * t56;
	t53 = t56 * t29;
	t52 = t34 * r_i_i_C(1) + t37 * r_i_i_C(2);
	t23 = -t39 * t35 - t38 * t55;
	t17 = -t23 * t29 + t33 * t62;
	t24 = -t35 * t55 + t39 * t38;
	t40 = t23 * t33 + t29 * t62;
	t7 = -t24 * t27 + t31 * t40;
	t47 = t17 * t28 + t32 * t7;
	t20 = -t30 * t38 * t29 + t56 * t33;
	t46 = (t31 * t53 + (-t27 * t35 + t31 * t58) * t30) * t32 + t20 * t28;
	t9 = t21 * t27 - t22 * t60;
	t44 = t22 * t64 + t32 * t9;
	t11 = -t23 * t27 - t24 * t60;
	t42 = t11 * t32 + t24 * t64;
	t14 = t30 * t35 * t31 + (t30 * t58 + t53) * t27;
	t12 = t23 * t31 - t24 * t65;
	t10 = -t21 * t31 - t22 * t65;
	t8 = t24 * t31 + t27 * t40;
	t2 = t47 * t34 + t8 * t37;
	t1 = -t8 * t34 + t47 * t37;
	t3 = [t73 * r_i_i_C(1) + t72 * r_i_i_C(2) + t6 * pkin(3) - t22 * pkin(2) - t36 * pkin(1) + pkin(10) * t61 + t15 * qJ(3) + t68 * (t15 * t32 - t5 * t28), (t12 * t37 + t42 * t34) * r_i_i_C(1) + (-t12 * t34 + t42 * t37) * r_i_i_C(2) + t12 * pkin(3) + t23 * pkin(2) + t24 * t57 + t68 * (-t11 * t28 + t24 * t63), t17, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0; t39 * pkin(1) + t24 * pkin(2) + t8 * pkin(3) + pkin(10) * t62 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t17 * qJ(3) + t68 * (t17 * t32 - t7 * t28), (t10 * t37 + t44 * t34) * r_i_i_C(1) + (-t10 * t34 + t44 * t37) * r_i_i_C(2) + t10 * pkin(3) - t21 * pkin(2) + t22 * t57 + t68 * (t22 * t63 - t9 * t28), -t15, -t72 * r_i_i_C(1) + t73 * r_i_i_C(2), 0, 0; 0, ((t37 * r_i_i_C(1) - t34 * r_i_i_C(2) + pkin(3)) * (-t27 * t59 + t31 * t38) + (-t68 * t28 + t52 * t32) * (-t27 * t38 - t31 * t59) + t38 * pkin(2) + (t52 * t28 + t68 * t32 + qJ(3)) * t35 * t29) * t30, t20, (-t14 * t34 + t46 * t37) * r_i_i_C(1) + (-t14 * t37 - t46 * t34) * r_i_i_C(2), 0, 0;];
	Ja_transl = t3;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:58
	% EndTime: 2019-10-10 11:10:59
	% DurationCPUTime: 0.95s
	% Computational Cost: add. (666->109), mult. (1904->195), div. (0->0), fcn. (2507->16), ass. (0->79)
	t58 = sin(qJ(2));
	t59 = sin(qJ(1));
	t62 = cos(qJ(2));
	t63 = cos(qJ(1));
	t79 = cos(pkin(6));
	t76 = t63 * t79;
	t44 = t58 * t76 + t59 * t62;
	t49 = sin(pkin(14));
	t53 = cos(pkin(14));
	t43 = t59 * t58 - t62 * t76;
	t51 = sin(pkin(7));
	t55 = cos(pkin(7));
	t52 = sin(pkin(6));
	t84 = t52 * t63;
	t68 = t43 * t55 + t51 * t84;
	t26 = t44 * t49 + t68 * t53;
	t37 = -t43 * t51 + t55 * t84;
	t50 = sin(pkin(8));
	t54 = cos(pkin(8));
	t15 = t26 * t50 - t37 * t54;
	t56 = sin(qJ(5));
	t27 = -t44 * t53 + t68 * t49;
	t57 = sin(qJ(4));
	t61 = cos(qJ(4));
	t73 = t26 * t54 + t37 * t50;
	t6 = t27 * t61 + t73 * t57;
	t60 = cos(qJ(5));
	t104 = t15 * t60 + t6 * t56;
	t103 = -t15 * t56 + t6 * t60;
	t100 = t27 * t57 - t73 * t61;
	t77 = t59 * t79;
	t45 = -t63 * t58 - t62 * t77;
	t85 = t52 * t59;
	t39 = -t45 * t51 + t55 * t85;
	t46 = -t58 * t77 + t63 * t62;
	t67 = t45 * t55 + t51 * t85;
	t65 = t46 * t49 - t67 * t53;
	t95 = -t39 * t50 + t65 * t54;
	t94 = r_i_i_C(3) + pkin(12);
	t82 = t55 * t58;
	t40 = (-t49 * t62 - t53 * t82) * t52;
	t92 = t40 * t50;
	t89 = t49 * t55;
	t88 = t50 * t51;
	t87 = t51 * t54;
	t86 = t52 * t58;
	t83 = t53 * t55;
	t81 = t55 * t62;
	t80 = t51 * qJ(3);
	t78 = t51 * t86;
	t75 = t79 * t51;
	t35 = t53 * t75 + (-t49 * t58 + t53 * t81) * t52;
	t42 = -t52 * t62 * t51 + t79 * t55;
	t72 = t35 * t54 + t42 * t50;
	t71 = t60 * r_i_i_C(1) - t56 * r_i_i_C(2) + pkin(4);
	t29 = t43 * t49 - t44 * t83;
	t20 = -t29 * t50 + t44 * t87;
	t70 = t29 * t54 + t44 * t88;
	t31 = -t45 * t49 - t46 * t83;
	t21 = -t31 * t50 + t46 * t87;
	t69 = t31 * t54 + t46 * t88;
	t66 = t40 * t54 + t50 * t78;
	t17 = t39 * t54 + t65 * t50;
	t41 = (-t49 * t82 + t53 * t62) * t52;
	t36 = t53 * t86 + (t52 * t81 + t75) * t49;
	t33 = t54 * t78 - t92;
	t32 = t45 * t53 - t46 * t89;
	t30 = -t43 * t53 - t44 * t89;
	t28 = t46 * t53 + t67 * t49;
	t23 = -t35 * t50 + t42 * t54;
	t19 = t41 * t61 + t66 * t57;
	t14 = t36 * t61 + t72 * t57;
	t12 = t32 * t61 + t69 * t57;
	t10 = t30 * t61 + t70 * t57;
	t8 = t28 * t61 - t95 * t57;
	t7 = t28 * t57 + t95 * t61;
	t2 = t17 * t56 + t8 * t60;
	t1 = t17 * t60 - t8 * t56;
	t3 = [-t59 * pkin(1) - t44 * pkin(2) + t27 * pkin(3) + t6 * pkin(4) + pkin(10) * t84 - t15 * pkin(11) + t103 * r_i_i_C(1) - t104 * r_i_i_C(2) + t37 * qJ(3) + t94 * t100, (t12 * t60 + t21 * t56) * r_i_i_C(1) + (-t12 * t56 + t21 * t60) * r_i_i_C(2) + t12 * pkin(4) + t32 * pkin(3) + t45 * pkin(2) + t46 * t80 + t94 * (t32 * t57 - t69 * t61) + t21 * pkin(11), t39, -t71 * t7 + t94 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t63 * pkin(1) + t46 * pkin(2) + t28 * pkin(3) + t8 * pkin(4) + pkin(10) * t85 + t17 * pkin(11) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t39 * qJ(3) + t94 * t7, (t10 * t60 + t20 * t56) * r_i_i_C(1) + (-t10 * t56 + t20 * t60) * r_i_i_C(2) + t10 * pkin(4) + t30 * pkin(3) - t43 * pkin(2) + t44 * t80 + t94 * (t30 * t57 - t70 * t61) + t20 * pkin(11), -t37, t71 * t100 - t6 * t94, t104 * r_i_i_C(1) + t103 * r_i_i_C(2), 0; 0, (t19 * t60 + t33 * t56) * r_i_i_C(1) + (-t19 * t56 + t33 * t60) * r_i_i_C(2) + t19 * pkin(4) + t41 * pkin(3) - pkin(11) * t92 + t94 * (t41 * t57 - t66 * t61) + (t62 * pkin(2) + (pkin(11) * t54 + qJ(3)) * t58 * t51) * t52, t42, t94 * t14 + t71 * (-t36 * t57 + t72 * t61), (-t14 * t56 + t23 * t60) * r_i_i_C(1) + (-t14 * t60 - t23 * t56) * r_i_i_C(2), 0;];
	Ja_transl = t3;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:58
	% EndTime: 2019-10-10 11:10:59
	% DurationCPUTime: 1.43s
	% Computational Cost: add. (1439->146), mult. (4127->251), div. (0->0), fcn. (5475->18), ass. (0->98)
	t118 = cos(qJ(4));
	t116 = sin(qJ(2));
	t117 = sin(qJ(1));
	t71 = cos(qJ(1));
	t107 = cos(pkin(6));
	t70 = cos(qJ(2));
	t97 = t70 * t107;
	t56 = t117 * t116 - t71 * t97;
	t92 = t107 * t116;
	t57 = t117 * t70 + t71 * t92;
	t60 = sin(pkin(14));
	t105 = cos(pkin(14));
	t106 = cos(pkin(7));
	t90 = t106 * t105;
	t62 = sin(pkin(7));
	t63 = sin(pkin(6));
	t98 = t63 * t105;
	t93 = t62 * t98;
	t40 = t56 * t90 + t57 * t60 + t71 * t93;
	t61 = sin(pkin(8));
	t64 = cos(pkin(8));
	t99 = t63 * t106;
	t87 = t56 * t62 - t71 * t99;
	t120 = t40 * t64 - t87 * t61;
	t109 = t63 * t71;
	t41 = (t106 * t56 + t62 * t109) * t60 - t57 * t105;
	t67 = sin(qJ(4));
	t18 = t41 * t118 + t120 * t67;
	t30 = t40 * t61 + t87 * t64;
	t66 = sin(qJ(5));
	t69 = cos(qJ(5));
	t6 = t18 * t69 - t30 * t66;
	t65 = sin(qJ(6));
	t129 = t6 * t65;
	t68 = cos(qJ(6));
	t128 = t6 * t68;
	t127 = t18 * t66 + t30 * t69;
	t124 = t41 * t67;
	t58 = -t117 * t92 + t71 * t70;
	t83 = t71 * t116 + t117 * t97;
	t80 = t83 * t105;
	t74 = t106 * t80 - t117 * t93 + t58 * t60;
	t78 = t117 * t99 + t83 * t62;
	t122 = -t78 * t61 + t74 * t64;
	t96 = t107 * t62;
	t77 = t105 * t96 + (-t116 * t60 + t70 * t90) * t63;
	t111 = t62 * t63;
	t84 = t107 * t106 - t70 * t111;
	t121 = t84 * t61 + t77 * t64;
	t119 = r_i_i_C(3) + pkin(13);
	t54 = (-t116 * t90 - t60 * t70) * t63;
	t113 = t54 * t61;
	t112 = t62 * t61;
	t110 = t62 * t64;
	t108 = t62 * qJ(3);
	t104 = t61 * t118;
	t103 = t63 * t117;
	t102 = t64 * t118;
	t100 = t60 * t106;
	t95 = t62 * t104;
	t94 = t116 * t111;
	t91 = t61 * t94;
	t89 = t68 * r_i_i_C(1) - t65 * r_i_i_C(2) + pkin(5);
	t88 = t65 * r_i_i_C(1) + t68 * r_i_i_C(2) + pkin(12);
	t43 = t56 * t60 - t57 * t90;
	t34 = t57 * t110 - t43 * t61;
	t45 = -t58 * t90 + t83 * t60;
	t35 = t58 * t110 - t45 * t61;
	t81 = -t119 * t66 - t89 * t69 - pkin(4);
	t72 = t74 * t61 + t78 * t64;
	t55 = (-t116 * t100 + t105 * t70) * t63;
	t51 = t116 * t98 + (t70 * t99 + t96) * t60;
	t47 = t64 * t94 - t113;
	t46 = -t58 * t100 - t80;
	t44 = -t57 * t100 - t56 * t105;
	t42 = t58 * t105 + (t62 * t103 - t83 * t106) * t60;
	t38 = -t77 * t61 + t84 * t64;
	t33 = t55 * t118 + (t54 * t64 + t91) * t67;
	t32 = -t54 * t102 - t118 * t91 + t55 * t67;
	t28 = t51 * t118 + t121 * t67;
	t27 = -t121 * t118 + t51 * t67;
	t26 = t33 * t69 + t47 * t66;
	t24 = t46 * t118 + (t58 * t112 + t45 * t64) * t67;
	t23 = -t45 * t102 + t46 * t67 - t58 * t95;
	t22 = t44 * t118 + (t57 * t112 + t43 * t64) * t67;
	t21 = -t43 * t102 + t44 * t67 - t57 * t95;
	t20 = t42 * t118 - t122 * t67;
	t19 = t122 * t118 + t42 * t67;
	t17 = -t40 * t102 + t104 * t87 + t124;
	t15 = t120 * t118 - t124;
	t14 = t28 * t69 + t38 * t66;
	t12 = t24 * t69 + t35 * t66;
	t10 = t22 * t69 + t34 * t66;
	t8 = t20 * t69 + t72 * t66;
	t7 = t20 * t66 - t72 * t69;
	t2 = t19 * t65 + t8 * t68;
	t1 = t19 * t68 - t8 * t65;
	t3 = [(t17 * t65 + t128) * r_i_i_C(1) + (t17 * t68 - t129) * r_i_i_C(2) + t6 * pkin(5) + t18 * pkin(4) + t17 * pkin(12) + t41 * pkin(3) - t57 * pkin(2) - t117 * pkin(1) + pkin(10) * t109 + t119 * t127 - t87 * qJ(3) - t30 * pkin(11), (t12 * t68 + t23 * t65) * r_i_i_C(1) + (-t12 * t65 + t23 * t68) * r_i_i_C(2) + t12 * pkin(5) + t24 * pkin(4) + t23 * pkin(12) + t46 * pkin(3) - t83 * pkin(2) + t58 * t108 + t119 * (t24 * t66 - t35 * t69) + t35 * pkin(11), t78, t81 * t19 + t88 * t20, t119 * t8 - t89 * t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t71 * pkin(1) + t58 * pkin(2) + t42 * pkin(3) + t20 * pkin(4) + t8 * pkin(5) + pkin(10) * t103 + t72 * pkin(11) + t19 * pkin(12) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t78 * qJ(3) + t119 * t7, (t10 * t68 + t21 * t65) * r_i_i_C(1) + (-t10 * t65 + t21 * t68) * r_i_i_C(2) + t10 * pkin(5) + t22 * pkin(4) + t21 * pkin(12) + t44 * pkin(3) - t56 * pkin(2) + t57 * t108 + t119 * (t22 * t66 - t34 * t69) + t34 * pkin(11), t87, t81 * t15 - t18 * t88, -t119 * t6 + t89 * t127, (t15 * t68 + t129) * r_i_i_C(1) + (-t15 * t65 + t128) * r_i_i_C(2); 0, (t26 * t68 + t32 * t65) * r_i_i_C(1) + (-t26 * t65 + t32 * t68) * r_i_i_C(2) + t26 * pkin(5) + t33 * pkin(4) + t32 * pkin(12) + t55 * pkin(3) - pkin(11) * t113 + t119 * (t33 * t66 - t47 * t69) + (t70 * pkin(2) + (pkin(11) * t64 + qJ(3)) * t62 * t116) * t63, t84, t81 * t27 + t88 * t28, t119 * t14 + t89 * (-t28 * t66 + t38 * t69), (-t14 * t65 + t27 * t68) * r_i_i_C(1) + (-t14 * t68 - t27 * t65) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end