% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:05
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRR7_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR7_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR7_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobia_transl_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (7->7), mult. (20->17), div. (0->0), fcn. (24->6), ass. (0->8)
	t4 = cos(pkin(6));
	t5 = sin(qJ(2));
	t8 = t4 * t5;
	t6 = cos(qJ(2));
	t7 = t4 * t6;
	t3 = cos(pkin(13));
	t1 = sin(pkin(13));
	t2 = [0, (-t1 * t7 - t3 * t5) * r_i_i_C(1) + (t1 * t8 - t3 * t6) * r_i_i_C(2), 0, 0, 0, 0; 0, (-t1 * t5 + t3 * t7) * r_i_i_C(1) + (-t1 * t6 - t3 * t8) * r_i_i_C(2), 0, 0, 0, 0; 1, (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * sin(pkin(6)), 0, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (37->15), mult. (103->31), div. (0->0), fcn. (131->10), ass. (0->18)
	t11 = cos(pkin(7));
	t5 = sin(pkin(14));
	t7 = sin(pkin(7));
	t9 = cos(pkin(14));
	t15 = (r_i_i_C(1) * t5 + r_i_i_C(2) * t9) * t11 - (r_i_i_C(3) + qJ(3)) * t7;
	t8 = sin(pkin(6));
	t21 = t11 * t8;
	t12 = cos(pkin(6));
	t13 = sin(qJ(2));
	t20 = t12 * t13;
	t14 = cos(qJ(2));
	t19 = t12 * t14;
	t16 = t9 * r_i_i_C(1) - t5 * r_i_i_C(2) + pkin(2);
	t10 = cos(pkin(13));
	t6 = sin(pkin(13));
	t3 = -t10 * t13 - t6 * t19;
	t1 = t10 * t19 - t6 * t13;
	t2 = [0, t16 * t3 + t15 * (-t10 * t14 + t6 * t20), t6 * t21 - t3 * t7, 0, 0, 0; 0, t16 * t1 + t15 * (-t10 * t20 - t14 * t6), -t1 * t7 - t10 * t21, 0, 0, 0; 1, (-t15 * t13 + t16 * t14) * t8, -t14 * t7 * t8 + t12 * t11, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (167->63), mult. (489->125), div. (0->0), fcn. (635->14), ass. (0->50)
	t57 = pkin(10) + r_i_i_C(3);
	t21 = sin(pkin(14));
	t28 = cos(pkin(7));
	t56 = t21 * t28;
	t22 = sin(pkin(8));
	t23 = sin(pkin(7));
	t55 = t22 * t23;
	t27 = cos(pkin(8));
	t54 = t23 * t27;
	t24 = sin(pkin(6));
	t26 = cos(pkin(13));
	t53 = t24 * t26;
	t25 = cos(pkin(14));
	t52 = t25 * t28;
	t30 = sin(qJ(2));
	t51 = t28 * t30;
	t32 = cos(qJ(2));
	t50 = t28 * t32;
	t49 = t23 * qJ(3);
	t48 = cos(pkin(6));
	t47 = sin(pkin(13));
	t46 = t24 * t47;
	t45 = t26 * t48;
	t44 = t48 * t23;
	t29 = sin(qJ(4));
	t31 = cos(qJ(4));
	t43 = t29 * r_i_i_C(1) + t31 * r_i_i_C(2);
	t16 = -t47 * t30 + t32 * t45;
	t11 = -t16 * t23 - t28 * t53;
	t17 = t30 * t45 + t47 * t32;
	t34 = t16 * t28 - t23 * t53;
	t42 = (-t17 * t21 + t34 * t25) * t27 + t11 * t22;
	t39 = t48 * t47;
	t18 = -t26 * t30 - t32 * t39;
	t12 = -t18 * t23 + t28 * t46;
	t19 = t26 * t32 - t30 * t39;
	t33 = t18 * t28 + t23 * t46;
	t41 = t12 * t22 + t27 * (-t19 * t21 + t33 * t25);
	t15 = -t24 * t32 * t23 + t48 * t28;
	t40 = t15 * t22 + t27 * (t25 * t44 + (-t21 * t30 + t25 * t50) * t24);
	t5 = -t16 * t21 - t17 * t52;
	t37 = t17 * t55 + t27 * t5;
	t7 = -t18 * t21 - t19 * t52;
	t35 = t19 * t55 + t27 * t7;
	t10 = t24 * t30 * t25 + (t24 * t50 + t44) * t21;
	t8 = t18 * t25 - t19 * t56;
	t6 = t16 * t25 - t17 * t56;
	t4 = t19 * t25 + t33 * t21;
	t2 = t17 * t25 + t34 * t21;
	t1 = [0, (t35 * t29 + t8 * t31) * r_i_i_C(1) + (-t8 * t29 + t35 * t31) * r_i_i_C(2) + t8 * pkin(3) + t18 * pkin(2) + t19 * t49 + t57 * (t19 * t54 - t7 * t22), t12, (-t4 * t29 + t41 * t31) * r_i_i_C(1) + (-t41 * t29 - t4 * t31) * r_i_i_C(2), 0, 0; 0, (t37 * t29 + t6 * t31) * r_i_i_C(1) + (-t6 * t29 + t37 * t31) * r_i_i_C(2) + t6 * pkin(3) + t16 * pkin(2) + t17 * t49 + t57 * (t17 * t54 - t5 * t22), t11, (-t2 * t29 + t42 * t31) * r_i_i_C(1) + (-t2 * t31 - t42 * t29) * r_i_i_C(2), 0, 0; 1, ((t31 * r_i_i_C(1) - t29 * r_i_i_C(2) + pkin(3)) * (-t21 * t51 + t25 * t32) + (-t57 * t22 + t43 * t27) * (-t21 * t32 - t25 * t51) + t32 * pkin(2) + (t43 * t22 + t57 * t27 + qJ(3)) * t30 * t23) * t24, t15, (-t10 * t29 + t40 * t31) * r_i_i_C(1) + (-t10 * t31 - t40 * t29) * r_i_i_C(2), 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.44s
	% Computational Cost: add. (466->90), mult. (1356->175), div. (0->0), fcn. (1783->16), ass. (0->73)
	t79 = r_i_i_C(3) + pkin(11);
	t39 = sin(pkin(14));
	t42 = sin(pkin(6));
	t43 = cos(pkin(14));
	t52 = cos(qJ(2));
	t46 = cos(pkin(7));
	t49 = sin(qJ(2));
	t71 = t46 * t49;
	t31 = (-t39 * t52 - t43 * t71) * t42;
	t40 = sin(pkin(8));
	t78 = t31 * t40;
	t77 = t39 * t46;
	t41 = sin(pkin(7));
	t76 = t40 * t41;
	t45 = cos(pkin(8));
	t75 = t41 * t45;
	t44 = cos(pkin(13));
	t74 = t42 * t44;
	t73 = t42 * t49;
	t72 = t43 * t46;
	t70 = t46 * t52;
	t69 = t41 * qJ(3);
	t68 = cos(pkin(6));
	t67 = sin(pkin(13));
	t66 = t41 * t73;
	t65 = t42 * t67;
	t64 = t44 * t68;
	t63 = t68 * t41;
	t35 = t49 * t64 + t67 * t52;
	t34 = -t67 * t49 + t52 * t64;
	t55 = t34 * t46 - t41 * t74;
	t18 = -t35 * t39 + t55 * t43;
	t29 = -t34 * t41 - t46 * t74;
	t62 = t18 * t45 + t29 * t40;
	t59 = t68 * t67;
	t37 = t44 * t52 - t49 * t59;
	t36 = -t44 * t49 - t52 * t59;
	t53 = t36 * t46 + t41 * t65;
	t20 = -t37 * t39 + t53 * t43;
	t30 = -t36 * t41 + t46 * t65;
	t61 = t20 * t45 + t30 * t40;
	t27 = t43 * t63 + (-t39 * t49 + t43 * t70) * t42;
	t33 = -t42 * t52 * t41 + t68 * t46;
	t60 = t27 * t45 + t33 * t40;
	t47 = sin(qJ(5));
	t50 = cos(qJ(5));
	t58 = t50 * r_i_i_C(1) - t47 * r_i_i_C(2) + pkin(4);
	t22 = -t34 * t39 - t35 * t72;
	t13 = -t22 * t40 + t35 * t75;
	t57 = t22 * t45 + t35 * t76;
	t24 = -t36 * t39 - t37 * t72;
	t14 = -t24 * t40 + t37 * t75;
	t56 = t24 * t45 + t37 * t76;
	t54 = t31 * t45 + t40 * t66;
	t51 = cos(qJ(4));
	t48 = sin(qJ(4));
	t32 = (-t39 * t71 + t43 * t52) * t42;
	t28 = t43 * t73 + (t42 * t70 + t63) * t39;
	t26 = t45 * t66 - t78;
	t25 = t36 * t43 - t37 * t77;
	t23 = t34 * t43 - t35 * t77;
	t21 = t37 * t43 + t53 * t39;
	t19 = t35 * t43 + t55 * t39;
	t17 = -t27 * t40 + t33 * t45;
	t16 = t32 * t51 + t54 * t48;
	t12 = -t20 * t40 + t30 * t45;
	t11 = -t18 * t40 + t29 * t45;
	t10 = t28 * t51 + t60 * t48;
	t8 = t25 * t51 + t56 * t48;
	t6 = t23 * t51 + t57 * t48;
	t4 = t21 * t51 + t61 * t48;
	t2 = t19 * t51 + t62 * t48;
	t1 = [0, (t14 * t47 + t8 * t50) * r_i_i_C(1) + (t14 * t50 - t8 * t47) * r_i_i_C(2) + t8 * pkin(4) + t25 * pkin(3) + t36 * pkin(2) + t37 * t69 + t79 * (t25 * t48 - t56 * t51) + t14 * pkin(10), t30, t79 * t4 + t58 * (-t21 * t48 + t61 * t51), (t12 * t50 - t4 * t47) * r_i_i_C(1) + (-t12 * t47 - t4 * t50) * r_i_i_C(2), 0; 0, (t13 * t47 + t6 * t50) * r_i_i_C(1) + (t13 * t50 - t6 * t47) * r_i_i_C(2) + t6 * pkin(4) + t23 * pkin(3) + t34 * pkin(2) + t35 * t69 + t79 * (t23 * t48 - t57 * t51) + t13 * pkin(10), t29, t79 * t2 + t58 * (-t19 * t48 + t62 * t51), (t11 * t50 - t2 * t47) * r_i_i_C(1) + (-t11 * t47 - t2 * t50) * r_i_i_C(2), 0; 1, (t16 * t50 + t26 * t47) * r_i_i_C(1) + (-t16 * t47 + t26 * t50) * r_i_i_C(2) + t16 * pkin(4) + t32 * pkin(3) - pkin(10) * t78 + t79 * (t32 * t48 - t54 * t51) + (t52 * pkin(2) + (pkin(10) * t45 + qJ(3)) * t49 * t41) * t42, t33, t79 * t10 + t58 * (-t28 * t48 + t60 * t51), (-t10 * t47 + t17 * t50) * r_i_i_C(1) + (-t10 * t50 - t17 * t47) * r_i_i_C(2), 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.72s
	% Computational Cost: add. (1087->120), mult. (3155->222), div. (0->0), fcn. (4179->18), ass. (0->89)
	t50 = sin(pkin(8));
	t54 = cos(pkin(8));
	t53 = cos(pkin(13));
	t58 = sin(qJ(2));
	t61 = cos(qJ(2));
	t93 = sin(pkin(13));
	t96 = cos(pkin(6));
	t82 = t96 * t93;
	t47 = t53 * t61 - t58 * t82;
	t49 = sin(pkin(14));
	t76 = t53 * t58 + t61 * t82;
	t94 = cos(pkin(14));
	t73 = t76 * t94;
	t51 = sin(pkin(7));
	t52 = sin(pkin(6));
	t87 = t52 * t93;
	t83 = t51 * t87;
	t95 = cos(pkin(7));
	t63 = t47 * t49 + t95 * t73 - t94 * t83;
	t69 = t76 * t51 + t95 * t87;
	t107 = -t69 * t50 + t63 * t54;
	t86 = t53 * t96;
	t46 = t58 * t86 + t93 * t61;
	t77 = t93 * t58 - t61 * t86;
	t72 = t77 * t94;
	t88 = t52 * t94;
	t65 = t53 * t51 * t88 + t46 * t49 + t95 * t72;
	t89 = t52 * t95;
	t71 = t77 * t51 - t53 * t89;
	t106 = -t71 * t50 + t65 * t54;
	t81 = t95 * t94;
	t85 = t96 * t51;
	t68 = t94 * t85 + (-t49 * t58 + t61 * t81) * t52;
	t100 = t51 * t52;
	t78 = -t61 * t100 + t96 * t95;
	t105 = t78 * t50 + t68 * t54;
	t104 = r_i_i_C(3) + pkin(12);
	t103 = cos(qJ(4));
	t44 = (-t49 * t61 - t58 * t81) * t52;
	t102 = t44 * t50;
	t101 = t50 * t51;
	t99 = t51 * t54;
	t98 = t58 * t51;
	t97 = t51 * qJ(3);
	t92 = t52 * t98;
	t91 = t54 * t103;
	t90 = t49 * t95;
	t84 = t103 * t101;
	t55 = sin(qJ(6));
	t59 = cos(qJ(6));
	t80 = t59 * r_i_i_C(1) - t55 * r_i_i_C(2) + pkin(5);
	t79 = t55 * r_i_i_C(1) + t59 * r_i_i_C(2) + pkin(11);
	t34 = -t46 * t81 + t77 * t49;
	t25 = -t34 * t50 + t46 * t99;
	t36 = -t47 * t81 + t76 * t49;
	t26 = -t36 * t50 + t47 * t99;
	t56 = sin(qJ(5));
	t60 = cos(qJ(5));
	t74 = -t104 * t56 - t80 * t60 - pkin(4);
	t57 = sin(qJ(4));
	t45 = (-t58 * t90 + t94 * t61) * t52;
	t42 = t58 * t88 + (t61 * t89 + t85) * t49;
	t38 = t54 * t92 - t102;
	t37 = -t47 * t90 - t73;
	t35 = -t46 * t90 - t72;
	t33 = t47 * t94 + (-t76 * t95 + t83) * t49;
	t32 = t46 * t94 + (-t53 * t100 - t77 * t95) * t49;
	t31 = -t68 * t50 + t78 * t54;
	t28 = t45 * t103 + (t44 * t54 + t50 * t92) * t57;
	t27 = -t52 * t58 * t84 - t44 * t91 + t45 * t57;
	t24 = t63 * t50 + t69 * t54;
	t23 = t65 * t50 + t71 * t54;
	t22 = t42 * t103 + t105 * t57;
	t21 = -t105 * t103 + t42 * t57;
	t20 = t28 * t60 + t38 * t56;
	t18 = t37 * t103 + (t47 * t101 + t36 * t54) * t57;
	t17 = -t36 * t91 + t37 * t57 - t47 * t84;
	t16 = t35 * t103 + (t46 * t101 + t34 * t54) * t57;
	t15 = -t34 * t91 + t35 * t57 - t46 * t84;
	t14 = t33 * t103 - t107 * t57;
	t13 = t107 * t103 + t33 * t57;
	t12 = t32 * t103 - t106 * t57;
	t11 = t106 * t103 + t32 * t57;
	t10 = t22 * t60 + t31 * t56;
	t8 = t18 * t60 + t26 * t56;
	t6 = t16 * t60 + t25 * t56;
	t4 = t14 * t60 + t24 * t56;
	t2 = t12 * t60 + t23 * t56;
	t1 = [0, (t17 * t55 + t8 * t59) * r_i_i_C(1) + (t17 * t59 - t8 * t55) * r_i_i_C(2) + t8 * pkin(5) + t18 * pkin(4) + t17 * pkin(11) + t37 * pkin(3) - t76 * pkin(2) + t47 * t97 + t104 * (t18 * t56 - t26 * t60) + t26 * pkin(10), t69, t74 * t13 + t79 * t14, t104 * t4 + t80 * (-t14 * t56 + t24 * t60), (t13 * t59 - t4 * t55) * r_i_i_C(1) + (-t13 * t55 - t4 * t59) * r_i_i_C(2); 0, (t15 * t55 + t6 * t59) * r_i_i_C(1) + (t15 * t59 - t6 * t55) * r_i_i_C(2) + t6 * pkin(5) + t16 * pkin(4) + t15 * pkin(11) + t35 * pkin(3) - t77 * pkin(2) + t46 * t97 + t104 * (t16 * t56 - t25 * t60) + t25 * pkin(10), t71, t74 * t11 + t79 * t12, t104 * t2 + t80 * (-t12 * t56 + t23 * t60), (t11 * t59 - t2 * t55) * r_i_i_C(1) + (-t11 * t55 - t2 * t59) * r_i_i_C(2); 1, (t20 * t59 + t27 * t55) * r_i_i_C(1) + (-t20 * t55 + t27 * t59) * r_i_i_C(2) + t20 * pkin(5) + t28 * pkin(4) + t27 * pkin(11) + t45 * pkin(3) - pkin(10) * t102 + t104 * (t28 * t56 - t38 * t60) + (t61 * pkin(2) + (pkin(10) * t54 + qJ(3)) * t98) * t52, t78, t74 * t21 + t79 * t22, t104 * t10 + t80 * (-t22 * t56 + t31 * t60), (-t10 * t55 + t21 * t59) * r_i_i_C(1) + (-t10 * t59 - t21 * t55) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end