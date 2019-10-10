% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPR12
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:50
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR12_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR12_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR12_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_jacobia_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:49:58
	% EndTime: 2019-10-10 12:49:58
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:49:58
	% EndTime: 2019-10-10 12:49:58
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
	% StartTime: 2019-10-10 12:49:58
	% EndTime: 2019-10-10 12:49:58
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
	t11 = (pkin(9) + r_i_i_C(3)) * t5;
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
	% StartTime: 2019-10-10 12:49:58
	% EndTime: 2019-10-10 12:49:58
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (93->46), mult. (252->84), div. (0->0), fcn. (316->10), ass. (0->33)
	t39 = pkin(10) + r_i_i_C(3);
	t12 = cos(pkin(7));
	t13 = sin(qJ(3));
	t16 = cos(qJ(3));
	t10 = sin(pkin(7));
	t11 = sin(pkin(6));
	t18 = cos(qJ(1));
	t35 = t11 * t18;
	t27 = t10 * t35;
	t14 = sin(qJ(2));
	t15 = sin(qJ(1));
	t17 = cos(qJ(2));
	t28 = cos(pkin(6));
	t24 = t18 * t28;
	t3 = t15 * t14 - t17 * t24;
	t4 = t14 * t24 + t15 * t17;
	t38 = (t12 * t3 + t27) * t16 + t4 * t13;
	t36 = t11 * t15;
	t34 = t12 * t13;
	t33 = t12 * t16;
	t32 = t13 * t14;
	t31 = t13 * t17;
	t30 = t14 * t16;
	t29 = t16 * t17;
	t26 = t10 * t39;
	t25 = t15 * t28;
	t5 = -t18 * t14 - t17 * t25;
	t22 = t10 * t36 + t12 * t5;
	t19 = t13 * t27 - t4 * t16 + t3 * t34;
	t6 = -t14 * t25 + t18 * t17;
	t2 = t22 * t13 + t6 * t16;
	t1 = -t6 * t13 + t22 * t16;
	t7 = [t19 * r_i_i_C(1) + t38 * r_i_i_C(2) - t4 * pkin(2) - t15 * pkin(1) + pkin(9) * t35 + t39 * (-t3 * t10 + t12 * t35), (t5 * t16 - t6 * t34) * r_i_i_C(1) + (-t5 * t13 - t6 * t33) * r_i_i_C(2) + t5 * pkin(2) + t6 * t26, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0; t18 * pkin(1) + t6 * pkin(2) + pkin(9) * t36 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t39 * (-t5 * t10 + t12 * t36), (-t3 * t16 - t4 * t34) * r_i_i_C(1) + (t3 * t13 - t4 * t33) * r_i_i_C(2) - t3 * pkin(2) + t4 * t26, -t38 * r_i_i_C(1) + t19 * r_i_i_C(2), 0, 0, 0; 0, ((-t12 * t32 + t29) * r_i_i_C(1) + (-t12 * t30 - t31) * r_i_i_C(2) + t17 * pkin(2) + t14 * t26) * t11, ((t12 * t29 - t32) * r_i_i_C(1) + (-t12 * t31 - t30) * r_i_i_C(2)) * t11 + (t16 * r_i_i_C(1) - t13 * r_i_i_C(2)) * t10 * t28, 0, 0, 0;];
	Ja_transl = t7;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:49:58
	% EndTime: 2019-10-10 12:49:58
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (242->73), mult. (663->129), div. (0->0), fcn. (854->12), ass. (0->49)
	t35 = sin(qJ(2));
	t36 = sin(qJ(1));
	t39 = cos(qJ(2));
	t40 = cos(qJ(1));
	t48 = cos(pkin(6));
	t44 = t40 * t48;
	t22 = t36 * t35 - t39 * t44;
	t30 = sin(pkin(7));
	t32 = cos(pkin(7));
	t31 = sin(pkin(6));
	t55 = t31 * t40;
	t15 = -t22 * t30 + t32 * t55;
	t33 = sin(qJ(4));
	t37 = cos(qJ(4));
	t23 = t35 * t44 + t36 * t39;
	t34 = sin(qJ(3));
	t38 = cos(qJ(3));
	t46 = t30 * t55;
	t54 = t32 * t34;
	t6 = t22 * t54 - t23 * t38 + t34 * t46;
	t62 = t15 * t37 - t6 * t33;
	t61 = t15 * t33 + t6 * t37;
	t60 = r_i_i_C(3) + pkin(11);
	t59 = t30 * pkin(10);
	t58 = t30 * t33;
	t57 = t30 * t37;
	t56 = t31 * t36;
	t53 = t32 * t38;
	t52 = t34 * t35;
	t51 = t34 * t39;
	t50 = t35 * t38;
	t49 = t38 * t39;
	t47 = t30 * t56;
	t45 = t36 * t48;
	t43 = t48 * t30;
	t42 = t37 * r_i_i_C(1) - t33 * r_i_i_C(2) + pkin(3);
	t24 = -t40 * t35 - t39 * t45;
	t17 = -t24 * t30 + t32 * t56;
	t41 = -t23 * t34 + (-t22 * t32 - t46) * t38;
	t25 = -t35 * t45 + t40 * t39;
	t21 = -t31 * t39 * t30 + t48 * t32;
	t14 = t34 * t43 + (t32 * t51 + t50) * t31;
	t12 = t24 * t38 - t25 * t54;
	t10 = -t22 * t38 - t23 * t54;
	t8 = t25 * t38 + (t24 * t32 + t47) * t34;
	t7 = -t24 * t53 + t25 * t34 - t38 * t47;
	t2 = t17 * t33 + t8 * t37;
	t1 = t17 * t37 - t8 * t33;
	t3 = [-t36 * pkin(1) - t23 * pkin(2) + t6 * pkin(3) + pkin(9) * t55 + t15 * pkin(10) + t61 * r_i_i_C(1) + t62 * r_i_i_C(2) + t60 * t41, (t12 * t37 + t25 * t58) * r_i_i_C(1) + (-t12 * t33 + t25 * t57) * r_i_i_C(2) + t12 * pkin(3) + t24 * pkin(2) + t25 * t59 + t60 * (t24 * t34 + t25 * t53), -t42 * t7 + t60 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0; t40 * pkin(1) + t25 * pkin(2) + t8 * pkin(3) + pkin(9) * t56 + t17 * pkin(10) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t60 * t7, (t10 * t37 + t23 * t58) * r_i_i_C(1) + (-t10 * t33 + t23 * t57) * r_i_i_C(2) + t10 * pkin(3) - t22 * pkin(2) + t23 * t59 + t60 * (-t22 * t34 + t23 * t53), t42 * t41 - t6 * t60, -t62 * r_i_i_C(1) + t61 * r_i_i_C(2), 0, 0; 0, (t42 * (-t32 * t52 + t49) + t39 * pkin(2) + (t33 * r_i_i_C(1) + t37 * r_i_i_C(2) + pkin(10)) * t35 * t30 + t60 * (t32 * t50 + t51)) * t31, t60 * t14 + t42 * (t38 * t43 + (t32 * t49 - t52) * t31), (-t14 * t33 + t21 * t37) * r_i_i_C(1) + (-t14 * t37 - t21 * t33) * r_i_i_C(2), 0, 0;];
	Ja_transl = t3;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:49:58
	% EndTime: 2019-10-10 12:49:58
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (338->74), mult. (792->125), div. (0->0), fcn. (1017->14), ass. (0->51)
	t41 = sin(qJ(4));
	t69 = t41 * pkin(4) + pkin(10);
	t45 = cos(qJ(4));
	t33 = t45 * pkin(4) + pkin(3);
	t36 = qJ(4) + pkin(13);
	t34 = sin(t36);
	t35 = cos(t36);
	t50 = t35 * r_i_i_C(1) - t34 * r_i_i_C(2) + t33;
	t68 = t34 * r_i_i_C(1) + t35 * r_i_i_C(2) + t69;
	t66 = r_i_i_C(3) + qJ(5) + pkin(11);
	t43 = sin(qJ(2));
	t44 = sin(qJ(1));
	t47 = cos(qJ(2));
	t48 = cos(qJ(1));
	t56 = cos(pkin(6));
	t52 = t48 * t56;
	t24 = t43 * t52 + t44 * t47;
	t42 = sin(qJ(3));
	t65 = t24 * t42;
	t38 = sin(pkin(6));
	t64 = t38 * t44;
	t63 = t38 * t48;
	t39 = cos(pkin(7));
	t62 = t39 * t42;
	t46 = cos(qJ(3));
	t61 = t39 * t46;
	t60 = t42 * t43;
	t59 = t42 * t47;
	t58 = t43 * t46;
	t57 = t46 * t47;
	t37 = sin(pkin(7));
	t55 = t37 * t64;
	t54 = t37 * t63;
	t53 = t44 * t56;
	t51 = t56 * t37;
	t23 = t44 * t43 - t47 * t52;
	t15 = -t23 * t37 + t39 * t63;
	t25 = -t48 * t43 - t47 * t53;
	t17 = -t25 * t37 + t39 * t64;
	t6 = t23 * t62 - t24 * t46 + t42 * t54;
	t49 = t68 * t37;
	t26 = -t43 * t53 + t48 * t47;
	t22 = -t38 * t47 * t37 + t56 * t39;
	t14 = t42 * t51 + (t39 * t59 + t58) * t38;
	t13 = -t46 * t51 + (-t39 * t57 + t60) * t38;
	t8 = t26 * t46 + (t25 * t39 + t55) * t42;
	t7 = -t25 * t61 + t26 * t42 - t46 * t55;
	t3 = t23 * t61 + t46 * t54 + t65;
	t2 = t17 * t34 + t8 * t35;
	t1 = t17 * t35 - t8 * t34;
	t4 = [-t24 * pkin(2) - t44 * pkin(1) + pkin(9) * t63 + t66 * (-t65 + (-t23 * t39 - t54) * t46) + t50 * t6 + t68 * t15, t25 * pkin(2) + t50 * (t25 * t46 - t26 * t62) + t26 * t49 + t66 * (t25 * t42 + t26 * t61), -t50 * t7 + t66 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2) + (t17 * t45 - t41 * t8) * pkin(4), t7, 0; t48 * pkin(1) + t26 * pkin(2) + pkin(9) * t64 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t69 * t17 + t8 * t33 + t66 * t7, -t23 * pkin(2) + t66 * (-t23 * t42 + t24 * t61) + t50 * (-t23 * t46 - t24 * t62) + t24 * t49, -t50 * t3 - t6 * t66, (-t15 * t35 + t6 * t34) * r_i_i_C(1) + (t15 * t34 + t6 * t35) * r_i_i_C(2) + (-t15 * t45 + t41 * t6) * pkin(4), t3, 0; 0, (t66 * (t39 * t58 + t59) + t50 * (-t39 * t60 + t57) + t47 * pkin(2) + t43 * t49) * t38, -t50 * t13 + t66 * t14, (-t14 * t34 + t22 * t35) * r_i_i_C(1) + (-t14 * t35 - t22 * t34) * r_i_i_C(2) + (-t14 * t41 + t22 * t45) * pkin(4), t13, 0;];
	Ja_transl = t4;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:49:58
	% EndTime: 2019-10-10 12:49:59
	% DurationCPUTime: 0.76s
	% Computational Cost: add. (700->114), mult. (1612->190), div. (0->0), fcn. (2102->16), ass. (0->75)
	t56 = sin(qJ(2));
	t59 = cos(qJ(2));
	t60 = cos(qJ(1));
	t78 = cos(pkin(6));
	t71 = t60 * t78;
	t86 = sin(qJ(1));
	t37 = t86 * t56 - t59 * t71;
	t38 = t56 * t71 + t86 * t59;
	t55 = sin(qJ(3));
	t51 = cos(pkin(7));
	t87 = cos(qJ(3));
	t74 = t51 * t87;
	t49 = sin(pkin(7));
	t50 = sin(pkin(6));
	t82 = t50 * t60;
	t76 = t49 * t82;
	t15 = t37 * t74 + t38 * t55 + t87 * t76;
	t81 = t51 * t55;
	t16 = -t37 * t81 + t38 * t87 - t55 * t76;
	t30 = -t37 * t49 + t51 * t82;
	t48 = qJ(4) + pkin(13);
	t46 = sin(t48);
	t47 = cos(t48);
	t4 = t16 * t47 - t30 * t46;
	t53 = sin(qJ(6));
	t57 = cos(qJ(6));
	t95 = -t15 * t57 + t4 * t53;
	t94 = -t15 * t53 - t4 * t57;
	t54 = sin(qJ(4));
	t91 = pkin(4) * t54 + pkin(10);
	t58 = cos(qJ(4));
	t45 = t58 * pkin(4) + pkin(3);
	t89 = r_i_i_C(3) + pkin(12);
	t90 = t89 * t46 + t45;
	t85 = t46 * t49;
	t84 = t47 * t49;
	t83 = t49 * t50;
	t80 = t55 * t56;
	t79 = t55 * t59;
	t77 = t56 * t83;
	t75 = t50 * t86;
	t73 = t87 * t56;
	t72 = t87 * t59;
	t70 = t78 * t49;
	t69 = t49 * t75;
	t68 = t91 * t49;
	t67 = t78 * t86;
	t66 = t57 * r_i_i_C(1) - t53 * r_i_i_C(2) + pkin(5);
	t52 = -qJ(5) - pkin(11);
	t65 = t53 * r_i_i_C(1) + t57 * r_i_i_C(2) - t52;
	t64 = t60 * t56 + t59 * t67;
	t63 = t64 * t87;
	t62 = -t66 * t47 - t90;
	t61 = t64 * t49 + t51 * t75;
	t39 = -t56 * t67 + t60 * t59;
	t36 = t78 * t51 - t59 * t83;
	t35 = (-t51 * t80 + t72) * t50;
	t34 = (t51 * t73 + t79) * t50;
	t29 = t55 * t70 + (t51 * t79 + t73) * t50;
	t28 = -t87 * t70 + (-t51 * t72 + t80) * t50;
	t26 = -t39 * t81 - t63;
	t25 = t39 * t74 - t64 * t55;
	t24 = -t37 * t87 - t38 * t81;
	t23 = -t37 * t55 + t38 * t74;
	t22 = t35 * t47 + t46 * t77;
	t20 = t39 * t87 + (-t64 * t51 + t69) * t55;
	t19 = t39 * t55 + t51 * t63 - t87 * t69;
	t14 = t29 * t47 + t36 * t46;
	t12 = t26 * t47 + t39 * t85;
	t10 = t24 * t47 + t38 * t85;
	t8 = t20 * t47 + t61 * t46;
	t7 = t20 * t46 - t61 * t47;
	t2 = t19 * t53 + t8 * t57;
	t1 = t19 * t57 - t8 * t53;
	t3 = [t94 * r_i_i_C(1) + t95 * r_i_i_C(2) - t4 * pkin(5) + t15 * t52 - t38 * pkin(2) - t86 * pkin(1) + pkin(9) * t82 - t90 * t16 + (-t89 * t47 + t91) * t30, (t12 * t57 + t25 * t53) * r_i_i_C(1) + (-t12 * t53 + t25 * t57) * r_i_i_C(2) + t12 * pkin(5) + t26 * t45 - t25 * t52 - t64 * pkin(2) + t39 * t68 + t89 * (t26 * t46 - t39 * t84), t62 * t19 + t65 * t20, t89 * t8 + (-t20 * t54 + t61 * t58) * pkin(4) - t66 * t7, t19, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t60 * pkin(1) + t39 * pkin(2) + t8 * pkin(5) + pkin(9) * t75 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t19 * t52 + t20 * t45 + t91 * t61 + t89 * t7, (t10 * t57 + t23 * t53) * r_i_i_C(1) + (-t10 * t53 + t23 * t57) * r_i_i_C(2) + t10 * pkin(5) + t24 * t45 - t23 * t52 - t37 * pkin(2) + t89 * (t24 * t46 - t38 * t84) + t38 * t68, t62 * t15 + t65 * t16, t89 * t4 + (-t16 * t54 - t30 * t58) * pkin(4) + t66 * (-t16 * t46 - t30 * t47), t15, -t95 * r_i_i_C(1) + t94 * r_i_i_C(2); 0, (t22 * t57 + t34 * t53) * r_i_i_C(1) + (-t22 * t53 + t34 * t57) * r_i_i_C(2) + t22 * pkin(5) + t35 * t45 - t34 * t52 + t89 * (t35 * t46 - t47 * t77) + (t59 * pkin(2) + t56 * t68) * t50, t62 * t28 + t65 * t29, t89 * t14 + (-t29 * t54 + t36 * t58) * pkin(4) + t66 * (-t29 * t46 + t36 * t47), t28, (-t14 * t53 + t28 * t57) * r_i_i_C(1) + (-t14 * t57 - t28 * t53) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end